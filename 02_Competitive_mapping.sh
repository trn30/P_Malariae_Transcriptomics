#!/bin/bash
# ==============================================================================
# Complete Pipeline: Competitive Mapping, Filtering, Extraction, and Plotting
# Description: Maps RNA-seq reads to a combined reference, filters by quality
# and size, extracts reads by species, and finally generates abundance plots
# ==============================================================================

# Stop the script if any command fails
set -e

# ==============================================================================
# 1. SETUP
# ==============================================================================
CORES=30

# Input and output directories
REF_DIR="refs_comp"
REF_FASTA="${REF_DIR}/Plasmodium_vector.fasta"
INDEX_PREFIX="${REF_DIR}/Plasmodium_vector_hisat2_index"
RAW_DIR="sortmerna_raw/output"

OUT_DIR="out_hisat2_Plamo5_vector"
FINAL_DIR="${OUT_DIR}/final_bam"
FINAL_FILT_DIR="${FINAL_DIR}/final_filt"
TMP_DIR="tmp_picard"

# Samples
SAMPLES=("SRR24471066" "SRR24471067" "SRR24471068" "SRR24471069" "SRR24471070")

# Annotation files per species (Associative array)
declare -A GFF_FILES=(
    ["Pm"]="PlasmoDB-68_PmalariaeUG01.coco_biotype.gtf"
    ["Pf"]="PlasmoDB-68_Pfalciparum3D7.gff"
    ["Pv"]="PlasmoDB-68_PvivaxP01.gff"
    ["Poc"]="PlasmoDB-68_PovalecurtisiGH01.gff"
    ["Pow"]="PlasmoDB-68_PovalewallikeriPowCR01.gff"
    ["Agam"]="VectorBase-68_AgambiaePEST.gff"
)

# Create necessary directories
mkdir -p ${FINAL_DIR} ${FINAL_FILT_DIR} ${TMP_DIR}

echo "=== STARTING PIPELINE ==="

# ==============================================================================
# 2. MAPPING AND FILTERING (HISAT2 + PICARD + SAMTOOLS)
# ==============================================================================
echo "Building HISAT2 index..."
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    hisat2-build -p ${CORES} ${REF_FASTA} ${INDEX_PREFIX}
else
    echo "Index already exists. Skipping step..."
fi

for srr in "${SAMPLES[@]}"; do
    echo "================ Processing sample: ${srr} ================"
    SAMPLE_OUT="${OUT_DIR}/${srr}_reads"
    mkdir -p ${SAMPLE_OUT}

    echo "[1/4] Competitive mapping with HISAT2..."
    hisat2 -t -p ${CORES} --very-sensitive \
        --summary-file ${SAMPLE_OUT}/${srr}_hisat2_mapping.log \
        --un-gz ${SAMPLE_OUT}/${srr}_un.fq.gz \
        --al-gz ${SAMPLE_OUT}/${srr}_al.fq.gz \
        --un-conc-gz ${SAMPLE_OUT}/${srr}_un.fq.gz \
        --al-conc-gz ${SAMPLE_OUT}/${srr}_al_conc.fq.gz \
        -x ${INDEX_PREFIX} \
        -1 ${RAW_DIR}/${srr}_sortmerna_R1.fastq.gz \
        -2 ${RAW_DIR}/${srr}_sortmerna_R2.fastq.gz | \
    samtools sort -@ ${CORES} -m 2G -l 9 --write-index \
        -o ${SAMPLE_OUT}/${srr}_hisat2.bam##idx##${SAMPLE_OUT}/${srr}_hisat2.bam.bai -

    echo "[2/4] Marking duplicates (Picard) and filtering by quality..."
    picard MarkDuplicates \
        I=${SAMPLE_OUT}/${srr}_hisat2.bam \
        O=/dev/stdout \
        M=${FINAL_DIR}/${srr}_dedup_metrics.txt \
        TMP_DIR=${TMP_DIR} \
        REMOVE_SEQUENCING_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate | \
    samtools view -@ ${CORES} -f 2 -F 1024 -q 60 -b | \
    samtools sort -@ ${CORES} --write-index \
        -o ${FINAL_DIR}/${srr}_tmp.bam##idx##${FINAL_DIR}/${srr}_tmp.bam.bai -

    echo "[3/4] Filtering alignments >= 100 bp..."
    samtools view -h ${FINAL_DIR}/${srr}_tmp.bam | \
    perl -nle '
        if (/^@/) { print }
        else {
            @ar=split(/\t/);
            if ($ar[5] =~ /(\d{3,})M/) { print }
        }
    ' | \
    samtools view -@ ${CORES} -b --write-index \
        -o ${FINAL_DIR}/${srr}_final.bam##idx##${FINAL_DIR}/${srr}_final.bam.bai -

    echo "[4/4] Fixing mates and final sorting..."
    samtools fixmate -@ ${CORES} -r -u <(samtools sort -n ${FINAL_DIR}/${srr}_final.bam) - | \
    samtools sort -@ ${CORES} - | \
    samtools view -@ ${CORES} -b -f 2 --write-index \
        -o ${FINAL_FILT_DIR}/${srr}_fixmate.bam##idx##${FINAL_FILT_DIR}/${srr}_fixmate.bam.bai -

    # Clean temporary files to save storage
    rm -f ${FINAL_DIR}/${srr}_tmp.bam ${FINAL_DIR}/${srr}_tmp.bam.bai
done

# ==============================================================================
# 3. SPECIES-SPECIFIC READ EXTRACTION
# ==============================================================================
echo "================ SPECIES EXTRACTION ================"
BAMS=$(ls ${FINAL_FILT_DIR}/*_fixmate.bam)

for species_id in "${!GFF_FILES[@]}"; do
    annot_file="${REF_DIR}/${GFF_FILES[$species_id]}"
    contig_file="${REF_DIR}/${species_id}_contigs.txt"

    echo "--- Extracting contigs and reads for: ${species_id} ---"

    # Extract unique contigs
    awk '{print $1}' "${annot_file}" | sort -u > "${contig_file}"

    for bam in ${BAMS}; do
        base_name=$(basename ${bam} _fixmate.bam)
        out_bam="${FINAL_FILT_DIR}/${base_name}_${species_id}_only.bam"

        # Filter reads mapped to specific contigs
        samtools view -@ ${CORES} -h "${bam}" | \
        grep -F -f "${contig_file}" | \
        samtools view -@ ${CORES} -b -o "${out_bam}"

        # Index BAM
        samtools index -@ ${CORES} "${out_bam}"

        # Save flagstats
        echo "### ${out_bam}" >> "${FINAL_FILT_DIR}/flagstat_${species_id}.txt"
        samtools flagstat -@ ${CORES} "${out_bam}" >> "${FINAL_FILT_DIR}/flagstat_${species_id}.txt"
        echo "" >> "${FINAL_FILT_DIR}/flagstat_${species_id}.txt"
    done
done

# ==============================================================================
# 4. PLOTTING
# ==============================================================================
echo "================ GENERATING PLOTS WITH R ================"
# Using <<'EOF' prevents Bash from evaluating variables ($) inside the R script block
Rscript - <<'EOF'
# Load libraries quietly

  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(paletteer)
  library(gridExtra)

# Create output folder
dir.create("Figures", showWarnings = FALSE)

# Data Preparation
df_parasites <- data.frame(
  Species = c("P. falciparum", "P. malariae", "P. ovale curtisi", "P. ovale wallikeri", "P. vivax"),
  SG_I3 = c(4, 103080, 1147, 602, 189),
  SG_I2 = c(4, 113385, 1035, 543, 135),
  SG_I1 = c(262, 298718, 1610, 1457, 1002),
  MG_I2 = c(6, 28832, 1001, 532, 187),
  MG_I1 = c(216, 213976, 1555, 1232, 1067)
)

# Long format & Percentage calculation
df_long <- df_parasites %>%
  pivot_longer(cols = -Species, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(Percentage = (Reads / sum(Reads)) * 100) %>%
  ungroup()

# Set Color Palette
color_palette <- paletteer_c("grDevices::Fall", 5)

# Plot 1: 100% Stacked Bar Chart
p1 <- ggplot(df_long, aes(x = Sample, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", color = "white", width = 0.7) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  labs(title = "Presence of Plasmodium species in samples (%)",
       y = "Percentage (%)", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("Figures/plasmodium_species_abundance_stacked.pdf", plot = p1, width = 12, height = 6)

# Plot 2: Log-Scaled Bar Chart
p2 <- ggplot(df_long, aes(x = Species, y = Reads, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_log10(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~Sample, ncol = 5) +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  labs(title = "Read Count of Plasmodium species across samples",
       y = "Read count (Log10)", x = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

ggsave("Figures/plasmodium_species_abundance_log.pdf", plot = p2, width = 12, height = 6)

message("Plots successfully generated in the 'Figures' folder.")
EOF
