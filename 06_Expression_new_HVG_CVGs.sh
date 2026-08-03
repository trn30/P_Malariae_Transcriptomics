#!/bin/bash
# ==============================================================================
# Pipeline: Identification and Analysis of Hypervariable Genes (HVGs)
# Description: Processes multigene families, uses genomic locations to find
# potential HVGs via BLAST, and validates expression using bulk RNA-Seq and
# Single-Cell RNA-Seq (Malaria Cell Atlas).
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status

# ==============================================================================
# 1. DIRECTORY SETUP
# ==============================================================================
DATA_DIR="data"                   # Put all your input tables and fasta here
MCA_DIR="data/mca"                # Put Malaria Cell Atlas files here
RESULTS_DIR="results"             # Output directory
TMP_DIR="tmp"                     # Temporary files for BED/BLAST

mkdir -p ${RESULTS_DIR} ${TMP_DIR}

echo "===  CLONALLY VARIANT GENES AND MULTIGENE FAMILIES (R) ==="

Rscript - <<'EOF'
# Load libraries
(library(dplyr))

# Define paths
DATA_DIR <- "data"
RESULTS_DIR <- "results"

# ------------------------------------------------------------------------------
# 1. Extraction of Multigene Families (Known Families)
# ------------------------------------------------------------------------------
genes_pmal <- read.delim(file.path(DATA_DIR, "Genes_info_Pmalariae_rel68.txt"))

# Extract specific gene families and filter out pseudogenes where necessary
pir_genes <- genes_pmal[grep("PIR protein", genes_pmal$Product.Description), ]
pir_filt <- pir_genes[-grep("PIR protein,", pir_genes$Product.Description), ]

famm_genes <- genes_pmal[grep("fam-m protein", genes_pmal$Product.Description), ]
famm_filt <- famm_genes[-grep("fam-m protein,", famm_genes$Product.Description), ]

faml_genes <- genes_pmal[grep("fam-l protein", genes_pmal$Product.Description), ]
faml_filt <- faml_genes[-grep("fam-l protein,", faml_genes$Product.Description), ]

stp1_genes <- genes_pmal[grep("STP1 protein", genes_pmal$Product.Description), ]
stp1_filt <- stp1_genes[-grep("STP1 protein,", stp1_genes$Product.Description), ]

cys_genes <- genes_pmal[grep("6-cysteine protein", genes_pmal$Product.Description), ]
cys_filt <- cys_genes[-grep("6-cysteine protein, pseudogene", cys_genes$Product.Description), ]

phist_genes <- genes_pmal[grep("(PHIST)", genes_pmal$Product.Description), ]
phist_filt <- phist_genes[-grep("pseudogene", phist_genes$Product.Description), ]
phist_filt <- phist_filt[-grep("PHISTc", phist_filt$Product.Description), ]

fikk_genes <- genes_pmal[grep("FIKK family", genes_pmal$Product.Description), ]

# Combine and save
multigene_families <- rbind(pir_filt, famm_filt, faml_filt, stp1_filt, cys_filt, phist_filt, fikk_genes)
write.table(multigene_families, file.path(RESULTS_DIR, "multigene_families_from_gff.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ------------------------------------------------------------------------------
# 2. Find Potential HVGs using Product Descriptions and Homology
# ------------------------------------------------------------------------------
pmal_info <- read.delim(file.path(DATA_DIR, "Pmalariae_plasmoDB_68_information.txt"), header=TRUE)

# Group by product description to find multicopy genes
pmal_grouped <- pmal_info %>% group_by(Product.Description) %>% summarise(n = n())
pmal_multicopy <- pmal_grouped[pmal_grouped$n != 1, ]

# Manual selection of indices known to be in tandem arrays (based on user criteria)
target_indices <- c(10,13,14,15,17,20,21,24,30,31,38,43,44,46,47,52,53,56,60,64,65,66,69,
                    74,77,78,80,81,82,83,84,85,87,88,89,93,94,95,96,97,100,101,103,104,
                    107,109,110,117,121,122,123,128,138,139,141,148,149,150,151,153,154,
                    155,158,159,164,165,166,168,170,174,176,177,181,182,185,186,187,189,
                    191,194,196,198,204,205,208,211,215,219,220,221,222,226,228,233,237,
                    238,243,249,254,259,260,261,269,270,275,277)

pmal_tandem <- pmal_multicopy[target_indices, ]
write.table(pmal_tandem, file.path(RESULTS_DIR, "Master_list_multicopy_all_genes_pseudo.txt"),
            sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

message(" R preprocessing finished.")
EOF

echo "===  GENOMIC LOCATION & BLAST ANALYSIS (BASH) ==="

# Clean formatting of the provided list
dos2unix ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt
awk -F'\t' 'BEGIN{OFS=FS} {gsub(" ", "_", $4); print}' ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt > ${TMP_DIR}/Master_list_multicopy_tmp.txt
mv ${TMP_DIR}/Master_list_multicopy_tmp.txt ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt

# Create BED file
cat ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo_genomic_location_transcript_plasmoDB.txt | \
    sed 's*,*_*g; s, ,_,g; s,__,_,g' | \
    sed '1d; s,:,\t,g; s/\./\t/g; s,(,\t,g; s,),,g; s,\t\t,\t,g' | \
    awk '{print $4"\t"$5"\t"$6"\t"$1"\t0\t"$7}' > ${TMP_DIR}/multicopy_genes.bed

# Extract sequences using BEDTools
bedtools getfasta -fi ${DATA_DIR}/PlasmoDB-68_PmalariaeUG01_Genome.fasta -bed ${TMP_DIR}/multicopy_genes.bed -s -name > ${TMP_DIR}/multicopy_genes.fasta

# Build BLAST DB and run BLASTN (Identity 80%)
makeblastdb -in ${TMP_DIR}/multicopy_genes.fasta -dbtype nucl -out ${TMP_DIR}/multicopy_genes -title multicopy_genes

blastn -query ${TMP_DIR}/multicopy_genes.fasta -db ${TMP_DIR}/multicopy_genes -num_threads 10 \
    -evalue 0.0001 -perc_identity 80.0 \
    -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps nident qseq sseq qlen slen sstrand" \
    -out ${TMP_DIR}/multicopy_genes_blast.txt

# Add headers to BLAST output
echo -e "Sequence_ID\tContig_ID\tPercentage_identical_matches\tAlignment_length\tNumber_mismatches\tNumber_gap_openings\tStart_alignment_query\tEnd_alignment_query\tStart_alignment_subject\tEnd_alignment_subject\tEvalue\tBit_score\tNumber_total_gaps\tNumber_identical_matches\tSequence_query\tSequence_subject\tSequence_query_length\tSequence_subject_length\tStrand" > ${RESULTS_DIR}/multicopy_genes_blast_final.txt
egrep -v "^#" ${TMP_DIR}/multicopy_genes_blast.txt >> ${RESULTS_DIR}/multicopy_genes_blast_final.txt

echo "BLAST process completed. Extracting hits..."

# Extract summary of hits matching family members
> ${RESULTS_DIR}/summary_of_blast_out_of_families.txt
for id in $(cut -f4 ${TMP_DIR}/multicopy_genes.bed | sort); do
    family=$(grep $id ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt | cut -f4 | sed 's/[[:cntrl:]]//g')
    total_size=$(cut -f4 ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt | sort | uniq -c | sort -k1,1nr | grep -w "$family" | head -1 | awk '{print $1}')
    rest_ids=$(grep "$family" ${DATA_DIR}/Master_list_multicopy_all_genes_pseudo.txt | grep -v $id | cut -f1 | tr '\n' '|' | sed 's,|$,,g' | sed 's/[[:cntrl:]]//g')

    if [ ! -z "$rest_ids" ]; then
        hits=$(awk -F '\t' '$1 ~ /'"$id"'/' ${RESULTS_DIR}/multicopy_genes_blast_final.txt | awk -F '\t' '$1 != $2' | egrep "$rest_ids" | grep -v -c '^$' || true)
        mean_ident=$(awk -F '\t' '$1 ~ /'"$id"'/' ${RESULTS_DIR}/multicopy_genes_blast_final.txt | awk -F '\t' '$1 != $2' | egrep "$rest_ids" | awk '{sum+=$3; count++} END {if (count > 0) print sum/count}' || echo "0")
        echo "$id has $hits blast hits out of $total_size with $mean_ident mean identity" >> ${RESULTS_DIR}/summary_of_blast_out_of_families.txt
    fi
done

grep -v "has 0 blast hits" ${RESULTS_DIR}/summary_of_blast_out_of_families.txt > ${RESULTS_DIR}/summary_of_blast_families_final.txt

echo "===  SINGLE-CELL ANALYSIS (MALARIA CELL ATLAS) (R) ==="

Rscript - <<'EOF'

    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(Hmisc)
    library(data.table)

MCA_DIR <- "data/mca"
RES_DIR <- "results"

# ------------------------------------------------------------------------------
# HELPER FUNCTION: Process Seurat Data & Extract HVGs
# Eliminates code duplication for processing different stages/species
# ------------------------------------------------------------------------------
analyze_sc_stage <- function(counts_mat, metadata, stage_col, stage_val, num_features, prefix) {
    message(paste("Processing Stage:", stage_val))

    # Create Seurat Object
    seu <- CreateSeuratObject(counts = counts_mat, meta.data = metadata)
    Idents(seu) <- seu[[stage_col]]

    # Subset by stage
    seu_sub <- subset(seu, idents = stage_val)
    seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst",
                                    mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf),
                                    nfeatures = num_features)

    hvg_seurat <- VariableFeatures(seu_sub)

    # Prepare dataframe for custom IQR calculation
    valid_cells <- rownames(metadata)[metadata[[stage_col]] == stage_val]
    # Replace dashes with underscores to match rownames if necessary
    hvg_query <- sub("-", "_", hvg_seurat)
    expr_df <- as.data.frame(counts_mat[hvg_query, valid_cells, drop=FALSE])
    expr_df$gene <- rownames(expr_df)

    df_long <- pivot_longer(expr_df, cols = -gene, names_to = "sample", values_to = "expression")

    # Fast IQR Calculation using Boxplot stats logic (equivalent to user's loop)
    df_iqr <- df_long %>%
        group_by(gene) %>%
        summarise(
            # boxplot stats[5] - stats[1] roughly equates to max non-outlier range
            IQR_val = diff(boxplot.stats(expression)$stats[c(1, 5)])
        ) %>%
        filter(IQR_val > 0) %>%
        arrange(desc(IQR_val))

    # Assign ranks
    df_iqr$rank <- as.character(cut2(df_iqr$IQR_val, g=5))
    rank_names <- sort(unique(df_iqr$rank))
    rank_map <- setNames(c("Very_Low", "Low", "Medium", "High", "Very_High"), rank_names)
    df_iqr$rank <- unname(rank_map[df_iqr$rank])

    # Extract High/Very High Genes (IQR Approach)
    hvg_iqr <- df_iqr$Gene_ID[df_iqr$rank %in% c("Very_High", "High")]

    # Generate Plot
    plot_df <- df_long %>% filter(gene %in% hvg_iqr)
    plot_df$gene <- factor(plot_df$gene, levels = rev(hvg_iqr))

    pdf(file.path(RES_DIR, paste0(prefix, "_top_iqr.pdf")))
    p <- ggplot(plot_df, aes(x = gene, y = expression, fill = sample)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
        labs(x = "Gene", y = "Expression", title = paste("HVGs in", stage_val)) +
        theme_minimal() + theme(legend.position = "none") + coord_flip()
    print(p)
    dev.off()

    # Return both lists to find intersection later
    return(list(
        iqr_hvgs = hvg_iqr,
        seurat_hvgs = hvg_seurat
    ))
}

# ------------------------------------------------------------------------------
# 1. P.falciparum MCA Analysis
# ------------------------------------------------------------------------------
pf_meta <- fread(file.path(MCA_DIR, "pf-ss2-set1-ss2-data.csv")) %>% select(1, 9, 10) %>% as.data.frame()
pf_exp <- fread(file.path(MCA_DIR, "pf-ss2-set1-ss2-exp.csv"))

pf_counts <- as.matrix(pf_exp[,-1])
rownames(pf_counts) <- pf_exp$V1
pf_meta$STAGE_HR <- gsub(" ", "_", pf_meta$STAGE_HR)
rownames(pf_meta) <- pf_meta$CELL_ID

# Run stages using function
pf_ookinete <- analyze_sc_stage(pf_counts, pf_meta, "STAGE_HR", "ookinete", 5720, "Pf_ookinete")
pf_sporo_ooc <- analyze_sc_stage(pf_counts, pf_meta, "STAGE_HR", "sporozoite_(oocyst)", 5720, "Pf_sporo_oocyst")
pf_sporo_sgl <- analyze_sc_stage(pf_counts, pf_meta, "STAGE_HR", "sporozoite_(salivary_gland)", 5720, "Pf_sporo_salivary")
pf_sporo_all <- analyze_sc_stage(pf_counts, pf_meta, "STAGE_LR", "sporozoite", 5720, "Pf_sporo_all")

# Combine results
pf_final_hvgs <- unique(c(
    pf_ookinete$iqr_hvgs, pf_ookinete$seurat_hvgs,
    pf_sporo_ooc$iqr_hvgs, pf_sporo_ooc$seurat_hvgs,
    pf_sporo_sgl$iqr_hvgs, pf_sporo_sgl$seurat_hvgs
))

# ------------------------------------------------------------------------------
# 2. P.berghei MCA Analysis
# ------------------------------------------------------------------------------
pb_meta <- fread(file.path(MCA_DIR, "pb-ss2-set1-ss2-data.csv")) %>% select(1, 9, 10) %>% as.data.frame()
pb_exp <- fread(file.path(MCA_DIR, "pb-ss2-set1-ss2-exp.csv"))

pb_counts <- as.matrix(pb_exp[,-1])
rownames(pb_counts) <- pb_exp$V1
pb_meta$STAGE_HR <- gsub(" ", "_", pb_meta$STAGE_HR)
rownames(pb_meta) <- pb_meta$CELL_ID

# Run stages using function
pb_ookinete <- analyze_sc_stage(pb_counts, pb_meta, "STAGE_LR", "ookinete", 5254, "Pb_ookinete")
pb_oocyst <- analyze_sc_stage(pb_counts, pb_meta, "STAGE_LR", "oocyst", 5254, "Pb_oocyst")
pb_sporo_sgl <- analyze_sc_stage(pb_counts, pb_meta, "STAGE_HR", "sporozoite_(salivary_gland)", 5254, "Pb_sporo_salivary")
pb_sporo_all <- analyze_sc_stage(pb_counts, pb_meta, "STAGE_LR", "sporozoite", 5254, "Pb_sporo_all")

message("MCA Single-Cell Feature Analysis Complete. Check 'results/' folder for PDFs.")
EOF
