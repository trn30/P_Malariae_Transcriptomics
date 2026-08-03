#!/bin/bash
# ==============================================================================
# Pipeline: Plasmodium malariae RNA-Seq Processing, Quantification & DESeq2
# Description: Decompresses raw data, filters taxonomy (KrakenUniq), removes
# rRNA (SortMeRNA), maps to reference (HISAT2), quantifies (CoCo), and
# performs Differential Expression and Ortholog Analysis using R.
# ==============================================================================

set -e

# ==============================================================================
# 1. SETUP
# ==============================================================================
CORES=60
export _JAVA_OPTIONS='-Xms5g -Xmx240g'

# Directories
DATA_DIR="data"
REF_DIR="data/references"
RES_DIR="results"
BAM_DIR="${RES_DIR}/bams"
COCO_DIR="${RES_DIR}/coco_out"
QC_DIR="${RES_DIR}/qc"
TMP_DIR="tmp"

mkdir -p ${BAM_DIR} ${COCO_DIR} ${QC_DIR} ${TMP_DIR}

# Map sample names to their specific raw file identifiers (e.g., 0058, 0059)
declare -A SAMPLES=(
    ["SG_I1"]="0058"
    ["MG_I1"]="0059"
    ["SG_I2"]="0060"
    ["MG_I2"]="0061"
    ["SG_I3"]="0062"
)

# Reference databases
KRAKEN_EUPATH="data/krakenuniq_eupathdb48"
KRAKEN_PM_AG="data/krakenuniq_Pm_Ag"
RRNA_DB="data/rRNA_databases/smr_v4.3_sensitive_db.fasta"
HISAT2_IDX="${REF_DIR}/PlasmoDB-59_PmalariaeUG01_hisat2_index"
GTF_FILE="${REF_DIR}/PlasmoDB-59_PmalariaeUG01_utr_biotype.gtf"

echo "=== STARTING RNA-SEQ PIPELINE ==="

# ==============================================================================
# 2. RAW DATA PROCESSING, FILTERING & MAPPING
# ==============================================================================
for sample in "${!SAMPLES[@]}"; do
    file_id="${SAMPLES[$sample]}"
    echo "================ Processing Sample: ${sample} (${file_id}) ================"

    # 2.1 Decompress .dsrc files to .fastq.gz
    # Note: Using half cores for dsrc decompression to avoid memory truncation issues
    echo "[1/7] Decompressing raw data..."
    dsrc d -t$((CORES/2)) -b512 -s $(find ${DATA_DIR}/raw* -name "*.fastq.dsrc" | grep "R1_001." | grep ${file_id} | head -1) | \
        cat - <(dsrc d -t$((CORES/2)) -b512 -s $(find ${DATA_DIR}/raw* -name "*.fastq.dsrc" | grep "R1_001." | grep ${file_id} | tail -1)) | \
        pigz -p ${CORES} --fast -c > ${TMP_DIR}/${sample}_R1.fastq.gz

    dsrc d -t$((CORES/2)) -b512 -s $(find ${DATA_DIR}/raw* -name "*.fastq.dsrc" | grep "R2_001." | grep ${file_id} | head -1) | \
        cat - <(dsrc d -t$((CORES/2)) -b512 -s $(find ${DATA_DIR}/raw* -name "*.fastq.dsrc" | grep "R2_001." | grep ${file_id} | tail -1)) | \
        pigz -p ${CORES} --fast -c > ${TMP_DIR}/${sample}_R2.fastq.gz

    fastqc -t ${CORES} -o ${QC_DIR} ${TMP_DIR}/${sample}_R1.fastq.gz ${TMP_DIR}/${sample}_R2.fastq.gz

    # 2.2 Taxonomy Filtering (KrakenUniq - 2 Passes)
    echo "[2/7] Taxonomy filtering with KrakenUniq..."
    # Pass 1: Keep only Eukaryotic pathogens (removes unmapped/bacteria)
    krakenuniq --db ${KRAKEN_EUPATH} --threads ${CORES} --output ${TMP_DIR}/${sample}_kraken1 --report-file ${QC_DIR}/${sample}_kraken1.tsv --preload-size 200G --paired ${TMP_DIR}/${sample}_R1.fastq.gz ${TMP_DIR}/${sample}_R2.fastq.gz
    extract_kraken_reads.py -k ${TMP_DIR}/${sample}_kraken1 -s1 ${TMP_DIR}/${sample}_R1.fastq.gz -s2 ${TMP_DIR}/${sample}_R2.fastq.gz -t 5858 --fastq-output -o ${TMP_DIR}/${sample}_k1_R1.fq -o2 ${TMP_DIR}/${sample}_k1_R2.fq

    # Pass 2: Separate P. malariae from Anopheles gambiae
    krakenuniq --db ${KRAKEN_PM_AG} --threads ${CORES} --output ${TMP_DIR}/${sample}_kraken2 --report-file ${QC_DIR}/${sample}_kraken2.tsv --preload-size 150G --paired ${TMP_DIR}/${sample}_k1_R1.fq ${TMP_DIR}/${sample}_k1_R2.fq
    extract_kraken_reads.py -k ${TMP_DIR}/${sample}_kraken2 -s1 ${TMP_DIR}/${sample}_k1_R1.fq -s2 ${TMP_DIR}/${sample}_k1_R2.fq -t 1000000002 --fastq-output -o ${TMP_DIR}/${sample}_Pm_R1.fq -o2 ${TMP_DIR}/${sample}_Pm_R2.fq

    # 2.3 Remove rRNA (SortMeRNA)
    echo "[3/7] Removing rRNA with SortMeRNA..."
    sortmerna --idx-dir data/sortmerna_idx --ref ${RRNA_DB} --reads ${TMP_DIR}/${sample}_Pm_R1.fq --reads ${TMP_DIR}/${sample}_Pm_R2.fq --workdir ${TMP_DIR}/${sample}_smr --fastx --threads ${CORES} --out2 --aligned ${TMP_DIR}/${sample}_rRNA --other ${TMP_DIR}/${sample}_no_rRNA -v

    # 2.4 Repair paired ends
    echo "[4/7] Repairing paired ends..."
    repair.sh in=${TMP_DIR}/${sample}_no_rRNA_fwd.fq in2=${TMP_DIR}/${sample}_no_rRNA_rev.fq out=${TMP_DIR}/${sample}_clean_R1.fq out2=${TMP_DIR}/${sample}_clean_R2.fq

    # 2.5 Mapping (HISAT2)
    echo "[5/7] Mapping reads with HISAT2..."
    hisat2 -t -p ${CORES} --very-sensitive -x ${HISAT2_IDX} \
        --summary-file ${QC_DIR}/${sample}_hisat2.log \
        -1 ${TMP_DIR}/${sample}_clean_R1.fq -2 ${TMP_DIR}/${sample}_clean_R2.fq | \
        samtools sort -l 9 -m 2G -@ ${CORES} --write-index -o ${TMP_DIR}/${sample}_raw.bam##idx##${TMP_DIR}/${sample}_raw.bam.bai -

    # 2.6 Filtering (MAPQ 60, Paired, Dedup, Length >= 100)
    echo "[6/7] Filtering BAM files..."
    samtools view -@ ${CORES} -bq 60 -f 0x2 ${TMP_DIR}/${sample}_raw.bam > ${TMP_DIR}/${sample}_q60.bam
    picard MarkDuplicates -I ${TMP_DIR}/${sample}_q60.bam -O ${TMP_DIR}/${sample}_dedup.bam -M ${QC_DIR}/${sample}_dedup.txt --REMOVE_SEQUENCING_DUPLICATES true

    # Filter reads mapping >= 100bp (removes short mapping noise)
    samtools view -@ ${CORES} -h ${TMP_DIR}/${sample}_dedup.bam | perl -nle 'if (/^@/){print} else {@ar=split(/\t/); if($ar[5] =~ /(\d{3,}M)/){print}}' | \
        samtools view -@ ${CORES} -b - > ${TMP_DIR}/${sample}_len.bam

    samtools fixmate -@ ${CORES} -u <(samtools sort -n ${TMP_DIR}/${sample}_len.bam) - | \
        samtools sort -@ ${CORES} --write-index -o ${BAM_DIR}/${sample}_final.bam##idx##${BAM_DIR}/${sample}_final.bam.bai -

    samtools flagstat -@ ${CORES} ${BAM_DIR}/${sample}_final.bam > ${QC_DIR}/${sample}_flagstat.txt

    # 2.7 Quantification (CoCo)
    echo "[7/7] Quantifying reads with CoCo..."
    coco correct_count -s 0 -p -t ${CORES} ${GTF_FILE} ${BAM_DIR}/${sample}_final.bam ${COCO_DIR}/${sample}_coco

    # Cleanup tmp
    rm -rf ${TMP_DIR}/${sample}_*
done

# ==============================================================================
# 3. PCA & CORRELATION
# ==============================================================================
echo "Running DeepTools PCA and Correlation..."
multiBamSummary bins -p ${CORES} --bamfiles ${BAM_DIR}/*_final.bam -out ${QC_DIR}/MultiBamSummary.npz
plotCorrelation --corData ${QC_DIR}/MultiBamSummary.npz --corMethod spearman --plotFile ${QC_DIR}/Correlation_Spearman.pdf --whatToPlot heatmap --skipZeros --plotNumbers
plotPCA --corData ${QC_DIR}/MultiBamSummary.npz --plotFile ${QC_DIR}/PCA.pdf --outFileNameData ${QC_DIR}/PCA.tab

# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION & ORTHOLOG ANALYSIS (R)
# ==============================================================================
echo "=== RUNNING R ANALYSIS (DESeq2, HEATMAPS, ORTHOLOGS) ==="

Rscript - <<'EOF'
# Load required libraries

    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(data.table)
    library(iheatmapr)
    library(RColorBrewer)
    library(corrplot)
    library(Hmisc)
    library(dplyr)

# Paths
COCO_DIR <- "results/coco_out"
RES_DIR <- "results"
REF_DIR <- "data/references"

# ------------------------------------------------------------------------------
# 4.1 Load & Merge CoCo Counts
# ------------------------------------------------------------------------------
samples <- c("MG_I1", "MG_I2", "SG_I1", "SG_I2", "SG_I3")
count_list <- lapply(samples, function(s) {
    df <- read.table(file.path(COCO_DIR, paste0(s, "_coco")), header=TRUE)
    df <- df[order(df$gene_id), c("gene_id", "count", "tpm")]
    colnames(df) <- c("Gene_ID", paste0(s, "_count"), paste0(s, "_tpm"))
    return(df)
})

counts_full <- Reduce(function(x, y) merge(x, y, by="Gene_ID"), count_list)
counts_tpm <- counts_full[, c("Gene_ID", grep("_tpm", colnames(counts_full), value=TRUE))]

# ------------------------------------------------------------------------------
# 4.2 Filtering & Expression Categorization
# ------------------------------------------------------------------------------
# Filter out non-expressed genes (Threshold determined dynamically via bimodality)
keep <- rowSums(counts_tpm[, -1]) >= 25
counts_tpm_filt <- counts_tpm[keep, ]
message(sprintf("Retained %d genes after TPM filtering.", sum(keep)))

# Categorize TPM into High, Medium, Low using cut2
categorize_tpm <- function(column) {
    cat <- as.character(cut2(column, g=3))
    levels <- sort(unique(cat))
    cat[cat == levels[1]] <- "Low"
    cat[cat == levels[2]] <- "Medium"
    cat[cat == levels[3]] <- "High"
    return(cat)
}

for(col in grep("_tpm", colnames(counts_tpm_filt), value=TRUE)) {
    new_col <- sub("_tpm", "_categ", col)
    counts_tpm_filt[[new_col]] <- categorize_tpm(counts_tpm_filt[[col]])
}

# Merge with PlasmoDB annotations
table_plasmodb <- fread(file.path(REF_DIR, "PlasmoDB59_Pmalariae_all_genes_info.txt"))
counts_tpm_annot <- merge(counts_tpm_filt, table_plasmodb[, c(1,15,3,4,11,22,23,25:38), with=FALSE],
                          by.x="Gene_ID", by.y="Gene ID", all.x=TRUE)
write.table(counts_tpm_annot, file.path(RES_DIR, "Pm_tpm_counts_filtered.txt"), row.names=FALSE, sep="\t", quote=FALSE)

# ------------------------------------------------------------------------------
# 4.3 DESeq2 Differential Expression (MG vs SG)
# ------------------------------------------------------------------------------
raw_counts <- counts_full[counts_full$Gene_ID %in% counts_tpm_filt$Gene_ID, grep("_count", colnames(counts_full))]
rownames(raw_counts) <- counts_tpm_filt$Gene_ID
colnames(raw_counts) <- sub("_count", "", colnames(raw_counts))

meta <- data.frame(
    row.names = colnames(raw_counts),
    condition = c("MG", "MG", "SG", "SG", "SG"),
    infection = c("I1", "I2", "I1", "I2", "I3")
)

dds <- DESeqDataSetFromMatrix(countData = round(raw_counts), colData = meta, design = ~ infection + condition)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "MG", "SG"))

# Add annotations to DESeq results
res_df <- as.data.frame(res)
res_df$Gene_ID <- rownames(res_df)
res_annot <- merge(res_df, counts_tpm_annot, by="Gene_ID", all.x=TRUE)
write.table(res_annot, file.path(RES_DIR, "Deseq_res_MGvsSG.txt"), row.names=FALSE, sep="\t", quote=FALSE)

# Volcano Plot
res_df$diff_exp <- "NO"
res_df$diff_exp[res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "SGs (Upregulated)"
res_df$diff_exp[res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "MGs (Upregulated)"

pdf(file.path(RES_DIR, "Volcano_plot.pdf"))
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=diff_exp)) +
    geom_point(alpha=0.8) +
    theme_minimal() +
    scale_color_manual(values=c("black", "red", "blue")) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="darkgray") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgray") +
    labs(title="MG vs SG Differential Expression")
dev.off()

# ------------------------------------------------------------------------------
# 4.4 Inter-Species Correlation (Orthologs Prep for CIBERSORTx)
# ------------------------------------------------------------------------------
orthologs <- fread(file.path(REF_DIR, "PlasmoDB59_Pmalariae_orthologs.txt"))

# Helper Function for Correlation
plot_species_corr <- function(ortho_tag, target_tpm_file, cols_to_keep, new_colnames, out_prefix) {
    # Filter syntenic orthologs
    ortho_filt <- orthologs[grep(ortho_tag, orthologs$Organism) & orthologs$`is syntenic` == "yes", ]
    ortho_filt <- ortho_filt[!duplicated(ortho_filt$Ortholog), c("Gene ID", "Ortholog")]

    # Merge with Pm TPMs
    pm_tpm_ortho <- merge(ortho_filt, counts_tpm_filt, by.x="Gene ID", by.y="Gene_ID")

    # Load target species TPM
    target_tpm <- fread(file.path(REF_DIR, target_tpm_file))
    target_data <- target_tpm[, cols_to_keep, with=FALSE]
    colnames(target_data) <- new_colnames

    # Merge & Correlate
    merged_data <- merge(pm_tpm_ortho[, 1:7], target_data, by.x="Ortholog", by.y="Gene_ID")
    merged_data <- merged_data[complete.cases(merged_data), ]

    if(nrow(merged_data) > 0) {
        M <- cor(merged_data[, -c(1:2)], method="pearson")
        pdf(file.path(RES_DIR, paste0(out_prefix, "_corrplot.pdf")))
        corrplot(M, method="circle", title=paste(out_prefix, "Correlation"), mar=c(0,0,1,0))
        dev.off()
        write.table(merged_data, file.path(RES_DIR, paste0(out_prefix, "_deconvolution_input.txt")), row.names=FALSE, sep="\t", quote=FALSE)
    }
}

# P. falciparum
plot_species_corr(
    ortho_tag = "falciparum 3D7",
    target_tpm_file = "PlasmoDB59_Pfalciparum3D7_transcriptomics.txt",
    cols_to_keep = c(1, 68, 69, 46:52, 4:6, 66, 67, 30, 31, 38, 39, 43:45),
    new_colnames = c("Gene_ID", "male_gamet", "female_gamet", "ring", "troph_early", "troph_late", "schz",
                     "gamet_ii", "gamet_v", "ook", "sporo", "ring2", "ooc", "ooc2", "sporo2", "ooc_sporo",
                     "sporo3", "sporo4", "sporo_cult", "ring3", "schz2", "troph2"),
    out_prefix = "Pm_Pf"
)

# P. berghei
plot_species_corr(
    ortho_tag = "ANKA",
    target_tpm_file = "PlasmoDB59_Pberghei_transcriptomics.txt",
    cols_to_keep = c(1, 6:27),
    new_colnames = c("Gene_ID", "f_gamet", "m_gamet", "eryth", "sporo", "liv_6", "liv_24", "liv_48", "liv_54", "liv_60",
                     "liv_detach", "gamet2", "ook", "schz", "ring_cl15", "troph_cl15", "schz_cl15", "gamet_cl15", "ook_cl15",
                     "troph_cl1", "ring_cl1", "schz_cl1", "gamet_cl1"),
    out_prefix = "Pm_Pb"
)

# ------------------------------------------------------------------------------
# 4.5 Genes lacking orthologs in Pf, Pb, Pv
# ------------------------------------------------------------------------------
gff <- read.delim(file.path(REF_DIR, "PlasmoDB-68_PmalariaeUG01.gff"), comment.char="#", header=FALSE)
gff$Gene_ID <- gsub(";.*", "", gsub(".*PmUG01", "PmUG01", gff$V9))

species_filter <- c("Plasmodium vivax P01", "Plasmodium falciparum 3D7", "Plasmodium berghei ANKA")
ortho_species <- orthologs[orthologs$Organism %in% species_filter, ]
no_ortho <- gff[!(gff$Gene_ID %in% ortho_species$`Gene ID`) & gff$V3 == "protein_coding_gene", ]

# Merge with our counts to see which unique genes are expressed
no_ortho_expr <- merge(no_ortho, counts_tpm_filt, by="Gene_ID")
write.table(unique(no_ortho_expr$Gene_ID), file.path(RES_DIR, "No_Orthologs_Pf_Pb_Pv_Expressed.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

message("=== R ANALYSIS COMPLETED SUCCESSFULLY ===")
EOF

echo "Pipeline finished. Outputs saved in ${RES_DIR}"
