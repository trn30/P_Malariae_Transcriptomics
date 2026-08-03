# Pipeline: ApiAP2 Family Expression and Promoter Motif Discovery (HOMER)

```bash
#!/bin/bash
# ==============================================================================
# Pipeline: ApiAP2 Family Expression and Promoter Motif Discovery (HOMER)
# Description: Evaluates the expression of ApiAP2 transcription factors across
# blood and mosquito stages, plots heatmaps/barplots, and performs motif
# discovery on the 5' UTRs of highly expressed genes.
# ==============================================================================

set -e

# --- Directories ---
DATA_DIR="data"
RES_DIR="results"
FIG_DIR="figures"
HOMER_DIR="${RES_DIR}/homer_motifs"
REF_DIR="data/references"

mkdir -p ${RES_DIR} ${FIG_DIR} ${HOMER_DIR}

echo "=== API-AP2 EXPRESSION ANALYSIS (R) ==="

Rscript - <<'EOF'
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(scales)
    library(iheatmapr)
    library(data.table)
})

# Paths
DATA_DIR <- "data"
FIG_DIR <- "figures"
RES_DIR <- "results"

# 1. Load ApiAP2 genes
ap2_data <- read.delim(file.path(DATA_DIR, "ApiAP2_Pm_29_orthology.txt"), header=TRUE)
colnames(ap2_data)[1] <- "Gene_ID"
ap2_genes <- unique(ap2_data$Gene_ID)

# Manually add specific genes of interest based on literature/orthology
ap2_genes <- c(ap2_genes, "PmUG01_10026800", "PmUG01_03014100", "PmUG01_13050300", "PmUG01_11058700")

# 2. Load Blood Transcriptome (Woodford et al. 2020)
blood <- read.delim(file.path(DATA_DIR, "Woodford_et_al_2020_Pm_transcriptome_blood_reanalyzed_tpm_ERR2580896_ERR2690580.txt"), header=TRUE)
colnames(blood) <- c("Gene_ID", "Blood_1", "Blood_2")
blood$blood_mean <- rowMeans(blood[, c("Blood_1", "Blood_2")])

# 3. Load our Mosquito Stages Transcriptome results
mosquito <- read.delim(file.path(DATA_DIR, "conteos_def_suppl_table.txt"), header=TRUE)
colnames(mosquito) <- c("Gene_ID", "I1_ooc", "I2_ooc", "I1_sporo", "I2_sporo", "I3_sporo")
mosquito$Sporo_mean <- rowMeans(mosquito[, c("I1_sporo", "I2_sporo", "I3_sporo")])
mosquito$Ooc_mean <- rowMeans(mosquito[, c("I1_ooc", "I2_ooc")])

# 4. Merge Expression Data and Filter for expressed AP2
express_merged <- merge(blood, mosquito, by="Gene_ID", all=TRUE)
express_filtered <- express_merged[, c("Gene_ID", "Ooc_mean", "Sporo_mean", "blood_mean")]
ap2_expr <- express_filtered[express_filtered$Gene_ID %in% ap2_genes, ]
rownames(ap2_expr) <- ap2_expr$Gene_ID

# 5. Load Log2FoldChange Data
logfc_data <- read.delim(file.path(DATA_DIR, "genes_Log2FC_def_suppl_table.txt"), header=TRUE)
# Assuming top 263 genes are DEGs based on original script logic
logfc_data$DEG <- "NO"
logfc_data$DEG[1:263] <- "YES"

# Merge AP2 expression with Log2FC
ap2_expr_log <- merge(ap2_expr, logfc_data, by="Gene_ID")
rownames(ap2_expr_log) <- ap2_expr_log$Gene_ID

# Prepare matrix for heatmap
log2fc_vals <- ap2_expr_log[, c("Gene_ID", "log2FoldChange")]

expr_matrix_df <- ap2_expr_log[, c("Gene_ID", "Ooc_mean", "Sporo_mean", "blood_mean")]
expr_matrix <- as.matrix(expr_matrix_df[, -1])
rownames(expr_matrix) <- expr_matrix_df$Gene_ID

# Sort by Log2FC
sort_order <- order(ap2_expr_log$log2FoldChange, decreasing = FALSE)
expr_matrix_sorted <- expr_matrix[sort_order, ]
log2fc_sorted <- log2fc_vals[sort_order, "log2FoldChange", drop=FALSE]

# Column normalization (Z-score)
expr_matrix_scaled <- scale(expr_matrix_sorted)

# 6. Generate Interactive Heatmap (iheatmapr)
my_palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(255)
ylgnbu_palette <- brewer.pal(7, "RdYlGn")

heatmap_ap2 <- iheatmap((expr_matrix_scaled[, 1:2] + 1),
                        name = "Expression AP2
(Mosquito Stages)",
                        colors = my_palette, cluster_rows = "none",
                        scale = "cols", scale_method = "center",
                        col_labels = TRUE, row_labels = TRUE)

# Add Blood Stage Signal
heatmap_ap2 <- add_row_signal(heatmap_ap2, signal = (expr_matrix_scaled[, 3] + 1),
                              name = "Expression AP2
(Blood Stage)",
                              colors = my_palette, size = 0.25)

# Add Log2FC Signal
heatmap_ap2 <- add_row_signal(heatmap_ap2, signal = log2fc_sorted$log2FoldChange,
                              name = "Log2FoldChange
(MG vs SG)",
                              colors = ylgnbu_palette, side = "left")

# Save Heatmap
save_iheatmap(heatmap_ap2, file.path(FIG_DIR, "Heatmap_ApiAP2.html"))

# 7. Generate Barplot of AP2 expression
ap2_barplot_data <- logfc_data[logfc_data$Gene_ID %in% ap2_genes, ]

p1 <- ggplot(ap2_barplot_data, aes(x = reorder(Gene_ID, -log2FoldChange), y = log2FoldChange, fill=DEG)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values = c("steelblue3", "coral3")) +
    labs(x = "AP2 Expressed Genes", y = "Log2 Fold Change", title = "ApiAP2 Differential Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 9))

pdf(file.path(FIG_DIR, "Barplot_ApiAP2_Log2FC.pdf"), width = 10, height = 6)
print(p1)
dev.off()

message("AP2 Expression Analysis Complete. Heatmaps and Barplots saved.")
EOF

echo "=== HOMER MOTIF DISCOVERY (BASH) ==="
CORES=40

# Define Reference Genomes and Annotations
REF_FASTA="${REF_DIR}/PlasmoDB-59_PmalariaeUG01_Genome.fasta"
REF_GTF="${REF_DIR}/PlasmoDB-59_PmalariaeUG01_utr_biotype.gtf"
REF_GFF_GENES="${REF_DIR}/PlasmoDB-59_PmalariaeUG01.gff.gtf_biotype.gtf"

# 1. Extract 5' UTR regions (1.5 Kb upstream)
# Subtract gene body coordinates from full regions containing UTRs
awk '{gsub(/gene_id|;|"/,"")}1' ${REF_GTF} | awk -v OFS="	" '{ print $1,$4,$5,$9,$6,$7 }' | sort | uniq > ${RES_DIR}/genes_including_utr.bed
awk '{gsub(/gene_id|;|"/,"")}1' ${REF_GFF_GENES} | awk -v OFS="	" '{ print $1,$4,$5,$9,$6,$7 }' | sort | uniq > ${RES_DIR}/genes_body_only.bed
bedtools subtract -a ${RES_DIR}/genes_including_utr.bed -b ${RES_DIR}/genes_body_only.bed -s > ${RES_DIR}/5_UTR_regions.bed

# Create BED for the 1005 Highly Expressed Genes
# (Assumes 1005_High_expr_5utrs.bed is generated by subsetting 5_UTR_regions.bed with your highly expressed genes list)
HIGH_EXPR_BED="${RES_DIR}/1005_High_expr_5utrs.bed"

# 2. Prepare Background (Shuffle)
bedtools shuffle -i ${HIGH_EXPR_BED} -g <(cat ${REF_FASTA}.fai | awk '{ print $1"	"$3}') > ${RES_DIR}/1005_High_expr_5utrs_shuffl.bed

# 3. Run HOMER Motif Discovery
echo "Running HOMER findMotifsGenome.pl (De Novo)..."
findMotifsGenome.pl ${HIGH_EXPR_BED} ${REF_FASTA} ${HOMER_DIR}/de_novo_motifs     -bg ${RES_DIR}/1005_High_expr_5utrs_shuffl.bed -p ${CORES} -len 8,10,12 -size given -mset all -dumpFasta

echo "Running HOMER findMotifsGenome.pl (Known Campbell et al. AP2 Motifs)..."
KNOWN_MOTIFS="${DATA_DIR}/homer_motifs_apiAP2_campbell.motif"
findMotifsGenome.pl ${HIGH_EXPR_BED} ${REF_FASTA} ${HOMER_DIR}/known_motifs     -bg ${RES_DIR}/1005_High_expr_5utrs_shuffl.bed -p ${CORES} -len 8,10,12 -size given -mset all -dumpFasta -find ${KNOWN_MOTIFS} > ${HOMER_DIR}/1005_High_expr_find_apiAP2.txt

# 4. Annotate Peaks for specific AP2-L and AP2-G motifs
echo "Running HOMER annotatePeaks.pl..."
# AP2-L Motifs (assuming Homer assigned them motif1, motif11, motif13 in de novo results)
annotatePeaks.pl ${HIGH_EXPR_BED} ${REF_FASTA} -m ${HOMER_DIR}/de_novo_motifs/homerResults/motif1.motif ${HOMER_DIR}/de_novo_motifs/homerResults/motif11.motif ${HOMER_DIR}/de_novo_motifs/homerResults/motif13.motif > ${HOMER_DIR}/1005_High_expr_AP2_L_instances.txt

# AP2-G Motifs (assuming motif10)
annotatePeaks.pl ${HIGH_EXPR_BED} ${REF_FASTA} -m ${HOMER_DIR}/de_novo_motifs/homerResults/motif10.motif > ${HOMER_DIR}/1005_High_expr_AP2_G_instances.txt


echo "=== INSPECTION HOMER OUTPUT (R) ==="

Rscript - <<'EOF'
library(data.table)
library(dplyr)

# 1. PATH DEFINITIONS AND GLOBAL VARIABLES
# ------------------------------------------------------------------------------
# Note for users: Change 'base_dir' to the directory where data is hosted
base_dir <- "path/to/your/data/" 
deg_file <- paste0(base_dir, "tableS5.txt")

# List of Transcription Factors and the column(s) where HOMER saved the motif
tf_info <- list(
  "AP2_L"  = 22, # HAS 3 MOTIFS BUT JUST 1 BETWEEN 10 FIRST
  "AP2_G5" = 22:24, # AP2_G5 has 3 motifs
  "AP2_P"  = 22,
  "AP2"    = 22,
  "AP2_G"  = 22
)

# 2. FUNCTION TO PROCESS HOMER FILES
# ------------------------------------------------------------------------------
process_homer_motifs <- function(tf_name, motif_cols, base_dir) {
  
  in_file  <- paste0(base_dir, "1005_High_expr_5utrs.bed_", tf_name, "_motifs_instances.txt")
  out_file <- paste0(base_dir, "1005_High_expr_5utrs.bed_", tf_name, "_motifs_instances_2.bed")
  
  a <- fread(in_file)
  setnames(a, 1, "Gene_ID")
  
  bed_list <- list()
  genes_list <- list()
  
  # Iterate over motif columns
  for (i in seq_along(motif_cols)) {
    col_idx <- motif_cols[i]
    
    # Extract coordinates and motif data
    b_temp <- a[, c(2, 3, 4, col_idx, 6, 5), with = FALSE]
    setnames(b_temp, c("chr", "start", "end", "ID", "score", "strand"))
    
    # Filter empty rows and format the ID
    b_temp <- b_temp[ID != "" & !is.na(ID)]
    b_temp[, ID := paste0(tf_name, "_m", i, "_", ID)]
    
    bed_list[[i]] <- b_temp
    genes_list[[i]] <- unique(a$Gene_ID[a[[col_idx]] != "" & !is.na(a[[col_idx]])])
  }
  
  # Save unified BED file
  b_final <- rbindlist(bed_list)
  write.table(b_final, file = out_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Return the unique list of genes with motifs
  unique_genes <- unique(unlist(genes_list))
  message(sprintf("Total genes with motifs for %s: %d", tf_name, length(unique_genes)))
  
  return(unique_genes)
}

# 3. EXTRACT MOTIFS FOR ALL TFs
# ------------------------------------------------------------------------------
# Run the function for each TF and store the gene lists
tf_genes <- lapply(names(tf_info), function(tf) {
  process_homer_motifs(tf_name = tf, motif_cols = tf_info[[tf]], base_dir = base_dir)
})
names(tf_genes) <- names(tf_info)


# 4. DIFFERENTIAL EXPRESSION ANALYSIS (RNA-seq)
# ------------------------------------------------------------------------------
tableS5 <- read.delim(deg_file)

# Filter significant DEGs (padj < 0.05) that contain the motif for each TF
deg_counts <- sapply(names(tf_genes), function(tf) {
  
  # Retrieve the gene list for the current TF
  gene_list <- tf_genes[[tf]]
  
  # Intersect with significant DEGs
  deg_with_motif <- tableS5 %>%
    filter(padj < 0.05 & Gene_ID %in% gene_list)
  
  # Return the row count (number of genes)
  return(nrow(deg_with_motif))
})

# Print the final count of DEGs per TF
print("--- Number of Significant DEGs with Motifs per TF ---")
print(deg_counts)

```
