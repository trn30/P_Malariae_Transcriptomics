# ==============================================================================
# Pipeline: Blood vs. Mosquito Stage Expression Comparison
# Description: Evaluates genes expressed in mosquito stages for co-expression
# in blood stages (using Woodford et al., 2020 dataset, TPM >= 25).
# Compares expression overlaps across DEGs, predefined HVG families,
# Single-Cell derived HVGs, and homology-derived HVGs.
# ==============================================================================

library(dplyr)

# Define Paths
DATA_DIR <- "data"
RES_DIR  <- "results/blood_comparison"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. BLOOD VS MOSQUITO
# ==============================================================================
# Load Blood Stage Expression (Threshold TPM >= 25 per Woodford et al. 2020)
blood <- read.delim(file.path(DATA_DIR, "Woodford_et_al_2020_Pm_transcriptome_blood_reanalyzed_tpm_ERR2580896_ERR2690580.txt"), header = TRUE)
blood$expre_ERR2580896 <- ifelse(blood$TPM_ERR2580896 >= 25, "YES", "NO")
blood$expre_ERR2690580 <- ifelse(blood$TPM_ERR2690580 >= 25, "YES", "NO")
blood$expressed_blood <- ifelse(blood$expre_ERR2580896 == "YES" | blood$expre_ERR2690580 == "YES", "YES", "NO")
blood$expressed_both_blood <- ifelse(blood$expre_ERR2580896 == "YES" & blood$expre_ERR2690580 == "YES", "YES", "NO")

# Load Mosquito Counts
norm_counts <- read.delim(file.path(DATA_DIR, "conteos_def_suppl_table.txt"), header = TRUE)
norm_counts$expressed_mosquito <- "YES"

# Merge and identify overlap
merged_counts <- merge(blood, norm_counts, by = "Gene_ID")
both_expressed <- merged_counts[merged_counts$expressed_blood == "YES" & merged_counts$expressed_mosquito == "YES", ]
both_expressed$expr_both_blood_mosquito <- ifelse(both_expressed$expressed_both_blood == "YES" & both_expressed$expressed_mosquito == "YES", "YES", "NO")

write.table(both_expressed, file.path(RES_DIR, "List_genes_expressed_blood_and_mosquito.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# ==============================================================================
# 2. OVERLAP WITH DEGs
# ==============================================================================
logFC <- read.delim(file.path(DATA_DIR, "genes_Log2FC_def_suppl_table.txt"), header = TRUE)
logFC$DEG <- "NO"
logFC$DEG[1:263] <- "YES" # Assumes pre-sorted top 263 are DEGs
deg_list <- logFC[logFC$DEG == "YES", ]

deg_merged <- merge(both_expressed, deg_list, by = "Gene_ID")
message(sprintf("DEGs expressed in both blood and mosquito stages: %d", sum(deg_merged$expr_both_blood_mosquito == "YES")))

# Categorize by expression levels
expr_levels <- read.delim(file.path(DATA_DIR, "gene_list_expression_logfc_levels.txt"), header = TRUE)[, c(1, 33:37)]
deg_levels <- merge(deg_merged, expr_levels, by = "Gene_ID")

high_mg <- deg_levels[deg_levels$I1_MG_categ == "High" & deg_levels$I2_MG_categ == "High", ]
high_sg <- deg_levels[deg_levels$I1_SG_categ == "High" & deg_levels$I2_SG_categ == "High" & deg_levels$I3_SG_categ == "High", ]

write.table(deg_levels, file.path(RES_DIR, "List_genes_expressed_blood_DEG_mosquito.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# ==============================================================================
# 3. OVERLAP WITH PREDEFINED CVG FAMILIES
# ==============================================================================
hvg <- read.delim(file.path(DATA_DIR, "List_ID_HVG_families_ALL_info.txt"), header = TRUE)
hvg_blood_merge <- merge(blood, hvg, by = "Gene_ID")
hvg$expre_blood_stages <- ifelse(hvg$Gene_ID %in% hvg_blood_merge$Gene_ID, "YES", "NO")

hvg_levels <- merge(hvg, expr_levels, by = "Gene_ID")
write.table(hvg_levels, file.path(RES_DIR, "List_genes_expressed_blood_mosquito_HVG_families.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# ==============================================================================
# 4. OVERLAP WITH SINGLE-CELL DERIVED HVGs (Pf & Pb)
# ==============================================================================
# Helper function to avoid code duplication
compare_sc_blood <- function(sc_filename, blood_df, out_filename) {
  sc_data <- read.delim(file.path(DATA_DIR, sc_filename), header = TRUE, sep = "\t")
  merged_data <- merge(blood_df, sc_data, by = "Gene_ID")

  sc_data$expre_blood_stages <- ifelse(sc_data$Gene_ID %in% merged_data$Gene_ID, "YES", "NO")
  write.table(sc_data, file.path(RES_DIR, out_filename), row.names = FALSE, quote = FALSE, sep = "\t")
  return(sc_data)
}

# Process P. falciparum files
pf_ook  <- compare_sc_blood("Pf_ookinete_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pf_ookinete_HVGs_Pm.txt")
pf_ooc  <- compare_sc_blood("Pf_sporozoite_ooc_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pf_sporozoite_ooc_HVGs_Pm.txt")
pf_sgl  <- compare_sc_blood("Pf_sporozoite_sglands_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pf_sporozoite_sglands_HVGs_Pm.txt")

# Process P. berghei files
pb_ook  <- compare_sc_blood("Pb_ookinete_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pb_ookinete_HVGs_Pm.txt")
pb_ooc  <- compare_sc_blood("Pb_oocyst_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pb_oocyst_HVGs_Pm.txt")
pb_sgl  <- compare_sc_blood("Pb_sporo_sgl_HVGs_Pm.txt", blood, "List_genes_expressed_blood_mosquito_Pb_sporo_sgl_HVGs_Pm.txt")

# ==============================================================================
# 5. OVERLAP WITH HOMOLOGY-DERIVED HVGs
# ==============================================================================
homology_hvg <- read.delim(file.path(DATA_DIR, "summary_of_blast_families_final_expressed_orthologs_GOterms.txt"), header = TRUE)
colnames(homology_hvg)[1] <- "Gene_ID"

hom_blood_merge <- merge(homology_hvg, blood, by = "Gene_ID")
homology_hvg$expre_blood_stages <- ifelse(homology_hvg$Gene_ID %in% hom_blood_merge$Gene_ID, "YES", "NO")

write.table(homology_hvg, file.path(RES_DIR, "List_genes_expressed_blood_mosquito_HVG_blast_families.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

message("=== PIPELINE COMPLETED SUCCESSFULLY ===")
