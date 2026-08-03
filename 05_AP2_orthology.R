# ==============================================================================
# Pipeline: Best Hit (BH) Ortholog Detection
# Description: Identifies orthologous proteins between Plasmodium malariae and
# P. falciparum and P. vivax using orthologr (BLAST-based BH).
# ==============================================================================

# 1. Load Libraries
  library(orthologr)

# 2. Define Paths
DATA_DIR <- "data/AP2"
RES_DIR  <- "results/AP2"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# Input FASTA files
query_pmal <- file.path(DATA_DIR, "29_AP2_Pmal_protein_0626.fasta")
target_pfal <- file.path(DATA_DIR, "28_AP2_Pfalc_protein_0626.fasta")
target_pviv <- file.path(DATA_DIR, "29_AP2_Pvivax_protein_0626.fasta")
target_pber <- file.path(DATA_DIR, "27_AP2_Pberghei_protein_0626.fasta")

# ==============================================================================
# 3. BH ORTHOLOG DETECTION
# ==============================================================================
CORES <- 1
E_VAL <- "1E-3"

message("Running BH for P. malariae vs P. falciparum...")
ortho_pfal <- orthologs(query_file = query_pmal, subject_files = target_pfal,
                        seq_type = "protein", eval = E_VAL, ortho_detection = "BH", comp_cores = CORES)
write.csv(ortho_pfal, file.path(RES_DIR, "Orthologr_BH_Results_Pfalc.csv"), row.names = FALSE)

message("Running BH for P. malariae vs P. vivax...")
ortho_pviv <- orthologs(query_file = query_pmal, subject_files = target_pviv,
                        seq_type = "protein", eval = E_VAL, ortho_detection = "BH", comp_cores = CORES)
write.csv(ortho_pviv, file.path(RES_DIR, "Orthologr_BH_Results_Pviv.csv"), row.names = FALSE)

message("Running BH for P. malariae vs P. berghei...")
ortho_pber <- orthologs(query_file = query_pmal, subject_files = target_pber,
                        seq_type = "protein", eval = E_VAL, ortho_detection = "BH", comp_cores = CORES)
write.csv(ortho_pber, file.path(RES_DIR, "Orthologr_BH_Results_Pber.csv"), row.names = FALSE)

# ==============================================================================
# 4. CLEANING AND MERGING
# ==============================================================================
message("Creating master gene list and merging tables...")

# Helper function to clean IDs (removes ".1" or other suffixes)
clean_r_id <- function(ids) {
  sapply(strsplit(as.character(ids), "\\."), `[`, 1)
}

# Extract ALL original IDs directly from the P. malariae FASTA
fasta_lines <- readLines(query_pmal)
headers <- fasta_lines[grepl("^>", fasta_lines)]
raw_ids <- sapply(strsplit(headers, " "), `[`, 1)
base_ids <- clean_r_id(gsub(">", "", raw_ids))

# Create the backbone of the final matrix (Ensures all 29 genes are present)
final_matrix <- data.frame(Pm_Gene_ID = base_ids, stringsAsFactors = FALSE)

# Clean and prepare individual target results
prep_ortho_df <- function(df, target_col_name) {
  colnames(df)[1:2] <- c("Pm_Gene_ID", target_col_name)
  df <- df[, c("Pm_Gene_ID", target_col_name)]
  df$Pm_Gene_ID <- clean_r_id(df$Pm_Gene_ID)
  df[[target_col_name]] <- clean_r_id(df[[target_col_name]])
  return(df)
}

df_pfal <- prep_ortho_df(ortho_pfal, "Pf_Ortholog_BH")
df_pviv <- prep_ortho_df(ortho_pviv, "Pv_Ortholog_BH")
df_pber <- prep_ortho_df(ortho_pber, "Pb_Ortholog_BH")

# Merge everything using the master P. malariae list as the base
final_matrix <- merge(final_matrix, df_pfal, by = "Pm_Gene_ID", all.x = TRUE)
final_matrix <- merge(final_matrix, df_pviv, by = "Pm_Gene_ID", all.x = TRUE)
final_matrix <- merge(final_matrix, df_pber, by = "Pm_Gene_ID", all.x = TRUE)

# Replace empty values (NA) with a hyphen "-"
final_matrix[is.na(final_matrix)] <- "-"

message("Preview of Final Master Ortholog Table:")
print(head(final_matrix))

write.csv(final_matrix, file.path(RES_DIR, "Orthologr_BH_Results_Clean_Master.csv"), row.names = FALSE)
message("Pipeline completed successfully.")
