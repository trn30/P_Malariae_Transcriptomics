# ==============================================================================
# Pipeline: Analysis of Genes with Unknown Function
# Description: Evaluates expression, differential expression (DEGs), and
# performs Gene Ontology (GO) enrichment analysis (TopGO & PANNZER2) for
# Plasmodium malariae genes lacking functional annotation.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------

  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(Hmisc)
  library(topGO)
  library(GO.db)
  library(wordcloud)
  library(RColorBrewer)

# Define paths
DATA_DIR <- "data/NO_function"
RES_DIR  <- "results/NO_function"
FIG_DIR  <- "figures"

dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. DATA PROCESSING: EXPRESSION CATEGORIZATION
# ------------------------------------------------------------------------------
# Load data: no funtion genes were obtained from PlasmoDB
no_func <- read.delim(file.path(DATA_DIR, "Master_list_No_function.txt"), header = TRUE)
norm_counts <- read.delim("data/conteos_def_suppl_table.txt", header = TRUE)

# Merge and calculate mean expression per tissue
no_func_expr <- merge(no_func, norm_counts, by = "Gene_ID")
no_func_expr$SG_mean <- rowMeans(no_func_expr[, c(6, 7, 8)])
no_func_expr$MG_mean <- rowMeans(no_func_expr[, c(4, 5)])

# Determine dominant tissue
no_func_expr$group <- ifelse(no_func_expr$SG_mean > no_func_expr$MG_mean, "SG", "MG")

# Separate by dominant tissue and recombine
expr_sg <- no_func_expr[no_func_expr$group == "SG", c(1, 2, 9)]
colnames(expr_sg) <- c("Gene", "Description", "Mean_exp")
expr_sg$group <- "SG"

expr_mg <- no_func_expr[no_func_expr$group == "MG", c(1, 2, 10)]
colnames(expr_mg) <- c("Gene", "Description", "Mean_exp")
expr_mg$group <- "MG"

expr_combined <- rbind(expr_sg, expr_mg)

# Assign labels for plotting highly expressed genes
expr_combined$delabel <- ifelse(expr_combined$Mean_exp >= 100, expr_combined$Description, NA)

# Categorize expression into High, Medium, Low using cut2
cats <- split(expr_combined, cut2(expr_combined$Mean_exp, g = 3))
expr_levels <- expr_combined

# Define thresholds based on the cut2 splits (adjust if data changes)
expr_levels$exp_categ <- "High"
expr_levels$exp_categ[expr_levels$Mean_exp < 147.6] <- "Medium"
expr_levels$exp_categ[expr_levels$Mean_exp < 36.49] <- "Low"

message("Expression Categories Distribution:")
print(table(expr_levels$exp_categ))

# Filter High Expression subset
high_exp <- subset(expr_levels, exp_categ == "High")
high_exp$Description <- gsub("%2C", ",", high_exp$Description)
high_exp$delabel <- ifelse(high_exp$Mean_exp >= 500, high_exp$Description, NA)

pdf(file.path(FIG_DIR, "High_Expression_Unknown_Function.pdf"), width = 10, height = 6)
ggplot(high_exp, aes(x = Gene, y = Mean_exp, fill = group)) +
  geom_bar(stat = "identity") +
  geom_label_repel(aes(label = delabel), nudge_x = 3, colour = "black", fill = "white", nudge_y = -2) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("powderblue", "olivedrab3")) +
  labs(title = "Highly Expressed Genes with Unknown Function", y = "Mean Expression", x = "Gene")
dev.off()

# ------------------------------------------------------------------------------
# 2. DIFFERENTIAL EXPRESSION (DEGs) OVERLAP
# ------------------------------------------------------------------------------
degs <- read.delim("data/List_DEGs_Pmal.txt", header = TRUE)
no_func_deg <- merge(no_func_expr, degs, by = "Gene_ID")

no_func_deg$Tissue <- ifelse(no_func_deg$log2FoldChange <= 0, "MG", "SG")
no_func_deg$Description <- gsub("%2C", ",", no_func_deg$Description)
no_func_deg$Label <- ifelse(no_func_deg$Gene.Name.or.Symbol != "N/A", no_func_deg$Gene.Name.or.Symbol, NA)

pdf(file.path(FIG_DIR, "DEGs_Unknown_Function.pdf"), width = 10, height = 6)
ggplot(no_func_deg, aes(x = Gene_ID, y = log2FoldChange, fill = Tissue)) +
  geom_bar(stat = "identity") +
  geom_text_repel(aes(label = Label), vjust = 5, colour = "black",
                  nudge_y = ifelse(no_func_deg$Tissue == "SG", 4, -2)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Differentially Expressed Genes with Unknown Function")
dev.off()

# Tag genes in the main dataset
deg_genes <- no_func_expr[no_func_expr$Gene_ID %in% no_func_deg$Gene_ID, "Gene_ID"]
no_func_all <- no_func_expr
no_func_all$DEG <- ifelse(no_func_all$Gene_ID %in% deg_genes, "YES", "NO")

# ------------------------------------------------------------------------------
# 3. HELPER FUNCTIONS: TopGO ENRICHMENT
# ------------------------------------------------------------------------------
run_topgo_enrichment <- function(candidate_genes, all_genes, gene2GO_map, output_prefix) {
  message("Running TopGO Enrichment for: ", output_prefix)

  geneList <- factor(as.integer(all_genes %in% candidate_genes))
  names(geneList) <- all_genes

  # Create TopGO Data Object
  GOdata <- new('topGOdata', ontology = 'BP', allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = gene2GO_map)

  # Run Tests
  weight_fisher_res <- runTest(GOdata, algorithm = 'weight01', statistic = 'fisher')

  # Generate Results Table
  allGO <- usedGO(GOdata)
  res_table <- GenTable(GOdata, weightFisher = weight_fisher_res,
                        orderBy = 'weightFisher', topNodes = length(allGO))

  # Clean and filter significant p-values
  res_table$weightFisher <- as.numeric(gsub("<", "", res_table$weightFisher))
  sig_res <- res_table[res_table$weightFisher < 0.05, ]

  # Map genes back to significant GO terms
  mygenes <- genesInTerm(GOdata, sig_res$GO.ID)
  mapped_genes <- data.frame(GO.ID = names(mygenes),
                             genes = sapply(mygenes, paste, collapse = ", "))

  final_res <- merge(mapped_genes, sig_res, by = "GO.ID")

  write.table(final_res, file.path(RES_DIR, paste0(output_prefix, "_enrichment.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  return(final_res)
}

plot_go_wordcloud <- function(go_res_df, output_filename, max_words = 100) {
  # Expand genes to count frequencies
  go_res_long <- go_res_df %>% separate_rows(genes, sep = ", ")
  go_counts <- go_res_long %>% group_by(GO.ID) %>% summarise(gene_count = n())
  go_names <- merge(go_counts, go_res_df, by = "GO.ID")

  pdf(file.path(FIG_DIR, output_filename), paper = "A4")
  wordcloud(words = go_names$Term, freq = go_names$gene_count,
            min.freq = 1, max.words = max_words, random.order = FALSE,
            rot.per = 0.1, colors = brewer.pal(8, "Dark2"))
  dev.off()
}

# ------------------------------------------------------------------------------
# 4. GO ENRICHMENT
# ------------------------------------------------------------------------------
# Prepare Background mapping
all_go <- read.delim(file.path(DATA_DIR, "GOterms_all_genome_release68.txt"), header = TRUE, comment.char = "#")
all_go[all_go == "N/A"] <- NA
all_go <- all_go %>% rowwise() %>%
  mutate(ALL_goterms = paste(na.omit(c(Computed.GO.Component.IDs, Computed.GO.Function.IDs,
                                       Computed.GO.Process.IDs, Curated.GO.Component.IDs,
                                       Curated.GO.Function.IDs, Curated.GO.Process.IDs)), collapse = ", ")) %>% ungroup()

all_go$ALL_goterms <- gsub(";", ", ", all_go$ALL_goterms)


# ------------------------------------------------------------------------------
# 5. PANNZER2 vs PLASMODB ANNOTATION COMPARISON
# ------------------------------------------------------------------------------
pannzer <- read.delim(file.path(DATA_DIR, "Pannzer2_GO_pred_955_def.txt"), header = TRUE)
pannzer <- pannzer[, c("qpid", "goid")]
colnames(pannzer) <- c("Gene.ID", "GO_PANNZER")

pannzer$Gene.ID <- gsub(".1-.*", "", pannzer$Gene.ID)
pannzer$GO_PANNZER <- paste0("GO:", sprintf("%07d", pannzer$GO_PANNZER))

# Base data with existing GO terms
base_info <- all_info[, c("Gene.ID", "GO_terms", "DEG")]

# Merge logic
comparison_df <- merge(base_info, pannzer, by = "Gene.ID", all = TRUE)

# Function to check if PANNZER GO term exists in the PlasmoDB/Fantasia GO terms list
check_coincidence <- function(go_list, target_go) {
  if (is.na(target_go) || is.na(go_list)) return(FALSE)

  # Clean string (remove newlines, extra spaces)
  go_list <- gsub("\n", " ", go_list)
  go_list <- gsub("\\s+", " ", go_list)

  # Split into array and check
  term_array <- trimws(unlist(strsplit(go_list, split = "[,\\s]+")))
  return(target_go %in% term_array)
}

# Vectorized apply
comparison_df$Match <- mapply(check_coincidence, comparison_df$GO_terms, comparison_df$GO_PANNZER)

message(sprintf("Total coincidences found between PANNZER2 and reference DBs: %d", sum(comparison_df$Match)))

# Save final comprehensive table
write.table(comparison_df, file.path(RES_DIR, "Final_Unknown_Function_Orthologs_PANNZER2.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")
