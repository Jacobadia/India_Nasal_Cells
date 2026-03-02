# options(repos = c(CRAN = "https://cloud.r-project.org"))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
suppressMessages(suppressWarnings(library(DESeq2)))
source("./parse_counts.R")
source("./filtration.R")

matrix_counts_file <- "../artifacts/gene_counts_corrected.tsv"
metadata_file <- "../artifacts/metadata.tsv"
gene_type_file <- "../artifacts/genetype_lookup.txt"

full_deg_results_file <- "../artifacts/deg_results_full.txt"
significant_deg_results_file <- "../artifacts/deg_results_significant.txt"

get_counts_and_condition_data_output <- get_counts_and_condition_data(matrix_counts_file, metadata_file)
counts <- get_counts_and_condition_data_output$counts
conditionData <- get_counts_and_condition_data_output$conditionData

genetype_lookup <- get_genetype_lookup(gene_type_file)
counts <- filter_total_counts(counts, 0)

counts <- filter_protein_coding_genes(counts, genetype_lookup)

counts <- filter_mean_counts(counts, 10)

run_deg_analysis <- function(counts, conditionData, design_formula) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = conditionData,
                                design = design_formula)

  dds <- DESeq(dds)
  res <- results(dds)

  res <- res[!is.na(res$padj), ]
  res_filtered <- res[base::order(res$padj), ]

  write.table(res_filtered, full_deg_results_file,
              sep = "\t", quote = FALSE, row.names = TRUE)

  res_filtered <- res_filtered[res_filtered$padj < 0.05, ]

  write.table(res_filtered, significant_deg_results_file,
              sep = "\t", quote = FALSE, row.names = TRUE)
}

run_deg_no_control <- function(counts, conditionData) {
  run_deg_analysis(counts, conditionData, ~ condition)
}

run_deg_control_sex <- function(counts, conditionData) {
  run_deg_analysis(counts, conditionData, ~ sex + condition)
}

run_deg_control_age <- function(counts, conditionData) {
  run_deg_analysis(counts, conditionData, ~ age + condition)
}

run_deg_control_sex_and_age <- function(counts, conditionData) {
  run_deg_analysis(counts, conditionData, ~ sex + age + condition)
}