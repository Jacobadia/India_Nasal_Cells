# options(repos = c(CRAN = "https://cloud.r-project.org"))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
suppressMessages(suppressWarnings(library(DESeq2)))

run_deg_analysis <- function(counts, conditionData, design_formula, 
full_deg_results_file, significant_deg_results_file) {
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

run_deg_control_nothing <- function(counts, conditionData, 
full_deg_results_file, significant_deg_results_file) {
  run_deg_analysis(counts, conditionData, ~ condition, 
  full_deg_results_file, significant_deg_results_file)
}

run_deg_control_sex <- function(counts, conditionData, 
full_deg_results_file, significant_deg_results_file) {
  run_deg_analysis(counts, conditionData, ~ sex + condition, 
  full_deg_results_file, significant_deg_results_file)
}

run_deg_control_age <- function(counts, conditionData, 
full_deg_results_file, significant_deg_results_file) {
  run_deg_analysis(counts, conditionData, ~ age + condition, 
  full_deg_results_file, significant_deg_results_file)
}

run_deg_control_sex_and_age <- function(counts, conditionData, 
full_deg_results_file, significant_deg_results_file) {
  run_deg_analysis(counts, conditionData, ~ sex + age + condition, 
  full_deg_results_file, significant_deg_results_file)
}

# --- limma-based analysis for pre-transformed (log2 CPM+1) data ---
suppressMessages(suppressWarnings(library(limma)))

run_limma_analysis <- function(counts, conditionData, design_formula,
full_deg_results_file, significant_deg_results_file) {
  design <- model.matrix(design_formula, data = conditionData)

  fit <- lmFit(as.matrix(counts), design)
  fit <- eBayes(fit)

  coef_name <- grep("conditionactive", colnames(design), value = TRUE)
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")

  # Rename columns to match DESeq2 output convention used elsewhere
  colnames(res)[colnames(res) == "P.Value"]   <- "pvalue"
  colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"

  write.table(res, full_deg_results_file, sep = "\t", quote = FALSE, row.names = TRUE)

  res_significant <- res[!is.na(res$padj) & res$padj < 0.05, ]
  write.table(res_significant, significant_deg_results_file,
              sep = "\t", quote = FALSE, row.names = TRUE)
}

run_limma_control_nothing <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~ condition,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_sex <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~ sex + condition,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_age <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~ age + condition,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_sex_and_age <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~ sex + age + condition,
  full_deg_results_file, significant_deg_results_file)
}