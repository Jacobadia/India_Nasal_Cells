# options(repos = c(CRAN = "https://cloud.r-project.org"))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install(c("DESeq2", "limma", "edgeR"), update = FALSE, ask = FALSE)

suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(edgeR)))

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
  # Ensure sample order matches counts matrix columns
  conditionData <- conditionData[colnames(counts), , drop = FALSE]

  # Create DGEList with group information from condition column
  group <- factor(conditionData$condition)
  d0 <- DGEList(counts = as.matrix(counts), group = group)

  # Filter lowly expressed genes
  keep <- filterByExpr(d0, group = group)
  d <- d0[keep, , keep.lib.sizes = FALSE]

  # Normalize using TMM normalization
  d <- calcNormFactors(d)

  # Create design matrix from the provided formula
  design <- model.matrix(design_formula, data = conditionData)

  # Voom transformation for variance stabilization
  voom.d <- voom(d, design, plot = FALSE)

  # Fit linear model
  fit <- lmFit(voom.d, design)

  # Match standalone limma workflow when design is ~ condition:
  # use no-intercept group design and an explicit active-vs-latent contrast.
  if (identical(as.character(design_formula), as.character(~ condition))) {
    group <- relevel(factor(conditionData$condition), ref = "latent")
    design_no_intercept <- model.matrix(~ 0 + group)
    voom_group <- voom(d, design_no_intercept, plot = FALSE)
    fit_group <- lmFit(voom_group, design_no_intercept)

    contrast <- makeContrasts(active_vs_latent = groupactive - grouplatent,
                              levels = colnames(design_no_intercept))
    fit_contrast <- contrasts.fit(fit_group, contrast)
    fit_contrast <- eBayes(fit_contrast)
    res <- topTable(fit_contrast, number = Inf, sort.by = "P")
  } else {
    # For adjusted models (sex/age covariates), use the condition coefficient
    fit <- eBayes(fit)
    coef_name <- grep("condition", colnames(design), value = TRUE)
    res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  }

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