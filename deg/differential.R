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

  # Create DGEList
  d0 <- DGEList(counts = as.matrix(counts))
  # technically I'm not supplying a filter here so its
  # using the default. We'll need to discuss that.
  keep <- filterByExpr(d0, group = conditionData$condition)
  d <- d0[keep, , keep.lib.sizes = FALSE]

  # perfrom a TMM normalization to account for some samples
  # having more counts than others.
  d <- calcNormFactors(d)

  # Create design matrix from the provided formula (should always be ~0 + ...)
  design <- model.matrix(design_formula, data = conditionData)
  # performs voom transformation to log CPM
  voom.d <- voom(d, design, plot = FALSE)
  fit <- lmFit(voom.d, design)
  fit <- eBayes(fit)
  # Always auto-generate contrast for 'active' vs 'latent'
  colnames_design <- colnames(design)
  active_col <- grep("active", colnames_design, value = TRUE)
  latent_col <- grep("latent", colnames_design, value = TRUE)
  if (length(active_col) == 1 && length(latent_col) == 1) {
    contrast_str <- sprintf("%s - %s", active_col, latent_col)
    print(paste("Using contrast:", contrast_str))
    contrast <- makeContrasts(contrasts = contrast_str, levels = colnames(design))
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    res <- topTable(fit2, number = Inf, sort.by = "P")
  } else {
    stop("Could not automatically determine 'active' and 'latent' columns in design matrix.")
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
  run_limma_analysis(counts, conditionData, ~0 + condition,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_sex <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~0 + condition + sex,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_age <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~0 + condition + age,
  full_deg_results_file, significant_deg_results_file)
}

run_limma_control_sex_and_age <- function(counts, conditionData,
full_deg_results_file, significant_deg_results_file) {
  run_limma_analysis(counts, conditionData, ~0 + condition + sex + age,
  full_deg_results_file, significant_deg_results_file)
}