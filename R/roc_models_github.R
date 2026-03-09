# ---- 1. Libraries ------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(caret)
library(ggplot2)
library(doParallel)
library(devEMF)
library(MLeval)
library(sva)
library(patchwork)
library(cowplot)
library(ranger)
library(glmnet)
library(kernlab)
library(pls)
library(readr)


# ---- 2. Theme ----------------------------------------------------------------

theme_SL2 <- function() {
  theme_bw() %+replace%
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_rect(colour = "black", fill = NA, size = 1),
      plot.background   = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key        = element_rect(fill = "transparent", colour = NA),
      legend.text       = element_text(size = 12, family = "Arial", face = "bold"),
      legend.title      = element_blank()
    )
}


# ---- 3. Load Data ------------------------------------------------------------

gene_types <- read_tsv(
  "R/genetype_lookup.txt",
  col_names = c("Geneid", "name", "type"),
  show_col_types = FALSE
)

gene_counts_raw <- read_tsv(
  "R/gene_counts_corrected.tsv",
  show_col_types = FALSE
)

metadata <- read_tsv(
  "R/metadata.tsv",
  show_col_types = FALSE
) %>%
  mutate(sample = `Nasal ID`)

# NOTE: One sample ID was corrected in the source data:
#   Original: 565-00103  →  Corrected: 656-00103


# ---- 4. Format Gene Counts ---------------------------------------------------

# Join gene type annotations
gene_counts <- left_join(gene_counts_raw, gene_types, by = "Geneid")

# Standardise sample column names (dots → dashes, strip extra info)
gene_counts <- gene_counts %>%
  rename_with(
    ~ str_extract(.x, "\\d{4}\\.\\d{5}\\.[A-Z]{2}"),
    matches("\\d{4}\\.\\d{5}\\.[A-Z]{2}")
  ) %>%
  rename_with(~ gsub("\\.", "-", .x), -name)

# Remove rows where gene name is still an Ensembl ID; drop genomic coordinate cols
gene_counts <- gene_counts %>%
  filter(!startsWith(name, "ENSG")) %>%
  select(-c(Geneid, Chr, Start, End, Strand, Length)) %>%
  select(name, everything())

# Collapse duplicate gene names by summing counts
gene_counts <- gene_counts %>%
  group_by(name) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


# ---- 5. Prepare Matrix & Metadata for DESeq2 / ComBat_seq -------------------

# Convert to genes-as-rows matrix
counts_data <- gene_counts %>%
  column_to_rownames("name") %>%
  as.matrix()

# Build colData for DESeq2 / ComBat
col_data <- metadata %>%
  select(sample, Age, Sex) %>%
  mutate(status = str_sub(sample, -1) == "A") %>%   # LA = index/TB case
  as.data.frame()

rownames(col_data) <- col_data$sample

# Ensure column order in counts matches row order in col_data
col_data <- col_data[colnames(counts_data), ]

stopifnot(all(colnames(counts_data) == rownames(col_data)))


# ---- 6. Batch Correction (ComBat_seq, correcting for Sex) --------------------

batch         <- factor(col_data$Sex)
sample_group  <- factor(col_data$status)

counts_corrected <- ComBat_seq(
  counts = counts_data,
  batch  = batch,
  group  = sample_group
)


# ---- 7. Define Feature Set ---------------------------------------------------

# Genes identified as significant from DESeq2 analysis
# feats <- c(
# "CST1", "CLC", "LGALS12", "ITLN1", "SPIB", "FETUB", "CEACAM21", "CD69",
# "PTGDR2", "VSTM1", "CRISP2", "CLEC12A", "CEBPE", "CCR3", "SHISA2", "TPSAB1",
# "TESPA1", "CISH", "CD1B", "ADAM19", "SIGLEC6", "PRSS33", "SLC9A3", "SORD",
# "OTOF", "CLEC12B", "GAS1", "BTNL8", "FHL3", "POSTN", "IFIT2", "ZBTB16",
#"FABP6", "RSAD2", "CCR7", "CASS4", "DYDC1", "GAPT", "IL31RA", "CNR2",
# "CDH26", "HS3ST3A1", "COL28A1", "HRH4" )

# NIS 4-gene signature
nis_signatures <- c("SPIB", "SHISA2", "TESPA1", "CD1B")


# ---- 8. Build Model Data Frames ----------------------------------------------

# Helper: transpose corrected counts, attach metadata, subset to feature set
make_model_data <- function(feature_set) {
  as.data.frame(t(counts_corrected)) %>%
    rownames_to_column("sample") %>%
    left_join(col_data, by = "sample") %>%
    select(all_of(feature_set), Age, Sex, status) %>%
    mutate(
      Sex    = factor(Sex),
      status = factor(
        ifelse(status, "case", "control"),
        levels = c("control", "case")
      )
    )
}

# model_data     <- make_model_data(feats)           # full gene set
model_data_nis <- make_model_data(nis_signatures)  # NIS 4-gene set


# ---- 9. Cross-Validation Control ---------------------------------------------

set.seed(42)
fit_control <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 10,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE,
  classProbs      = TRUE,
  verboseIter     = TRUE,
  search          = "random"
)


# ---- 10. ML Model Helper Functions -------------------------------------------

# Each function accepts a data frame and fit_control; returns the evalm object.

run_random_forest <- function(data, fit_control) {
  set.seed(42)
  fit <- caret::train(
    status ~ .,
    data       = data,
    method     = "ranger",
    metric     = "ROC",
    tuneLength = 15,
    trControl  = fit_control,
    importance = "permutation"
  )
  fm  <- evalm(fit, gnames = "random forest", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit  # return the fit for importance extraction
}

run_elasticnet <- function(data, fit_control) {
  set.seed(42)
  fit <- train(
    status ~ .,
    data       = data,
    method     = "glmnet",
    metric     = "ROC",
    trControl  = fit_control,
    preProcess = c("center", "scale")
  )
  fm  <- evalm(fit, gnames = "elastic net", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}

run_svm_radial <- function(data, fit_control) {
  set.seed(42)
  fit <- train(
    status ~ .,
    data       = data,
    method     = "svmRadial",
    metric     = "ROC",
    trControl  = fit_control,
    tuneLength = 15,
    preProcess = c("center", "scale")
  )
  fm  <- evalm(fit, gnames = "SVM (radial)", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}

run_svm_linear <- function(data, fit_control) {
  set.seed(42)
  fit <- train(
    status ~ .,
    data       = data,
    method     = "svmLinear",
    metric     = "ROC",
    trControl  = fit_control,
    tuneLength = 15,
    preProcess = c("center", "scale")
  )
  fm  <- evalm(fit, gnames = "SVM (linear)", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}

run_knn <- function(data, fit_control) {
  set.seed(42)
  fit <- train(
    status ~ .,
    data       = data,
    method     = "knn",
    metric     = "ROC",
    trControl  = fit_control,
    tuneLength = 15,
    preProcess = c("center", "scale")
  )
  fm  <- evalm(fit, gnames = "kNN", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}

run_pls <- function(data, fit_control) {
  set.seed(42)
  fit <- train(
    status ~ .,
    data       = data,
    method     = "pls",
    metric     = "ROC",
    trControl  = fit_control,
    preProcess = c("center", "scale")
  )
  fm  <- evalm(fit, gnames = "PLS", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}


# ---- 11. Run Models — Full Feature Set ---------------------------------------

# r_fit       <- run_random_forest(model_data, fit_control)
# glmnet_fit  <- run_elasticnet(model_data, fit_control)
# svmr_fit    <- run_svm_radial(model_data, fit_control)
# svm_fit     <- run_svm_linear(model_data, fit_control)
# knn_fit     <- run_knn(model_data, fit_control)
# pls_fit     <- run_pls(model_data, fit_control)



# ---- 13. Run Models — NIS 4-Gene Subset --------------------------------------

r_fit_nis      <- run_random_forest(model_data_nis, fit_control)
glmnet_fit_nis <- run_elasticnet(model_data_nis, fit_control)
svmr_fit_nis   <- run_svm_radial(model_data_nis, fit_control)
svm_fit_nis    <- run_svm_linear(model_data_nis, fit_control)
knn_fit_nis    <- run_knn(model_data_nis, fit_control)
pls_fit_nis    <- run_pls(model_data_nis, fit_control)


#### Figures — Individual ROC plots per model ####

make_roc_plot <- function(fit, model_name) {
  evalm(fit, gnames = model_name, plots = "r", fsize = 11)$roc +
    theme_SL2() +
    theme(legend.position = "bottom") +
    ggtitle(model_name) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

roc_rf      <- make_roc_plot(r_fit_nis,      "Random Forest")
roc_enet    <- make_roc_plot(glmnet_fit_nis,  "Elastic Net")
roc_svmr    <- make_roc_plot(svmr_fit_nis,    "SVM (radial)")
roc_svml    <- make_roc_plot(svm_fit_nis,     "SVM (linear)")
roc_knn     <- make_roc_plot(knn_fit_nis,     "kNN")
roc_pls     <- make_roc_plot(pls_fit_nis,     "PLS")

# Assemble into a 2x3 grid
(roc_rf | roc_enet | roc_svmr) /
  (roc_svml | roc_knn | roc_pls)
