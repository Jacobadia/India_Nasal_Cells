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


### Theme

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


### Load Data
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


### Format Gene Counts

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


### Prepare Matrix & Metadata for DESeq2 / ComBat_seq 

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


### Batch Correction (ComBat_seq, correcting for Sex)

batch         <- factor(col_data$Sex)
sample_group  <- factor(col_data$status)

counts_corrected <- ComBat_seq(
  counts = counts_data,
  batch  = batch,
  group  = sample_group
)


### Define Feature Set

# Genes identified as significant from DESeq2 analysis
signatures <- c("RNU6-289P","RN7SKP270","FAM90A11","RNU6ATAC36P","HAVCR1P1", "ENSG00000278215")

# NIS 4-gene signature
# signatures <- c("SPIB", "SHISA2", "TESPA1", "CD1B")

#profiler genes
#signatures <- c("FCGR1A", "GBP5", "GBP6", "C1QB", "FCGR1B", "SEPT4", "ANDKRD22")

valid_genes <- intersect(
  signatures,
  colnames(as.data.frame(t(counts_corrected)))
)


### Build Model Data Frames 

# Helper: transpose corrected counts, attach metadata, subset to feature set
make_model_data <- function(feature_set) {
  as.data.frame(t(counts_corrected)) %>%
    rownames_to_column("sample") %>%
    left_join(col_data, by = "sample") %>%
    select(all_of(feature_set), status) %>%
    mutate(
      status = factor(
        ifelse(status, "case", "control"),
        levels = c("control", "case")
      )
    )
}

model_data <- make_model_data(valid_genes)


### Cross-Validation Control

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


### ML Model Helper Functions 
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


### Run Models

r_fit       <- run_random_forest(model_data, fit_control)
glmnet_fit  <- run_elasticnet(model_data, fit_control)
svmr_fit    <- run_svm_radial(model_data, fit_control)
svm_fit     <- run_svm_linear(model_data, fit_control)
knn_fit     <- run_knn(model_data, fit_control)
pls_fit     <- run_pls(model_data, fit_control)


### Figures

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




