library(tidyverse)
library(DESeq2)
library(caret)
library(ggplot2)
library(MLeval)
library(sva)
library(patchwork)
library(ranger)
library(glmnet)
library(kernlab)
library(pls)
library(pROC)


### Theme

theme_SL2 <- function() {
  theme_bw(base_family = "sans") %+replace%
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.background   = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key        = element_rect(fill = "transparent", colour = NA),
      legend.text       = element_text(size = 12, face = "bold"),
      legend.title      = element_blank()
    )
}


### Load Data
gene_types <- read_tsv(
  "data/genetype_lookup.txt",
  col_names = c("Geneid", "name", "type"),
  show_col_types = FALSE
)

gene_counts_raw <- read_tsv(
  "data/gene_counts.tsv",
  show_col_types = FALSE
)

metadata <- read_tsv(
  "data/metadata.tsv",
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

# counts_corrected <- counts_data

stopifnot(all(colnames(counts_corrected) == rownames(col_data)))

### Define Feature Set

# Genes identified as significant from DESeq2 analysis
# signatures <- c("RNU6-289P","RN7SKP270","FAM90A11","RNU6ATAC36P","HAVCR1P1", "ENSG00000278215")

#top 2 genes identified by Lima
signatures <-  c("ENSG00000156738", "ENSG00000167157")

#top 100 from limma
# signatures <- c("MS4A1","PRRX2","TXLNGY","CAMK4","H2AC12",
#                "SPC24","PIK3R6","NOX1","MTNR1A","RPS4Y1",
#                "DDX3Y","EPOR","KDM5D","USP9Y","ASPM",
#                "PTPN20","ZFP57","NLGN4Y","GPR27","SLC26A4-AS1",
#                "PREX2","GPR68","ZFY","C3orf70","ZNF683",
#                "CPED1","RHBDL3","NCAPG","XIST","GPR143","ESCO2",
#                "PRKY","TTTY14","LINC00278","AKAP5","ENSG00000300770",
#                "CLCNKA","ENSG00000307688","H2BC14","PLCL1",
#                "UTY","SAA1","SGO1","PADI2","SYT8","GCSAM","HENMT1",
#                "ICOS","ENSG00000294508","EIF1AY","TLR10",
#                "SAA2", "FCRL5","ATP6AP1-DT","NLRP2","XK",
#                "CCDC141","HOXB3","ITGAD","CD180",
#                "ENSG00000288049","MARCHF1","ARHGAP15","ENSG00000295911",
#                "ENSG00000248242","ENSG00000273906","SLC26A4",
#                "ENSG00000262714","F11-AS1",
#                "KBTBD8","SH3D21","CD3D","CD8B","LINC02649","DPP4",
#                "RAPGEF6","ENSG00000307289","TDRP","LAX1","TSPAN2",
#                "CRELD2","CENPE","SCAND3","ADGRG5","PDE11A",
#                "ENSG00000307241","HPS1-AS1","TBL1Y","TRPC6",
#                "GASK1B","ACTL10","FCRL4","ENSG00000267317","MEF2C","TOX",
#                "TMEM156","SHC4","P2RY10","TBC1D22A-DT","CIB2"
# )

##top 50
# signatures <- c("MS4A1","PRRX2","TXLNGY","CAMK4","H2AC12",
#                 "SPC24","PIK3R6","NOX1","MTNR1A","RPS4Y1",
#                 "DDX3Y","EPOR","KDM5D","USP9Y","ASPM",
#                 "PTPN20","ZFP57","NLGN4Y","GPR27","SLC26A4-AS1",
#                 "PREX2","GPR68","ZFY","C3orf70","ZNF683",
#                 "CPED1","RHBDL3","NCAPG","XIST","GPR143",
#                 "ESCO2","PRKY","TTTY14","LINC00278","AKAP5",
#                 "ENSG00000300770","CLCNKA","ENSG00000307688","H2BC14","PLCL1",
#                 "UTY","SAA1","SGO1","PADI2","SYT8",
#                 "GCSAM","HENMT1","ICOS","ENSG00000294508","EIF1AY"
# )

# NIS 4-gene signature
# signatures <- c("SPIB", "SHISA2", "TESPA1", "CD1B")

#profiler genes
#signatures <- c("FCGR1A", "GBP5", "GBP6", "C1QB", "FCGR1B", "SEPT4", "ANDKRD22")

# remove version numbers from Geneid
gene_types <- gene_types %>%
  mutate(Geneid_clean = str_remove(Geneid, "\\..*"))

# convert ENS IDs to gene names if present
lookup <- setNames(gene_types$name, gene_types$Geneid_clean)
signatures <- ifelse(signatures %in% names(lookup),
                     lookup[signatures],
                     signatures)

valid_genes <- intersect(signatures, rownames(counts_corrected))

if (length(valid_genes) == 0) {
  stop("No valid genes found in dataset.")
}

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

### Train/Test Split

set.seed(42)

train_index <- createDataPartition(model_data$status, p = 0.80, list = FALSE)

train_data <- model_data[train_index, ]
test_data  <- model_data[-train_index, ]


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
  max_comp <- min(ncol(data) - 1, 20) 
  fit <- train(
    status ~ .,
    data       = data,
    method     = "pls",
    metric     = "ROC",
    trControl  = fit_control,
    preProcess = c("center", "scale"),
    tuneGrid = expand.grid(ncomp = 1:max_comp)
  )
  fm  <- evalm(fit, gnames = "PLS", plots = "r", fsize = 11)
  fm$roc + theme_SL2() + theme(legend.position = "bottom")
  fit
}


### Run Models

r_fit       <- run_random_forest(train_data, fit_control)
glmnet_fit  <- run_elasticnet(train_data, fit_control)
svmr_fit    <- run_svm_radial(train_data, fit_control)
svm_fit     <- run_svm_linear(train_data, fit_control)
knn_fit     <- run_knn(train_data, fit_control)
pls_fit     <- run_pls(train_data, fit_control)


### Figures

make_roc_plot_test <- function(fit, test_data, model_name) {
  probs <- predict(fit, newdata = test_data, type = "prob")
  
  roc_obj <- pROC::roc(
    response  = test_data$status,
    predictor = probs$case,
    levels    = c("control", "case"),
    direction = "<"
  )
  
  auc_val <- round(pROC::auc(roc_obj), 3)
  
  roc_df <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  
  ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(colour = "red", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    annotate("text", x = 0.75, y = 0.05,
             label = paste0("AUC = ", auc_val),
             size = 4, fontface = "bold") +
    labs(title = model_name, x = "1 - Specificity", y = "Sensitivity") +
    theme_SL2()
}

roc_rf   <- make_roc_plot_test(r_fit,      test_data, "Random Forest")
roc_enet <- make_roc_plot_test(glmnet_fit, test_data, "Elastic Net")
roc_svmr <- make_roc_plot_test(svmr_fit,   test_data, "SVM (radial)")
roc_svml <- make_roc_plot_test(svm_fit,    test_data, "SVM (linear)")
roc_knn  <- make_roc_plot_test(knn_fit,    test_data, "kNN")
roc_pls  <- make_roc_plot_test(pls_fit,    test_data, "PLS")

# Assemble into a 2x3 grid

roc_grid <- (roc_rf | roc_enet | roc_svmr) /
  (roc_svml | roc_knn | roc_pls)

ggsave("models/roc_plots_figure_5_Internal_Validate.pdf", plot = roc_grid, width = 12, height = 8)

roc_grid
