# ROC and models

####   Libraries   ####
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

### extract gene_type_lookup.txt ###
gene_types = read_tsv("R/genetype_lookup.txt")
head(gene_types)

### Prepping gene_counts data ###
gene_counts_raw = read_tsv("R/gene_counts.tsv")
head(gene_counts_raw)

metadata = read_tsv("R/metadata.tsv") 
head(metadata)


print(metadata, n=80)
metadata = metadata %>% mutate(sample = `Nasal ID`)




# joining the two tibbles
gene_counts = left_join(gene_counts_raw, gene_types)


gene_counts <- gene_counts %>%
  rename_with(~ str_extract(.x, "\\d{4}\\.\\d{5}\\.[A-Z]{2}"), 
              matches("\\d{4}\\.\\d{5}\\.[A-Z]{2}"))

gene_counts <- gene_counts %>%
  rename_with(~ gsub("\\.", "-", .x), -name)

gene_counts = gene_counts %>%
  filter(!startsWith(name, "ENSG")) %>%
  select(-c(Geneid,Chr,Start,End,Strand,Length)) %>%
  select(name, everything())

gene_counts

print(gene_counts, n=80)

nis_signatures = c("SPIB", "SHISA2", "TESPA1", "CD1B") 

nis_gene_counts = gene_counts %>%
  filter(name %in% nis_signatures) %>%
  select(-type)

nis_gene_counts


nis_gene_counts <- nis_gene_counts %>%
  pivot_longer(
    cols = -name,
    names_to = "sample",
    values_to = "count"
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = count
  )

print(nis_gene_counts, n=80)

nis_gene_counts = left_join(nis_gene_counts, metadata) %>% select(-`Nasal ID`)

nis_gene_counts = mutate(nis_gene_counts, status = (str_sub(sample, -1) == "A"))

nis_train_counts = select(nis_gene_counts, -sample)
nis_train_counts = nis_train_counts %>% 
  mutate(Sex = factor(Sex)) %>%
  mutate(status = factor(status))

nis_train_counts


# All ID's ending with LA are index cases i.e participants with pulmonary tuberculosis. 
# All ID's ending with LB are Household contacts i.e participants that are IGRA positive. 
# The ID number is common as they are  from the same family. 
# Also, there is an error in one of the ID numbers. I have changed it. 
# The ID was 565-00103 but the correct ID is 656-00103.



print(nis_train_counts, n=80)


## theme ###

### modling data frame ###

prepare_nasal_data <- function(counts_path, coldata_path) {

  counts_data <- read.csv(counts_path, header = TRUE, row.names = 1)
  coldata <- read.csv(coldata_path, header = TRUE, row.names = 1)

  coldata <- coldata %>%
    dplyr::mutate(
      status = factor(status, levels = c("control", "case")),
      sex = factor(sex, levels = c("male", "female"))
    )

  counts_data <- as.matrix(counts_data)

  batch <- factor(coldata$sex)
  sample_group <- factor(coldata$status)

  counts_corrected <- ComBat_seq(
    counts_data,
    batch = batch,
    group = sample_group
  )

  return(list(
    counts_raw = counts_data,
    counts_corrected = counts_corrected,
    coldata = coldata
  ))
}


# prep_nasal_data = prepare_nasal_data(gene_counts_raw, NULL)

counts_data = prep_nasal_data$counts_raw
counts_corrected = prep_nasal_data$counts_corrected
col_data = prep_nasal_data$col_data

# Convert to matrix

counts_data <- as.matrix(counts_data)
all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))






batch <- factor(coldata$sex)
sample_group <- factor(coldata$status)
counts_corrected <- ComBat_seq(as.matrix(counts_data),
    batch=batch,
    group=sample_group)

# DESeqq2, matrix -> DESeqq object for machine learing

#Get and filter for significant genes
# Placeholder for genes
feats <- c(
  "CST1", "CLC", "LGALS12", "ITLN1", "SPIB", "FETUB", "CEACAM21", "CD69",
  "PTGDR2", "VSTM1", "CRISP2", "CLEC12A", "CEBPE", "CCR3", "SHISA2", "TPSAB1",
  "TESPA1", "CISH", "CD1B", "ADAM19", "SIGLEC6", "PRSS33", "SLC9A3", "SORD",
  "OTOF", "CLEC12B", "GAS1", "BTNL8", "FHL3", "POSTN", "IFIT2", "ZBTB16",
  "FABP6", "RSAD2", "CCR7", "CASS4", "DYDC1", "GAPT", "IL31RA", "CNR2",
  "CDH26", "HS3ST3A1", "COL28A1", "HRH4"
)


### machine learnig models ###

### prepping cross fold validation ###
set.seed(42)
fit_control <- trainControl(method="repeatedcv",
                            number=10,
                            repeats=10, 
                            summaryFunction=twoClassSummary,
                            savePredictions = TRUE,
                            classProbs = TRUE,
                            verboseIter = TRUE,
                            search = 'random')


#### RANGER ####
set.seed(42)
r_fit <- caret::train(status~.,
                      data = nis_train_counts,
                      method = 'ranger',
                      metric = 'ROC',
                      tuneLength = 15,
                      trControl = fit_control,
                      importance = 'permutation')

fm_model_r <- evalm(r_fit, gnames='random forest', plots="r", fsize=11)
ggr <- fm_model_r$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_r


#### ELASTIC NET ####

elasticnet = function(feats, fit_control) {
  set.seed(42)
  glmnet_fit <- train(
    status ~ .,
    data = feats,
    method = "glmnet",
    metric = "ROC",
    trControl = fit_control,
    preProcess = c("center","scale")
  )

  fm_model_e <- evalm(glmnet_fit, gnames='elastic net', plots="r", fsize=11)
  gge <- fm_model_e$roc + theme_SL2() + theme(legend.position = "bottom")
  dev.off()
  fm_model_e
}

elasticnet(feats, fit_control)


#### SVM radial ####

svmradial = function(feats, fit_control) {
  set.seed(42)
  svmr_fit <- train(
    status ~ .,
    data = feats,
    method = "svmRadial",
    metric = "ROC",
    trControl = fit_control,
    tuneLength = 15,
    preProcess = c("center","scale")
  )
  
  fm_model_svmr <- evalm(svmr_fit, gnames='SVM (radial)', plots="r", fsize=11)
  ggsvmr <- fm_model_svmr$roc + theme_SL2() + theme(legend.position = "bottom")
  dev.off()
  fm_model_svmr
}

svmradial(feats, fit_control)

#### SVM linear ####

svmlinear = function(feats, fit_control) {
  set.seed(42)
  svm_fit <- train(
    status ~ ., data = feats,
    method = "svmLinear",
    metric = "ROC",
    trControl = fit_control,
    tuneLength = 15,
    preProcess = c("center","scale")
  )
  
  fm_model_svm <- evalm(svm_fit, gnames='SVM (linear)', plots="r", fsize=11)
  ggsvm <- fm_model_svm$roc + theme_SL2() + theme(legend.position = "bottom")
  dev.off()
  fm_model_svm
  
}

svmlinear(feats, fit_control)

#### KNN ####

knn = function(feats, fit_control) {
  
set.seed(42)
knn_fit <- train(
  status ~ ., data = feats,
  method = "knn",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale"))


fm_model_knn <- evalm(knn_fit, gnames='kNN', plots="r", fsize=11)
ggk <- fm_model_knn$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_knn
}

knn(feats, fit_control)

#### PLS####

pls = function(feats, fit_control) {
  set.seed(42)
  pls_fit <- train(
    status ~ ., data = feats,
    method   = "pls",
    metric   = "ROC",
    trControl = fit_control,
    preProcess = c("center","scale"))
  
  
  fm_model_pls <- evalm(pls_fit, gnames='PLS', plots="r", fsize=11)
  ggp <- fm_model_pls$roc + theme_SL2() + theme(legend.position = "bottom")
  dev.off()
  fm_model_pls
}

pls(feats, fit_control)


### Variable importance

plot_importance <- function(fit, model_name, top = 10) {
  imp <- varImp(fit)
  ggplot(imp, top = top) +
    theme_SL2() +
    ggtitle(paste("Top", top, "Genes by Importance in", model_name)) +
    xlab("Genes") +
    ylab("Variable Importance") +
    theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1))
}

#Generate importance plots

#Find overlaps
#overlaps
importance_knn <- varImp(knn_fit,)
importance_svmr <- varImp(svmr_fit)
importance_svml <- varImp(svm_fit)
importance_pls <- varImp(pls_fit)
importance_glm <- varImp(glmnet_fit)
importance_rf <- varImp(r_fit)


importance_rf_df <- as.data.frame(importance_rf$importance)
sorted_importance_rf_df <- importance_rf_df[order(importance_rf_df$Overall,
    decreasing = TRUE), , drop = FALSE]


importance_knn_df <- as.data.frame(importance_knn$importance)
sorted_importance_knn_df <- importance_knn_df[order(importance_knn_df$control,
    decreasing = TRUE), , drop = FALSE]


importance_pls_df <- as.data.frame(importance_pls$importance)
sorted_importance_pls_df <- importance_pls_df[order(importance_pls_df$Overall,
    decreasing = TRUE), , drop = FALSE]


importance_glm_df <- as.data.frame(importance_glm$importance)
sorted_importance_glm_df <- importance_glm_df[order(importance_glm_df$Overall,
    decreasing = TRUE), , drop = FALSE]


# which overlap between machine learning methods ####
# top 10 important
top_predictive_genes_knn <- rownames(sorted_importance_knn_df)[1:10]
top_predictive_genes_rf <- rownames(sorted_importance_rf_df)[1:10]
top_predictive_genes_pls <- rownames(sorted_importance_pls_df)[1:10]
top_predictive_genes_glm <- rownames(sorted_importance_glm_df)[1:10]

overlapping_genes_ml <- Reduce(intersect, list(
  top_predictive_genes_knn,
  top_predictive_genes_rf,
  top_predictive_genes_pls,
  top_predictive_genes_glm
))
overlapping_genes_ml

### refit on overlapping genes ###

### machine learnig models ###

#### Figures  ####



##notes
    # LA = Diseased pulmonary (case), LB = infection (control)
    # Three pipelines 1. 80 all 2. 80 feature 65 train, 14 test 3. 14 feature, 50 train, 15 test
    # function for each model ethan
    # modling data frame function jacob
    # find overlap function ethan
    #variable importance function jacob
