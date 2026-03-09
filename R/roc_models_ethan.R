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

theme_SL2 <- function () {
  theme_bw() %+replace%
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_blank(),
          legend.background = element_rect(fill="transparent", colour=NA),
          legend.key = element_rect(fill="transparent", colour=NA),
          legend.text=element_text(size=12, family="Arial", face="bold"),
          legend.title=element_blank())
}

### extract gene_type_lookup.txt ###
gene_types = read_tsv(
  "R/genetype_lookup.txt",
  col_names = c("Geneid", "name", "type")
)
head(gene_types)

### Prepping gene_counts data ###
gene_counts_raw = read_tsv("R/gene_counts_corrected.tsv")
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


# Collapse duplicate gene names
gene_counts <- gene_counts %>%
  group_by(name) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

gene_counts

# Convert to matrix with genes as rows
counts_data <- gene_counts %>%
  column_to_rownames("name") %>%
  as.matrix()

counts_data

### Prepare metadata

col_data <- metadata %>%
  select(sample, Age, Sex) %>%
  mutate(status = str_sub(sample, -1) == "A") %>%
  as.data.frame()

col_data

# rownames must match count matrix columns
rownames(col_data) <- col_data$sample

### Check alignment

# ensure metadata matches counts
all(colnames(counts_data) %in% rownames(col_data))

# reorder metadata if needed
col_data <- col_data[colnames(counts_data), ]

# confirm perfect match
all(colnames(counts_data) == rownames(col_data))

### Batch correction using ComBat_seq


batch <- factor(col_data$Sex)
sample_group <- factor(col_data$status)

counts_corrected <- ComBat_seq(
  counts = counts_data,
  batch = batch,
  group = sample_group
)

counts_corrected <- as.data.frame(counts_corrected)
counts_corrected$name <- rownames(counts_corrected)
counts_corrected


# nis_signatures = c("SPIB", "SHISA2", "TESPA1", "CD1B") 

nis_signatures = c("RNU6-289P","RN7SKP270","FAM90A11","RNU6ATAC36P","HAVCR1P1")


nis_gene_counts = counts_corrected %>%
  filter(name %in% nis_signatures)


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


nis_train_counts

### Prepare counts matrix

# Collapse duplicate gene names
gene_counts <- gene_counts %>%
  group_by(name) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

gene_counts

# Convert to matrix with genes as rows
counts_data <- gene_counts %>%
  column_to_rownames("name") %>%
  as.matrix()

counts_data

### Prepare metadata

col_data <- metadata %>%
  select(sample, Age, Sex) %>%
  mutate(status = str_sub(sample, -1) == "A") %>%
  as.data.frame()

col_data

# rownames must match count matrix columns
rownames(col_data) <- col_data$sample

### Check alignment

# ensure metadata matches counts
all(colnames(counts_data) %in% rownames(col_data))

# reorder metadata if needed
col_data <- col_data[colnames(counts_data), ]

# confirm perfect match
all(colnames(counts_data) == rownames(col_data))

### Batch correction using ComBat_seq


batch <- factor(col_data$Sex)
sample_group <- factor(col_data$status)

counts_corrected <- ComBat_seq(
  counts = counts_data,
  batch = batch,
  group = sample_group
)

counts_corrected <- as.data.frame(counts_corrected)
counts_corrected$gene <- rownames(counts_corrected)
counts_corrected


counts_corrected_df <- as.data.frame(t(counts_corrected))

model_data <- counts_corrected_df %>%
  rownames_to_column("sample") %>%
  left_join(col_data, by = "sample") %>%
  select(all_of(feats), status)

model_data$status <- factor(
  ifelse(model_data$status, "case", "control"),
  levels = c("control","case")
)

model_data = tibble(model_data)
#Get and filter for significant genes
# Placeholder for genes


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

nis_train_counts = nis_train_counts %>% select(-Age)
nis_train_counts = nis_train_counts %>% select(-Sex)


nis_train_counts$status <- factor(
  nis_train_counts$status,
  levels = c(FALSE, TRUE),
  labels = c("control", "case")
)

nis_train_counts


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
fm_model_r$roc

colnames(model_data)



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
  return(fm_model_e)
}

enet = elasticnet(nis_train_counts, fit_control)


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
  return(fm_model_svmr)
}

svmr = svmradial(nis_train_counts, fit_control)

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
  return(fm_model_svm)
}

svml = svmlinear(nis_train_counts, fit_control)

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
  returnValue(fm_model_knn)
}

kneigh = knn(nis_train_counts, fit_control)

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

pl = pls(nis_train_counts, fit_control)


fm_model_r$roc
enet$roc
svmr$roc
svml$roc
kneigh$roc
pl$roc




library(patchwork)

roc_plot <- (fm_model_r$roc +
               enet$roc +
               svmr$roc +
               svml$roc +
               kneigh$roc +
               pl$roc) +
              plot_layout(ncol = 2)

ggsave("roc_models_byu.pdf", roc_plot, width = 14, height = 10)



##notes
# LA = Diseased pulmonary (case), LB = infection (control)
# Three pipelines 1. 80 all 2. 80 feature 65 train, 14 test 3. 14 feature, 50 train, 15 test
# function for each model ethan
# modling data frame function jacob
# find overlap function ethan
#variable importance function jacob
#variable importance function jacob