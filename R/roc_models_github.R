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

### theme ###

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


# Convert to matrix

# DESeqq2, matrix -> DESeqq object for machine learing

#Get and filter for significant genes

### machine learnig models ###

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

### Refit nasal models based on overlaps

### machine learnig models ###

#### Figures  ####



##notes
    #Three pipelines 1. 80 all 2. 80 feature 65 train, 14 test 3. 14 feature, 50 train, 15 test
    #function for each model ethan
    #modling data frame function jacob
    #find overlap function ethan
    #variable importance function jacob