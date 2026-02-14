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

# Convert to matrix

# DESeqq2, matrix -> DESeqq object for machine learing

#Get and filter for significant genes

### machine learnig models ###

### Variable importance

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