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

### Prepare Matrix & Metadata for DESeq2 / ComBat_seq 

# Remove featureCounts annotation columns and convert to matrix
counts_data <- gene_counts_raw %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# Build metadata dataframe
col_data <- metadata %>%
  select(sample, Age, Sex) %>%
  mutate(status = str_sub(sample, -1) == "A") %>%   # LA = index/TB case
  as.data.frame()

rownames(col_data) <- col_data$sample

# Align metadata with count matrix
col_data <- col_data[colnames(counts_data), ]
stopifnot(all(colnames(counts_data) == rownames(col_data)))

### Batch Correction (ComBat_seq correcting for Sex)

batch        <- factor(col_data$Sex)
sample_group <- factor(col_data$status)

batch_corrected <- ComBat_seq(
  counts = counts_data,
  batch  = batch,
  group  = sample_group
)

# Convert corrected matrix back to dataframe
batch_corrected_df <- as.data.frame(batch_corrected) %>%
  rownames_to_column("Geneid")

# Save batch corrected counts
write_tsv(batch_corrected_df, "pathways/counts_sex_corrected.tsv")
