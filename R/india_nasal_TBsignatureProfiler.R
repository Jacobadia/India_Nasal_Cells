library(tidyverse)
library(ggplot2)
library(TBSignatureProfiler)
library(sva)
library(readr)
library(SummarizedExperiment)
library(cowplot)


#"PRRX2"
#"MS4A1"
# PRRX2 shared no commonality with any of the TB signatures
#shared_ms4_tb_signatures = c('Bloom_RES_558', 'Blankley_380', 'Esmail_OD_893')

set.seed(42)
gene_types = read_tsv(
  "../data/genetype_lookup.txt",
  col_names = c("Geneid", "name", "type")
)

gene_counts_raw = read_tsv("../data/gene_counts.tsv")

metadata = read_tsv("../data/metadata.tsv") 

metadata = metadata %>% mutate(sample = `Nasal ID`)

gene_counts <- gene_counts_raw %>%
  rename_with(~ str_extract(.x, "\\d{4}\\.\\d{5}\\.[A-Z]{2}"), 
              matches("\\d{4}\\.\\d{5}\\.[A-Z]{2}"))

gene_counts = left_join(gene_counts_raw, gene_types)

gene_counts = gene_counts %>%
#   #  filter(!startsWith(name, "ENSG")) %>%
#   select(-c(Chr,Start,End,Strand,Length)) %>%
  select(name, everything()) 

# Collapse duplicate gene names
gene_counts <- gene_counts %>%
  group_by(name) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


# Convert to matrix with genes as rows
counts_data <- gene_counts %>%
  column_to_rownames("name") %>%
  as.matrix()

### Prepare metadata

col_data <- metadata %>%
  select(sample, Age, Sex) %>%
  mutate(status = str_sub(sample, -1) == "A") %>%
  as.data.frame()

col_data = tibble(col_data)

# rownames must match count matrix columns
rownames(col_data) <- col_data$sample

# ensure metadata matches counts
all(colnames(counts_data) %in% rownames(col_data))

# reorder metadata if needed
col_data <- col_data[colnames(counts_data), ]

# confirm perfect match
all(colnames(counts_data) == rownames(col_data))

## Batch correction using ComBat_seq
batch <- factor(col_data$Sex)
sample_group <- factor(col_data$status)

counts_corrected <- ComBat_seq(
  counts = counts_data,
  batch = batch,
  group = sample_group
)


counts_corrected <- as.data.frame(counts_corrected)
counts_corrected$name <- rownames(counts_corrected)
# [1] "Nasal ID" "Age" "Sex" "sample"  # or maybe only the first 3
counts_corrected <- counts_corrected[, colnames(counts_corrected) != "name", drop = FALSE]

counts_corrected[] <- lapply(counts_corrected, as.numeric)
counts_corrected <- as.matrix(counts_corrected)

india_tb <- SummarizedExperiment(
  assays  = list(counts = counts_corrected),
  colData = col_data
)

india_tb <- mkAssay(india_tb, log = TRUE, counts_to_CPM = TRUE)
siglist_indiatb <- names(TBsignatures)

## Run the TBSignatureProfiler to score the signatures in the data
ssgsea_result <- runTBsigProfiler(input = india_tb,
                                  useAssay = "counts",
                                  signatures = TBsignatures,
                                  algorithm = "ssGSEA",
                                  combineSigAndAlgorithm = TRUE,
                                  parallel.sz = 1)


cd <- as.data.frame(colData(ssgsea_result))
good_sigs = colnames(cd[5:82])

booted <- bootstrapAUC(SE_scored = ssgsea_result, annotationColName = "status",
                       signatureColNames = good_sigs, num.boot = 30)
boot_auc <- as.data.frame(booted$`Boot AUC Values`)

boot_results <- tibble(
  signature   = colnames(boot_auc),
  auc_nonboot = booted$`Non-Boot AUC Values`,
  auc_boot    = colMeans(boot_auc),
  ci_lower    = booted$`pROC Lower`,
  ci_upper    = booted$`pROC Upper`,
  p_value     = booted$`P-values`
)

readr::write_csv(boot_results, "../data/bootstrapAUC_results_india_nasal.csv")
significant_signatures = boot_results %>% 
    filter(p_value <= 0.05) %>%
    pull(signature) %>%
    print()

boot_results <- boot_results %>%
  mutate(
    neg_log_p = -log(p_value),
    significant = p_value < 0.05
  )

plot = ggplot(boot_results, aes(signature, neg_log_p)) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
  geom_text(
    data = subset(boot_results, significant),
    aes(label = signature),
    vjust = -0.3,
    hjust = 0,
    angle = 30,   # tilt labels
    size = 3,
    color = "red"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),  
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("Signatures from TBSignatureProfiler") +
  ylab("Negative Log P-value") + 
  labs(color = "P value below 0.05", title = "P-Value Significance of Signatures from TBsignatureProfiler") 
  
ggsave2("../data/bootstrap_results.pdf", plot = plot, width = 8, height = 10)
