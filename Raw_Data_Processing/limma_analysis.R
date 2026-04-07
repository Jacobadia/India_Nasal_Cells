library(edgeR)
library(limma)
library(readr)
library(tidyverse)
library(DESeq2)
library(sva)
library(org.Hs.eg.db)

setwd("/India_Nasal_Cells/data") # change this to the working directory where your gene counts table is

# Read and format the featureCounts table properly

counts <- read_tsv("gene_counts.tsv")
column_to_rownames("Geneid")

setwd("/India_Nasal_Cells/Part_1") # change this to working directory with metadata folder

# Extract groups and create the DGEList

cols <- colnames(counts)
group <- factor(substr(cols, 12, 13))

d0 <- DGEList(counts = counts, group = group)

metadata <- read_tsv("metadata.tsv")
sex <- metadata %>%
  dplyr::select(Sex) %>%
  pull() %>%
  factor()

# Filter lowly expressed genes

keep <- filterByExpr(d0, group = group)
d <- d0[keep, , keep.lib.sizes = FALSE]

# Normalize

d <- calcNormFactors(d)

plotMDS(d, col = as.numeric(d$samples$group))

# Voom

mm <- model.matrix(~0 + group)

# Heatmap of model matrix (optional)
# library(pheatmap)
# pheatmap(mm, cluster_rows = FALSE, cluster_cols = FALSE)

voom.y.d <- voom(d, mm, plot = TRUE)

# Fit linear model and Contrasts

fit <- lmFit(voom.y.d, mm)

contr <- makeContrasts(LA_vs_LB = groupLA - groupLB, levels = colnames(mm))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Extract Results

top.table <- topTable(tmp, sort.by = "P", n = Inf)

# View DEGs

DEGs <- top.table %>% 
  arrange(adj.P.Val) %>% 
  filter(adj.P.Val < 0.5)

print(paste("Significant DEGs found:", nrow(DEGs)))
head(DEGs, 10)

clean_ids <- gsub("\\..*", "", rownames(DEGs))

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = clean_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

DEGs$Gene_Symbol <- gene_symbols

DEGs_final <- DEGs %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::select(Gene_Symbol, Ensembl_ID, everything())

head(DEGs_final)