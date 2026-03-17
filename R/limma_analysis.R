library(edgeR)
library(limma)
library(readr)
library(tidyverse)

setwd("/grphome/grp_tb/processing_scripts/results/quality_test_featureCounts")

# Read and format the featureCounts table properly
counts <- read_tsv("gene_counts_corrected.tsv")

counts <- counts %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid") # Moves Geneid to row names so the matrix is strictly numeric

# Extract groups and create the DGEList
cols <- colnames(counts)
group <- factor(substr(cols, 12, 13))

d0 <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(d0, group = group)
d <- d0[keep, , keep.lib.sizes = FALSE]

# Normalize
d <- calcNormFactors(d)

# Check how many genes are left
#print(paste("Genes left after filtering:", nrow(d)))

plotMDS(d, col = as.numeric(d$samples$group))

# Voom
mm <- model.matrix(~0 + group)

# Heatmap of model matrix 
#library(pheatmap)
#pheatmap(mm, cluster_rows = FALSE, cluster_cols = FALSE)

voom.y.d <- voom(d, mm, plot = TRUE)

# Fit linear model and Contrasts
fit <- lmFit(voom.y.d, mm)

contr <- makeContrasts(LA_vs_LB = groupLA - groupLB, levels = colnames(mm))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

#Extract Results
top.table <- topTable(tmp, sort.by = "P", n = Inf)

# View DEGs
DEGs <- top.table %>% 
  arrange(adj.P.Val) #%>% 
  #filter(adj.P.Val < 0.5)

print(paste("Significant DEGs found:", nrow(DEGs)))
head(DEGs, 5)