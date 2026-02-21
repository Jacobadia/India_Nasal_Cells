# options(repos = c(CRAN = "https://cloud.r-project.org"))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
library(DESeq2)

counts <- read.delim(
  "gene_counts_matrix.txt",
  row.names = 1,
  comment.char = "#"
)

counts <- counts[, 6:ncol(counts)]

sample_names <- colnames(counts)

condition <- ifelse(grepl("LB", sample_names), "latent",
                    ifelse(grepl("LA", sample_names), "active", NA))
condition <- factor(condition, levels = c("latent", "active"))

conditionData <- data.frame(condition = condition)
rownames(conditionData) <- sample_names

# Filter out genes with no counts
keep_genes <- rowSums(counts) > 0
throw_away_genes <- rowSums(counts) <= 0
print(paste("Number of genes with no counts:", sum(throw_away_genes)))
print(paste("Number of genes with counts:", sum(keep_genes)))
counts <- counts[keep_genes, ]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conditionData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

res <- res[!is.na(res$padj), ]
res_filtered <- res[base::order(res$padj), ]

write.table(res_filtered, "deg_results_full.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)

res_filtered <- res_filtered[res_filtered$padj < 0.05, ]

write.table(res_filtered, "deg_results_significant.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)