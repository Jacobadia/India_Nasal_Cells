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
non_zero_genes <- rowSums(counts) > 0
zero_genes <- rowSums(counts) <= 0
print(paste("Number of genes with no counts:", sum(zero_genes)))
print(paste("Number of genes with counts:", sum(non_zero_genes)))
counts <- counts[non_zero_genes, ]

# Filter out genes that are not protein coding

genetype_lookup <- read.delim("genetype_lookup.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(genetype_lookup) <- c("gene_id", "gene_name", "gene_type")

# Keep only protein coding genes
protein_coding_genes <- genetype_lookup$gene_id[genetype_lookup$gene_type == "protein_coding"]
print(paste("Number of protein coding genes:", length(protein_coding_genes)))
non_protein_coding_genes <- length(unique(rownames(counts))) - length(protein_coding_genes)
print(paste("Number of non-protein coding genes:", non_protein_coding_genes))
counts <- counts[rownames(counts) %in% protein_coding_genes, ]

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