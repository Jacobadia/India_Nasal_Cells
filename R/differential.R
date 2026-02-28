# options(repos = c(CRAN = "https://cloud.r-project.org"))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
library(DESeq2)

matrix_counts_file <- "../artifacts/gene_counts_corrected.tsv"
metadata_file <- "../artifacts/metadata.tsv"
gene_type_file <- "../artifacts/genetype_lookup.txt"

full_deg_results_file <- "../artifacts/deg_results_full.txt"
significant_deg_results_file <- "../artifacts/deg_results_significant.txt"

counts <- read.delim(
  matrix_counts_file,
  row.names = 1,
  comment.char = "#",
  check.names = FALSE
)

counts <- counts[, 6:ncol(counts)]

sample_names <- colnames(counts)

condition <- ifelse(grepl("LB", sample_names), "latent",
                    ifelse(grepl("LA", sample_names), "active", NA))
condition <- factor(condition, levels = c("latent", "active"))

# Load the tsv of Sex and Age data
metadata <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)

ordering <- match(sample_names, metadata$Nasal.ID)

sex <- metadata$Sex[ordering]
sex <- factor(sex, levels = c("Male", "Female"))
age_raw <- metadata$Age[ordering]
age_scaled <- scale(age_raw)

conditionData <- data.frame(condition = condition, sex = sex, age = age_scaled)
rownames(conditionData) <- sample_names

# Filter out genes with no counts
non_zero_genes <- rowSums(counts) > 0
zero_genes <- rowSums(counts) <= 0
print(paste("Number of genes with no counts:", sum(zero_genes)))
print(paste("Number of genes with counts:", sum(non_zero_genes)))
counts <- counts[non_zero_genes, ]

# Filter out genes that are not protein coding

genetype_lookup <- read.delim(gene_type_file, header = FALSE, stringsAsFactors = FALSE)
colnames(genetype_lookup) <- c("gene_id", "gene_name", "gene_type")

# Keep only protein coding genes
protein_coding_genes <- genetype_lookup$gene_id[genetype_lookup$gene_type == "protein_coding"]
print(paste("Number of protein coding genes:", length(protein_coding_genes)))
non_protein_coding_genes <- length(unique(rownames(counts))) - length(protein_coding_genes)
print(paste("Number of non-protein coding genes:", non_protein_coding_genes))
counts <- counts[rownames(counts) %in% protein_coding_genes, ]


# keep only that have an average count of at least 10 across all samples
average_counts <- rowMeans(counts)
counts <- counts[average_counts >= 10, ]


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conditionData,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

res <- res[!is.na(res$padj), ]
res_filtered <- res[base::order(res$padj), ]

write.table(res_filtered, full_deg_results_file,
            sep = "\t", quote = FALSE, row.names = TRUE)

res_filtered <- res_filtered[res_filtered$padj < 0.05, ]

write.table(res_filtered, significant_deg_results_file,
            sep = "\t", quote = FALSE, row.names = TRUE)