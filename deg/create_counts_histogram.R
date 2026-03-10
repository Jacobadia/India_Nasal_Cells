library(ggplot2)

counts_input_file <- "../artifacts/gene_counts_corrected.tsv"
metadata_file <- "../artifacts/metadata.tsv"
gene_type_file <- "../artifacts/genetype_lookup.txt"

output_plot_name <- "gene_counts_histogram"
output_histogram_file <- paste0("../artifacts/", output_plot_name, ".png")


# Load the gene counts data

counts <- read.delim(
  counts_input_file,
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

# Create a histogram of average gene counts across all samples
average_counts <- rowMeans(counts)

histogram_plot <- ggplot(data.frame(average_counts), aes(x = average_counts)) +
  geom_histogram(bins = 50, fill = "blue", color = "white") +
  scale_x_log10(breaks = 10^(-2:6), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title = "Histogram of Average Gene Counts Across All Samples",
       x = "Average Gene Count (log10 scale)",
       y = "Frequency") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave(output_histogram_file, width = 8, height = 6, plot = histogram_plot)