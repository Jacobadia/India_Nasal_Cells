filter_total_counts <- function(counts, threshold = 0) {
  # Filter out genes with counts less than or equal
  # to the threshold across all samples
  passing_genes <- rowSums(counts) > threshold
  failing_genes <- rowSums(counts) <= threshold
  print(paste("Number of genes with total counts >", threshold, ":", sum(passing_genes)))
  print(paste("Number of genes  with total counts <=", threshold, ":", sum(failing_genes)))
  counts <- counts[passing_genes, ]
  return (counts)
}

filter_protein_coding_genes <- function(counts, genetype_lookup) {
  # Keep only protein coding genes
  protein_coding_genes <- genetype_lookup$gene_id[genetype_lookup$gene_type == "protein_coding"]
  print(paste("Number of protein coding genes:", length(protein_coding_genes)))
  non_protein_coding_genes <- length(unique(rownames(counts))) - length(protein_coding_genes)
  print(paste("Number of non-protein coding genes:", non_protein_coding_genes))
  counts <- counts[rownames(counts) %in% protein_coding_genes, ]
  return (counts)
}

filter_mean_counts <- function(counts, threshold) {
  # keep only that have an average count of at least n across all samples
  average_counts <- rowMeans(counts)
  passing_genes <- average_counts >= threshold
  failing_genes <- average_counts < threshold
  print(paste("Number of genes with average counts < ", threshold, ":", sum(failing_genes)))
  print(paste("Number of genes with average counts >= ", threshold, ":", sum(passing_genes)))
  counts <- counts[average_counts >= threshold, ]
  return(counts)
}

filter_protein_mean_counts <- function(counts, genetype_lookup, threshold) {
  counts <- filter_total_counts(counts, 0)
  counts <- filter_protein_coding_genes(counts, genetype_lookup)
  counts <- filter_mean_counts(counts, threshold)
  return (counts)
}

filter_out_hemoglobin_genes <- function(counts, hemoglobin_lookup) {
  # Remove hemoglobin genes
  # Strip version suffix (e.g. ENSG00000244734.3 -> ENSG00000244734) before matching
  base_ids <- gsub("\\.[0-9]+$", "", rownames(counts))
  hemoglobin_gene_ids <- hemoglobin_lookup$Ensembl.gene.ID
  passing_genes <- !base_ids %in% hemoglobin_gene_ids
  failing_genes <- base_ids %in% hemoglobin_gene_ids
  print(paste("Number of hemoglobin genes:", sum(failing_genes)))
  print(paste("Number of non-hemoglobin genes:", sum(passing_genes)))
  counts <- counts[passing_genes, ]
  return(counts)
}

filter_protein_hemoglobin <- function(counts, genetype_lookup, hemoglobin_lookup) {
  counts <- filter_total_counts(counts, 0)
  counts <- filter_protein_coding_genes(counts, genetype_lookup)
  counts <- filter_out_hemoglobin_genes(counts, hemoglobin_lookup)
  return (counts)
}

filter_protein_hemoglobin_mean_counts <- function(counts, genetype_lookup, hemoglobin_lookup, threshold) {
  counts <- filter_protein_hemoglobin(counts, genetype_lookup, hemoglobin_lookup)
  counts <- filter_mean_counts(counts, threshold)
  return (counts)
}

# Log2(CPM + 1) transform: normalize each sample to counts per million, then apply log2(x + 1)
lpm_cpm_transform <- function(counts) {
  library_sizes <- colSums(counts)
  cpm <- sweep(counts, 2, library_sizes, FUN = "/") * 1e6
  return(log2(cpm + 1))
}

# Complementary functions that apply log per million transform before filtering
lpm_filter_pure <- function(counts) {
  counts <- lpm_cpm_transform(counts)
  return(counts)
}

lpm_filter_protein_coding <- function(counts, genetype_lookup) {
  counts <- filter_total_counts(counts, 0)
  counts <- filter_protein_coding_genes(counts, genetype_lookup)
  counts <- lpm_cpm_transform(counts)
  return(counts)
}

lpm_filter_protein_mean_counts <- function(counts, genetype_lookup, threshold) {
  counts <- filter_total_counts(counts, 0)
  counts <- filter_protein_coding_genes(counts, genetype_lookup)
  counts <- lpm_cpm_transform(counts)
  counts <- filter_mean_counts(counts, threshold)
  return(counts)
}

lpm_filter_protein_hemoglobin <- function(counts, genetype_lookup, hemoglobin_lookup) {
  counts <- filter_total_counts(counts, 0)
  counts <- filter_protein_coding_genes(counts, genetype_lookup)
  counts <- filter_out_hemoglobin_genes(counts, hemoglobin_lookup)
  counts <- lpm_cpm_transform(counts)
  return(counts)
}

lpm_filter_protein_hemoglobin_mean_counts <- function(counts, genetype_lookup, hemoglobin_lookup, threshold) {
  counts <- filter_protein_hemoglobin(counts, genetype_lookup, hemoglobin_lookup)
  counts <- lpm_cpm_transform(counts)
  counts <- filter_mean_counts(counts, threshold)
  return(counts)
}