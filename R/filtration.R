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