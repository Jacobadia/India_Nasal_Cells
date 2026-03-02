filter_total_counts <- function(counts, threshold = 0) {
  # Filter out genes with counts less than or equal
  # to the threshold across all samples
  passing_genes <- rowSums(counts) > threshold
  failing_genes <- rowSums(counts) <= threshold
  print(paste("Number of genes  with total counts <=", threshold, ":", sum(failing_genes)))
  print(paste("Number of genes with total counts >", threshold, ":", sum(passing_genes)))
  counts <- counts[passing_genes, ]
  return (counts)
}