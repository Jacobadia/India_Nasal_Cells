raw_gene_counts_file <- "../artifacts/gene_counts_corrected.tsv"

counts <- read.delim(
    raw_gene_counts_file,
    row.names = 1,
    comment.char = "#",
    check.names = FALSE
  )

counts <- counts[, 6:ncol(counts)]

write.table(counts, "gene_counts_matrix_cut.txt", sep = "\t", quote = FALSE, row.names = TRUE)