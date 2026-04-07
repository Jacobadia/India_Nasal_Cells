args <- commandArgs(trailingOnly = TRUE)

raw_gene_counts_file <- args[1]
outfile_name <- args[2]

counts <- read.delim(
    raw_gene_counts_file,
    row.names = 1,
    comment.char = "#",
    check.names = FALSE
  )

counts <- counts[, 6:ncol(counts)]

write.table(counts, outfile_name, sep = "\t", quote = FALSE, row.names = TRUE)