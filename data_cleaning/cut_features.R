args <- commandArgs(trailingOnly = TRUE)

raw_gene_counts_file <- args[1]
outfile_name <- args[2]

counts <- read.delim(
    raw_gene_counts_file,
    comment.char = "#",
    check.names = FALSE
  )

counts <- counts[, c(1, 7:ncol(counts))]

write.table(counts, outfile_name, sep = "\t", quote = FALSE, row.names = FALSE)