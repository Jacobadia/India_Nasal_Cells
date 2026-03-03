get_counts_data <- function(matrix_counts_file) {
  print(paste("Loading counts data from:", matrix_counts_file))
  counts <- read.delim(
    matrix_counts_file,
    row.names = 1,
    comment.char = "#",
    check.names = FALSE
  )

  counts <- counts[, 6:ncol(counts)]
  return (counts)
}

get_condition_data <- function(counts, metadata_file) {
  print(paste("Loading metadata from:", metadata_file))
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
  return(conditionData)
}

get_genetype_lookup <- function(gene_type_file) {
  # parse gene type information
  print(paste("Loading gene type lookup from:", gene_type_file))
  genetype_lookup <- read.delim(gene_type_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(genetype_lookup) <- c("gene_id", "gene_name", "gene_type")
  return(genetype_lookup)
}

get_hemoglobin_lookup <- function(hemoglobin_file) {
  print(paste("Loading hemoglobin gene lookup from:", hemoglobin_file))
  hemoglobin_genes <- read.delim(hemoglobin_file, header = TRUE, stringsAsFactors = FALSE)
  return(hemoglobin_genes)
}