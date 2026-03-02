get_counts_and_condition_data <- function(matrix_counts_file, metadata_file) {
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
  return (list(counts = counts, conditionData = conditionData))
}