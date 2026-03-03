source("./parse_counts.R")
source("./filtration.R")
source("./differential.R")
source("./create_pvalue_histogram.R")

artifacts_dir <- "../artifacts/"
matrix_counts_file <- paste0(artifacts_dir, "gene_counts_corrected.tsv")
metadata_file <- paste0(artifacts_dir, "metadata.tsv")
gene_type_file <- paste0(artifacts_dir, "genetype_lookup.txt")

# counts <- get_counts_data(matrix_counts_file)
# conditionData <- get_condition_data(counts, metadata_file)

# genetype_lookup <- get_genetype_lookup(gene_type_file)
# counts <- filter_total_counts(counts, 0)

# counts <- filter_protein_coding_genes(counts, genetype_lookup)

# counts <- filter_mean_counts(counts, 10)

pure_no_control <- function(matrix_counts_file, metadata_file, gene_type_file) {
    dir_name <- "pure_no_control"
    full_dir_path <- paste0(artifacts_dir, dir_name)

    dir.create(full_dir_path, showWarnings = FALSE)

    full_deg_results_file <- paste0(full_dir_path, "/deg_results_full.txt")
    significant_deg_results_file <- paste0(full_dir_path, "/deg_results_significant.txt")
    counts <- get_counts_data(matrix_counts_file)
    conditionData <- get_condition_data(counts, metadata_file)
    run_deg_control_nothing(counts, conditionData, 
    full_deg_results_file, significant_deg_results_file)
    create_pvalue_histogram(dir_name, full_deg_results_file)
}

get_protein_coding_non_zero <- function(matrix_counts_file, gene_type_file) {
    counts <- get_counts_data(matrix_counts_file)
    genetype_lookup <- get_genetype_lookup(gene_type_file)
    counts <- filter_total_counts(counts, 0)
    counts <- filter_protein_coding_genes(counts, genetype_lookup)
    return (counts)
}

protein_mean_10 <- function(matrix_counts_file, metadata_file, gene_type_file) {
    dir_name <- "protein_mean_10"
    full_dir_path <- paste0(artifacts_dir, dir_name)
    dir.create(full_dir_path, showWarnings = FALSE)
    full_deg_results_file <- paste0(full_dir_path, "/deg_results_full.txt")
    significant_deg_results_file <- paste0(full_dir_path, "/deg_results_significant.txt")
    counts <- get_protein_coding_non_zero(matrix_counts_file, gene_type_file)
    counts <- filter_mean_counts(counts, 10)
    conditionData <- get_condition_data(counts, metadata_file)
    run_deg_control_nothing(counts, conditionData, 
    full_deg_results_file, significant_deg_results_file)
    create_pvalue_histogram(dir_name, full_deg_results_file)
}

pure_no_control(matrix_counts_file, metadata_file, gene_type_file)