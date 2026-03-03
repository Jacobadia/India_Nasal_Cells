source("./parse_counts.R")
source("./filtration.R")
source("./differential.R")
source("./create_pvalue_histogram.R")

artifacts_dir <- "../artifacts/"
matrix_counts_file <- paste0(artifacts_dir, "gene_counts_corrected.tsv")
metadata_file <- paste0(artifacts_dir, "metadata.tsv")
gene_type_file <- paste0(artifacts_dir, "genetype_lookup.txt")

create_test_func <- function(matrix_counts_file, metadata_file, gene_type_file, dir_name, deg_function, filter_function) {
    return (function() {
        full_dir_path <- paste0(artifacts_dir, dir_name)
        dir.create(full_dir_path, showWarnings = FALSE)
        full_deg_results_file <- paste0(full_dir_path, "/deg_results_full.txt")
        significant_deg_results_file <- paste0(full_dir_path, "/deg_results_significant.txt")
        counts <- get_counts_data(matrix_counts_file)
        counts <- filter_function(counts, gene_type_file)
        conditionData <- get_condition_data(counts, metadata_file)
        deg_function(counts, conditionData, 
        full_deg_results_file, significant_deg_results_file)
        create_pvalue_histogram(dir_name, full_deg_results_file)
    })
}

pure_no_control <- create_test_func(matrix_counts_file, metadata_file, gene_type_file, 
"pure_no_control", run_deg_control_nothing, function(counts, gene_type_file) {
    return (counts)
})

pure_control_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"pure_control_sex", run_deg_control_sex, function(counts, gene_type_file) {
    return (counts)
})

pure_control_age <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"pure_control_age", run_deg_control_age, function(counts, gene_type_file) {
    return (counts)
})

pure_control_age_and_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"pure_control_age_and_sex", run_deg_control_sex_and_age, function(counts, gene_type_file) {
    return (counts)
})

protein_mean_10_control_nothing <- create_test_func(matrix_counts_file, metadata_file, gene_type_file, 
"protein_mean_10_control_nothing", run_deg_control_nothing, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 10))
})

protein_mean_10_control_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file, 
"protein_mean_10_control_sex", run_deg_control_sex, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 10))
})

protein_mean_10_control_age <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_10_control_age", run_deg_control_age, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 10))
})

protein_mean_10_control_age_and_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_10_control_age_and_sex", run_deg_control_sex_and_age, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 10))
})

protein_mean_100_control_nothing <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_100_control_nothing", run_deg_control_nothing, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 100))
})

protein_mean_100_control_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_100_control_sex", run_deg_control_sex, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 100))
})

protein_mean_100_control_age <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_100_control_age", run_deg_control_age, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 100))
})

protein_mean_100_control_age_and_sex <- create_test_func(matrix_counts_file, metadata_file, gene_type_file,
"protein_mean_100_control_age_and_sex", run_deg_control_sex_and_age, function(counts, gene_type_file) {
    return (filter_protein_mean_counts(counts, get_genetype_lookup(gene_type_file), 100))
})


pure_no_control()
pure_control_sex()
pure_control_age()
pure_control_age_and_sex()
protein_mean_10_control_nothing()
protein_mean_10_control_sex()
protein_mean_10_control_age()
protein_mean_10_control_age_and_sex()
protein_mean_100_control_nothing()
protein_mean_100_control_sex()
protein_mean_100_control_age()
protein_mean_100_control_age_and_sex()