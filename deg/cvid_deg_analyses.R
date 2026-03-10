source("./parse_counts.R")
source("./filtration.R")
source("./differential.R")
source("./create_pvalue_histogram.R")
source("./create_volcano_plot.R")

artifacts_dir <- "../cvid_artifacts/"
matrix_counts_file <- paste0(artifacts_dir, "features_combined.txt")
metadata_file <- paste0(artifacts_dir, "meta_data.txt")

create_test_func <- function(matrix_counts_file, metadata_file, dir_name, deg_function) {
    return (function() {
        print("Running this DEG analysis")
        full_dir_path <- paste0(artifacts_dir, dir_name)
        dir.create(full_dir_path, showWarnings = FALSE)
        full_deg_results_file <- paste0(full_dir_path, "/deg_results_full.txt")
        significant_deg_results_file <- paste0(full_dir_path, "/deg_results_significant.txt")
        counts <- get_counts_data_cvim(matrix_counts_file)
        conditionData <- get_condition_data_cvim(counts, metadata_file)
        deg_function(counts, conditionData,
        full_deg_results_file, significant_deg_results_file)
        create_pvalue_histogram(dir_name, full_deg_results_file, artifacts_dir)
        create_volcano_plot(dir_name, full_deg_results_file, artifacts_dir)
    })
}

run_deg_control_nothing_test <- create_test_func(matrix_counts_file, metadata_file, "deg_control_nothing", run_deg_control_nothing)
run_deg_control_nothing_test()