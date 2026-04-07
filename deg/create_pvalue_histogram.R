library(ggplot2)

create_pvalue_histogram <- function(plot_dir, deg_results_file, artifacts_dir = "../data/") {

  print(paste("Creating p-value histogram for:", deg_results_file))
  # Load the full DEG results
  full_deg_results_file <- deg_results_file
  output_histogram_file <- paste0(artifacts_dir, plot_dir, "/pvalue_histogram.png")

  deg_results <- read.delim(full_deg_results_file, header = TRUE, stringsAsFactors = FALSE)

  # Create a histogram of p-values
  histogram_plot <- ggplot(deg_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.05, fill = "blue", color = "white", boundary = 0) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = paste0("Histogram of p-values from DESeq2 analysis for ", plot_dir),
        x = "p-value",
        y = "Frequency") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
    
  ggsave(output_histogram_file, width = 8, height = 6, plot = histogram_plot, create.dir = TRUE)
}