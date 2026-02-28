library(ggplot2)

output_plot_name <- "pvalue_histogram_controlling_for_sex"

# Load the full DEG results
full_deg_results_file <- "../artifacts/deg_results_full.txt"
output_histogram_file <- paste0("../artifacts/", output_plot_name, ".png")

deg_results <- read.delim(full_deg_results_file, header = TRUE, stringsAsFactors = FALSE)

# Create a histogram of p-values
histogram_plot <- ggplot(deg_results, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "white", boundary = 0) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Histogram of p-values from DESeq2 analysis Controlling for Sex",
       x = "p-value",
       y = "Frequency") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
  
ggsave(output_histogram_file, width = 8, height = 6, plot = histogram_plot)