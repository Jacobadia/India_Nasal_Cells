library(ggplot2)

create_volcano_plot <- function(plot_dir, deg_results_file, artifacts_dir = "../artifacts/") {
  print(paste("Creating volcano plot for:", deg_results_file))

  output_volcano_file <- paste0(artifacts_dir, plot_dir, "/volcano_plot.png")

  deg_results <- read.delim(deg_results_file, header = TRUE, stringsAsFactors = FALSE)

  # Support both DESeq2 (log2FoldChange) and limma (logFC) column naming
  if ("log2FoldChange" %in% colnames(deg_results)) {
    deg_results$logFC <- deg_results$log2FoldChange
  }

  deg_results$significant <- factor(
    !is.na(deg_results$padj) & deg_results$padj < 0.05,
    levels = c(FALSE, TRUE)
  )
  deg_results$neg_log10_pvalue <- -log10(deg_results$pvalue)

  volcano_plot <- ggplot(deg_results, aes(x = logFC, y = neg_log10_pvalue, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                       labels = c("FALSE" = "Not significant (padj \u2265 0.05)", "TRUE" = "Significant (padj < 0.05)"),
                       name = "Significance",
                       drop = FALSE) +
    labs(title = paste0("Volcano plot for ", plot_dir),
         x = "log2 Fold Change",
         y = expression(-log[10](p-value))) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(output_volcano_file, width = 8, height = 6, plot = volcano_plot)
}
