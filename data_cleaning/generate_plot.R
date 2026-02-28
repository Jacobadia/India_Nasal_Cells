library(ggplot2)
library(dplyr)
library(stringr)

# ── 1. Read & parse the data ──────────────────────────────────────────────────
input_file <- "fastqc_per_sequence_quality_scores_plot.txt"
output_file <- "fastqc_per_sequence_quality_scores.png"
raw <- readLines(input_file)

# Drop the header line
raw <- raw[-1]
raw <- raw[nzchar(trimws(raw))]

parse_row <- function(line) {
  # Sample name is everything before the first tab
  sample <- str_extract(line, "^[^\t]+")
  
  # Extract ALL (score, count) pairs from the entire line at once
  pairs  <- str_extract_all(line, "\\(([0-9.]+),\\s*([0-9.]+)\\)")[[1]]
  
  if (length(pairs) == 0) return(NULL)
  
  scores <- as.numeric(str_match(pairs, "\\(([0-9.]+),")[, 2])
  counts <- as.numeric(str_match(pairs, ",\\s*([0-9.]+)\\)")[, 2])
  
  data.frame(sample = sample, score = scores, count = counts,
             stringsAsFactors = FALSE)
}

df <- bind_rows(lapply(raw, parse_row))

cat("Parsed", nrow(df), "rows across", length(unique(df$sample)), "samples\n")
cat("Score range:", range(df$score), "\n")

# ── 2. Normalize to percentage within each sample ────────────────────────────

df <- df %>%
  group_by(sample) %>%
  mutate(pct = count / sum(count) * 100) %>%
  ungroup()

# ── 3. Plot ───────────────────────────────────────────────────────────────────

p <- ggplot(df, aes(x = score, y = count, group = sample)) +
  annotate("rect", xmin =  0, xmax = 20, ymin = -Inf, ymax = Inf, fill = "#ff0000", alpha = 0.25) +
  annotate("rect", xmin = 20, xmax = 27, ymin = -Inf, ymax = Inf, fill = "#ffbf00", alpha = 0.25) +
  annotate("rect", xmin = 27, xmax = 40, ymin = -Inf, ymax = Inf, fill = "#00ff08", alpha = 0.25) +
  geom_line(linewidth = 0.6, alpha = 0.6, color = "#2196F3") +
  scale_x_continuous(
    name   = "Accuracy Score (-10Log(Error Probability))",
    breaks = seq(0, 40, by = 5),
    limits = c(0, 40),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name   = "Number of Sequences (millions)",
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) x / 1e6
  ) +
  labs(
    title    = "RNA Sequence Quality Score Distribution",
    subtitle = sprintf("n = %d samples", length(unique(df$sample)))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(output_file, plot = p,
       width = 10, height = 5, dpi = 300)

message("Plot saved to", output_file)