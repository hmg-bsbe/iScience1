# ============================================================
# Whole-transcriptome scatterplots and Pearson correlation
# ------------------------------------------------------------
# Input:
#   - 'vsd' : a DESeq2 variance-stabilized (VST) object created upstream

# Tested with:
#   R 4.5.1 on Ubuntu 22.04.5 LTS (see sessionInfo file)
# ============================================================

library(DESeq2)
library(ggplot2)

# -------------------------------
# 1) Load VST matrix from DESeq2 object
# -------------------------------
vst_obj <- vsd
mat <- assay(vst_obj)

# -------------------------------
# 2) Detect group labels from sample names
#    (expects sample names like "W8L8_1", "W8L8_2", etc.)
# -------------------------------
sample_names <- colnames(mat)
group <- sub("_[0-9]+$", "", sample_names)

cat("Group counts:\n")
print(table(group))

groups <- sort(unique(group))
cat("Detected groups:", paste(groups, collapse = ", "), "\n")

# -------------------------------
# 3) Keep expressed genes
#    (basic filter: mean VST > 0 across all samples)
# -------------------------------
keep <- rowMeans(mat, na.rm = TRUE) > 0
mat <- mat[keep, , drop = FALSE]
cat("Genes retained after expression filter:", nrow(mat), "\n")

# -------------------------------
# 4) Create output directory
# -------------------------------
outdir <- file.path(getwd(), "plots_allgenes")
dir.create(outdir, showWarnings = FALSE)
cat("Saving plots to:", outdir, "\n")

# -------------------------------
# 5) scatterplot for a single pair of groups
# -------------------------------
make_scatter <- function(g1, g2) {

  x <- rowMeans(mat[, group == g1, drop = FALSE], na.rm = TRUE)
  y <- rowMeans(mat[, group == g2, drop = FALSE], na.rm = TRUE)

  df <- data.frame(x = x, y = y)

  r <- cor(df$x, df$y, method = "pearson", use = "complete.obs")

  p <- ggplot(df, aes(x = x, y = y)) +

    # Points: light grey fill with black outline
    geom_point(
      shape = 21,
      fill  = "grey80",
      color = "black",
      size  = 1.5,
      stroke = 0.4,
      alpha = 0.9
    ) +

    # Regression line (for visualization only)
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +

    theme_classic(base_size = 14) +

    labs(
      title = paste0("All expressed genes: ", g1, " vs ", g2,
                     " | Pearson r = ", round(r, 3)),
      x = paste0("Mean VST expression (", g1, ")"),
      y = paste0("Mean VST expression (", g2, ")")
    )

  print(p)

  outfile <- file.path(outdir, paste0("Scatter_AllGenes_", g1, "_vs_", g2, ".png"))
  ggsave(filename = outfile, plot = p, width = 6, height = 5, dpi = 300)

  cat("Saved:", outfile, "| r =", round(r, 3), "\n")

  invisible(r)
}

# -------------------------------
# 6) Run all pairwise group comparisons
# -------------------------------
pairs <- combn(groups, 2, simplify = FALSE)
r_values <- sapply(pairs, function(pr) make_scatter(pr[1], pr[2]))

cat("\nPairwise Pearson r summary:\n")
for (i in seq_along(pairs)) {
  cat(pairs[[i]][1], "vs", pairs[[i]][2], ": r =", round(r_values[i], 3), "\n")
}

cat("\nGenerated plot files:\n")
print(list.files(outdir, pattern = "\\.png$", full.names = TRUE))

# -------------------------------
# 7) Reproducibility: save session information
# -------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo_WholeTranscriptomeScatterPlot.txt"))
