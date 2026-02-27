# ============================================================
# Whole-transcriptome heatmap, correlations and PCA 
# Top variable genes heatmap: Top 10,000
# ============================================================

# -------------------------------
# 1) Packages required
# -------------------------------
req_pkgs <- c("tximport", "DESeq2", "GenomicFeatures", "txdbmaker",
              "edgeR", "vegan", "pheatmap", "matrixStats")

missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nPlease install them before running this script.")
}

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(GenomicFeatures)
  library(txdbmaker)
  library(edgeR)
  library(vegan)
  library(pheatmap)
  library(matrixStats)
})

cat("Working directory:", getwd(), "\n")

# -------------------------------
# 2) Locate GTF
# -------------------------------
gtf_file <- list.files(pattern = "\\.gtf(\\.gz)?$", ignore.case = TRUE)
if (length(gtf_file) == 0) {
  stop("No GTF found in the working directory. Please place a *.gtf or *.gtf.gz file here.")
}
gtf_file <- gtf_file[1]
cat("Using GTF:", gtf_file, "\n")

# -------------------------------
# 3) Locate Salmon quant files (generic)
#    A) recursive quant.sf (typical)
#    B) fallback: flat *.sf in current directory
# -------------------------------
quant_sf_paths <- list.files(pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
flat_sf_paths  <- list.files(pattern = "\\.sf$", recursive = FALSE, full.names = TRUE)

if (length(quant_sf_paths) > 0) {
  files_vec  <- quant_sf_paths
  sample_ids <- basename(dirname(files_vec))  # sample = folder name
  cat("Detected quant.sf files:", length(files_vec), "\n")
} else if (length(flat_sf_paths) > 0) {
  files_vec  <- flat_sf_paths
  sample_ids <- sub("\\.sf$", "", basename(files_vec))
  cat("Detected flat .sf files:", length(files_vec), "\n")
} else {
  stop("No Salmon quantification files found. Expected quant.sf (recursive) or *.sf in current directory.")
}

names(files_vec) <- sample_ids

# -------------------------------
# 4) Sample table (Condition inferred from sample IDs)
#    Expected sample IDs like W8_1, W22_2, M8_3 -> Condition = W8, W22, M8
# -------------------------------
Condition <- sub("_[0-9]+$", "", sample_ids)

samples <- data.frame(
  sample = sample_ids,
  Condition = factor(Condition),
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

cat("Detected conditions:\n")
print(table(samples$Condition))



# -------------------------------
# 5) Build tx2gene from GTF
# -------------------------------
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
tx2gene$GENEID <- sub("\\..*", "", tx2gene$GENEID)

# -------------------------------
# 6) tximport (gene-level)
# -------------------------------
txi <- tximport(
  files_vec,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE,
  ignoreAfterBar = TRUE
)

# -------------------------------
# 7) DESeq2 + VST
# -------------------------------
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
save(dds, vsd, file = "dds_vsd.RData")
cat("Saved: dds_vsd.RData\n")

# -------------------------------
# 8) edgeR TMM normalization -> log2CPM matrix (M)
# -------------------------------
cts <- counts(dds)
y <- DGEList(counts = cts)
y <- calcNormFactors(y, method = "TMM")
M <- cpm(y, log = TRUE, prior.count = 1)  # log2CPM (can include negatives)

# match sample order
M <- M[, rownames(samples), drop = FALSE]

ann <- data.frame(Condition = samples$Condition)
rownames(ann) <- rownames(samples)

# -------------------------------
# 9) PERMANOVA (whole transcriptome)
# -------------------------------
dist_all <- dist(t(M))
adon_all <- adonis2(dist_all ~ Condition, data = ann, permutations = 999)

sink("WholeTx_PERMANOVA_TMMlog2CPM.txt")
print(adon_all)
sink()
cat("Saved: WholeTx_PERMANOVA_TMMlog2CPM.txt\n")

# Pairwise PERMANOVA 
conds <- levels(ann$Condition)
if (length(conds) >= 2) {
  sink("WholeTx_PERMANOVA_pairwise.txt")
  for (i in 1:(length(conds) - 1)) {
    for (j in (i + 1):length(conds)) {
      g1 <- conds[i]; g2 <- conds[j]
      keep <- ann$Condition %in% c(g1, g2)
      da <- dist(t(M[, keep, drop = FALSE]))
      a2 <- adonis2(da ~ Condition, data = ann[keep, , drop = FALSE], permutations = 999)
      cat(sprintf("%s vs %s\n", g1, g2))
      print(a2); cat("\n")
    }
  }
  sink()
  cat("Saved: WholeTx_PERMANOVA_pairwise.txt\n")
}

# -------------------------------
# 10) Mean cross-group correlations
# -------------------------------
cor_mat <- cor(M, use = "pairwise.complete.obs")

by_group <- function(g1, g2) {
  s1 <- rownames(ann)[ann$Condition == g1]
  s2 <- rownames(ann)[ann$Condition == g2]
  mean(cor_mat[s1, s2], na.rm = TRUE)
}

if (length(conds) >= 2) {
  out <- data.frame()
  for (i in 1:(length(conds) - 1)) {
    for (j in (i + 1):length(conds)) {
      out <- rbind(out, data.frame(
        Group1 = conds[i],
        Group2 = conds[j],
        MeanCorrelation = by_group(conds[i], conds[j])
      ))
    }
  }
  write.csv(out, "group_mean_correlations.csv", row.names = FALSE)
  save(out, file = "group_mean_correlations.RData")
  cat("Saved: group_mean_correlations.csv and group_mean_correlations.RData\n")
}

# -------------------------------
# 11) Top variable genes heatmap (Top 10,000)
# -------------------------------
n_top <- 10000
rv <- rowVars(M)
topVar <- head(order(rv, decreasing = TRUE), min(n_top, nrow(M)))

pdf("WholeTx_topVar10000_rowsCLUSTER_colsFIXED_TMMlog2CPM.pdf", width = 8, height = 8)
pheatmap(
  M[topVar, , drop = FALSE],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = ann,
  show_rownames = FALSE,
  color = colorRampPalette(c("navy","white","firebrick3"))(101),
  main = "Top variable genes (Top 10,000; rows clustered; columns fixed)"
)
dev.off()

cat("Saved: WholeTx_topVar10000_rowsCLUSTER_colsFIXED_TMMlog2CPM.pdf\n")

library(ggplot2)

# ----------------------------------------------------
# PCA on ALL genes
# ---------------------------------------------------
pca <- prcomp(t(M), scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Condition = ann$Condition
)

percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    title = "PCA – Whole transcriptome (TMM log2CPM)"
  )

ggsave("WholeTranscriptome_PCA.pdf", width = 6, height = 5)

# -------------------------------
# 14) Reproducibility
# -------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo_WholeTxPipeline.txt")
cat("Saved: sessionInfo_WholeTxPipeline.txt\n")

cat("\nDone.\n")
