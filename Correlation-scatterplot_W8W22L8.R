# ================================
# W8W22 vs W8l8: Merge + All-genes plot + SigBoth coding/noncoding plots
# Tested with: R 4.5.1 on Ubuntu 22.04.5 LTS (see sessionInfo_*.txt output)

# ================================
# W8W22 vs W8l8: Merge + All-genes plot + SigBoth coding/noncoding plots
# ================================

library(readxl)
library(dplyr)
library(ggplot2)

# -------------------------------
# 1) Read Excel Sheets (ALL genes)
# -------------------------------
w8w22 <- read_excel("W8W22.xlsx")
w8l8  <- read_excel("W8L8.xlsx")

# -------------------------------
# 2) Select required columns
#    (Include gene_biotype as Type if present)
# -------------------------------
# If your file has gene_biotype, keep it. If not, we'll set Type = NA.

if ("gene_biotype" %in% colnames(w8w22)) {
  w8w22_sel <- w8w22 %>% select(gene_name, log2FoldChange, pvalue, gene_biotype)
  colnames(w8w22_sel) <- c("Gene_name", "W8W22_Log2FoldChange", "P-Value [W8W22]", "Type")
} else {
  w8w22_sel <- w8w22 %>% select(gene_name, log2FoldChange, pvalue)
  colnames(w8w22_sel) <- c("Gene_name", "W8W22_Log2FoldChange", "P-Value [W8W22]")
  w8w22_sel$Type <- NA_character_
}

w8l8_sel <- w8l8 %>% select(gene_name, log2FoldChange, pvalue)
colnames(w8l8_sel) <- c("Gene_name", "W8L8_Log2FoldChange", "P-Value [W8L8]")

# -------------------------------
# 3) Merge (ALL genes present in BOTH sheets)
#    inner_join = intersection by gene_name
# -------------------------------
merged_all <- inner_join(w8w22_sel, w8l8_sel, by = "Gene_name")

# -------------------------------
# 4) Coerce to numeric and remove missing log2FC (for correlation/plot)
# -------------------------------
merged_all <- merged_all %>%
  mutate(
    W8W22_Log2FoldChange = as.numeric(W8W22_Log2FoldChange),
    W8L8_Log2FoldChange  = as.numeric(W8L8_Log2FoldChange)
  )

n_total_merged <- nrow(merged_all)

merged_all_clean <- merged_all %>%
  filter(!is.na(W8W22_Log2FoldChange) & !is.na(W8L8_Log2FoldChange))

n_used_all <- nrow(merged_all_clean)
n_excluded_all <- n_total_merged - n_used_all

cat("\n=== ALL GENES (no significance filter) ===\n")
cat("Merged genes (present in both sheets):", n_total_merged, "\n")
cat("Used for correlation/plot (both log2FC present):", n_used_all, "\n")
cat("Excluded due to missing/non-numeric log2FC:", n_excluded_all, "\n")

# Pearson for all genes
cor_allgenes <- cor(merged_all_clean$W8W22_Log2FoldChange,
                    merged_all_clean$W8L8_Log2FoldChange,
                    method = "pearson",
                    use = "complete.obs")
cat("Pearson r (all genes):", round(cor_allgenes, 3), "\n")

# -------------------------------
# 5) Add Significant_Both column (still computed on full merged table)
# -------------------------------
merged_all_clean <- merged_all_clean %>%
  mutate(Significant_Both = ifelse(`P-Value [W8W22]` < 0.05 & `P-Value [W8L8]` < 0.05,
                                  "YES", "NO"))

# -------------------------------
# 6) Plot 1: ALL genes scatter (no filtering)
# -------------------------------
p_all <- ggplot(merged_all_clean, aes(x = W8W22_Log2FoldChange, y = W8L8_Log2FoldChange)) +
  geom_point(color = "steelblue", size = 1.6, alpha = 0.5) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("All genes | r = ", round(cor_allgenes,3),
                   " | excluded = ", n_excluded_all),
    x = "Log2FC (W8W22)",
    y = "Log2FC (W8L8)"
  )

print(p_all)
ggsave("Scatter_AllGenes_W8W22_vs_W8L8.png", p_all, width = 6, height = 5, dpi = 300)

# -------------------------------
# 7) Significant in BOTH (p<0.05) then split coding/noncoding
# -------------------------------
sig_both <- merged_all_clean %>% filter(Significant_Both == "YES")

coding <- sig_both %>% filter(Type == "protein_coding")
noncoding <- sig_both %>% filter(is.na(Type) | Type != "protein_coding")

# Pearson for coding/noncoding (guard against too-few rows)
cor_safe <- function(df) {
  if (nrow(df) < 3) return(NA_real_)
  cor(df$W8W22_Log2FoldChange, df$W8L8_Log2FoldChange, method = "pearson", use = "complete.obs")
}

cor_coding <- cor_safe(coding)
cor_noncoding <- cor_safe(noncoding)

cat("\n=== SIGNIFICANT IN BOTH (p<0.05) ===\n")
cat("Total significant-both genes:", nrow(sig_both), "\n")
cat("Protein-coding:", nrow(coding), " | r =", ifelse(is.na(cor_coding), "NA", round(cor_coding,3)), "\n")
cat("Non-coding/other:", nrow(noncoding), " | r =", ifelse(is.na(cor_noncoding), "NA", round(cor_noncoding,3)), "\n")

# -------------------------------
# 8) Plot 2: Protein-coding (significant in both)
# -------------------------------
p1 <- ggplot(coding, aes(x = W8W22_Log2FoldChange, y = W8L8_Log2FoldChange)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("Protein coding (sig in both) | r = ",
                   ifelse(is.na(cor_coding), "NA", round(cor_coding,3))),
    x = "Log2FC (W8W22)",
    y = "Log2FC (W8L8)"
  )

print(p1)
ggsave("Scatter_SigBoth_ProteinCoding_W8W22_vs_W8L8.png", p1, width = 6, height = 5, dpi = 300)

# -------------------------------
# 9) Plot 3: Non-coding/other (significant in both)
# -------------------------------
p2 <- ggplot(noncoding, aes(x = W8W22_Log2FoldChange, y = W8L8_Log2FoldChange)) +
  geom_point(color = "tomato", size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("Non-coding/other (sig in both) | r = ",
                   ifelse(is.na(cor_noncoding), "NA", round(cor_noncoding,3))),
    x = "Log2FC (W8W22)",
    y = "Log2FC (W8L8)"
  )

print(p2)
ggsave("Scatter_SigBoth_NonCoding_W8W22_vs_W8L8.png", p2, width = 6, height = 5, dpi = 300)

# -------------------------------
# 10) Optional: export tables
# -------------------------------
write.csv(merged_all_clean, "Merged_AllGenes_Clean_W8W22_vs_W8L8.csv", row.names = FALSE)
write.csv(sig_both, "Merged_SigBoth_W8W22_vs_W8L8.csv", row.names = FALSE)
write.csv(coding, "Merged_SigBoth_ProteinCoding_W8W22_vs_W8L8.csv", row.names = FALSE)
write.csv(noncoding, "Merged_SigBoth_NonCoding_W8W22_vs_W8L8.csv", row.names = FALSE)

# -------------------------------
# Reproducibility: save session information
# -------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo_CombinedW22L8M8.txt")
