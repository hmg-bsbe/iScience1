# ================================
# Combine Two DESeq2 Results
# Tested with: R 4.5.1 on Ubuntu 22.04.5 LTS (see sessionInfo_*.txt output)
# ================================


library(readxl)
library(dplyr)
library(writexl)

# -------------------------------
# 1. Read Excel Sheets
# -------------------------------

w2w10 <- read_excel("W2W10_HDC.xlsx")
w2h2  <- read_excel("W2H2_HD.xlsx")

# -------------------------------
# 2. Select Required Columns
# -------------------------------

w2w10_sel <- w2w10 %>%
  select(gene_name, log2FoldChange, pvalue)

w2h2_sel <- w2h2 %>%
  select(gene_name, log2FoldChange, pvalue)

# -------------------------------
# 3. Rename Columns
# -------------------------------

colnames(w2w10_sel) <- c("Gene_name",
                          "W2W10_Log2FoldChange",
                          "P-Value [W2W10]")

colnames(w2h2_sel) <- c("Gene_name",
                         "W2H2_Log2FoldChange",
                         "P-Value [W2H2]")

# -------------------------------
# 4. Merge Tables By Gene Name
# -------------------------------

merged <- inner_join(w2w10_sel, w2h2_sel, by = "Gene_name")

# -------------------------------
# 5. Compute Significant_Both
# -------------------------------

merged <- merged %>%
  mutate(Significant_Both = ifelse(`P-Value [W2W10]` < 0.05 &
                                     `P-Value [W2H2]` < 0.05,
                                   "YES", "NO"))

# -------------------------------
# 6. Reorder Columns
# -------------------------------

final <- merged %>%
  select(Gene_name,
         W2W10_Log2FoldChange, `P-Value [W2W10]`,
         W2H2_Log2FoldChange,  `P-Value [W2H2]`,
         Significant_Both)

# -------------------------------
# 7. Export Final Excel File
# -------------------------------

write_xlsx(final, "Final_Combined_ResultsW2W10H2.xlsx")

# ================================
# END
# ================================

# -------------------------------
# Reproducibility: save session information
# -------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo_Dseq-filtering-W2W10H2filtering.txt")
