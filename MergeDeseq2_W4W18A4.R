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

w4w18 <- read_excel("W4W18_AlzC.xlsx")
w4a4  <- read_excel("W4A4_Alz.xlsx")

# -------------------------------
# 2. Select Required Columns
# -------------------------------

w4w18_sel <- w4w18 %>%
  select(gene_name, log2FoldChange, pvalue)

w4a4_sel <- w4a4 %>%
  select(gene_name, log2FoldChange, pvalue)

# -------------------------------
# 3. Rename Columns
# -------------------------------

colnames(w4w18_sel) <- c("Gene_name",
                          "W4W18_Log2FoldChange",
                          "P-Value [W4W18]")

colnames(w4a4_sel) <- c("Gene_name",
                         "W4A4_Log2FoldChange",
                         "P-Value [W4A4]")

# -------------------------------
# 4. Merge Tables By Gene Name
# -------------------------------

merged <- inner_join(w4w18_sel, w4a4_sel, by = "Gene_name")

# -------------------------------
# 5. Compute Significant_Both
# -------------------------------

merged <- merged %>%
  mutate(Significant_Both = ifelse(`P-Value [W4W18]` < 0.05 &
                                     `P-Value [W4A4]` < 0.05,
                                   "YES", "NO"))

# -------------------------------
# 6. Reorder Columns
# -------------------------------

final <- merged %>%
  select(Gene_name,
         W4W18_Log2FoldChange, `P-Value [W4W18]`,
         W4A4_Log2FoldChange,  `P-Value [W4A4]`,
         Significant_Both)

# -------------------------------
# 7. Export Final Excel File
# -------------------------------

write_xlsx(final, "Final_Combined_ResultsW4W18A4.xlsx")

# ================================
# END
# ================================

# -------------------------------
# Reproducibility: save session information
# -------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo_Dseq-filtering-W4W18A4filtering.txt")
