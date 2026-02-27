# ================================
# Combine Two DESeq2 Results
# ================================

library(readxl)
library(dplyr)
library(writexl)

# -------------------------------
# 1. Read Excel Sheets
# -------------------------------

w8w22 <- read_excel("W8W22_Dseq_linux.xlsx")
w8m8  <- read_excel("W8M8_Dseq_linux.xlsx")

# -------------------------------
# 2. Select Required Columns
# -------------------------------

w8w22_sel <- w8w22 %>%
  select(gene_name, log2FoldChange, pvalue, gene_biotype)

w8m8_sel <- w8m8 %>%
  select(gene_name, log2FoldChange, pvalue)

# -------------------------------
# 3. Rename Columns
# -------------------------------

colnames(w8w22_sel) <- c("Gene_name",
                          "W8W22_Log2FoldChange",
                          "P-Value [W8W22]",
                          "Type")

colnames(w8m8_sel) <- c("Gene_name",
                         "W8M8_Log2FoldChange",
                         "P-Value [W8M8]")

# -------------------------------
# 4. Merge Tables By Gene Name
# -------------------------------

merged <- inner_join(w8w22_sel, w8m8_sel, by = "Gene_name")

# -------------------------------
# 5. Compute Significant_Both
# -------------------------------

merged <- merged %>%
  mutate(Significant_Both = ifelse(`P-Value [W8W22]` < 0.05 &
                                     `P-Value [W8M8]` < 0.05,
                                   "YES", "NO"))

# -------------------------------
# 6. Reorder Columns
# -------------------------------

final <- merged %>%
  select(Gene_name,
         W8W22_Log2FoldChange, `P-Value [W8W22]`,
         W8M8_Log2FoldChange,  `P-Value [W8M8]`,
         Significant_Both,
         Type)

# -------------------------------
# 7. Export Final Excel File
# -------------------------------

write_xlsx(final, "Final_CombinedW8W22M8_Results.xlsx")

# ================================
# END
# ================================
