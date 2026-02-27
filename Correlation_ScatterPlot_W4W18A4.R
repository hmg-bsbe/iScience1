# ================================
# Tested with: R 4.5.1 on Ubuntu 22.04.5 LTS (see sessionInfo_*.txt output)

# Correlation + Scatter Plot
# ================================

library(readxl)
library(dplyr)
library(ggplot2)

# -------------------------------
# 1. Load Excel File
# -------------------------------

data <- read_excel("Final_Combined_ResultsW4W18A4.xlsx")

# -------------------------------
# 2. Convert Log2FC to numeric
#    (important for Linux + Excel imports)
# -------------------------------

data$W4W18_Log2FoldChange <- as.numeric(data$W4W18_Log2FoldChange)
data$W4A4_Log2FoldChange  <- as.numeric(data$W4A4_Log2FoldChange)

# -------------------------------
# 3. Filter Significant + Valid LogFC
# -------------------------------

sig <- data %>%
  filter(`P-Value [W4W18]` < 0.05 &
         `P-Value [W4A4]`  < 0.05 &
         !is.na(W4W18_Log2FoldChange) &
         !is.na(W4A4_Log2FoldChange))

# -------------------------------
# 4. Pearson Correlation
# -------------------------------

cor_all <- cor(sig$W4W18_Log2FoldChange,
               sig$W4A4_Log2FoldChange,
               method = "pearson",
               use = "complete.obs")

cat("\nPearson correlation (All significant genes):",
    round(cor_all,3), "\n")

# -------------------------------
# 5. Scatter Plot
# -------------------------------

p <- ggplot(sig, aes(x = W4W18_Log2FoldChange,
                     y = W4A4_Log2FoldChange)) +
  
  geom_point(color = "steelblue", size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) +
  
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = paste("Pearson r =", round(cor_all,3)),
    x = "Log2 Fold Change (W4W18)",
    y = "Log2 Fold Change (W4A4)"
  )

print(p)

# -------------------------------
# 6. Export Cleaned Data
# -------------------------------

write.csv(sig, "Significant_Both_Clean_W4W18A4.csv",
          row.names = FALSE)

# --------------------------------
# 7.Quality Control - to check how many genes were excluded:
#----------------------------------
cat("Total genes:", nrow(data), "\n")
cat("Significant + valid genes:", nrow(sig), "\n")
# ================================
# END
# ================================

# -------------------------------
# Reproducibility: save session information
# -------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo_Significant_coding_noncodingW4W18A4.txt")
