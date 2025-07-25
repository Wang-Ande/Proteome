#' Analysis of Avariance, ANOVA 
#'

# Packages ----
library(readr)
library(openxlsx)
library(ggplot2)

# Data input ----
data <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
group <- read.xlsx("./01_Data/IC50_group.xlsx")
head(data)
head(group)
rownames(group) <- group$id
data_target <- data[,grep("OCI", colnames(data))]
group_target <- group[grep("OCI", group$id), ]

# 创建分组因子
group_target$group <- factor(group_target$group, levels =  c("Ctrl", "Low", "High"))  

# 转置表达矩阵：行为样本，列为蛋白
expr_matrix <- t(log2(data_target))        # log2 标准化处理

# 对比样本名和分组，需确保完全对应
identical(rownames(expr_matrix),rownames(group_target)) 

if(T){
# ----------------------------- ANOVA ------------------------------------------
anova_results <- apply(expr_matrix, 2, function(x) {
  fit <- aov(x ~ group_target$group)
  summary(fit)[[1]][["Pr(>F)"]][1]  # 提取 p 值
})
# ------------------------- Kruskal-Wallis -------------------------------------
kruskal_results <- apply(expr_matrix, 2, function(x) {
  kruskal.test(x ~ group_target$group)$p.value
})

# 整理结果表格
stats_table <- data.frame(
  Protein_ID = colnames(expr_matrix),
  ANOVA_p = anova_results,
  KW_p = kruskal_results
)

# 添加 FDR 校正
stats_table$ANOVA_FDR <- p.adjust(stats_table$ANOVA_p, method = "BH")
stats_table$KW_FDR <- p.adjust(stats_table$KW_p, method = "BH")

# 筛选显著蛋白（如 ANOVA 或 Kruskal-Wallis）
sig_anova <- subset(stats_table, ANOVA_p < 0.05)
sig_kw    <- subset(stats_table, KW_p < 0.05)
}

# Res output ----
write.csv(stats_table , "./03_Result/ANOVA/OCI_AML2/Result_table.csv")
