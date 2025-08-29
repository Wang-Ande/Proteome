#' Analysis of Avariance, ANOVA 
#'

# Packages ----
library(readr)
library(openxlsx)
library(ggplot2)

# Data input ----
data <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
group <- read.xlsx("./01_Data/IC50_group.xlsx")
data <- data[,-grep("OCI_4W",colnames(data))]
group <- group[-grep("OCI_4W",group$id),]
head(data)
head(group)
rownames(group) <- group$id
data_target <- data[,grep("MV4", colnames(data))]
group_target <- group[grep("MV4", group$id), ]
# 创建分组因子
group_target$group <- factor(group_target$group, levels =  c("WT", "Low", "High"))  

# 转置表达矩阵：行为样本，列为蛋白
expr_matrix <- t(log2(data_target+1))        # log2 标准化处理

# 对比样本名和分组，需确保完全对应
identical(rownames(expr_matrix),rownames(group_target)) 

if(T){
# ----------------------------- ANOVA ------------------------------------------
anova_results <- apply(expr_matrix, 2, function(x) {
  fit <- aov(x ~ group_target$group)
  summary(fit)[[1]][["Pr(>F)"]][1]  # 提取 p 值
})
# ------------------------- Welch ANOVA -------------------------------------
welch_results <- apply(expr_matrix, 2, function(x) {
  oneway.test(x ~ group_target$group, var.equal = FALSE)$p.value
})
# ------------------------- Kruskal-Wallis -------------------------------------
kruskal_results <- apply(expr_matrix, 2, function(x) {
  kruskal.test(x ~ group_target$group)$p.value
})
# 汇总结果
stats_table <- data.frame(
  Protein_ID = colnames(expr_matrix),
  ANOVA_p = anova_results,
  Welch_p = welch_results,
  KW_p = kruskal_results
)

# 添加 FDR 校正
stats_table$ANOVA_FDR <- p.adjust(stats_table$ANOVA_p, method = "BH")
stats_table$Welch_FDR <- p.adjust(stats_table$Welch_p, method = "BH")
stats_table$KW_FDR <- p.adjust(stats_table$KW_p, method = "BH")

# 筛选显著蛋白（如 ANOVA 或 Kruskal-Wallis）
sig_anova <- subset(stats_table, ANOVA_p < 0.05)
sig_welch <- subset(stats_table, Welch_p < 0.05)
sig_kw    <- subset(stats_table, KW_p < 0.05)
}

# Res output ----
write.csv(stats_table , "./03_Result/ANOVA/MV4_11/Result_table.csv")
