# 必要包
library(car)
library(dplyr)


# data：行为蛋白，列为样本的 log2 表达矩阵（data.frame or matrix）
# group：因子变量，表示每列样本的组别，长度 == ncol(data)
data <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
group <- read.xlsx("./01_Data/IC50_group.xlsx")
head(data)
head(group)
data_target <- data[,grep("MV4_11", colnames(data))]
group_target <- group[grep("MV4_11", group$id),]

# 创建分组因子
group_target <- factor(group_target$group, levels =  c("Ctrl", "Low", "High"))  

# 转置表达矩阵：行为样本，列为蛋白
expr_matrix <- t(log2(data_target))        # log2 标准化处理

# 对比样本名和分组，需确保完全对应
data.frame(
  Sample = rownames(expr_matrix),
  Group = group_target
)


# 准备结果存储表
check_results <- data.frame(
  Protein_ID = colnames(expr_matrix),
  normality_pass = NA,
  variance_homogeneity_pass = NA,
  recommended_method = NA,
  stringsAsFactors = FALSE
)

# 循环检测每个蛋白
for (i in seq_len(ncol(expr_matrix))) {
  protein_expr <- expr_matrix[, i]
  
  # 分组表达
  group_expr <- split(protein_expr, group_target)
  
  # 正态性检验（每组都要满足）
  normality_p <- sapply(group_expr, function(x) {
    if (length(x) < 3) return(NA)  # Shapiro 最少3个值
    shapiro.test(x)$p.value
  })
  
  normality_pass <- all(normality_p > 0.05, na.rm = TRUE)
  
  # 方差齐性检验（Levene's test）
  levene_p <- tryCatch({
    leveneTest(protein_expr ~ group_target)$`Pr(>F)`[1]
  }, error = function(e) NA)
  
  variance_pass <- !is.na(levene_p) && levene_p > 0.05
  
  # 推荐方法判断
  method <- if (normality_pass && variance_pass) "ANOVA" else "Kruskal-Wallis"
  
  # 填入结果
  check_results[i, c("normality_pass", "variance_homogeneity_pass", "recommended_method")] <-
    list(normality_pass, variance_pass, method)
}

# 输出结果
table(check_results$recommended_method)
table(check_results$normality_pass)            # 单独看正态性检验
table(check_results$variance_homogeneity_pass) # 单独看方差齐性检验

# 筛选出可用于ANOVA的蛋白 （满足正态性和方差齐性）
anova_proteins <- check_results %>% filter(recommended_method == "ANOVA")

# 筛选出建议用Kruskal-Wallis的蛋白
nonparametric_proteins <- check_results %>% filter(recommended_method == "Kruskal-Wallis")

# 保存为csv（可选）
write.csv(check_results, "./03_Result/ANOVA/MV4_11/anova_check_results.csv", row.names = FALSE)
