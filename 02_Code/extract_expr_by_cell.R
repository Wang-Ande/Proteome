process_cell_line <- function(x) {
  # 选取表达矩阵的列
  expr_subset <- expr_log2[, grep(x, colnames(expr_log2))]
  
  # 选取样本分组信息的行
  group_subset <- sample_group[grep(x, rownames(sample_group)), ]
  
  # 修改表达矩阵的列名
  colnames(expr_subset) <- group_subset$group
  
  # 赋值到全局变量
  assign(paste0("expr_", x), expr_subset, envir = .GlobalEnv)
}
