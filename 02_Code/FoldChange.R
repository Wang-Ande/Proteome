# Foldchange计算
# single sample foldchange
# foldchange定义为对照组和实验组中基因平均表达值的比值，单样本直接计算比值即可
library(openxlsx)
library(readr)
library(ggplot2)
library(ggrepel)
library(corrplot)
# 一、single sample logfc ----
## 1. data input ----
expr <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
rownames(expr) <- expr$X
expr <- expr[,-1]

# filter lowexpr genes
expr <- log2(expr+1)
expr_filter <- expr[apply(expr,1,mean)>1,]        # 根据histogram设定阈值

# 去除样本名多余信息减少绘图占用空间
colnames(expr_filter) <- gsub("cas9","",colnames(expr_filter)) 
colnames(expr_filter) <- gsub("w","W",colnames(expr_filter))

# 确保所有的实验组和对照组配平
expr_filter$MOLM13_2W_1 <- 1
expr_filter <- expr_filter[,-grep("4W",colnames(expr_filter))]
colnames(expr_filter)

## 2. subset symbol names ----
y <- rownames(expr_filter)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
expr_filter$gene <- gene
expr_filter <- expr_filter[, order(names(expr_filter))]

# 若进行log2，先转置回原格式
expr_filter <- 2^expr_filter[,2:28]-1
expr_filter <- expr_filter+1
expr <- expr_filter

## 3. creat fc df ----
# 定义细胞系和对应的时间点
cell_lines <- c("MOLM13", "MV4_11", "OCI_AML2")
time_points <- c("WT", "2W", "6W")

# 创建一个空的 fold_change 数据框
fold_change <- data.frame(Gene = gene)

# 计算单样本 fold change
for (cell_line in cell_lines) {
  # 获取对应细胞系的列名
  wt_samples <- grep(paste0(cell_line, "_WT"), colnames(expr), value = TRUE)
  two_w_samples <- grep(paste0(cell_line, "_2W"), colnames(expr), value = TRUE)
  six_w_samples <- grep(paste0(cell_line, "_6W"), colnames(expr), value = TRUE)
  
  # 单样本对照组 (WT) 和实验组 (2W, 6W) 之间的 Fold Change
  for (i in 1:length(two_w_samples)) {
    fold_change[paste0(cell_line, "_2W_", i, "_vs_WT_", i)] <- expr[[two_w_samples[i]]] / expr[[wt_samples[i]]]
  }
  
  for (i in 1:length(six_w_samples)) {
    fold_change[paste0(cell_line, "_6W_", i, "_vs_WT_", i)] <- expr[[six_w_samples[i]]] / expr[[wt_samples[i]]]
  }
}

fold_change[,c(2:19)] <- log2(fold_change[,c(2:19)]) 

## 4. res output ----
write.xlsx(fold_change,file = "./01_Data/Cpm_data/Foldchange.xlsx")

# foldchange比对
# transcriptome
transcriptome_fc <- read.xlsx("./01_data/Cpm_data/Foldchange.xlsx")

# proteome
proteome_fc <- read.xlsx("../Proteome/01_Data/Foldchange.xlsx")

# common set
common_set <- intersect(transcriptome_fc$Gene, proteome_fc$Gene)

# 定义数据处理函数
process_fc_data <- function(fc_data, common_set) {
  fc_data <- fc_data[fc_data$Gene %in% common_set, ]  # 筛选共同基因
  fc_data <- fc_data[!duplicated(fc_data$Gene), ]    # 去重
  rownames(fc_data) <- fc_data$Gene                  # 设置行名
  fc_data <- fc_data[, -1]                           # 去除 Gene 列
  fc_data[order(rownames(fc_data)), ]                # 按基因名排序
}

# 处理数据
transcriptome_fc_select <- process_fc_data(transcriptome_fc, common_set)
proteome_fc_select <- process_fc_data(proteome_fc, common_set)

# 合并
colnames(merged_df) <- c(paste0("Transcript_", seq_len(ncol(transcriptome_fc_select))), 
                         paste0("Protein_", seq_len(ncol(proteome_fc_select))))

## 5. fc cor analysis ----
correlation_matrix <- cor(transcriptome_fc_select, proteome_fc_select)

# plot
pdf("../FC/fc_cor_cpm.pdf", width = 8,height = 6)
colnames(correlation_matrix) <- c(paste0("Transcript_", c(1:17)))
rownames(correlation_matrix) <- c(paste0("Protein_", c(1:17)))
cor <- corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.7, tl.cex = 0.7)
dev.off()

## 6. scatter plot ----
library(ggplot2)
start_column <- 1
end_column <- start_column + 1     # 因为循环3次，所以结束列为 start_column + 2
dir_scatter <- "../FC/fc_scatter_cpm/"


for(i in start_column:end_column)  {
# 提取对应fc
scatter_df <- merged_df[,c(i, i+17)]                      # 每次取 i 和 i+17 作为一对
colnames(scatter_df) <- c("Transcript","Protein")
scatter_df$Transcript <- scale(scatter_df$Transcript)
scatter_df$Protein <- scale(scatter_df$Protein)
scatter_df <- as.data.frame(scatter_df)
scatter_df <- scatter_df[order(scatter_df$Transcript), ] 
cor_value <- cor(scatter_df$Transcript, scatter_df$Protein)

# 绘制散点图
# 创建 PDF 文件名和标题
pdf_filename <- paste0(dir_scatter,"fc_scatter_", i, ".pdf")  # 文件名：fc_scatter_1.pdf...
plot_title <- paste0("MOLM13_2W", i - start_column + 1, 
                     ".vs.WT_", i - start_column + 1)         # 图标题：MOLM13_6W_1.vs.WT_1

ggplot_plot <- ggplot(scatter_df, aes(x = Transcript, y = Protein)) +
  geom_point(shape = 1, size = 3) +  # x 点颜色
  labs(x = "Fc of Transcript", y = "Fc of Protein", title = plot_title) +
  theme_minimal() +
  theme(legend.position = "none")+  # 隐藏图例
  annotate("text", x = min(scatter_df$Transcript), y = max(scatter_df$Protein), 
           label = paste("Correlation: ", round(cor_value, 2)), hjust = 0, vjust = 1, size = 5, color = "black")
ggsave(filename = pdf_filename, plot = ggplot_plot, width = 8, height = 6 )
}

# 二、mean sample logfc ----
## 1. data_input ----
transcriptome_expr_select <- read.csv()
proteome_expr <- read.csv("../FC/Proteome_selected_expr.csv")

# 统一列名
colnames(proteome_expr)[1] <- "gene"
colnames(transcriptome_expr_select)[1] <- "gene"

# 获取共同基因集
common_set <- intersect(transcriptome_expr_select$gene, proteome_expr$gene)

# 定义数据处理函数
process_expr_data <- function(expr_data, common_set) {
  expr_data <- expr_data[expr_data$gene %in% common_set, ]  # 筛选共同基因
  expr_data <- expr_data[!duplicated(expr_data$gene), ]     # 去重
  rownames(expr_data) <- expr_data$gene                     # 设定行名
  expr_data <- expr_data[, -1]                              # 去掉 gene 列
  expr_data <- expr_data[order(rownames(expr_data)), ]      # 按基因排序
  expr_data[, order(colnames(expr_data))]                   # 按列名排序
}

# 处理蛋白组和转录组数据
proteome_expr_select <- process_expr_data(proteome_expr, common_set)
transcriptome_expr_select <- process_expr_data(transcriptome_expr_select, common_set)

# 函数：计算logFC
calculate_logfc <- function(expr_data, group1_pattern, group2_pattern) {
  
  # log2
  expr_data <- log2(expr_data)
  
  # 获取group1（如2W）和group2（如WT）的列名
  group1_samples <- grep(group1_pattern, colnames(expr_data), value = TRUE)
  group2_samples <- grep(group2_pattern, colnames(expr_data), value = TRUE)
  
  # 计算组内每个样本的均值
  group1_mean <- rowMeans(expr_data[, group1_samples], na.rm = TRUE)
  group2_mean <- rowMeans(expr_data[, group2_samples], na.rm = TRUE)
  
  # 计算log2 Fold Change (logFC)
  logfc <- group1_mean - group2_mean
  
  return(logfc)
}

cell_lines <- c("MOLM13", "MV4_11", "OCI_AML2")
time_points <- c("2W", "6W")

# 结果保存的列表
logfc_results <- list()

# 迭代处理每个细胞系和时间点

for (cell in cell_lines) {
  for (time_point in time_points) {
    # 构造2W和WT样本名称的模式
    group1_pattern <- paste0(cell, "_", time_point)
    group2_pattern <- paste0(cell, "_WT")
    
    # 提取和计算logFC
    logfc_results[[paste0(cell, "_", time_point,".vs.WT")]] <- calculate_logfc(transcriptome_expr_select, group1_pattern, group2_pattern)
  }
}

##  2. res output ----
transcriptome_logfc_results <- as.data.frame(logfc_results)
proteome_logfc_results <- as.data.frame(logfc_results)

# 合并
merged_expr_df <- cbind(transcriptome_logfc_results, proteome_logfc_results)
colnames(merged_expr_df) <- c(paste0("Trans_", colnames(transcriptome_logfc_results)), 
                         paste0("Prote_", colnames(proteome_logfc_results)))

## 3. cor analysis ----
library(corrplot)
expr_correlation_matrix <- cor(transcriptome_logfc_results, proteome_logfc_results)

# plot
pdf("../FC/fc_all_cpm_cor.pdf", width = 8,height = 6)
colnames(expr_correlation_matrix) <- c(paste0("Trans_", colnames(transcriptome_logfc_results)))
rownames(expr_correlation_matrix) <- c(paste0("Prote_", colnames(proteome_logfc_results)))
cor <- corrplot(expr_correlation_matrix, method = "color", type = "upper", 
                tl.col = "black", tl.srt = 45, addCoef.col = "black",
                number.cex = 0.7, tl.cex = 0.7)
dev.off()

# 三、expr cor ----
## 1. data_input ----
prote_expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
trans_expr <- read.csv("../Transcriptome/01_data/Cpm_data/IC50_merged_adjusted.csv",row.names = 1)
#colnames(trans_expr) <- gsub("cas9","",colnames(trans_expr))
prote_anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
prote_anno <- prote_anno[rownames(prote_anno)%in%rownames(prote_expr),]

# subset genenames
y <- rownames(trans_expr)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
trans_expr$gene <- gene 

# 将rownames由protein.group改为gene
rownames(trans_expr) <- trans_expr$gene    # gene做为行名重复，进行去重处理

tmp <- trans_expr
if(T){
  tmp <- tmp[!is.na(tmp$gene), ]                       # remove NA rows
  tmp <- aggregate(. ~ gene, data = tmp, FUN = max)    # remove repeat genes
  rownames(tmp) <- tmp$gene                            # set rownames
  tmp <- tmp[,-1]
}
trans_expr <- tmp

# 取交集
common_set <- intersect(rownames(trans_expr), rownames(prote_expr))

# 定义数据处理函数
process_expr_data <- function(expr_data, common_set) {
  expr_data <- expr_data[rownames(expr_data) %in% common_set, ]  # 筛选共同基因
  expr_data <- expr_data[order(rownames(expr_data)), ]           # 按基因排序
  expr_data[, order(colnames(expr_data))]                        # 按列名排序
}

# 处理蛋白组和转录组数据
prote_expr_select <- process_expr_data(prote_expr, common_set)
trans_expr_select <- process_expr_data(trans_expr, common_set)
# write.csv(prote_expr_select, file = "../FC/Prote_selected_expr.csv")
prote_log2 <- log2(prote_expr_select + 1) 
trans_log2 <- log2(trans_expr_select + 1)

# scale 标准化
prote_scaled <- scale(prote_log2)
trans_scaled <- scale(trans_log2)

## 2. cor analysis ---- 
colnames(prote_log2) <- c(paste0("P_", colnames(prote_log2)))
colnames(trans_log2) <- c(paste0("T_", colnames(trans_log2)))
merge_expr_select <- merge(prote_log2,trans_log2,by = 0)
rownames(merge_expr_select) <- merge_expr_select$Row.names
merge_expr_select <- subset(merge_expr_select,select = -c(Row.names))

# cor_matrix <- cor(trans_log2,prote_log2)
corr_matrix <- cor(prote_log2,trans_log2)

## 3. plot ----
# two visualization methods
# methods 1
pdf("../FC/P&T_log2_cor.pdf", width = 15,height = 12)
cor <- corrplot(corr_matrix, method = "color", type = "upper", 
                tl.col = "black", tl.srt = 45, addCoef.col = "black",
                number.cex = 0.7, tl.cex = 0.7)
dev.off()

# methods 2
# 转换相关性矩阵为长格式
library(reshape2)
library(ggplot2)
corr_long <- melt(corr_matrix)

# 使用层次聚类重新排序列名
hc <- hclust(dist(corr_matrix))
ordered_names <- rownames(corr_matrix)[hc$order]

# 重新调整数据框的因子顺序
corr_long$Var1 <- factor(corr_long$Var1, levels = ordered_names)
corr_long$Var2 <- factor(corr_long$Var2, levels = ordered_names)

# plot
ggplot(corr_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white","#AED4E5","#81B5D5","#5795C7","#3371B3","#345D82","#1E4C9C"), 
                       limits = c(0, 1),     # 设定颜色映射范围
                       name = expression(R^2)) +            # 更改图例标题为 R²
  labs(title = "Pearson correlation") +     # 添加标题
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
        axis.title = element_blank(),                       # 隐藏 X 和 Y 轴标题
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中、加粗
        legend.position = "right")                          # 保持图例在右侧
ggsave("../FC/P&T_log2_cor_1.pdf", width = 8, height = 7)
