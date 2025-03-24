# PCA_Track
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(dplyr)

#  样本信息
sample_info <- read.xlsx("./01_Data/IC50_group.xlsx")
rownames(sample_info) <- sample_info$ID
target_sample <- sample_info[grep("MOLM13",rownames(sample_info)),]

# 表达矩阵 data_matrix (行为基因，列为样本)
data_matrix <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
data_matrix <- log2(data_matrix)

# 提取目标分组样本
target_matrix <- data_matrix[,rownames(target_sample)]

#  获取 PCA 结果
pca_res <- prcomp(t(target_matrix), scale. = TRUE)

# 计算主成分贡献值（方差百分比）
var_percent <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x)
pca_df$ID <- rownames(pca_df)

#  结合样本信息
pca_df <- merge(pca_df, target_sample, by.y = "ID")

#  提取样本分组信息（颜色、形状）
pca_df$Sample_Pair <- gsub("_(WT|2W|6W)", "", pca_df$ID)  # 只保留 MOLM13_x 作为配对标识

color_palette <- c("MOLM13_1" = "#66C1A4", 
                   "MOLM13_2" = "#8C9FCA", 
                   "MOLM13_3" = "#FB8C63")
shape_values <- c("High" = 15, 
                  "Low" = 17, 
                  "Ctrl" = 16)

#  绘制 PCA 图
dir_pca <- "./03_Result/QC/MOLM13/"
# 确保数据按照 Sample_Pair 和 Group (Ctrl → Low → High) 排序
pca_df$Group <- factor(pca_df$Group, levels = c("Ctrl", "Low", "High"))
pca_df <- pca_df %>%
  arrange(Sample_Pair,Group)

p <- ggplot(pca_df, aes(x = PC1, y = PC2 )) +
  geom_point(aes(color = Sample_Pair, shape = Group), size = 3.5) +  # 绘制散点
  geom_text_repel(aes(label = ID), force_pull = 30,size = 4) +     # 避免标签重叠
  # 添加有方向的虚线
  geom_path(aes(group = Sample_Pair, color = Sample_Pair), 
            arrow = arrow(type = "closed", length = unit(0.6, "cm"),angle = 30), 
            linetype = "dotdash", linewidth = 1.4, alpha = 0.4) +  
  #stat_ellipse(aes(color = Group), level = 0.95)+
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_values)  +  # 形状
  labs(title = "PCA Plot", x = paste0("PC1 (", var_percent[1], "%)"),
       y = paste0("PC2 (", var_percent[2], "%)"), color = "Sample Pair", shape = "Group") +
  theme_bw() 
print(p)
ggsave(
  filename = paste0(dir_pca,"Pca_Track.pdf"),    
  plot = p,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 20,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 16,         # 图像高度
  units = "cm")          # 尺寸单位c("in", "cm", "mm", "px")

