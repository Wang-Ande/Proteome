# PCA_Track
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(dplyr)
library(openxlsx)
library(readr)

dir_pca <- "./03_Result/QC/OCI_AML2/"
# Group Input ----
sample_info <- read.xlsx("./01_Data/IC50_group.xlsx")
rownames(sample_info) <- sample_info$ID
target_sample <- sample_info[grep("OCI",rownames(sample_info)),]

# Data Input ----
# 表达矩阵 data_matrix (行为基因，列为样本)
data_matrix <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
data_matrix <- log2(data_matrix)

# 提取目标分组样本
target_matrix <- data_matrix[,rownames(target_sample)]

# Caculate PCA ----
pca_res <- prcomp(t(target_matrix), scale. = TRUE)

# 计算主成分贡献值（方差百分比）
var_percent <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x)
pca_df$ID <- rownames(pca_df)

#  结合样本信息
pca_df <- merge(pca_df, target_sample, by.y = "ID")

#  提取样本分组信息（颜色、形状）
pca_df$Sample_Pair <- gsub("_(WT|2W|6W)", "", pca_df$ID)  # 只保留 MOLM13_x 作为配对标识

# 按样本配对映射颜色 -----------------------------------------------------------
color_palette <- c("MOLM13_1" = "#66C1A4", 
                   "MOLM13_2" = "#8C9FCA", 
                   "MOLM13_3" = "#FB8C63")
shape_values <- c("High" = 15, 
                  "Low" = 17, 
                  "Ctrl" = 16)

#  绘制 PCA 图
# 确保数据按照 Sample_Pair 和 Group (Ctrl → Low → High) 排序
pca_df$Group <- factor(pca_df$Group, levels = c("Ctrl", "Low", "High"))
pca_df <- pca_df %>%
  arrange(Sample_Pair,Group)
 
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2 )) +
  geom_point(aes(color = Sample_Pair, shape = Group), size = 3.5) +  # 绘制散点
  geom_text_repel(aes(label = ID), force_pull = 30,size = 4) +     # 避免标签重叠
  # 添加有方向的虚线
  geom_path(aes(group = Sample_Pair, color = Sample_Pair), 
            arrow = arrow(type = "closed", length = unit(0.6, "cm"),angle = 30), 
            linetype = "dotdash", linewidth = 1.4, alpha = 0.4,lineend = "square",  # 优化连线末端形状
            linejoin = "round") +  #优化连线拐角形状
  #stat_ellipse(aes(color = Group), level = 0.95)+
  scale_color_manual(values = color_palette,guide = guide_legend(ncol = 1, order = 2)) +
  scale_shape_manual(values = shape_values,guide = guide_legend(ncol = 1, order = 1))  +  # 形状
  labs(title = "PCA Plot", x = paste0("PC1 (", var_percent[1], "%)"),
       y = paste0("PC2 (", var_percent[2], "%)"), color = "Sample Pair", shape = "Group") +
  theme_bw()  +
  theme(legend.box = "vertical") # 图例垂直排列
print(p1)
ggsave(
  filename = paste0(dir_pca,"Pca_sample_Track.pdf"),    
  plot = p1,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 20,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 16,         # 图像高度
  units = "cm")          # 尺寸单位c("in", "cm", "mm", "px")

# 按样本分段映射颜色 -----------------------------------------------------------
pca_df$Group <- factor(pca_df$Group, levels = c("Ctrl", "Low", "High"))
shape_values <- c(Ctrl = 16, Low = 17, High = 15)  # 示例形状映射

# 生成每个阶段的起点和终点坐标
pca_phase <- pca_df %>%
  group_by(Sample_Pair) %>%
  arrange(Group) %>%
  mutate(
    next_Group = lead(Group),
    Phase = case_when(
      Group == "Ctrl" & next_Group == "Low" ~ "Ctrl→Low",  # 直接使用箭头符号
      Group == "Ctrl" & next_Group == "High" ~ "Ctrl→High",
      Group == "Low" & next_Group == "High" ~ "Low→High",
      TRUE ~ NA_character_
    ),
    next_PC1 = lead(PC1),
    next_PC2 = lead(PC2)
  ) %>%
  filter(!is.na(Phase)) %>%
  ungroup()

# 定义阶段颜色映射 
phase_colors <- c(
  "Ctrl→Low" = "#1B9E77",    
  "Ctrl→High" = "#7570B3",
  "Low→High" = "#D95F02"
)


p2 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # 绘制阶段线段（颜色按 Phase 映射）
  geom_segment(
    data = pca_phase,
    aes(
      x = PC1, y = PC2,
      xend = next_PC1, yend = next_PC2,
      color = Phase
    ),
    linewidth = 0.9,
    alpha = 0.8,              # 设置透明度（0=全透明，1=不透明）
    linetype = "solid",       # 线段类型（可选 "solid", "dashed", "dotted", "dotdash" 等）
    arrow = arrow(
      type = "closed",          # 箭头类型（"open" 或 "closed"）
      length = unit(0.5, "cm"),
      angle = 30,             # 箭头头部角度（30~90，越小头部越尖）
      ends = "last"           # 箭头位置（"last" 线段末端，"first" 起点，"both" 两端）
    )
  ) +
  # 绘制样本点（保留原有颜色和形状）
  geom_point(aes(color = Sample_Pair, shape = Group), size = 3.5) +
  # 标签防重叠
  geom_text_repel(
    aes(label = ID),
    size = 4.5,                  # 字体大小
    color = "#333333",           # 字体颜色（深灰色）
    fontface = "plain",           # "plain", "bold", "italic"
    family = "sans",             # 字体类型（如Arial）
    box.padding = 0.3,           # 标签与点的间距
    point.padding = 0.5,         # 标签之间的最小间距
    min.segment.length = 0.2,    # 短距离标签显示线段的最小长度
    segment.color = "#666666",   # 连接线颜色（中灰色）
    segment.alpha = 0.6,         # 连接线透明度
    segment.linetype = "dotted", # 连接线类型
    force = 0.5,                 # 避让算法的力度（值越大，标签越分散）
    max.overlaps = 20,           # 允许的最大重叠次数
    max.time = 1,                # 计算避让的最大时间（秒）
    max.iter = 1e4,              # 迭代次数上限
    direction = "both"           # 避让方向（"x"/"y"/"both"）
  )  +
  # 颜色和形状比例尺
  scale_color_manual(
    name = "Phase",
    values = c(phase_colors), 
    guide = guide_legend(ncol = 1, order = 2),
  ) +
  scale_shape_manual(values = shape_values) +
  labs(
    x = paste0("PC1 (", var_percent[1], "%)"),
    y = paste0("PC2 (", var_percent[2], "%)")
  ) +
  guides(
    shape = guide_legend(ncol = 1, order = 1)
  ) + 
  theme_bw() +
  theme_minimal(base_size = 12) +  # 基础字体大小
  theme(
      axis.title = element_text(size = 12, face = "plain"),  # 轴标题加粗，字号14
    # 图例标题 
    legend.title = element_text(
      family = "serif",        # 字体类型（如serif衬线体）
      size = 12,               # 字体大小
      face = "bold",           # 加粗
      color = "#333333",       # 深灰色
      margin = margin(b = 5)   # 标题与标签的间距
    ),
    
    # 图例标签 
    legend.text = element_text(
      family = "sans",         # 无衬线字体（如Arial）
      size = 11,               
      color = "#666666",       # 中灰色
      margin = margin(r = 10)  # 标签右侧间距
    ),
    
    # 整体布局 
    legend.position = "right",
    legend.spacing.y = unit(0.3, "cm"),  # 图例项垂直间距
    legend.key.size = unit(0.8, "cm"),   # 图例符号大小
    legend.background = element_rect(fill = "white", color = NA)  # 背景优化
  )

# 输出图形
print(p2)
ggsave(
  filename = paste0(dir_pca,"Pca_phase_Track.pdf"),    
  plot = p2,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 7.91,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 6.58,         # 图像高度
  units = "in")          # 尺寸单位c("in", "cm", "mm", "px")

