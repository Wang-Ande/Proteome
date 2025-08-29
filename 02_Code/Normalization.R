# 0. Package&Function ----
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(openxlsx)
library(corrplot)
source("./02_Code/QC_PCA.R")
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
source("./02_Code/run_enrichment_analysis.R")

# 1. Data input ----
## 1.1 Group input ----
# 导入分组信息
data_group <- read_excel("./01_Data/IC50_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
data_group <- data_group[grep("OCI",data_group$id),c("id","group")]
table(data_group$group)
# 配色设置
value_colour <- c("Ctrl" = "#E64B35FF",
                  "Low" = "#4DBBD5FA",
                  "High" = "#F2A200")
rownames(data_group) <- data_group$id

## 1.2 DIA matrix input ----
data_input <- read.csv("./01_Data/report.pg_matrix.csv")
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Protein.Group

# 保留前四列注释
data_anno <- data_input[,1:4]
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
data_input <- data_input[,-1:-4]
write.xlsx(data_anno,file = "./01_Data/data_anno.xlsx")

## 1.3 NA_guideR ----
# 注意NA guide官网上也有log2和median，如果此时选择，则在官网上不要选择
# 是否进行log2运算？
data_input <- log2(data_input)

# 是否进行median normalization运算？

# 计算每列的中位数
median_values <- apply(data_input, 2, median, na.rm = TRUE)

# 进行中位数标准化
data_input <- data_input %>%                 
  mutate(across(everything() , ~ ./median(., na.rm = TRUE)))

# 使用此输出进行NA填充，填充网站 https://www.omicsolution.org/wukong/NAguideR/#
write.csv(data_input,file = "./01_Data/report.pg_matrix_fill_before.csv")

# 将填充后的数据导入
data_fill <- read_csv("./01_Data/OCI_report.pg_matrix_fill.csv")
data_fill <- as.data.frame(data_fill)
rownames(data_fill) <- data_fill$...1
data_fill <- subset(data_fill,select = -c(`...1`))

# 如果进行了median normalization运算，运行下面的函数
data_fill <- sweep(data_fill, 2, median_values, `*`)

# 如果进行log2处理，运行下面的函数
data_fill <- 2 ^ data_fill

## 1.4 10minimum fill ----
# 如果已经使用网站填充NA,直接进行normalization ！！！！！！！！！！！！！！！！！！！
intensity_all_protein <- data.frame(Protein.Group = rownames(data_input) ,
                        intensity = rowMeans(data_input,na.rm = T))
intensity = rowMeans(data_input,na.rm = T)
intensity <- na.omit(intensity)
intensity <- sort(intensity)
intensity_10 <- quantile(intensity, 0.1)
intensity_10
log2(intensity_10)
ggplot(intensity, aes(x = log2(intensity))) + 
  geom_density(color = "black", fill = "gray")

# 将第6列到最后一列的所有NA值替换为0
data_input[is.na(data_input)] <- 0

# 将第6列到最后一列的所有值加100000
data_input <- data_input + 100000
data_fill <- data_input

# 2. Normalization -----------------------------------------------------------
# 设置输出目录
#dir.create("./03_Result/")
dir <- "./03_Result/QC/OCI_AML2_single_fill/"

## 2.1 Intensity normalization ----
data_before <- log2(data_fill)

# 计算校正前各样本的intensity median
column_medians <- apply(data_before, 2, median, na.rm = TRUE)
column_medians

# 目标中位数
target_median <- max(column_medians)
data_after <- data_before
for (col in names(data_after)) {
  if (is.numeric(data_after[[col]])) {
    # 防止除以零的情况
    if (column_medians[col] != 0) {
      data_after[[col]] <- data_after[[col]] / column_medians[col] * target_median
    }
  }
}

# 验证校正后各样本intensity median是否一致
column_medians_2 <- apply(data_after, 2, median, na.rm = TRUE)
column_medians_2

# 返回log2之前的数据
data_fill_normalization <- 2 ^ data_after
write.csv(data_fill_normalization,file = "./01_Data/OCI_report.pg_matrix_fill_norma.csv")

# 3. QC --------------------------------------------------------------------------

## 3.1 Boxplot -----------------------------------------------------------------
# 函数内有log2（）
pdf(file = paste0(dir,"QC_boxplot_non-normalized.pdf"),
    width = 6,
    height = 4)

QC_boxplot(2 ^ data_before,data_group = data_group,
           value_colour = value_colour,
           title = "Non-normalized data")
dev.off()

pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_fill_normalization,data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
## 3.2 Heatmap -----------------------------------------------------------------
# 函数内有log2（）
pdf(file = paste0(dir,"QC_heatmap_non-normalized.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = 2 ^ data_before,data_group = data_group,
           value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_fill_normalization,data_group = data_group,
           value_colour = value_colour)
dev.off()

## 3.3 PCA ---------------------------------------------------------------------
# 函数内没有log2（）
pdf(file = paste0(dir,"QC_pca_non-normalized.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = 2 ^ data_before,
       data_group = data_group,
       value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = log2(data_fill_normalization),
       data_group = data_group,
       value_colour = value_colour)
dev.off()

## 3.4 Cor ---------------------------------------------------------------------
#计算样本之间的相关性
dir_cor <- "./03_Result/QC/ALL/"
prote_expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
sorted_colnames <- sort(colnames(prote_expr))       # 对列名进行排序
prote_expr <- prote_expr[, sorted_colnames]         # 按照排序后的列名重新排列

# control 
prote_WT <- prote_expr[,grep("WT",colnames(prote_expr))]

# treatment
prote_trt <- prote_expr[,grep("2w",colnames(prote_expr))]

# Cell line
prote_cell <- prote_expr[,grep("MV4_11",colnames(prote_expr))]

# log2(expr+1)
prote_log2 <-log2(prote_expr)  

# cor analysis
corr_matrix <- cor(prote_log2, method = 'pearson')         # cor函数计算两两样本（列与列）之间的相关系数
View(corr_matrix)	                            

# two visualization methods
# methods 1
if(T){
  cairo_pdf(paste0(dir_cor,'MV4_11_sample_cor.pdf'), width = 11, height = 11)	
  breaks <- seq(0.9, 1, length.out = 100)  # 强制颜色范围放大 0.9~1
  color_palette <- colorRampPalette(c("white", "#2166AC"))(99)
  corrplot(corr_matrix, type = 'upper',   # type='upper'：只显示右上角
           method = "color",      # ("circle", "square", "ellipse", "number", "shade", "color")
           col = color_palette,
           col.lim = c(0.9, 1),
           tl.col = 'black',       # tl.col='black'：字体颜色黑色
           order = 'hclust',       # order='hclust'：使用层次聚类算法
           tl.srt = 45,            # tl.srt = 45：x轴标签倾斜45度
           number.cex = 0.75,       # 相关性系数字体大小
           addCoef.col = 'white')	 # addCoef.col='white'：添加相关系数数值，颜色白色
  dev.off() 
}

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
  scale_fill_gradientn(colors = c("#D1E5F0", "#2166AC"), 
                       limits = c(min(corr_matrix), 1),                  # 设定颜色映射范围
                       name = expression(R)) +            # 更改图例标题为 R²
  labs(title = "Pearson correlation between samples") +     # 添加标题
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
        axis.title = element_blank(),                       # 隐藏 X 和 Y 轴标题
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中、加粗
        legend.position = "right")                          # 保持图例在右侧
ggsave("./03_Result/QC/ALL/All_sample_cor.pdf", width = 8, height = 7)
