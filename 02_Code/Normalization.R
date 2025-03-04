# 0. Package&Function ----
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(openxlsx)
source("./02_Code/QC_PCA.R")
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
source("./02_Code/run_enrichment_analysis.R")

# 1. Data input ----
## 1.1 Group input ----
# 导入分组信息
data_group <- read_excel("./01_Data/All_group_info.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
#data_group <- data_group[-grep("4W",data_group$id),] # 删除4w样本
data_group$group <- gsub("_WT","",data_group$group)
table(data_group$group)
# 配色设置
# 配色设置
value_colour <- c("MOLM13" = "#E64B35FF",
                  "MV4_11" = "#4DBBD5FA",
                  "OCI" = "#F2A200")
rownames(data_group) <- data_group$id

## 1.2 DIA matrix input ----
data_input <- read.csv("./01_Data/report.pg_matrix.csv")
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Protein.Group
data_input <- data_input[,-grep("4W",colnames(data_input))] # 删除4w样本

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
data_fill <- read_csv("./01_Data/report.pg_matrix_fill.csv")
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
dir.create("./03_Result/QC/OCI_AML2")
dir <- "./03_Result/QC/ALL/"

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
write.csv(data_fill_normalization,file = "./01_Data/report.pg_matrix_fill_norma.csv")

# 3. QC --------------------------------------------------------------------------

## 3.1 Boxplot -----------------------------------------------------------------
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
QC_PCA(data = data_fill_normalization,
       data_group = data_group,
       value_colour = value_colour)
dev.off()
