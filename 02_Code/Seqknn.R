# install.packages("multiUS")
library(multiUS)

# 1. Data input ----
data_input <- read.csv("./01_Data/report.pg_matrix_fill_before.csv",row.names = 1)
data_input <- data_input[,grep("MV4",colnames(data_input))]
data_input <- data_input[,-grep("MV4_11_6W_1",colnames(data_input))]
# 是否进行log2运算？
data_input <- log2(data_input)

# 是否进行median normalization运算？
# 计算每列的中位数
median_values <- apply(data_input, 2, median, na.rm = TRUE)

# 进行中位数标准化
data_input <- data_input %>%                 
  mutate(across(everything() , ~ ./median(., na.rm = TRUE)))

# 计算NA值比例
NA_ratio <- colSums(is.na(data_input))/dim(data_input)[1]
NA_ratio <- as.data.frame(NA_ratio)
NA_ratio$Samples <- rownames(NA_ratio)

# 设定NA值比例分类
NA_ratio[NA_ratio$NA_ratio < 0.2,"group"] <- "Good"
NA_ratio[NA_ratio$NA_ratio < 0.5&NA_ratio$NA_ratio > 0.2,"group"] <- "OK"
NA_ratio[NA_ratio$NA_ratio < 0.8&NA_ratio$NA_ratio > 0.5,"group"] <- "Bad"
NA_ratio[NA_ratio$NA_ratio > 0.8,"group"] <- "Remove"
NA_ratio$group <- factor(NA_ratio$group,levels = c("Good","OK","Bad","Remove"))

# 绘制NA值比例分布图
library(ggplot2)
p <- ggplot(NA_ratio,aes(x = NA_ratio, y = reorder(Samples,NA_ratio,decreasing = T), fill = group)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = sprintf("%.2f", NA_ratio)), hjust = -0.5) + 
  scale_fill_manual(values = c("Good" = "#62c882",
                               "OK" = "#b5e281",
                               "Bad" = "#febe80",
                               "Remove" = "#bb2022")) + 
  theme_bw() + 
  labs( x = "Number of missing rows",
        y = "Samples")
ggsave(filename = "NA ratio.pdf",device = pdf,plot = p,
       path = "./03_Result/QC/MV4_11_single_fill/",dpi = 300)

# 计算每个蛋白质的缺失比例
cutoff_NA_ratio <- 0.5 # defult

NA_ratio_protein <- as.data.frame(rowSums(is.na(data_input))/dim(data_input)[2])
colnames(NA_ratio_protein) <- "NA_ratio_protein"
NA_ratio_protein$Protein.Group <- rownames(NA_ratio_protein)
NA_ratio_protein <- NA_ratio_protein[NA_ratio_protein$NA_ratio_protein < cutoff_NA_ratio,]

# 2. multiUS::seqKNNimp ----
data <- data_input[rownames(data_input)%in%rownames(NA_ratio_protein),]
data_fill <- multiUS::seqKNNimp(data = data,k = 8)
min(data_fill)

write.csv(data_fill,file = "./01_Data/MV4_11_report.pg_matrix_fill.csv")

# 3. Normalization -----------------------------------------------------------
# 
# 将填充后的数据导入
data_fill <- read_csv("./01_Data/MV4_11_report.pg_matrix_fill.csv")
data_fill <- as.data.frame(data_fill)
rownames(data_fill) <- data_fill$...1
data_fill <- subset(data_fill,select = -c(`...1`))

# 如果进行了median normalization运算，运行下面的函数
data_fill <- sweep(data_fill, 2, median_values, `*`)

# 如果进行log2处理，运行下面的函数
data_fill <- 2 ^ data_fill

## 3.1 Intensity normalization ----
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
# 4. Res output ----
write.csv(data_fill_normalization,file = "./01_Data/MV4_11_report.pg_matrix_fill_norma.csv")
