# install.packages("multiUS")
library(multiUS)

# 1. Data input ----
data_input <- read.csv("./01_Data/report.pg_matrix_fill_before.csv",row.names = 1)
data_input <- data_input[,grep("OCI",colnames(data_input))]
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
       path = "./03_Result/QC/OCI_AML2/",dpi = 300)

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

# 3. Res output ----
write.csv(data_fill,file = "./01_Data/OCI_report.pg_matrix_fill.csv")
