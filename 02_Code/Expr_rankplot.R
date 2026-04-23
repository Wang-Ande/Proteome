# Expr rankplot ----
library(readxl)
library(readr)
library(dplyr)
library(openxlsx)
library(ggplot2)

## data input ----
raw_data <- read.csv("01_data/report.pg_matrix.csv", row.names = 1)
raw_data <- raw_data[,-c(1:2,4)]
raw_data <- raw_data[!is.na(raw_data$Genes),]
y <- raw_data$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
raw_data$Genes <- gene
colnames(raw_data)[1] <- "Protein"
# 求平均值排序
raw_data1 <- raw_data[,grep("MV4_11_WT|Protein", colnames(raw_data))]
raw_data1 <- as.data.frame(raw_data1)
raw_data1 <- raw_data1[rowSums(is.na(raw_data1[,2:ncol(raw_data1)])) != 3, ]
raw_data1[,2:ncol(raw_data1)] <- lapply(raw_data1[,2:ncol(raw_data1)], as.numeric)
raw_data1$aver_exp <- rowMeans(raw_data1[, 2:ncol(raw_data1)], na.rm = TRUE)
raw_data1 <- as.data.frame(raw_data1)
raw_data1 <- raw_data1[order(raw_data1$aver_exp),] # 从小到大排列
raw_data1$rank <- seq(1,dim(raw_data1)[1])
View(raw_data1)
write.csv(raw_data1, file = "./01_Data/MV4_11_WT_Rank.csv", row.names = FALSE)
## plot ---
## select ----
library(ggplot2)
library(ggrepel)
# 标记genelist
genelist <- read.xlsx("./01_data/AML surface protein.xlsx", sheet = 1)
genelist <- genelist$official_gene
# genelist <- "SLC27A4"
# 添加颜色分组列
raw_data1$color_group <- ifelse(raw_data1$official_gene %in% genelist, "highlight", "normal")

# 计算后20%的表达量阈值（即第20百分位数）
threshold_value <- quantile(raw_data1$aver_exp, probs = 0.2, na.rm = TRUE)
log2_threshold <- log2(threshold_value)

# 筛选需要添加标签的点（高亮且表达量高于后20%）
label_data <- subset(raw_data1, color_group == "highlight")

#raw_data$rank <- 
ggplot(raw_data1, aes(x = rank, y = log2(aver_exp))) +
  # 底层灰点
  geom_point(data = subset(raw_data1, color_group == "normal"),
             color = "#E2E2E2", fill = "#E2E2E2", alpha = 0.5,
             shape = 21,size = 3,   # 0~1，越小越透明
             stroke = 0.2) +
  # 高亮红点，稍大一点 + 黑边
  geom_point(data = subset(raw_data1, color_group == "highlight"),
             color = "#99494B", fill = "#99494B",
             size = 3, shape = 21, stroke = 0.4) +
  # 添加后20%指示线
  geom_hline(yintercept = log2_threshold, 
             color = "#99494B", linetype = "dashed", linewidth = 0.8) +
  # 可选：添加标签说明
  annotate("text", x = max(raw_data1$rank), y = log2_threshold, 
           label = "Bottom 20%", color = "#99494B", hjust = 1, vjust = -0.5, size = 4) +
  # 添加标签（仅高亮且表达量高于后20%的点）
  geom_text_repel(
    data = label_data[order(label_data$aver_exp, decreasing = TRUE)[1:20], ],  # 只标前20个,
    aes(label = official_gene),
    size = 5,
    color = "#99494B",
    max.overlaps = 50,  # 允许更多标签重叠
    box.padding = 1,   # 调整标签间距
    segment.color = "#E2E2E2"  # 标签连接线颜色
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Rank",
       y = "Aver log2 abundance")

ggsave(plot = last_plot(),
       height = 6,width = 6,
       filename = "03_result/01_abundance rank/sheet1/Average Expression Rank.pdf")
dev.off()

## all ----
file <- "D:/R/R-Project/Database/MassIVE/PXD032180/01_data/AML surface protein.xlsx"
sheets <- excel_sheets(file)
all_data <- lapply(sheets, function(x){
  read_excel(file, sheet = x) %>%
    dplyr::select(1) %>%
    mutate(sheet = x)
})
all_data <- bind_rows(all_data)
library(dplyr)
table(all_data$sheet)
all_data <- all_data %>%
  mutate(group = case_when(
    sheet == "免疫细胞表面" ~ "CD",
    sheet == "整合素家族" ~ "Integrin",
    sheet == "免疫球蛋白超家族" ~ "IgSF",
    sheet == "选择素" ~ "Selectin",
    sheet == "Toll样受体" ~ "TLR",
    sheet == "酶类表面蛋白" ~ "Enzyme",
    sheet == "粘附分子" ~ "CAM",
    sheet == "其他" ~ "Others",
    TRUE ~ "Other"
  ))
genelist <- all_data

# 对齐基因名
library(limma)
library(org.Hs.eg.db) # 假设是人类数据

# 1. 处理蛋白组基因名
# alias2SymbolTable 会查找别名并返回官方 Symbol
genelist$official_gene <- alias2SymbolTable(genelist$gene, species = "Hs")
# 2. 处理转录组基因名
raw_data1$official_gene <- alias2SymbolTable(raw_data1$Protein, species = "Hs")
# 3. 再次取交集（记得处理转换后可能产生的 NA）
common_genes1 <- intersect(genelist$gene, raw_data1$Protein)
common_genes2 <- intersect(genelist$official_gene, raw_data1$official_gene)
#diff_genes <- setdiff(official_gene1,official_gene2)
# 没有官方名的NA转为原用名
genelist$official_gene <- ifelse(is.na(genelist$official_gene), genelist$gene, genelist$official_gene)
raw_data1$official_gene <- ifelse(is.na(raw_data1$official_gene), raw_data1$Protein, raw_data1$official_gene)

# 添加颜色分组列
raw_data1$color_group <- genelist$group[
  match(raw_data1$official_gene, genelist$official_gene)
]
raw_data1$color_group[is.na(raw_data1$color_group)] <- "normal"

# 计算后20%的表达量阈值（即第20百分位数）
threshold_value <- quantile(raw_data1$aver_exp, probs = 0.2, na.rm = TRUE)
log2_threshold <- log2(threshold_value)

# 筛选需要添加标签的点（高亮且表达量高于后20%）
label_data <- subset(raw_data1, color_group %in%genelist$group)
table(label_data$color_group)
highlight_colors <- c(
  "CAM"       = "#E64B35",
  "CD"        = "#4DBBD5",
  "Enzyme"    = "#00A087",
  "IgSF"      = "#3C5488",
  "TLR"       = "#91D1C2",
  "Integrin"  = "#F39B7F",
  "Selectin"  = "#8491B4",
  "Others"    = "#DAA520"
)

ggplot(raw_data1, aes(x = rank, y = log2(aver_exp))) +
  # 底层灰点，固定颜色，不需要图例
  geom_point(data = subset(raw_data1, color_group == "normal"),
             color = "#E2E2E2", fill = "#E2E2E2", alpha = 0.5,
             shape = 21, size = 3, stroke = 0.2) +
  # 高亮点，颜色映射到 color_group → 会生成图例
  geom_point(data = subset(raw_data1, color_group != "normal"),
             aes(color = color_group, fill = color_group),
             shape = 21, size = 3, stroke = 0.4) +
  # 添加后20%指示线
  geom_hline(yintercept = log2_threshold, 
             color = "#99494B", linetype = "dashed", linewidth = 0.8) +
  # 可选：添加标签说明
  annotate("text", x = max(raw_data1$rank), y = log2_threshold, 
           label = "Bottom 20%", color = "#99494B", hjust = 1, vjust = -0.5, size = 4) +
  # 添加标签（仅高亮且表达量高于后20%的点）
  geom_text_repel(
    data = label_data[order(label_data$aver_exp, decreasing = TRUE)[1:32], ],  # 只标前20个,
    aes(label = official_gene, colour = color_group),
    size = 5,
    #color = "#99494B",
    max.overlaps = 50,  # 允许更多标签重叠
    box.padding = 1,  # 调整标签间距
    segment.size =  0.25,
    segment.alpha = 0.5,
    segment.color = "#E2E2E2",  # 标签连接线颜色
    show.legend = FALSE   # 不显示字体图例，即去掉字母图例图形中的字母a
  ) +
  scale_color_manual(values = highlight_colors, guide = "none") +
  scale_fill_manual(values = highlight_colors, name = "Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Rank",
       y = "Aver log2 abundance")

ggsave(plot = last_plot(),
       height = 6,width = 6,
       filename = "03_result/Rank/OCI_AML2/Average Expression Rank color group.pdf")
dev.off()

# 定义颜色向量
highlight_colors <- c(
  "CAM"       = "#E64B35",
  "CD"        = "#4DBBD5",
  "Enzyme"    = "#00A087",
  "IgSF"      = "#3C5488",
  "Integrin"  = "#F39B7F",
  "Selectin"  = "#8491B4",
  "Others"    = "#91D1C2",
  "Optional"  = "#DAA520"  # 如果有第8组
)
highlight_colors <- c(
  "CAM"       = "#E64B35FF",
  "CD"        = "#4DBBD5FF",
  "Enzyme"    = "#00A087FF",
  "IgSF"      = "#3C5488FF",
  "Integrin"  = "#F39B7FFF",
  "Selectin"  = "#8491B4FF",
  "Others"    = "#91D1C2FF"
)

# abundance rank with logfc ----
library(openxlsx)
library(ggplot2)
library(ggrepel)

## data input ----
raw_data <- read.xlsx("./01_data/Supplemental Table 3.xlsx", sheet = 3)
raw_data <- raw_data[, -c(49:ncol(raw_data))]
raw_data <- raw_data[!is.na(raw_data$Protein), ]
raw_data[, 2:ncol(raw_data)][raw_data[, 2:ncol(raw_data)] == 0] <- NA

# 求平均值排序
raw_data <- as.data.frame(raw_data)
raw_data <- raw_data[rowSums(is.na(raw_data[, 2:45])) != 44, ]
raw_data <- raw_data[rowSums(is.na(raw_data[, 46:48])) != 3, ]
raw_data[, 2:ncol(raw_data)] <- lapply(raw_data[, 2:ncol(raw_data)], as.numeric)
raw_data$aver_exp_patient <- rowMeans(raw_data[, 2:45], na.rm = TRUE)
raw_data$aver_exp_healthy <- rowMeans(raw_data[, 46:48], na.rm = TRUE)
raw_data$foldchange <- raw_data$aver_exp_patient/raw_data$aver_exp_healthy
raw_data <- raw_data[order(raw_data$aver_exp_patient), ]
raw_data$rank <- seq_len(nrow(raw_data))

## select label genes ----
file <- "./01_data/AML surface protein.xlsx"
sheets <- excel_sheets(file)
all_data <- lapply(sheets, function(x){
  read_excel(file, sheet = x) %>%
    dplyr::select(1) %>%
    mutate(sheet = x)
})
genelist <- bind_rows(all_data)
genelist <- genelist$official_gene

# 这里假设 raw_data 中已经有 logFC 列
# 例如 raw_data$logFC 代表 AML vs HD CD34p 的 log2 fold change
# 如果你的列名不是 logFC，请替换成实际列名
raw_data$foldchange <- as.numeric(raw_data$foldchange)
raw_data$logfc <- log2(raw_data$foldchange)
# 保留标签分组，只用于挑选要标注的蛋白，不再用于点颜色
raw_data$label_group <- ifelse(raw_data$official_gene %in% genelist, "highlight", "normal")

# 后20%丰度阈值
threshold_value <- quantile(raw_data$aver_exp_patient, probs = 0.2, na.rm = TRUE)
log2_threshold <- log2(threshold_value)

# 只给感兴趣基因加标签
label_data <- subset(raw_data, label_group == "highlight")

# 可选：限制颜色范围，避免极端值把中间颜色压缩
fc_limit_up <- max(subset(raw_data, label_group == "highlight")$logfc, na.rm = TRUE)
fc_limit_down <- min(subset(raw_data, label_group == "highlight")$logfc, na.rm = TRUE)

p <- ggplot(raw_data, aes(x = rank, y = log2(aver_exp_patient))) +
  
  # ① 普通点：固定灰色（不参与映射）
  geom_point(
    data = subset(raw_data, label_group == "normal"),
    color = "#E2E2E2",
    size = 3,
    alpha = 0.6
  ) +
  
  # ② 高亮点：用 logFC 上色
  geom_point(
    data = subset(raw_data, label_group == "highlight"),
    aes(color = logfc),
    size = 3.2
  ) +
  
  geom_hline(
    yintercept = log2_threshold,
    color = "#e06663",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = max(raw_data$rank),
    y = log2_threshold,
    label = "Bottom 20%",
    color = "#e06663",
    hjust = 1,
    vjust = 1.2,
    size = 4
  ) +
  
  geom_text_repel(
    data = label_data[order(label_data$aver_exp_patient, decreasing = TRUE)[1:20], ],
    aes(label = official_gene),
    size = 5,
    color = "#e06663",
    max.overlaps = 50,
    box.padding = 1,
    segment.color = "#D0D0D0",
    show.legend = FALSE
  ) +
  
  scale_color_gradient2(
    low = "#327eba",
    mid = "#D9D9D9",
    high = "#e06663",
    midpoint = 0,
    limits = c(fc_limit_down, fc_limit_up),
    name = "log2FC\n(vs HD CD34p)"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  ) +
  
  labs(
    x = "Rank",
    y = "Aver log2 abundance"
  )
print(p)

ggsave(
  plot = p,
  height = 6,
  width = 6.5,
  filename = "03_result/01_abundance rank/all/Average_Expression_Rank_logFC.pdf"
)
dev.off()

# abundance rank with logfc NA filled ----
library(openxlsx)
library(ggplot2)
library(ggrepel)

## data input ----
raw_data <- read.xlsx("./01_data/Supplemental Table 3.xlsx", sheet = 3)
raw_data <- raw_data[, -c(49:ncol(raw_data))]
raw_data <- raw_data[!is.na(raw_data$Protein), ]
raw_data[, 2:ncol(raw_data)][raw_data[, 2:ncol(raw_data)] == 0] <- NA

# 求平均值排序
raw_data <- as.data.frame(raw_data)
raw_data <- raw_data[rowSums(is.na(raw_data[, 2:45])) != 44, ]
raw_data <- raw_data[rowSums(is.na(raw_data[, 46:48])) != 3, ]
raw_data[, 2:ncol(raw_data)] <- lapply(raw_data[, 2:ncol(raw_data)], as.numeric)
raw_data$aver_exp_patient <- rowMeans(raw_data[, 2:45], na.rm = TRUE)
raw_data$aver_exp_healthy <- rowMeans(raw_data[, 46:48], na.rm = TRUE)
raw_data$foldchange <- raw_data$aver_exp_patient/raw_data$aver_exp_healthy
raw_data <- raw_data[order(raw_data$aver_exp_patient), ]
raw_data$rank <- seq_len(nrow(raw_data))

## select label genes ----
file <- "./01_data/AML surface protein.xlsx"
sheets <- excel_sheets(file)
all_data <- lapply(sheets, function(x){
  read_excel(file, sheet = x) %>%
    dplyr::select(1) %>%
    mutate(sheet = x)
})
genelist <- bind_rows(all_data)
genelist <- genelist$official_gene

# 这里假设 raw_data 中已经有 logFC 列
# 例如 raw_data$logFC 代表 AML vs HD CD34p 的 log2 fold change
# 如果你的列名不是 logFC，请替换成实际列名
raw_data$foldchange <- as.numeric(raw_data$foldchange)
raw_data$logfc <- log2(raw_data$foldchange)
# 保留标签分组，只用于挑选要标注的蛋白，不再用于点颜色
raw_data$label_group <- ifelse(raw_data$official_gene %in% genelist, "highlight", "normal")

# 后20%丰度阈值
threshold_value <- quantile(raw_data$aver_exp_patient, probs = 0.2, na.rm = TRUE)
log2_threshold <- log2(threshold_value)

# 只给感兴趣基因加标签
label_data <- subset(raw_data, label_group == "highlight")

# 可选：限制颜色范围，避免极端值把中间颜色压缩
fc_limit_up <- max(subset(raw_data, label_group == "highlight")$logfc, na.rm = TRUE)
fc_limit_down <- min(subset(raw_data, label_group == "highlight")$logfc, na.rm = TRUE)

p <- ggplot(raw_data, aes(x = rank, y = log2(aver_exp_patient))) +
  
  # ① 普通点：固定灰色（不参与映射）
  geom_point(
    data = subset(raw_data, label_group == "normal"),
    color = "#E2E2E2",
    size = 3,
    alpha = 0.6
  ) +
  
  # ② 高亮点：用 logFC 上色
  geom_point(
    data = subset(raw_data, label_group == "highlight"),
    aes(color = logfc),
    size = 3.2
  ) +
  
  geom_hline(
    yintercept = log2_threshold,
    color = "#e06663",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = max(raw_data$rank),
    y = log2_threshold,
    label = "Bottom 20%",
    color = "#e06663",
    hjust = 1,
    vjust = 1.2,
    size = 4
  ) +
  
  geom_text_repel(
    data = label_data[order(label_data$aver_exp_patient, decreasing = TRUE)[1:20], ],
    aes(label = official_gene),
    size = 5,
    color = "#e06663",
    max.overlaps = 50,
    box.padding = 1,
    segment.color = "#D0D0D0",
    show.legend = FALSE
  ) +
  
  scale_color_gradient2(
    low = "#327eba",
    mid = "#D9D9D9",
    high = "#e06663",
    midpoint = 0,
    limits = c(fc_limit_down, fc_limit_up),
    name = "log2FC\n(vs HD CD34p)"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  ) +
  
  labs(
    x = "Rank",
    y = "Aver log2 abundance"
  )
print(p)

ggsave(
  plot = p,
  height = 6,
  width = 6.5,
  filename = "03_result/01_abundance rank/all/Average_Expression_Rank_logFC.pdf"
)
dev.off()
