# kegg通路图 ----
library(openxlsx)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)

enrich_res <- read.csv("./03_Result/GO&KEGG/MOLM13/High_vs_Con/KEGG_down.csv")
pathway_targeted <- enrich_res[grep("reactive",enrich_res$Description),]
gene_ids <- pathway_targeted$geneID
gene_list <- unlist(strsplit(gene_ids, "/"))
gene <- clusterProfiler::bitr(gene_list, fromType = 'ENTREZID', 
                              toType = 'SYMBOL', OrgDb = org.Hs.eg.db)

write.xlsx(gene, file = "./03_Result/GO&KEGG/MOLM13/High_vs_Con/hsa05208_ROS_genes_up.xlsx")

# logFC
DE_gene <- read.csv("./03_Result/DEP/MOLM13/High_vs_Ctrl/result_DE.csv",row.names = 1)

y <- DE_gene$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
DE_gene$gene <- gene1

gene_logfc <- merge(gene, DE_gene, by.x = "SYMBOL", by.y = "gene", all.x = TRUE)

# 输入：一个命名向量，名称是Entrez ID，数值是logFC或任意表达量
enriched_genes  <- gene_logfc$ENTREZID
gene_logfc$logFC <- ifelse(gene_logfc$logFC < 0, -1, 1)
gene_data <- gene_logfc$logFC
names(gene_data) <- enriched_genes

# 绘制并高亮基因
pathview(gene.data = gene_data, 
         pathway.id = "05208",       # 通路ID
         species = "hsa",            # 人类
         gene.idtype = "entrez",     # ID类型
         limit = list(gene=1, cpd=1),
         bins = list(gene=10, cpd=10),
         low = list(gene="#93A5CB", cpd="#66F1A9"),      # 深蓝（下调）
         mid = list(gene="#F7F7F7", cpd="#F0F0F0"),      # 白灰色（中性）
         high = list(gene="#C6133B", cpd="#FC4E2A")     # 深红（上调）
         )    

# 箱线图+点图 ----
library(ggplot2)
library(ggpubr)
library(dplyr)
library(openxlsx)

# set output category
dir <- "./03_Result/DEP/MV4_11/Expr_boxplot/"

# load data
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
expr <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)
group <- read.xlsx("./01_Data/IC50_group.xlsx")

# 定义变量
gene <- "ELOVL5"           # 目标基因, 必须定义！！！！！！！！！！！！！！！！！！
sample <- "MV4_11"       # 箱线图标题
mycol <- c("#6388B4","#FFAE34","#EF6F6A")  

# process data
targeted_gene <- expr[grep(gene,expr$Genes),]
expr_targeted <- targeted_gene[1,grep(sample,colnames(targeted_gene))]
expr_targeted <- t(log2(expr_targeted))
colnames(expr_targeted)[1] <- "expr"
expr_targeted <- as.data.frame(expr_targeted)
group_targeted <- group[grep(sample,group$id),c(1,2)]
expr_targeted$group <- group_targeted$group


if(T){
# 计算样本数（假设 group1 是分组列）
dt <- expr_targeted
num_normal <- sum(dt$group == "Ctrl")
num_low_resistance <- sum(dt$group == "Low")
num_high_resistance <- sum(dt$group == "High")

# 确保分组列是因子，并指定顺序
dt$group <- factor(dt$group, levels = c("Ctrl", "Low", "High"))
rownames(dt)
# dt$pair <- c(1,1,1,2,3,2,3,2,3)
# table(dt$pair)

# 绘图
p1 <- ggplot(data = dt,
             aes(x = group, 
                 y = expr,  # 假设表达量列名为 expression
                 color = group, 
                 fill = group)) +
  #geom_line(aes(group = pair), color = "gray50", linewidth = 0.5, alpha = 0.6) +  # 加入配对线
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +  # 隐藏箱线图异常值
  geom_jitter(width = 0.2, size = 2.8, shape = 21, 
              stroke = 0.8,aes(fill = group), 
              alpha = 0.7) +
  theme_classic() +  scale_color_manual(name = "Group", values = mycol) +
  scale_fill_manual(name = "Group" ,values = mycol) +
  labs(
    title = paste(gene, "expression in ",sample ),
    x = "Resistance state",
    y = "Log2(expr + 1)"
  ) + 
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        aspect.ratio = 1.2,                      # 设置纵横比，调整为更高
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
        legend.background = element_rect(fill = "white", color = NA))  # 背景优化  
# 添加统计学检验（可选）
#p1 <- p1 + stat_compare_means(comparisons = list(c("Ctrl", "Low"),
                                                 #c("Low", "High"),
                                                 #c("Ctrl", "High")),
                              #method = "t.test", paired = TRUE,  # 非参数检验 wilcox.test, 参数检验 t.test 
                              #label = "p.signif")  # 显示显著性标记（***, ​**, *）

ggsave(filename = paste0(dir,"/",gene,"_expr_boxplot.pdf"),
       plot = p1, device = cairo_pdf, 
       width = 6, height = 5)
}

# 多蛋白表达热图 ----
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
expr_anno <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)

# 转换基因名 
y <- expr_anno$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
expr_anno$gene <- gene1

# 选择目标蛋白
targeted_prote <- expr_anno[grep("PRDX4|ECH1|ELOVL5|ACAT1", expr_anno$gene),]
rownames(targeted_prote) <- targeted_prote$gene
targeted_prote <- targeted_prote[,grep("MOLM13|MV4_11|OCI", colnames(targeted_prote))]

rownames(targeted_prote)
# 去除同源蛋白（可选）
targeted_prote <- targeted_prote[-grep("TP53I11|TP53BP1|TP53BP2|TP53I3|TP53RK|CDK11B|CDK13|CDK10|CDK19|CDK12|CDK11A|CDK2P1|URB1|ARRB1|ADARB1|GRB1|CCDC25|LILRB1|ZCRB1|RB1CC1|SCARB1",rownames(targeted_prote)),]
targeted_prote <- targeted_prote[,order(colnames(targeted_prote))]
colnames(targeted_prote)
Sample_order <- c("MOLM13_WT_1", "MOLM13_WT_2", "MOLM13_WT_3", 
                  "MOLM13_2W_2", "MOLM13_2W_3", 
                  "MOLM13_6W_1", "MOLM13_6W_2", "MOLM13_6W_3",
                  "MV4_11_WT_1", "MV4_11_WT_2", "MV4_11_WT_3",
                  "MV4_11_2W_1", "MV4_11_2W_2", "MV4_11_2W_3", 
                  "MV4_11_6W_1", "MV4_11_6W_2", "MV4_11_6W_3",
                  "OCI_WT_1", "OCI_WT_2", "OCI_WT_3",
                  "OCI_2W_1", "OCI_2W_2", "OCI_2W_3", 
                  "OCI_4W_2","OCI_6W_1", "OCI_6W_2", "OCI_6W_3")
# Sample_order <- gsub("OCI","MV4_11",Sample_order)
targeted_prote <- targeted_prote[,Sample_order]
colnames(targeted_prote)

library(circlize)
library(ComplexHeatmap)
expr_matrix <- targeted_prote
expr_matrix <- targeted_prote[,-grep("4W",colnames(targeted_prote))]

# 定义 Cell line 分组
cell_lines <- case_when(
  grepl("MOLM13", colnames(expr_matrix)) ~ "MOLM13",
  grepl("MV4_11", colnames(expr_matrix)) ~ "MV4_11",
  grepl("OCI", colnames(expr_matrix)) ~ "OCI_AML2"
)
cell_lines <- factor(cell_lines, levels = c("MOLM13", "MV4_11", "OCI_AML2"))

# 为 Cell line 设置颜色
cell_line_colors <- c(
  MOLM13 = "#66C2A5",
  MV4_11 = "#8DA0CB",
  OCI_AML2 = "#FC8D62"
)

# 创建分组信息
sample_groups <- case_when(
  grepl("WT", colnames(expr_matrix)) ~ "WT",
  grepl("2W|4W", colnames(expr_matrix)) ~ "Low",
  grepl("6W", colnames(expr_matrix)) ~ "High"
)

# 转换为因子，设置顺序
sample_groups <- factor(sample_groups, levels = c("WT", "Low", "High"))

# 创建颜色映射
group_colors <- c(WT = "#6388B4", Low = "#FFAE34", High = "#EF6F6A")

# approach1 同时显示细胞系和分组！！！！！！！！！！！！！！！！！！！！！！！！
ha_col <- HeatmapAnnotation(
  Cell_line = cell_lines,
  Group = sample_groups,
  col = list(
    Cell_line = cell_line_colors, #（Cell line 在上面）
    Group = group_colors
  ),
  annotation_name_side = "right",
  annotation_legend_param = list(
    Cell_line = list(title = "Cell line", title_gp = gpar(fontface = "bold", fontsize = 9), labels_gp = gpar(fontsize = 8)),
    Group = list(title = "Group", title_gp = gpar(fontface = "bold", fontsize = 9), labels_gp = gpar(fontsize = 8))
  ),
  annotation_height = unit(c(5, 5), "mm")  # 控制两个注释行的高度
)

# approach2 只显示分组！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# 创建 group 的列注释对象
ha_col <- HeatmapAnnotation(
  Group = sample_groups,
  col = list(Group = group_colors),
  annotation_name_side = "right",
  annotation_legend_param = list(title = "Group", 
                                 title_gp = gpar(fontface = "bold",fontsize = 10),
                                 labels_gp = gpar(fontsize = 8))
)

# scale
log_expr <- log2(expr_matrix + 1)
scaled_expr <- t(scale(t(log_expr)))  # 每行（每个基因）标准化

# plot
ht <- Heatmap(
  scaled_expr,
  name = "Z-score",
  top_annotation = ha_col,           # 加上分组注释
  col = circlize::colorRamp2(c(-2, 0, 2), colors = c("#1E90FF", "white", "#FF4500")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title = "Z-score", 
                              title_gp = gpar(fontface = "bold",fontsize = 10),
                              legend_height = unit(3, "cm"),   # 控制图例高度
                              grid_width = unit(0.35, "cm"),    # 色块宽度
                              labels_gp = gpar(fontsize = 8)   # 标签字体
  ))
ht <- draw(ht,
           heatmap_legend_side = "right",
           annotation_legend_side = "right",
           merge_legend = TRUE,
           padding = unit(c(2, 12, 2, 2), "mm"),  # 调整边距
           legend_gap = unit(6, "mm")) # 图例边距
# res output
cairo_pdf("./03_Result/DEP/Targeted_proteins_expr_heatmap.pdf", width = 8, height = 2)
ht
dev.off()

# NA fill ----
# 导入分组信息
data_group <- read_excel("./01_Data/IC50_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
# 选择分组样本
data_group <- data_group[grep("OCI",data_group$id),c("id","group")]
table(data_group$group)
# 导入intensity
data_input <- read.csv("./01_Data/report.pg_matrix.csv", row.names = 1)
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Protein.Group
# 保留前四列注释
data_anno <- data_input[,1:4]
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
data_input <- data_input[,-1:-4]
# 选择分组样本
data_input <- data_input[,data_group$id]

# 是否进行log2运算？
data_input <- log2(data_input)

# 是否进行median normalization运算？
# 计算每列的中位数
median_values <- apply(data_input, 2, median, na.rm = TRUE)

# 进行中位数标准化
data_input <- data_input %>%                 
  mutate(across(everything() , ~ ./median(., na.rm = TRUE)))

# 计算样本NA值比例
NA_ratio <- colSums(is.na(data_input))/dim(data_input)[1]
NA_ratio <- as.data.frame(NA_ratio)
NA_ratio$Samples <- rownames(NA_ratio)

# 设定样本NA值比例分类
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
        y = "Samples") +
  # 扩大横坐标范围
  coord_cartesian(xlim = c(0, max(NA_ratio$NA_ratio) * 1.1))
print(p)
ggsave(filename = "NA ratio.pdf", device = pdf,
       plot = p, width = 6, height = 5,
       path = "./03_Result/QC/OCI_AML2_single_fill/",dpi = 300)

# 筛除缺失比例过高的蛋白质
cutoff_NA_ratio <- 0.5 # defult

NA_ratio_protein <- as.data.frame(rowSums(is.na(data_input))/dim(data_input)[2])
colnames(NA_ratio_protein) <- "NA_ratio_protein"
NA_ratio_protein$Protein.Group <- rownames(NA_ratio_protein)
NA_ratio_protein <- NA_ratio_protein[NA_ratio_protein$NA_ratio_protein < cutoff_NA_ratio,]

#  multiUS::seqKNNimp 
data <- data_input[rownames(data_input)%in%rownames(NA_ratio_protein),]
sum(rowSums(is.na(data)) > 0)  # 查看包含NA的蛋白有多少行 
data_fill <- multiUS::seqKNNimp(data = data,k = 10) # 返回填充了多少行
min(data_fill)

#  Res output 
# 如果进行了median normalization运算，运行下面的函数
data_fill <- sweep(data_fill, 2, median_values, `*`)
# 如果进行log2处理，运行下面的函数
data_fill <- 2 ^ data_fill
write.csv(data_fill,file = "./01_Data/OCI_report.pg_matrix_fill.csv")

# 设置输出目录
dir.create("./03_Result/QC/OCI_AML2_single_fill/")
dir <- "./03_Result/QC/OCI_AML2_single_fill/"

# Intensity normalization 
data_fill <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
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

# QC --------------------------------------------------------------------------
# 导入分组信息
data_group <- read_excel("./01_Data/IC50_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
data_group <- data_group[grep("OCI",data_group$id),c("id","group")]
table(data_group$group)
# 设置因子排列顺序
data_group$group <- factor(data_group$group, levels = c("Ctrl", "Low", "High"))
# 配色设置
value_colour <- c("Ctrl" = "#E64B35FF",
                  "Low" = "#4DBBD5FA",
                  "High" = "#F2A200")
rownames(data_group) <- data_group$id

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
           title = "Normalized data")
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
                       name = expression(R^2)) +            # 更改图例标题为 R²
  labs(title = "Pearson correlation between samples") +     # 添加标题
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
        axis.title = element_blank(),                       # 隐藏 X 和 Y 轴标题
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中、加粗
        legend.position = "right")                          # 保持图例在右侧
ggsave("./03_Result/QC/ALL/All_sample_cor.pdf", width = 8, height = 7)

# DE ------------------------------------------------------------------------
# expr input
data_fill_norm <- read.csv("./01_Data/MV4_11_report.pg_matrix_fill_norma.csv", row.names = 1)

# anno input
# 注意，data和data_anno的行名应一致
data_anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
data_anno <- data_anno[rownames(data_anno)%in%rownames(data_fill_norm),]
data_fill_norm <- data_fill_norm[,order(colnames(data_fill_norm))]

## 4.1 Set output catagory ----
dir_DE <- "./03_Result/DEP/MV4_11_single_fill/"

## 4.2 Set group ----
data_group <- read.xlsx("./01_Data/IC50_group.xlsx")
rownames(data_group) <- data_group$id
data_group <- data_group[order(rownames(data_group)),]
table(data_group$group)
targeted_group <- data_group[grep("MV4_11",data_group$id),]

targeted_data <- data_fill_norm[,rownames(targeted_group)]

## 4.3 Res output ----
source("./02_Code/run_DE.R")
group_1 <- "Ctrl"         # group 1为Wild type
group_2 <- "Low"         # group 2为Treatment
result_merge <- run_DE(data = targeted_data,
                       data_group = targeted_group,
                       data_anno = data_anno,
                       group_1 = group_1,
                       group_2 = group_2,
                       log2 = TRUE,
                       logfc_threshold = 0.268,         # 对应fc为1.2倍
                       pvalue_threshold = 0.05,
                       paired = FALSE,
                       pair_col = "pair_id",
                       dir = dir_DE)
# 统计上下调基因数量
table(result_merge$result_merge$Sig)

## 4.4 Annotated volcano plot ----

# 5. GO&KEGG ----
## 5.1 Set output catagory----
dir_enrich <- "./03_Result/GO&KEGG/OCI_AML2_single_fill/Low_vs_Con/"

## 5.2 DE_res input ----
DP_result <- read.csv('./03_Result/DEP/OCI_AML2_single_fill/Low_vs_Ctrl/result_DE.csv')

## 5.3 set P.Value ----
GeneSymbol <- subset(DP_result, P.Value < 0.05)

## 5.4 set cutoff值 ----
cutoff <- 0.263                  # 对应fc约为1.25

# 转换基因名 
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
GeneSymbol$gene <- gene

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

if(T){
  ## 5.5 down genes ----
  down_genes <- subset(GeneSymbol, logFC < -cutoff)
  
  # gene ID转换 
  gene <- clusterProfiler::bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  ### 5.5.1 GO ----
  # GO富集分析
  go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                  OrgDb = GO_database, 
                                  keyType = "ENTREZID", 
                                  ont = "ALL",           # (ALL,BP,CC,MF）
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1,
                                  readable = T)    
  
  ### 5.5.2 KEGG ----
  # KEGG富集分析
  kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                      keyType = "kegg",
                                      organism = KEGG_database,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 1)
  
  ## GO、KEGG结果整合 
  result <- list(enrichGO = go, enrichKEGG = kegg)
  
  # 结果标记为下调 
  result_down <- result
  GO_down <- result_down$enrichGO
  KEGG_down <- result_down$enrichKEGG
  
  ### 5.5.3 down res_output ----
  # 导出下调enrichGO 
  write.csv(GO_down@result, file = paste0(dir_enrich, "/GO_down.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/GO_down.pdf"), width = 6, height = 7.5) # 如果用重叠 width = 6.5, height = 8.5 
  p1 <- dotplot(GO_down, showCategory = 5, split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
    theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # 控制每行最多显示40个字符 
  print(p1)
  dev.off()
  
  # 导出下调enrichKEGG
  write.csv(KEGG_down@result, file = paste0(dir_enrich, "/KEGG_down.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/KEGG_down.pdf"), width = 6, height = 5)
  p2 <- dotplot(KEGG_down,showCategory = 10)
  print(p2)
  dev.off()
}

# 若提示存在'select()' returned 1:many mapping between keys and columns，可能原因是symbol名同名了，
# 比如symbol"HBD",既表示hypophosphatemic bone disease，也可以表示hemoglobin subunit delta，因此为出现1：many
dup_genes <- gene$SYMBOL[duplicated(gene$SYMBOL)]

if(T){
  ## 5.6 up genes ----
  up_genes <- subset(GeneSymbol, logFC > cutoff)
  
  ### 5.6.1 GO-up ----
  # gene ID转换 
  gene <- clusterProfiler::bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  # GO富集分析
  go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                  OrgDb = GO_database, 
                                  keyType = "ENTREZID", 
                                  ont = "ALL", 
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1, 
                                  readable = T)
  ### 5.6.2 KEGG-up ----
  # KEGG富集分析
  kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                      keyType = "kegg",
                                      organism = KEGG_database,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 1)
  
  # GO、KEGG结果整合
  result <- list(enrichGO = go, enrichKEGG = kegg)
  
  # 结果标记为上调
  result_up <- result
  GO_up <- result_up$enrichGO
  KEGG_up <- result_up$enrichKEGG
  
  ### 5.6.3 up res_output ----
  
  # 导出上调enrichGO
  write.csv(GO_up@result, file = paste0(dir_enrich, "/GO_up.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/GO_up.pdf"), width = 6, height = 7.5) # 如果用重叠 width = 6.5, height = 8.5 
  p3 <- dotplot(GO_up, showCategory = 5, split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
    theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35))  # 控制每行最多显示40个字符
  print(p3)
  dev.off()
  
  # 导出上调enrichKEGG
  write.csv(KEGG_up@result, file = paste0(dir_enrich, "/KEGG_up.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/KEGG_up.pdf"), width = 6, height = 5)
  p4 <- dotplot(KEGG_up,showCategory = 10)
  print(p4)
  dev.off()
}

if(T){
  ## 5.7 All DE genes ----
  All_genes <- subset(GeneSymbol, abs(logFC) > cutoff)
  
  ### 5.7.1 GO ----
  # gene ID转换 
  gene <- clusterProfiler::bitr(All_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  # GO富集分析
  go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                  OrgDb = GO_database, 
                                  keyType = "ENTREZID", 
                                  ont = "ALL", 
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1, 
                                  readable = T)
  ### 5.7.2 KEGG ----
  # KEGG富集分析
  kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                      keyType = "kegg",
                                      organism = KEGG_database,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 1)
  
  # GO、KEGG结果整合
  result <- list(enrichGO = go, enrichKEGG = kegg)
  
  # 结果标记为上调
  result_up <- result
  GO <- result_up$enrichGO
  KEGG <- result_up$enrichKEGG
  
  ### 5.7.3 res_output ----
  
  # 导出上调enrichGO
  write.csv(GO@result, file = paste0(dir_enrich, "/GO_all_DE.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/GO_all_DE.pdf"), width = 6, height = 7.5) # 如果用重叠 width = 6.5, height = 8.5 
  p5 <- dotplot(GO, showCategory = 5, split = "ONTOLOGY") + 
    facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
    theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35))  # 控制每行最多显示40个字符
  print(p5)
  dev.off()
  
  # 导出上调enrichKEGG
  write.csv(KEGG@result, file = paste0(dir_enrich, "/KEGG_all_DE.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_enrich, "/KEGG_all_DE.pdf"), width = 6, height = 5)
  p6 <- dotplot(KEGG,showCategory = 10)
  print(p6)
  dev.off()
}

## 5.8 The number of up and down pathways ----
# 下调GO 
table(GO_down@result$p.adjust<0.05)
# 下调KEGG 
table(KEGG_down@result$p.adjust<0.05)
# 上调GO 
table(GO_up@result$p.adjust<0.05)
# 上调KEGG 
table(KEGG_up@result$p.adjust<0.05)

# 6. Selected pathway ----

kegg <- kk_up@result[kk_up@result$pvalue<0.05,]
p2 <- ggplot(kegg,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color= -log10(pvalue)))+
  theme_bw()+labs(y="",x="GeneRatio")+ 
  scale_color_gradient(low="blue",high="red")+
  #scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # 调整Y轴标签角度

ggsave(plot = p2,filename = "./03_Result/GO&KEGG/MOLM13/6W_vs_WT/downkegg.pdf")


# 6.heatmap------
library(reshape2)
library(ComplexHeatmap)
library(openxlsx)
library(gtable)
library(grid)
## df structure
# gene    module
# PCDH9   tan
# ALX1    tan
# CGREF1  darkorange
# 
df <- data.frame(gene = colnames(cleanedData), module = mergedColors)
df$gene <- gsub(".*\\|", "", df$gene)
kegg_list <- list()
i <- 0
all_paths <- c()
for (module in unique(df$module)) {
  i <- i + 1
  genes <- df[df$module == module,]$gene
  trans_id <- bitr(genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  kegg <- enrichKEGG(gene = trans_id$ENTREZID, organism = "hsa", 
                     keyType = "kegg", 
                     qvalueCutoff = 0.05)
  result <- kegg@result
  result$module <- module
  kegg_list[[i]] <- result
  significant_paths <- result[result$pvalue < 0.05, "Description"]
  all_paths <- c(all_paths, significant_paths)
  #all_paths <- c(all_paths, kegg@result$Description[c(1,2)])
  #all_paths <- c(all_paths, kegg@result$Description[1])
}

paths <- c("PI3K-Akt signaling pathway", "MAPK signaling pathway")
all_paths <- unique(c(all_paths, paths))
kegg_all <- do.call(rbind, kegg_list)
kegg_all <- kegg_all[kegg_all$Description %in% all_paths,]
plot_data <- dcast(kegg_all[,c("Description", "pvalue", "module")], Description ~ module, value.var = "pvalue", fill = 1)

write.xlsx(plot_data, file = "03_result/WGCNA/KEGGplotdata.xlsx",sep=",", quote=F)
plot_data <- read_excel("03_result/WGCNA/KEGGplotdataselected.xlsx")
plot_data<- as.data.frame(plot_data)
rownames(plot_data) <- plot_data$Description
plot_data$Description <- NULL
plot_data <- -log(plot_data)
plot_data[plot_data>5] <- 5
#

colnames(plot_data) <- paste0("Cluster", seq(ncol(plot_data)))
plot_data <- as.matrix(plot_data)
pdf(file =  "03_result/WGCNA/10.pathwayheatmap.pdf", width = 10, height = 10)
p<- pheatmap(plot_data,name =" -log(pvalue)",
             cluster_cols = T,
             cluster_rows = F,
             color = colorRampPalette(c("white", "blue"))(100))
