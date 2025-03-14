# ClusterGVis

# 1. Packages ----
# Note: please update your ComplexHeatmap to the latest version!
library(devtools)
devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(grid)
library(ComplexHeatmap)
library(clusterProfiler)
library(TCseq)
library(e1071)
library(tkWidgets)
library(Mfuzz)
library(circlize)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(openxlsx)

# 2. Data input ----
# load data(standardised)
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)

# anno input
# 注意，data和data_anno的行名应一致
data_anno <- read.xlsx("./01_Data/data_anno.xlsx")
data_anno <- as.data.frame(data_anno)
colnames(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
data_anno <- data_anno[rownames(data_anno)%in%rownames(expr),]
y <- data_anno$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
expr$gene <- gene 

# 将rownames由protein.group改为gene
rownames(expr) <- expr$gene    # gene做为行名重复，进行去重处理

# 1. 处理 NA 值
expr <- expr[!is.na(expr$gene), ]  # 删除 NA 行

# 2. 处理重复基因
expr <- aggregate(. ~ gene, data = expr, FUN = max)

# 3. 设定 rownames
rownames(expr) <- expr$gene
expr <- expr[,-1]


# data group
sample_group <- read.xlsx("./01_Data/IC50_group.xlsx")
rownames(sample_group) <- sample_group$id
table(sample_group$group)

# check
rownames(sample_group) == colnames(expr)

## 2.1 Data preprocess ----
# 要求输入归一化之后的数据（已进行中位值归一化）
expr_log2 <- log2(expr)
colnames(expr_log2)

# 修改表达矩阵的列名
colnames(expr_log2) <- sample_group$group
colnames(expr_log2)
#[1] "Medium" "Low"    "High"   "Low"    "High"   "Low"    "Medium" "High"   "High"   "High"  
#[11] "Low"    "Low"    "Medium" "Medium" "Medium" "High"   "High"   "Low"    "Low"    "Medium"
#[21] "Medium" "Medium" "Medium" "High"   "High"   "Low"    "Low"   

library(limma)
# limma::avereps：这是来自 limma 包的函数，avereps 用于对重复数据进行平均值计算。
# avereps 会根据指定的 ID 进行分组，并对相同 ID 的数据取平均值
avereps_df <- t(limma::avereps(t(expr_log2) , ID = colnames(expr_log2))) #对相同时间序列的表达值取平均
avereps_df[1:3,1:3]
avereps_df <- avereps_df[,c("Ctrl","Low","High")]
save(avereps_df,file = './03_Result/C_means_cluster/IC50 group/avereps_df.Rdata')

# 3. Cluster ----
## 3.1 Set output path ----
dir_cl <- "./03_Result/C_means_cluster/IC50 group/"
exps <- as.matrix(avereps_df) 

# check optimal cluster numbers
getClusters(obj = exps)
# choose 6

## 3.2 mfuzz for clustering ----
cm <- clusterData(obj = exps,
                  scaleData = TRUE,
                  seed = 123,
                  cluster.method = "mfuzz",
                  cluster.num = 6)
save(cm,file = paste0(dir_cl,"C-Means_res.Rdata"))

## 3.3 kmeans for clustering ----
km <- clusterData(obj = exps,
                  scaleData = TRUE,
                  cluster.method = "kmeans",
                  cluster.num = 6,
                  seed = 123)


# 4. Plot ----
## 4.1 plot line only ----
visCluster(object = cm,
           plot.type = "line")

# change color
visCluster(object = cm,
           plot.type = "line",
           ms.col = c("green","orange","red"))

# remove meadian line
pdf(paste0(dir_cl,"Cluster_lineplot.pdf"),width = 8,height = 5)
visCluster(object = cm,
           plot.type = "line",
           ms.col = c("green","orange","red"),
           add.mline = FALSE)
dev.off()

## 4.2 plot heatmap only ----
pdf(paste0(dir_cl,'Cluster_heatmap.pdf'),height = 10,width = 6)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45)
dev.off()

# 5. Enrich analysis ----
library(org.Hs.eg.db)

# enrich for clusters
enrich <- enrichCluster(object = cm,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        id.trans = TRUE,
                        fromType = "SYMBOL",
                        toType = c("ENTREZID"),
                        readable = TRUE,
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
save(enrich,file = (paste0(dir_cl,"enrich.Rdata")))
# check
head(enrich,3)

# cluster num
cl_num <- 6   # ggsci::pal_d3()(cl_num) 

# plot
pdf(paste0(dir_cl,'Term_enrichplot.pdf'),height = 10,width = 11,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = enrich,
           line.side = "left",
           go.col = rep(ggsci::pal_d3()(cl_num),each = 5), 
           go.size = "pval")
dev.off()
