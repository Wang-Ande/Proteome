# 1. Packages ----
rm(list = ls())
library(grid)
library(ComplexHeatmap)
library(clusterProfiler)
library(TCseq)
library(e1071)
library(tkWidgets)
library(Mfuzz)
library(circlize)
library(ClusterGVis)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(openxlsx)

# 2. Data input ----
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv")
rownames(expr) <- expr$X
expr <- expr[,-1]

# anno input
# 注意，data和data_anno的行名应一致
data_anno <- read.xlsx("./01_Data/data_anno.xlsx")
data_anno <- as.data.frame(data_anno)
colnames(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
data_anno <- data_anno[rownames(data_anno)%in%rownames(expr),]

## 2.1 Data preprocess ----
# 要求输入归一化之后的数据（已进行中位值归一化）
expr_log2 <- log2(expr)
colnames(expr_log2) <- gsub("WT","0W",colnames(expr_log2))
colnames(expr_log2)
dat <- expr_log2[,grep("MV4_11",colnames(expr_log2))]


# 取平均值(dat同一样本不同生物学重复列名需要相同)
colnames(dat) <- substr(colnames(dat), 1, nchar(colnames(dat)) - 2)     # 去除编号后缀
avereps_df  <- t(limma::avereps(t(dat) , ID = colnames(dat))) 

# 按时间序列排序
colnames(avereps_df)
colnames(avereps_df) <- gsub("MV4_11_","",colnames(avereps_df))
avereps_df = avereps_df[,c( "0W", "2W","6W" )]

## 2.2 ExpressionSet ----
# Mfuzz聚类时要求是一个ExpressionSet类型的对象
eset <- new("ExpressionSet",exprs = as.matrix(avereps_df))

# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset,min.std=0)

## 2.3 Standardized ----
# 聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 2.4 Set m value ----
m <- mestimate(eset)

## 2.5 Set c value ----
# cselection 函数，repeats 为重复聚类，visu=TRUE 意味着输出可视化
tmp  <- cselection(eset, m=m, crange=seq(5,40,5), 
                   repeats=5, visu=TRUE) # 此处为出现空聚类的示例，由此可确定 c 的范围

# Dmin 函数确定最优 c, D.min 在达到最优 c 后下降得更慢
# 也可以通过观察cluster图如果相似趋势过多，则适量减少cluster的数量，
#如果一些线变化趋势特别复杂，则考虑增加cluster数量
tmp  <- Dmin(eset, m=m, crange=seq(4,40,4), repeats=5, visu=TRUE)
c <- 12


# 3. Cluster ----
## 3.1 Set output path ----
dir_cl <- "./03_Result/Cmeas_cluster/"

## 3.2 Mfuzz ----
#需要设定随机数种子，以避免再次运行时获得不同的结果
set.seed(123)
cl <- mfuzz(eset, c=19, m=m)

# 查看结果
cl$size                           # 查看每个cluster中的基因个数

head(cl$cluster[cl$cluster == 2]) # 提取某个cluster下的基因

## 3.3 Res output ----
cluster_info <- data.frame(
  Gene = names(cl$cluster),         # 基因名
  Cluster = cl$cluster,             # 每个基因的聚类编号
  Membership = cl$membership        # 每个基因的隶属度值
)
write.xlsx(cluster_info, paste0(dir_cl, "cluster_info.xlsx")) 

#隶属度值也可以表示向量之间的相似性。

## 提取形成软聚类簇α核心的基因
clusters_genes <- acore(eset,cl,min.acore=0.5)         # 设置隶属值标准
head(clusters_genes[[1]])
saveRDS(clusters_genes, paste0(dir_cl, "clusters_coregenes.rds"))

# 查看集群之间的耦合情况
overlap_clusters <- overlap(cl)
pdf(paste0(dir_cl, "mfuzz_overlap_plot.pdf"),height = 4,width = 5)
p_overlaps <- overlap.plot(cl, over = overlap_clusters, thres = 0.05)
p_overlaps
dev.off()

## 3.4 Cluster plot ----
library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)
pdf(paste0(dir_cl, "mfuzz_clusters_plot.pdf"), height = 7,width = 12)
mfuzz.plot(eset,cl,mfrow=c(3,3),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)
dev.off()

pdf(paste0(dir_cl, "mfuzz_clusters_plot01.pdf"), height = 7,width = 12)
mfuzz.plot2(eset, cl, mfrow = c(3, 3),
            centre = T, 
            x11 = F, 
            centre.lwd = 0.2)
dev.off()

# 4. Enrich analysis ----
## 4.1 Set output path ----
dir_enrich <- "./03_Result/Cmeas_cluster/MV4_11/"

## 4.2 Gene prepare ----
#cluster_info <- read.xlsx("./03_Result/Cmeas_cluster/MOLM13/cluster_info.xlsx")

# subset target cluster gene
target_genes <- names(cl$cluster[cl$cluster == 2])

# 根据需要来决定是否筛选隶属度
if(T){
  # 获取每个基因在不同 cluster 中的隶属度
  membership_values <- cl$membership[target_genes, ]
  
  # 选择隶属度大于 0.7 的基因
  selected_genes <- target_genes[membership_values[, 9] > 0.5]  # 假设2是cluster 2
  
  # 查看筛选出的基因
  selected_genes
}

target_anno <- data_anno[rownames(data_anno)%in%target_genes,]
# 提取基因名 
y <- target_anno$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
target_anno$gene <- gene


# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

# gene ID转换 
gene <- bitr(target_anno$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

## 4.3 GO ----
# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)# 设定q值阈值

## 4.4 KEGG ----
# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

#GO、KEGG结果整合 
result <- list(enrichGO = kkd, enrichKEGG = kk)
GO_res <- result$enrichGO
KEGG_res <- result$enrichKEGG

## 4.5 Res output ----
dir_enrich <- "./03_Result/Cmeas_cluster/MV4_11/"

# 导出enrichGO 
write.xlsx(GO_res@result, file = paste0(dir_enrich, "/GO_down_membership_0.5.xlsx"))
pdf(file = paste0(dir_enrich, "/GO_down_membership_0.5.pdf"), width = 6, height = 8)
p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # 控制每行最多显示40个字符
print(p1)
dev.off()

# _membership_0.5
# 导出enrichKEGG
write.xlsx(KEGG_res@result, file = paste0(dir_enrich, "/KEGG_down_membership_0.5.xlsx"))
pdf(file = paste0(dir_enrich, "/KEGG_down_membership_0.5.pdf"), width = 6, height = 5)
p2 <- dotplot(KEGG_res,showCategory = 10)
print(p2)
dev.off()


# dotplot
# ggplot2画气泡图，scale_color_gradient设置蓝红配色
library(ggplot2)
library(stringr)

# Enrich_res input 
Enrich_res <- read.xlsx("./03_Result/Cmeas_cluster/OCI/KEGG_up_all_selected.xlsx")
Enrich_res <- Enrich_res[Enrich_res$pvalue<0.05,]
# 按 p 值升序排序后，取前 20 行（最显著的 20 个结果）
Enrich_res <- head(Enrich_res, 20)                    # 取前 20 行
Enrich_res <- Enrich_res[order(-Enrich_res$pvalue), ]  # 按 p 值排序

Enrich_res$Description <- factor(Enrich_res$Description, 
                                 levels = unique(Enrich_res$Description))
#Enrich_res <- Enrich_res[order(Enrich_res$Count),]

pdf(file ="./03_Result/Cmeas_cluster/OCI/KEGG_up_all_selected.pdf",
    width = 5, height = 6 )
p2 <- ggplot(Enrich_res,aes(x=Count,y=Description))+
  geom_point(aes(size=Count,color= -log10(pvalue)))+
  theme_bw()+labs(y="",x="Count")+ 
  scale_color_gradient(low = "lightblue", high = "darkblue")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +  # 动态换行
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  
print(p2)      
dev.off()
#scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
#调整Y轴标签角度
