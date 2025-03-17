# ClusterGVis

# 1. Packages ----
# Note: please update your ComplexHeatmap to the latest version!
library(pak)
pak::pak("junjunlab/ClusterGVis")
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
source("./02_Code/extract_expr_by_cell.R")
# 2. Data input ----
# load data(standardised)
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)

# anno input
# Attention! The rownames of data and data_anno should be consistent
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
expr <- expr[!is.na(expr$gene), ]  # remove NA rows

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

# 单个细胞系分析（根据目的是否需要这一步）
process_cell_line("OCI")  # 生成 expr_cellline

## 2.2 caculate average ----
### 2.2.1 approach 1 -----------------------------------------------------------
library(limma)
# limma::avereps：这是来自 limma 包的函数，avereps 用于对重复数据进行平均值计算。
# avereps 会根据指定的 ID 进行分组，并对相同 ID 的数据取平均值
avereps_df <- t(limma::avereps(t(expr_OCI) , ID = colnames(expr_OCI))) #对相同时间序列的表达值取平均
avereps_df[1:3,1:3]
avereps_df <- avereps_df[,c("Ctrl","Low","High")]
save(avereps_df,file = './03_Result/C_means_cluster/IC50 group/OCI-M2/avereps_df.Rdata')

### 2.2.2 approach 2 -----------------------------------------------------------
# 处理列名，去掉 `.1`, `.2` 这种后缀
colnames(expr_molm13) <- sub("\\..*", "", colnames(expr_molm13))
selected_mean <- apply(expr_molm13[, c(2, 7, 8)], 1, mean, na.rm = TRUE)
selected_mean <- as.data.frame(selected_mean)
expr_molm13[,"34"] <- selected_mean$selected_mean
expr_molm13 <- expr_molm13[,-c(7,8)]

# 按列名数值大小排序
expr_molm13 <- expr_molm13[, order(as.numeric(colnames(expr_molm13)))]
avereps_df <- expr_molm13

# 3. Cluster ----
## 3.1 Set output path ----
dir_cl <- "./03_Result/C_means_cluster/IC50 group/MOLM13/"
exps <- as.matrix(avereps_df) 

# check optimal cluster numbers
pdf(paste0(dir_cl,"Elbow_plot.pdf"),width = 6,height = 5)
getClusters(obj = exps)
dev.off()
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
exp <- filter.std(exps,min.std ,visu = FALSE)
exp <- exps
cluster.num = 6
cluster.method = "kmeans"
# whether zsocre data
if(TRUE){
  hclust_matrix <- exp %>% t() %>% scale() %>% t()
}else{
  hclust_matrix <- exp
}
# add kmeans func n stats
km <- stats::kmeans(x = hclust_matrix,centers = cluster.num,nstart = 10)

od.res <- data.frame(od = match(names(km$cluster),rownames(hclust_matrix)),
                     id = as.numeric(km$cluster),
                     check.names = FALSE)


cl.info <- data.frame(table(od.res$id),check.names = FALSE)

# reorder matrix
m <- hclust_matrix[od.res$od,]

# add cluster and gene.name
wide.r <- m %>%
  data.frame(check.names = FALSE) %>%
  dplyr::mutate(gene = rownames(.),
                cluster = od.res$id) %>%
  dplyr::arrange(cluster)

# whether subset clusters
if(!is.null(subcluster)){
  wide.r <- wide.r %>% dplyr::filter(cluster %in% subcluster)
}

# wide to long
df <- reshape2::melt(wide.r,
                     id.vars = c('cluster','gene'),
                     variable.name = 'cell_type',
                     value.name = 'norm_value')

# add cluster name
df$cluster_name <- paste('cluster ',df$cluster,sep = '')

# add gene number
cltn <- table(wide.r$cluster)
purrr::map_df(unique(df$cluster_name),function(x){
  tmp <- df %>%
    dplyr::filter(cluster_name == x)
  
  cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])
  
  tmp %>%
    dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
}) -> df

# cluster order
df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                         " (",cl.info$Freq,")",sep = ''))
# cluster list
wide <- wide.r
wide$cluster <- paste("C",wide$cluster,sep = "")
cluster.list <- split(wide$gene,wide$cluster)

# return
km <- list(wide.res = wide.r,
            long.res = df,
            cluster.list = cluster.list,
            type = cluster.method,
            geneMode = "none",
            geneType = "none")


# 4. Plot ----
## 4.1 plot line only ----
visCluster(object = cm,
           plot.type = "line")

# change color
visCluster(object = km,
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

##  5.1 enrich for clusters ----
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

## 5.2 plot ----
cairo_pdf(paste0(dir_cl,'Cluster_enrichplot.pdf'),
          height = 10,width = 11,onefile = F)
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

## 5.3 all Enrichpathway ----
# The above fun() only shows the top 5 pathways, so we need the next fun() to display all pathways.
# Gene prepare
source("./02_Code/Extract_genes.R")
load("./03_Result/C_means_cluster/IC50 group/OCI-M2/C-Means_res/C-Means_res.Rdata")

# subset target cluster gene
targeted_genes <- extract_genes(cm, cl = 4, membership_threshold = 0)

# set database
GO_database <- 'org.Hs.eg.db'  # GO is org.Hs.eg.db database
KEGG_database <- 'hsa'         # KEGG is hsa database

# gene ID转换 
gene <- clusterProfiler::bitr(targeted_genes, fromType = 'SYMBOL', 
                              toType = 'ENTREZID', OrgDb = GO_database)

# GO 
# GO富集分析
go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                 OrgDb = GO_database, 
                                 keyType = "ENTREZID", 
                                 ont = "ALL",          # (ALL,BP,CC,MF）
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 1)   

# KEGG 
# KEGG富集分析
kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                  keyType = "kegg",
                                  organism = KEGG_database,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1)

#GO、KEGG结果整合 
result <- list(enrichGO = go, enrichKEGG = kegg)
GO_res <- result$enrichGO
KEGG_res <- result$enrichKEGG

# Res output
dir_enrich <- "./03_Result/C_means_cluster/IC50 group/OCI-M2/Cluster_4/"
if(T){
# output enrichGO 
write.xlsx(GO_res@result, file = paste0(dir_enrich, "/GO_membership_0.xlsx"))
pdf(file = paste0(dir_enrich, "/GO_membership_0.pdf"), width = 6, height = 8)
p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # Show a maximum of 40 characters 
print(p1)
dev.off()

# output enrichKEGG
write.xlsx(KEGG_res@result, file = paste0(dir_enrich, "/KEGG_membership_0.xlsx"))
pdf(file = paste0(dir_enrich, "/KEGG_membership_0.pdf"), width = 6, height = 5)
p2 <- dotplot(KEGG_res,showCategory = 10)
print(p2)
dev.off()
}

# dotplot
# ggplot2画气泡图，scale_color_gradient设置蓝红配色
library(ggplot2)
library(stringr)

# Enrich_res input 
Enrich_res <- read.xlsx("./03_Result/C_means_cluster/")
Enrich_res <- Enrich_res[Enrich_res$pvalue<0.05,]
# 按 p 值升序排序后，取前 20 行（最显著的 20 个结果）
Enrich_res <- head(Enrich_res, 20)                    # 取前 20 行
Enrich_res <- Enrich_res[order(-Enrich_res$pvalue), ]  # 按 p 值排序

Enrich_res$Description <- factor(Enrich_res$Description, 
                                 levels = unique(Enrich_res$Description))
#Enrich_res <- Enrich_res[order(Enrich_res$Count),]

pdf(file =paste0(dir_enrich,"KEGG_dw_ms_0_selected.pdf") ,
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
