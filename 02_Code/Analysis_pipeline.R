# 0. Package&Function ----------------------------------------------------------
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
source("./02_Code/QC_PCA.R")
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
source("./02_Code/run_enrichment_analysis.R")

# 1. Data input ----------------------------------------------------------------
## 1.1 Group input ----
# 导入分组信息
data_group <- read.excel("./01_Data/IC50_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)

#table(data_group$group)

# 配色设置
# 配色设置
value_colour <- c("P53_WT" = "#E64B35FF",
                  "Ctrl" = "#4DBBD5FA")
                  #"Con" = "#F2A200"
rownames(data_group) <- data_group$id

## 1.2 DIA matrix input ----
fill_norm <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv")
colnames(fill_norm)
rownames(fill_norm) <- fill_norm$X
fill_norm <- fill_norm[,-1]
OCI_fill_norm <- fill_norm[,grep("OCI",colnames(fill_norm))]

# remove abnormal sample
data_group <- data_group[-grep("4W",data_group$id),]
fill_norm <- fill_norm[,-grep("4W",colnames(fill_norm))]

# Protein anno input
data_anno <- read.xlsx("./01_Data/data_anno.xlsx")
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group

# 2. Set output category -------------------------------------------------------------
#dir.create("./03_Result/GO&KEGG/MOLM13")
dir_QC <- "./03_Result/QC/P53/P53_wt/"
QC_data <- data_fill_norm
QC_group <- targeted_group
# 3. QC ------------------------------------------------------------------------
## 3.1 Boxplot -----------------------------------------------------------------
pdf(file = paste0(dir_QC,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(QC_data,data_group = QC_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
## 3.2 Heatmap -----------------------------------------------------------------
pdf(file = paste0(dir_QC,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(QC_data,data_group = QC_group,
           value_colour = value_colour)
dev.off()

## 3.3 PCA ---------------------------------------------------------------------
pdf(file = paste0(dir_QC,"QC_pca_normalization.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = log2(QC_data+1),
       data_group = QC_group,
       value_colour = value_colour)
dev.off()

# 4. DE ------------------------------------------------------------------------
# expr input
data_fill_norm <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)

# anno input
# 注意，data和data_anno的行名应一致
data_anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
data_anno <- data_anno[rownames(data_anno)%in%rownames(data_fill_norm),]
data_fill_norm <- data_fill_norm[,order(colnames(data_fill_norm))]

## 4.1 Set output catagory ----
dir_DE <- "./03_Result/Diff_Prote/P53/P53_wt_vs_Ctrl/"

## 4.2 Set group ----
data_group <- read.xlsx("./01_Data/IC50_group.xlsx")
rownames(data_group) <- data_group$id
data_group <- data_group[order(rownames(data_group)),]
table(data_group$group)
targeted_group <- data_group[c("MOLM13_6W_3","MOLM13_WT_3","MV4_11_6W_2","MV4_11_6W_3","MV4_11_WT_2","MV4_11_WT_3"),]
targeted_group$group <- gsub("High","P53_WT",targeted_group$group)

## 4.3 Res output ----
group_1 <- "Ctrl"         # group 1为Wild type
group_2 <- "High"         # group 2为Treatment
source("./02_Code/run_DE.R")

result_merge <- run_DE(data = data_fill_norm,
                       data_group = targeted_group,
                       data_anno = data_anno,
                       group_1 = group_1,
                       group_2 = group_2,
                       log2 = TRUE,
                       logfc_threshold = 0.25,         # 对应fc为1.25倍
                       pvalue_threshold = 0.05,
                       paired = TRUE,
                       pair_col = "pair_id",
                       dir = "./03_Result/")
# 统计上下调基因数量
table(result_merge$result_merge$Sig)

## 4.4 Annotated volcano plot ----

# 5. GO&KEGG ----
## 5.1 Set output catagory----
dir_enrich <- "./03_Result/GO&KEGG/P53/P53_WT/"

## 5.2 DE_res input ----
DP_result <- read.csv('./03_Result/Diff_Prote/P53/P53_wt_vs_Ctrl/P53_WT_vs_Ctrl/result_DE.csv')

## 5.3 set P.Value ----
GeneSymbol <- subset(DP_result, P.Value < 0.05)

## 5.4 set cutoff值 ----
cutoff <- 0.25                  # 对应fc约为1.25

# 转换基因名 
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
GeneSymbol$gene <- gene

if(T){
## 5.5 down genes ----
down_genes <- subset(GeneSymbol, logFC < -cutoff)

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

# gene ID转换 
gene <- clusterProfiler::bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### 5.5.1 GO ----
# GO富集分析
go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                 OrgDb = GO_database, 
                                 keyType = "ENTREZID", 
                                 ont = "ALL",           #(ALL,BP,CC,MF）
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 1)    

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
pdf(file = paste0(dir_enrich, "/GO_down.pdf"), width = 6, height = 7)
p1 <- dotplot(GO_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
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
pdf(file = paste0(dir_enrich, "/GO_up.pdf"), width = 6, height = 7)
p3 <- dotplot(GO_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
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

## 5.7 The number of up and down pathways ----
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
             
