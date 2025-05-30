#1.packages and function ----
library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pathview)
source("02_code/QC_PCA.R")
source("02_code/QC_boxplot.R")
source("02_code/QC_heatmap.R")
source("02_code/run_DE.R")
source("02_code/run_enrichment_analysis.R")
#2. data input ----
data_input <- read_delim("01_data/report.pg_matrix.tsv",
                         delim = "\t", escape_double = FALSE,
                         trim_ws = TRUE)
data_input <- data_input[,c(-1:-2,-4)]
colnames(data_input) <- c("id","C2.1", "C2.2","C2.3","C2.4","Csc.1", "Csc.2","Csc.3","Csc.4",
                          "G2.1", "G2.2","G2.3","G2.4","Gsc.1", "Gsc.2","Gsc.3","Gsc.4",
                          "V2.1", "V2.2","V2.3","V2.4","Vsc.1", "Vsc.2","Vsc.3","Vsc.4")
data_group <- read_xlsx("01_data/pg_group.xlsx")
value_colour <- c("sgCHD7" = "#E64B35FF",# Experimental group
                  "sgsc" = "#4DBBD5FF"# other group1
)
data_input <- aggregate(data_input, by = list(data_input$id), mean)
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Group.1
data_input <- subset(data_input[,-1:-2])
# 是否进行log2运算？
data_input <- log2(data_input)
# 是否进行median normalization运算？
median_values <- apply(data_input, 2, median, na.rm = TRUE)
data_input <- data_input %>%#除以每列中位数
  mutate(across(everything() , ~ ./median(., na.rm = TRUE)))
# 使用此输出进行NA填充，填充网站 https://www.omicsolution.org/wukong/NAguideR/#
write.csv(data_input,file = "01_data/fill_before.csv")

## 2.2 data fill----
data_fill <- read_csv("01_data/Pro_filtered.csv")
data_fill <- as.data.frame(data_fill)
col1 <- data_fill[,1]
rownames(data_fill) <- col1
data_fill <- data_fill[,-1]

# 如果进行了median normalization运算，运行下面的函数
data_fill <- sweep(data_fill, 2, median_values, `*`)
# 如果进行log2处理，运行下面的函数
data_fill <- 2 ^ data_fill

## 2.3 normalization ----
## intensity normalization 
data_before <- log2(data_fill + 1)
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
data_fill_normalization <- 2 ^ data_after -1
write.csv(data_fill_normalization,file = "01_data/data_fill_normalization.csv")

##2.4 QC ------
dir <- "03_result/DE/" #设置输出目录和输出PDF
pdf(file = paste0(dir,"QC_boxplot_before.pdf"),
    width = 6,
    height = 4)

QC_boxplot(data_before,data_group = data_group,
           value_colour = value_colour,
           title = "unormalized data") 
dev.off()
pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(log2(data_fill_normalization+1),data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
# heatmap
pdf(file = paste0(dir,"QC_heatmap_before.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_before,data_group = data_group,
           value_colour = value_colour) 
dev.off()
pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = log2(data_fill_normalization+1),data_group = data_group,
           value_colour = value_colour)
dev.off()

# PCA 
pdf(file = paste0(dir,"QC_pca_before.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = data_before,
       data_group = data_group,
       value_colour = value_colour)

dev.off()
pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = log2(data_fill_normalization+1),
       data_group = data_group,
       value_colour = value_colour)
dev.off()
# 3. C----
## 3.1 QC DE -----
data_fill_normalization <- read.csv("01_data/data_fill_normalization.csv")
row.names(data_fill_normalization) <- data_fill_normalization$id
data_fill_normalization <- data_fill_normalization[,-1]
group <- subset(data_group[1:8,])
data <- subset(data_fill_normalization[,1:8])
pdf(file = paste0(dir,"QC_pca_C.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = log2(data+1),
       data_group = group,
       value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_heatmap_C.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = log2(data+1),data_group = group,
           value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_boxplot_C.pdf"),
    width = 6,
    height = 4)
QC_boxplot(log2(data+1),data_group = group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
group_1 <- "sgCHD7"
group_2 <- "sgsc"
result_merge <- run_DE(data = data,
                       data_group = group,
                       log2 = T,
                      
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/CsgCHD_vs_Csg/")

res_data <- result_merge
mapkgenes <- read.csv("03_result/limma/mapkgenes in rna-seq.csv")
pi3kgenes <- read.csv("03_result/limma/pi3Kgenes in rna-seq.csv")
pi3k <- res_data[row.names(res_data) %in% pi3kgenes$SYMBOL,]
mapk <- res_data[row.names(res_data) %in% mapkgenes$SYMBOL,]

## 3.2 volcano plot ------------------------------------------------------------
data <- res_data
data$gene <- row.names(data)

# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3 sd=0.32
data$sig[data$P.Value >= 0.05 | abs(data$logFC) <= 0.6] <- "Not"

data$sig[data$P.Value < 0.05 & data$logFC > 0.6] <- "Up"

data$sig[data$P.Value < 0.05 & data$logFC < -0.6] <- "Down"

data_up <- subset(data, sig == "Up")
data_up_sorted <- data_up[order(data_up$logFC, decreasing = TRUE), ]
top10_up <- head(data_up_sorted, 10)
data_down <- subset(data, sig == "Down")
data_down_sorted <- data_down[order(data_down$logFC, decreasing = F), ]
top10_down <- head(data_down_sorted, 10)
select <- rbind(top10_up,top10_down)


volc <- ggplot(data = data, aes(x = logFC,
                                y = -log10(P.Value),
                                color = sig)) +
  scale_x_continuous(limits = c(-3,3))+
  geom_point(alpha = 0.7,size=1.5) +  theme_bw() +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#25377f","grey80","#ca0e12")) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_vline(xintercept = c(-0.6,0.6),lty=4,col="black",lwd=0.6) +
  labs(title = "CsgCHD7-Csgsc",size = 20)

p <- volc +
  geom_text_repel(data = select, 
                  aes(x = logFC, y = -log10(P.Value), label = gene),
                  show.legend = FALSE, size = 4,
                  nudge_x = 0.2,    
                  nudge_y = 0.4,      
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.linetype = 2,
                  segment.alpha = 0.5,
                  color = "black") 


ggsave("03_result/DE/volcano_C.pdf", plot = p, width = 8, height = 8)

## 3.3 KEGG GO -----------------------------------------------------------------
GeneSymbol <- subset(result_merge,P.Value< 0.05&abs(logFC)>0.3) #1192

y1 <- row.names(GeneSymbol)
GeneSymbol <- cbind(y1,GeneSymbol)
colnames(GeneSymbol)[1] <- "id"
## down-regulated genes
down_genes <- subset(GeneSymbol, logFC < 0) #641
gene.symbol.eg <- id2eg(ids = down_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_down <- result$enrichGO
kk_down <- result$enrichKEGG
dir <- "D:/2024/project/ydd/03_result/DE/"
library(openxlsx)
write.xlsx(kkd_down@result, file = paste0(dir, "/CsgCHD_vs_Csg/GO_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/CsgCHD_vs_Csg/GO_down.pdf"), width = 8, height = 8)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p1)
dev.off()

write.xlsx(kk_down@result, file = paste0(dir, "/CsgCHD_vs_Csg/KEGG_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/CsgCHD_vs_Csg/KEGG_down.pdf"), width = 6, height = 4.5)
p2 <- dotplot(kk_down,color= "pvalue")
p2
dev.off()

## up-regulated genes
up_genes <- subset(GeneSymbol, logFC > 0) #551
gene.symbol.eg <- id2eg(ids = up_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_up <- result$enrichGO
kk_up <- result$enrichKEGG

write.xlsx(kkd_up@result, file = paste0(dir, "/CsgCHD_vs_Csg/GO_up.xlsx"), sep =",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/CsgCHD_vs_Csg/GO_up.pdf"), width = 8, height = 8)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off()

write.xlsx(kk_up@result, file = paste0(dir, "/CsgCHD_vs_Csg/KEGG_up.xlsx"), sep=",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/CsgCHD_vs_Csg/KEGG_up.pdf"), width = 6, height = 4.5)
p4 <- dotplot(kk_up)
print(p4)
dev.off()
# 3.4 GSEA-----
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
dir <- "D:/2024/project/ydd/03_result/DE/"
df_C <- read_csv("03_result/DE/CsgCHD_vs_Csg/sgsc_vs_sgCHD7/result_DE.csv") #7469
df_C <- as.data.frame(df_C)
df_C <- subset(df_C[,c(1,10)])

rownames(df_C) <- df_C[,1]
SYMBOL <- row.names(df_C)
df_C <- cbind(SYMBOL,df_C)
df_C <- df_C[,-2]
df_id<-bitr(df_C$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_id <- dplyr::distinct(df_id,SYMBOL,.keep_all = T)
df_all<-merge(df_C,df_id,by="SYMBOL",all=F)
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]
gene_fc = df_all_sort$logFC 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
gseaKEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1)
gseaKEGG <- setReadable(gseaKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
sortKEGG<-gseaKEGG[order(gseaKEGG$enrichmentScore, decreasing = T),]
write.xlsx(sortKEGG, file = paste0(dir, "CsgCHD_vs_Csg/GSEA_KEGG_C.xlsx"), sep=",",quote = F, rowNames = F)

# pi3k
p <- gseaplot2(gseaKEGG, 
               geneSetID = "hsa04151", 
               title = "PI3K-Akt signaling pathway", 
               pvalue_table = TRUE)
p
ggsave("/2024/project/ydd/03_result/DE/CsgCHD_vs_Csg/GSEA_PI3K_Pathway_C.pdf", width = 8, height = 6)

# 4. G----
## 4.1 QC DE -----
group <- subset(data_group[9:16,])
data <- subset(data_fill_normalization[,9:16])
pdf(file = paste0(dir,"QC_pca_G.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = log2(data+1),
       data_group = group,
       value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_heatmap_G.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = log2(data+1),data_group = group,
           value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_boxplot_G.pdf"),
    width = 6,
    height = 4)
QC_boxplot(log2(data+1),data_group = group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()

group_1 <- "sgCHD7"
group_2 <- "sgsc"
result_merge <- run_DE(data = data,
                       data_group = group,
                       log2 = T,
                       
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/GsgCHD_vs_Gsg/")
res_data <- result_merge
mapkgenes <- read.csv("03_result/limma/mapkgenes in rna-seq.csv")
pi3kgenes <- read.csv("03_result/limma/pi3Kgenes in rna-seq.csv")
pi3k <- res_data[row.names(res_data) %in% pi3kgenes$SYMBOL,]
mapk <- res_data[row.names(res_data) %in% mapkgenes$SYMBOL,]
## 4.2 volcano plot ------------------------------------------------------------
data <- res_data
data$gene <- row.names(data)



# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
data$sig[data$P.Value >= 0.05 | abs(data$logFC) <= 0.6] <- "Not"

data$sig[data$P.Value < 0.05 & data$logFC > 0.6] <- "Up"

data$sig[data$P.Value < 0.05 & data$logFC < -0.6] <- "Down"

data_up <- subset(data, sig == "Up")
data_up_sorted <- data_up[order(data_up$logFC, decreasing = TRUE), ]
top10_up <- head(data_up_sorted, 10)
data_down <- subset(data, sig == "Down")
data_down_sorted <- data_down[order(data_down$logFC, decreasing = F), ]
top10_down <- head(data_down_sorted, 10)
select <- rbind(top10_up,top10_down)


volc <- ggplot(data = data, aes(x = logFC,
                                y = -log10(P.Value),
                                color = sig)) +
  scale_x_continuous(limits = c(-3.5,3.5))+
  geom_point(alpha = 0.7,size=1.5) +  theme_bw() +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#25377f","grey80","#ca0e12")) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_vline(xintercept = c(-0.6,0.6),lty=4,col="black",lwd=0.6) +
  labs(title = "GsgCHD7-Gsgsc",size = 20)

p <- volc +
  geom_text_repel(data = select, 
                  aes(x = logFC, y = -log10(P.Value), label = gene),
                  show.legend = FALSE, size = 4,
                  nudge_x = 0.2,    
                  nudge_y = 0.4,      
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.linetype = 2,
                  segment.alpha = 0.5,
                  color = "black") 


ggsave("03_result/DE/volcano_G.pdf", plot = p, width = 8, height = 8)

## 4.3 KEGG GO -----------------------------------------------------------------
GeneSymbol <- subset(result_merge,P.Value< 0.05&abs(logFC)>0.3) #1768

y1 <- row.names(GeneSymbol)
GeneSymbol <- cbind(y1,GeneSymbol)
colnames(GeneSymbol)[1] <- "id"
## down-regulated genes
down_genes <- subset(GeneSymbol, logFC < 0) #962
gene.symbol.eg <- id2eg(ids = down_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_down <- result$enrichGO
kk_down <- result$enrichKEGG
dir <- "D:/2024/project/ydd/03_result/DE/"
library(openxlsx)
write.xlsx(kkd_down@result, file = paste0(dir, "/GsgCHD_vs_Gsg/GO_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/GsgCHD_vs_Gsg/GO_down.pdf"), width = 8, height = 8)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p1)
dev.off()

write.xlsx(kk_down@result, file = paste0(dir, "/GsgCHD_vs_Gsg/KEGG_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/GsgCHD_vs_Gsg/KEGG_down.pdf"), width = 6, height = 4.5)
p2 <- dotplot(kk_down,color= "pvalue")
p2
dev.off()

## up-regulated genes
up_genes <- subset(GeneSymbol, logFC > 0) #806
gene.symbol.eg <- id2eg(ids = up_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_up <- result$enrichGO
kk_up <- result$enrichKEGG

write.xlsx(kkd_up@result, file = paste0(dir, "/GsgCHD_vs_Gsg/GO_up.xlsx"), sep =",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/GsgCHD_vs_Gsg/GO_up.pdf"), width = 8, height = 8)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off()

write.xlsx(kk_up@result, file = paste0(dir, "/GsgCHD_vs_Gsg/KEGG_up.xlsx"), sep=",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/GsgCHD_vs_Gsg/KEGG_up.pdf"), width = 6, height = 4.5)
p4 <- dotplot(kk_up)
print(p4)
dev.off()
## 4.4 GSEA-----
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
dir <- "D:/2024/project/ydd/03_result/DE/"
df_C <- read_csv("03_result/DE/GsgCHD_vs_Gsg/sgCHD7_vs_sgsc/result_DE.csv") #7469
df_C <- as.data.frame(df_C)
df_C <- subset(df_C[,c(1,10)])

rownames(df_C) <- df_C[,1]
SYMBOL <- row.names(df_C)
df_C <- cbind(SYMBOL,df_C)
df_C <- df_C[,-2]
df_id<-bitr(df_C$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_id <- dplyr::distinct(df_id,SYMBOL,.keep_all = T)
df_all<-merge(df_C,df_id,by="SYMBOL",all=F)
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]
gene_fc = df_all_sort$logFC 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
gseaKEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1)
gseaKEGG <- setReadable(gseaKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
sortKEGG<-gseaKEGG[order(gseaKEGG$enrichmentScore, decreasing = T),]
write.xlsx(sortKEGG, file = paste0(dir, "GsgCHD_vs_Gsg/GSEA_KEGG_G.xlsx"), sep=",",quote = F, rowNames = F)

# pi3k
p <- gseaplot2(gseaKEGG, 
               geneSetID = "hsa04151", 
               title = "PI3K-Akt signaling pathway", 
               pvalue_table = TRUE)
p
ggsave("/2024/project/ydd/03_result/DE/GsgCHD_vs_Gsg/GSEA_PI3K_Pathway_G.pdf", width = 8, height = 6)

# 5. V----
## 5.1 QC DE -----
group <- subset(data_group[17:24,])
data <- subset(data_fill_normalization[,17:24])
pdf(file = paste0(dir,"QC_pca_V.pdf"),
    width = 5,
    height = 5)
QC_PCA(data = log2(data+1),
       data_group = group,
       value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_heatmap_V.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = log2(data+1),data_group = group,
           value_colour = value_colour)
dev.off()
pdf(file = paste0(dir,"QC_boxplot_V.pdf"),
    width = 6,
    height = 4)
QC_boxplot(log2(data+1),data_group = group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()

group_1 <- "sgCHD7"
group_2 <- "sgsc"
result_merge <- run_DE(data = data,
                       data_group = group,
                       log2 = T,
                       
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/VsgCHD_vs_Vsg/")
res_data <- result_merge
pi3k <- res_data[row.names(res_data) %in% pi3kgenes$SYMBOL,]
mapk <- res_data[row.names(res_data) %in% mapkgenes$SYMBOL,]
## 5.2 volcano plot ------------------------------------------------------------
data <- res_data
data$gene <- row.names(data)

# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
data$sig[data$P.Value >= 0.05 | abs(data$logFC) <= 0.6] <- "Not"

data$sig[data$P.Value < 0.05 & data$logFC > 0.6] <- "Up"

data$sig[data$P.Value < 0.05 & data$logFC < -0.6] <- "Down"

data_up <- subset(data, sig == "Up")
data_up_sorted <- data_up[order(data_up$logFC, decreasing = TRUE), ]
top10_up <- head(data_up_sorted, 10)
data_down <- subset(data, sig == "Down")
data_down_sorted <- data_down[order(data_down$logFC, decreasing = F), ]
top10_down <- head(data_down_sorted, 10)
select <- rbind(top10_up,top10_down)


volc <- ggplot(data = data, aes(x = logFC,
                                y = -log10(P.Value),
                                color = sig)) +
  scale_x_continuous(limits = c(-4,4))+
  geom_point(alpha = 0.7,size=1.5) +  theme_bw() +
  theme(strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#25377f","grey80","#ca0e12")) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_vline(xintercept = c(-0.6,0.6),lty=4,col="black",lwd=0.6) +
  labs(title = "VsgCHD7-Vsgsc",size = 20)

p <- volc +
  geom_text_repel(data = select, 
                  aes(x = logFC, y = -log10(P.Value), label = gene),
                  show.legend = FALSE, size = 4,
                  nudge_x = 0.2,    
                  nudge_y = 0.4,      
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.linetype = 2,
                  segment.alpha = 0.5,
                  color = "black") 


ggsave("03_result/DE/volcano_V.pdf", plot = p, width = 8, height = 8)

## 5.3 KEGG GO -----------------------------------------------------------------
GeneSymbol <- subset(result_merge,P.Value< 0.05&abs(logFC)>0.3) #2838&

y1 <- row.names(GeneSymbol)
GeneSymbol <- cbind(y1,GeneSymbol)
colnames(GeneSymbol)[1] <- "id"
## down-regulated genes
down_genes <- subset(GeneSymbol, logFC < 0) #1354
gene.symbol.eg <- id2eg(ids = down_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_down <- result$enrichGO
kk_down <- result$enrichKEGG
dir <- "D:/2024/project/ydd/03_result/DE/"
library(openxlsx)
write.xlsx(kkd_down@result, file = paste0(dir, "/VsgCHD_vs_Vsg/GO_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/VsgCHD_vs_Vsg/GO_down.pdf"), width = 8, height = 8)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p1)
dev.off()

write.xlsx(kk_down@result, file = paste0(dir, "/VsgCHD_vs_Vsg/KEGG_down.xlsx"), sep=",", quote=F,rownames = F)
pdf(file = paste0(dir, "/VsgCHD_vs_Vsg/KEGG_down.pdf"), width = 6, height = 4.5)
p2 <- dotplot(kk_down,color= "pvalue")
p2
dev.off()

## up-regulated genes
up_genes <- subset(GeneSymbol, logFC > 0) #642
gene.symbol.eg <- id2eg(ids = up_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_up <- result$enrichGO
kk_up <- result$enrichKEGG

write.xlsx(kkd_up@result, file = paste0(dir, "/VsgCHD_vs_Vsg/GO_up.xlsx"), sep =",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/VsgCHD_vs_Vsg/GO_up.pdf"), width = 8, height = 8)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY",label_format=60) + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off()

write.xlsx(kk_up@result, file = paste0(dir, "/VsgCHD_vs_Vsg/KEGG_up.xlsx"), sep=",",quote = F, rowNames = F)
pdf(file = paste0(dir, "/VsgCHD_vs_Vsg/KEGG_up.pdf"), width = 6, height = 4.5)
p4 <- dotplot(kk_up)
print(p4)
dev.off()
## 4.4 GSEA-----
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
dir <- "D:/2024/project/ydd/03_result/DE/"
df_C <- read_csv("03_result/DE/VsgCHD_vs_Vsg/sgsc_vs_sgCHD7/result_DE.csv") #7469
df_C <- as.data.frame(df_C)
df_C <- subset(df_C[,c(1,10)])

rownames(df_C) <- df_C[,1]
SYMBOL <- row.names(df_C)
df_C <- cbind(SYMBOL,df_C)
df_C <- df_C[,-2]
df_id<-bitr(df_C$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_id <- dplyr::distinct(df_id,SYMBOL,.keep_all = T)
df_all<-merge(df_C,df_id,by="SYMBOL",all=F)
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]
gene_fc = df_all_sort$logFC 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
gseaKEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1)
gseaKEGG <- setReadable(gseaKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
sortKEGG<-gseaKEGG[order(gseaKEGG$enrichmentScore, decreasing = T),]
write.xlsx(sortKEGG, file = paste0(dir, "VsgCHD_vs_Vsg/GSEA_KEGG_V.xlsx"), sep=",",quote = F, rowNames = F)

# pi3k
p <- gseaplot2(gseaKEGG, 
               geneSetID = "hsa04151", 
               title = "PI3K-Akt signaling pathway", 
               pvalue_table = TRUE)
p
ggsave("/2024/project/ydd/03_result/DE/VsgCHD_vs_Vsg/GSEA_PI3K_Pathway_V.pdf", width = 8, height = 6)

# 5. venn -----
## 5.1 venn----
library(VennDiagram)
res_C1 <- read.csv("03_result/DE/CsgCHD_vs_Csg/sgCHD7_vs_sgsc/result_DE.csv")
res_G1 <- read.csv("03_result/DE/GsgCHD_vs_Gsg/sgCHD7_vs_sgsc/result_DE.csv")
res_V1 <- read.csv("03_result/DE/VsgCHD_vs_Vsg/sgCHD7_vs_sgsc/result_DE.csv")
dif_logFC_up_V <- (subset(res_V1,res_V1$logFC > 0))
dif_logFC_up_V <- dif_logFC_up_V[order(dif_logFC_up_V$logFC,decreasing = T),]
dif_logFC_up_C <- (subset(res_C1,res_C1$logFC > 0))
dif_logFC_up_C <- dif_logFC_up_C[order(dif_logFC_up_C$logFC,decreasing = T),]
dif_logFC_up_G <- (subset(res_G1,res_G1$logFC > 0))
dif_logFC_up_G <- dif_logFC_up_G[order(dif_logFC_up_G$logFC,decreasing = T),]

dif_adj_up_C <- subset(dif_logFC_up_C,P.Value< 0.05&logFC>0)
colnames(dif_adj_up_C)[1] <- "id"

dif_adj_up_G <- subset(dif_logFC_up_G,P.Value< 0.05&logFC>0)
colnames(dif_adj_up_G)[1] <- "id"

dif_adj_up_V <- subset(dif_logFC_up_V,P.Value< 0.05&logFC>0)
colnames(dif_adj_up_V)[1] <- "id"

venn_list_up <- list(group1 = dif_adj_up_C$id,
                     group2 = dif_adj_up_G$id, 
                     group3 = dif_adj_up_V$id)

pdf("03_result/DE/venn_up.pdf", width = 8, height = 8)
venn_up <- venn.diagram(venn_list_up, filename = NULL,
                        fill = c( '#CED1F0','#C0DEE2', '#E5C8BC'), alpha = 0.50, 
                        scaled = FALSE,
                        cat.col = "black", 
                        cat.cex = 1, cat.fontfamily = 'serif',
                        col = c( '#CED1F0','#C0DEE2', '#E5C8BC'), 
                        cex = 1.5, fontfamily = 'serif',
                        cat.dist = 0.05,
                        category.names = c(  "CsgCHD7-Csgsc","GsgCHD7-Gsgsc","VsgCHD7-Vsgsc"))
grid.draw(venn_up)
dev.off()

inter <- get.venn.partitions(venn_list_up)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.xlsx(inter, "03_result/DE/venn_intersections_up.xlsx", rowNames = FALSE)
## 5.2 commongenes-----
DE_C <- res_C1
DE_G <- res_G1
DE_V <- res_V1
colnames(DE_C)[1] <- "id"
colnames(DE_G)[1] <- "id"
colnames(DE_V)[1] <- "id"

genestring <- "CEACAM6, PTMS, CD48, SENP6, NCF1C, NCF1, MYOF, S100A8, S100A9, IFITM3, TES, HELZ2, THBD, OAS3, FNBP1, ME1, CD300A, CCR1, OCC1, NT5DC2, CDKN2C, NQO1, NDUFV3, CCL5, IFIT5, RUNX2, LTC4S, RNASE6, OAS2, PRG2, NRIP1, MMP14, AFF3, NIBAN2, TP53, SLC2A3, GPR84, MT-ND5, GPAT3, GLIPR1, NMRAL1, TXNIP, CD59, ALOX5AP, MMRN1, INPP4B, CCDC57, ACAD8, SOCS2, FCER1G, CST7, PTPRE, SH3PXD2B, AMPD3, TYMP, AKNA, MT-CO3, DNAJC12, MT-CYB, CRIP1, PDLIM2, NCF2, HIP1, ITGA2B, LMNA, MT-CO2, CYTH4, EXTL2, HOXA10, ISG15, PRDX3, AHR, CD14, PLPP1, NDUFA7, TMEM154, STAP1, BMI1, MT-ND4, CTSG, COX7C, SUN1, LGALS1, SDC2, ACSF2, CCND3, ADGRE2, ACOT4, AK1, LZTFL1, CTSV, ATN1, SPRY4, MT2A, MB21D2, CLEC2A, ASL, PLAUR, CHST12, VSIR, COX7A2, NDUFS4, FCGR1A, NDUFA12, MAP1LC3A, FHOD1, MAFF, EMILIN1, TNFAIP3, IFI16, MARCKSL1, FGR, AFF1, MEFV, ARHGEF12, RPS6KA5, COA5, FHL1, SNCA, CALCOCO2, CNK3/IPCEF1, ENOSF1, NOTCH1, CKLF, CD69, BNIP3L, NDUFB2, NDUFS1, QTRT2, ITGB7, SLC2A1, ANPEP, ZNF521, HOXA5, MYO1C, TTC21B, CAMK1, TPK1, PCBD1, ABI3, RIGI, NDUFS5, DOHH, PARP12, LDLR, SLC43A3, INAVA, FGD4, DNAJC1, SPRY2, SLC38A2, GSN, NUCB1, GPT2, CYLD, MX2, NDUFA4, BCL2L11, MOCS2, LIN7C, OGFRL1, PIK3R5, HCK, NFKB2, TRPS1, LCLAT1, CD44, NCS1, PLSCR1, LACTB, KCTD12, COX5B, NDUFB7, C1orf162, PRR20A;PRR20B;PRR20C;PRR20D;PRR20E, PBK, TP53I3, NDUFB8, EVL, CD53, HSCB, SORT1, ARHGAP26, NDUFA10, HERPUD1, BAIAP2, IREB2, CPPED1, NDUFV1, RALB, NCF4, DYNC2H1, SPAG5, GLTP, SLC20A1, PIK3C2A, PAK4, MALT1, CNDP2, SATB1, TFRC, GOT1, RBM23, MOB3A, NOTUM, CNN3, KRIT1, HPCAL1, UFSP2, EHD1, TPM4, ARMC5, MEF2C, PTPN6, GNL1, KIF3B, RECQL5, CRYBG1, PPFIBP1, HSH2D, ACTN1, BCAT1, AGGF1, RNASEL, ARID2, PECAM1, WBP2, ST3GAL4, TGFB1, AGMAT, CUL1, FADS2, NDUFA2, SLC1A4, PCBD2, IFT81, ALDH4A1, COA4, LCP2, SHCBP1, LYRM1, OSBPL5, KIF23, USP11, KCTD15, MYO1E, ELL, RANBP9, YRDC, SPATA33, PLEKHO2, RBM38, MYO1G, GLUL, ACSL4, NFATC2, FEZ1, RNF213, REEP5, SQSTM1, SP110, WDR54, CMC2, PARP10, NAP1L5, CST3, SFXN4, SGPL1, MLLT11, NF1, NAMPT, RRM2, LACTB2, SLC12A6, ZDHHC13, IFI30, ADI1, INF2, CYB5A, NCOA3, TOPORS, RFLNB, GCSH, UQCC1, CHP1, CRK, PFKP, PARP14, DDA1, ECD, ARRB1, NAB1, GSTM2, MPRIP, GOLGA6L4, CNN2, OGA, RASA1, HSPH1, CD276, INPPL1, LRRFIP2, DCP1A, EIF2AK2, CD151, WNK1, BAG3, KRR1, TYSND1, BSG, PDXP, CD55, ZFTRAF1, PNPO, NFIC, ODR4, H1-10, MINDY3, SLC1A5, PDLIM1, RASA2, GSDME, DLGAP5, DPCD, MLLT1, CMPK1, ECHS1, SLC16A1, PPFIA1, MCCC2, PTCD2, NDUFB3, SEPTIN9, MRPL37, APOBEC3A, VAV2, CMPK2, LARP4B, SAPCD2, CELF1, ITGB1BP1, FLNA, NCEH1, ENGASE, SLC38A1, NDUFB9, TNFAIP2, UHRF2, INPP5D, OVCA2, PSMB5, VDAC2, ASMTL, NDUFC2, ACSS1, LAP3, KIF11, PSME2, TRUB1, SPR, INPP5B, UBE2F, DDX21, KIF3A, NCAPD3, ATP6AP1, LGALS9, CTDSP1, SERPINB6, RASSF5, PNISR, PROSER2, TLE3, MTIF2, AIP, C7orf50, NAA16, ACACA, USP34, CBX8, FAM111A, AKAP10, CSNK2A2, IQGAP1, WDR6, AKAP8L, KIF14, FTSJ1, SH3BP2, SETD3, ITPR1, TRIP12, BRD7, TNIP1, GCC2, TRIM26, PCBP3, ANK2, HMCES, MOV10, ZHX2, AGL, DNAJC7, WRNIP1, AACS, COPG1, DHRS7B, PML, FKBP4, PTPN9, MID1IP1, PPP2R5D, NKAPD1, DNAJA4, CCPG1, NDE1, IRF2, HK2, LGALS8, EHD4, SLC30A1, ARFIP1, GINS2, SAR1B, CDCA4, SCYL2, GFPT1, PET100, SLC16A3, PPM1G, RABGAP1, TPD52L2, DNM2, UBE2C, GUF1, DCK, VASP, AGAP2, EXOSC9, DERPC, DTD1, RPP30, ESF1, STAM2, CBFB, TRAFD1, PDLIM5, MMP21, CD109, PBRM1, ETF1, CARS1, SUPV3L1, B4GALT3, RMDN3, AFG2A, IRGQ, NDUFA9, CTNND1, SH3BP1, CHEK2, SLC7A5, DUS1L, ATG13, DOCK10, MARK2, JMJD1C, TELO2, AP4E1, CTU2, RACGAP1, LPCAT1, PDLIM7, SDF2L1, CEP97, TNPO1, RABGGTA, ACYP1, HSD17B10, TIMM44, EMC8, GPD2, TTF1, IRF2BP2, NOSIP, PPP1R18, UBR1, HMMR, CALCOCO1, NSFL1C, EMILIN2, STAU1, THUMPD1, PTPN1, AHSA1, NIBAN1, MYL12B, LETM1, NAPG, EML4, NFKB1, TWF1, MRPL9, KIF5B, MRPL11, ABCC1, SLC3A2, UBFD1, NDUFB1, CT45A10, TOP2A, CDC37, PCM1, GTF2F2, HDGFL2, SNAP23, CMC1, MKI67, YARS1, EEF2, TOP1, CSRP1, RAF1, TRIP13, CYB5R3, LANCL2, G3BP2, TOMM34, BCL7A, EXOC6, DCUN1D5, NEDD1, EIF4G2, COA3, SH3GL1, UNC45A, IPO4, SRM, GBF1, AEBP2, SBDS, ATL3, REXO4, HSPA4, CTU1, BRD2, TFDP1, HECTD1, WWP2, TIMM8B, CBX4, EFHD2, SMAD2, GOLGA2, PRMT1, FERMT1, PSMC4, UBE2I, ARL2BP, TRIR, RPL8, MEF2D, DYNC1LI2, RBM5, KIF13B, WDHD1, GTF2F1, SMG9, TIGAR, GLRX3, LIMA1, CEP131, C9orf78, DDX52, SART3, GINS4, LMNB2, RAVER1, PAF1, TCP1, UBE2O, RBM25, PFDN1, EEF1D, CWC22, MAEA, MTOR, VCP, SBNO1, DNAJC2, PRPF40A, ANKFY1, HEATR1, UBE4A, EIF3J, MTA2, PHB2, EIF5B, EIF4B, KHSRP
"
genes <- unlist(strsplit(genestring, ", "))
gene_common <- data.frame( id = genes)
common_group1 <- left_join(gene_common,DE_C,by = "id")
common_group2 <- left_join(gene_common,DE_G,by = "id")
common_group3 <- left_join(gene_common,DE_V,by = "id")

wb <- createWorkbook()
# 添加工作表
addWorksheet(wb, "DE_C")  # 第一个工作表
addWorksheet(wb, "DE_G")  # 第二个工作表
addWorksheet(wb, "DE_V")  # 第三个工作表
addWorksheet(wb, "DE_CGV")  # 第三个工作表

# 将数据框写入工作表
writeData(wb, "DE_C", common_group1)
writeData(wb, "DE_G", common_group2)
writeData(wb, "DE_V", common_group3)
writeData(wb, "DE_CGV", c(common_group1,common_group2,common_group3))
# 保存 Excel 文件
saveWorkbook(wb, "03_result/DE/commongenes_total.xlsx", overwrite = TRUE)

## 5.3 KEGG-----
common_genes <- read.xlsx("D:/2024/project/ydd/03_result/DE/commongenes_total.xlsx",sheet = 4,rowNames = T) #485
gene.symbol.eg <- id2eg(ids = common_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_up <- result$enrichGO
kk_up <- result$enrichKEGG

write.xlsx(kk_up@result, file = paste0(dir, "KEGG_up_venn.xlsx"), sep=",",quote = F, rowNames = F)
pdf(file = paste0(dir, "KEGG_up_venn.pdf"), width = 6, height = 4.5)
p4 <- dotplot(kk_up)
print(p4)
dev.off()

## 5.4 GSEA -----
df_C <- subset(common_genes[,c(1,44)])

rownames(df_C) <- df_C[,1]
SYMBOL <- row.names(df_C)
df_C <- cbind(SYMBOL,df_C)
df_C <- df_C[,-2]
df_id<-bitr(df_C$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_id <- dplyr::distinct(df_id,SYMBOL,.keep_all = T)
df_all<-merge(df_C,df_id,by="SYMBOL",all=F)
df_all_sort <- df_all[order(df_all$mean, decreasing = T),]
gene_fc = df_all_sort$mean 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
gseaKEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1,scoreType = "pos")
gseaKEGG <- setReadable(gseaKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
sortKEGG<-gseaKEGG[order(gseaKEGG$enrichmentScore, decreasing = T),]
write.xlsx(sortKEGG, file = paste0(dir, "GSEA_KEGG_venn.xlsx"), sep=",",quote = F, rowNames = F)

# 6. 象限-----
## 6.1 象限图-----
pro_df <- read_csv("03_result/DE/VsgCHD_vs_Vsg/sgCHD7_vs_sgsc/result_DE.csv")
pro_df <- pro_df[,c(1,10,13)]
colnames(pro_df)[1] <- "id"
rna_df <- read_csv("03_result/limma/VsgCHD_vs_Vsg/DEG_results_V.csv")
rna_df <- rna_df[,c(1:2,5)]
colnames(rna_df) [1]<- "id"
data <- merge(rna_df, pro_df, by = "id", suffixes = c("_RNA","_Protein") ,all = FALSE)
#对数据进行分组；
#生成显著上下调数据标签；
data$part <- case_when(abs(data$logFC_RNA) >= 1 & abs(data$logFC_Protein) >= 0.3 ~ "part1379",
                       abs(data$logFC_RNA) < 1 & abs(data$logFC_Protein) > 0.3 ~ "part28",
                       abs(data$logFC_RNA) > 1 & abs(data$logFC_Protein) < 0.3 ~ "part46",
                       abs(data$logFC_RNA) < 1 & abs(data$logFC_Protein) < 0.3 ~ "part5")
head(data)

#开始尝试绘图；
p0 <-ggplot(data,aes(logFC_RNA,logFC_Protein,color=part))
#添加散点；
p1 <- p0+geom_point(size=1.2)+guides(color="none")
p1
#自定义半透明颜色；
mycolor <- c("#FF9999","#99CC00","#c77cff","gray80")
p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p2
#添加辅助线；
p3 <- p2+theme_bw()+geom_hline(yintercept = c(-0.3,0.3),
                    size = 0.5,
                    color = "grey40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 0.5,
             color = "grey40",
             lty = "dashed")
p3

#如果考虑差异的显著性，则需要进一步分组；
#生成至少在一个组学显著上下调的数据标签；
data$sig <- case_when(data$P.Value_RNA < 0.05 & data$P.Value_Protein <0.05 ~ "sig",
                      data$P.Value_RNA >= 0.05 | data$P.Value_Protein >=0.05 ~ "no")
part3 <- subset(data,logFC_RNA>1&logFC_Protein>0.3&sig=="sig")
write.table(part3, "03_result/DE/Vpart3.csv",            
            row.names=FALSE,col.names=TRUE,sep=",") 
head(data)
#将作图数据表格拆分成显著和不显著两部分；
sig <- filter(data,sig == "sig")
non <- filter(data,sig == "no")

#重新进行绘图；
p6 <-ggplot(data,aes(logFC_RNA,logFC_Protein))+
  geom_point(data=non,aes(logFC_RNA,logFC_Protein),size=1.2,color="gray90")+
  geom_point(data=sig,aes(logFC_RNA,logFC_Protein,color=part),size=1.5)+
  guides(color="none")
p6
#自定义半透明颜色；
mycolor <- c("#99CC00","#FF9999","#c77cff","gray90")
p7 <- p6 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p7
#添加辅助线并修改坐标轴范围；
correlation <- cor(data$logFC_RNA, data$logFC_Protein, method = "pearson")
print(correlation)

comdata <-  subset(data,data$logFC_RNA> 0 &data$logFC_Protein>0)
comdata$logFC_mean <- rowMeans(comdata[, c("logFC_RNA", "logFC_Protein")], na.rm = TRUE)
p8 <- p7+theme_bw()+geom_hline(yintercept = c(-0.3,0.3),
                    size = 0.5,
                    color = "gray40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 0.5,
             color = "gray40",
             lty = "dashed")+
  scale_y_continuous(expand=expansion(add = c(0.5, 0.5)))+
  scale_x_continuous(expand=expansion(add = c(0.5, 0.5)))+
  geom_text(x = -3, y = 2.5, size=4,label = "correlation=0.22", color="gray40")
p8
ggsave("03_result/Pro_RNA_quadrantV.pdf",width = 6,height = 6)

# 计算皮尔逊相关系数

## 6.2 GSEA----
df_C <- subset(comdata[,c(1,4)])

rownames(df_C) <- df_C[,1]
SYMBOL <- row.names(df_C)
df_C <- cbind(SYMBOL,df_C)
df_C <- df_C[,-2]
df_id<-bitr(df_C$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_id <- dplyr::distinct(df_id,SYMBOL,.keep_all = T)
df_all<-merge(df_C,df_id,by="SYMBOL",all=F)
df_all_sort <- df_all[order(df_all$logFC_Protein, decreasing = T),]
gene_fc = df_all_sort$logFC_Protein 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
gseaKEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1,scoreType = "pos")
gseaKEGG <- setReadable(gseaKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
sortKEGG<-gseaKEGG[order(gseaKEGG$enrichmentScore, decreasing = T),]
write.xlsx(sortKEGG, file = paste0(dir, "GSEA_KEGG_PR.xlsx"), sep=",",quote = F, rowNames = F)

p <- gseaplot2(gseaKEGG, 
               geneSetID = "hsa04151", 
               title = "PI3K-Akt signaling pathway", 
               pvalue_table = TRUE)
p
ggsave("/2024/project/ydd/03_result/GSEA_PI3K_Pathway_PR.pdf", width = 8, height = 6)

## 5.3 KEGG-----
common_genes <- subset(comdata[,c(1,4,5)])
common_genes <- subset(common_genes,P.Value_Protein< 0.05&abs(logFC_Protein)>0.3)
gene.symbol.eg <- id2eg(ids = common_genes$id, category = 'SYMBOL', org = 'Hs', na.rm = F)
gene.symbol.eg <- as.data.frame(gene.symbol.eg)
gene.symbol.eg <- na.omit(gene.symbol.eg)
kkd <- enrichGO(gene = gene.symbol.eg$ENTREZID,
                OrgDb = "org.Hs.eg.db",
                ont = "ALL",
                qvalueCutoff = 0.05,
                readable = TRUE)
kk <- enrichKEGG(gene = gene.symbol.eg$ENTREZID, organism = "hsa", 
                 keyType = "kegg", 
                 qvalueCutoff = 0.05)


result <- list(enrichGO = kkd, enrichKEGG = kk)
kkd_up <- result$enrichGO
kk_up <- result$enrichKEGG

write.xlsx(kk_up@result, file = paste0(dir, "KEGG_up_PR.xlsx"), sep=",",quote = F, rowNames = F)
pdf(file = paste0(dir, "KEGG_up_PR.pdf"), width = 6, height = 4.5)
p4 <- dotplot(kk_up)
print(p4)
dev.off()

# 6. heatmap-----
group <- data_group
data <- data_fill_normalization
new_order <- c("C2.1","C2.2","C2.3","C2.4",
               "G2.1","G2.2","G2.3","G2.4",
               "V2.1","V2.2","V2.3","V2.4",
               "Csc.1","Csc.2","Csc.3","Csc.4",
               "Gsc.1","Gsc.2","Gsc.3","Gsc.4",
               "Vsc.1","Vsc.2","Vsc.3","Vsc.4")

# 重新排列数据框
data <- data[, new_order]
group_1 <- "sgCHD7"
group_2 <- "sgsc"
result_merge <- run_DE(data = data,
                       data_group = group,
                       log2 = T,
                       
                       group_1 = group_1,group_2 = group_2,
                       dir = "03_result/DE/sgCHD7_vs_sg/")
res_data <- result_merge
filtered_genes <- rownames(common_genes)
protein_data <- data_fill_normalization[filtered_genes, ]

genestring <- "CSF1/NTRK1/ANGPT1/PTK2/ITGA2B/BDNF/C8orf44-SGK3/IL2RB/LAMC3/ITGB5/PDGFB/PDGFA/VEGFC/FGF8/TEK/FGF6/FLT3LG/VTN/CCND1/ITGAV/PIK3R6/NR4A1/KITLG/FLT3/EREG/IL6R/RASGRP3/RPS6KA5/CACNB1/RELB/TGFB1/IL1B/TNF
"
targetgenes <- unlist(strsplit(genestring, "/"))
scaled_data <- t(scale(t(protein_data)))
row_order <- c(setdiff(rownames(scaled_data), targetgenes), targetgenes)
library(KEGGREST)
library(ggplot2)
library(pheatmap)
mapk_pathway <- keggGet("hsa04010")

#查找所有基因 
mapkgenes<-unlist(lapply(mapk_pathway[[1]]$GENE,function(x) strsplit(x,';'))) 
mapkgenes <- mapkgenes[1:length(mapkgenes)%%3 ==2] 
mapkgenes <- data.frame(mapkgenes)  
#把结果写入表格中 
write.table(mapkgenes, "03_result/DE/hsa04010.csv",            
            row.names=FALSE,col.names=TRUE,sep=",") 

# 获取PI3K-Akt信号通路基因
pi3k_pathway <- keggGet("hsa04151")
#查找所有基因 
genes<-unlist(lapply(pi3k_pathway[[1]]$GENE,function(x) strsplit(x,';'))) 
pi3kgenes <- genes[1:length(genes)%%3 ==2] 
pi3kgenes <- data.frame(pi3kgenes)  
#把结果写入表格中 
write.table(pi3kgenes, "03_result/DE/hsa04151.csv",            
            row.names=FALSE,col.names=TRUE,sep=",") 

mapk_data <- protein_data[rownames(protein_data) %in% mapkgenes$mapkgenes, ]
pi3k_data <- protein_data[rownames(protein_data) %in% pi3kgenes$pi3kgenes, ]
new_order <- c("Csc.1","Csc.2","Csc.3","Csc.4","C2.1","C2.2","C2.3","C2.4",
               "Gsc.1","Gsc.2","Gsc.3","Gsc.4","G2.1","G2.2","G2.3","G2.4",
                "Vsc.1","Vsc.2","Vsc.3","Vsc.4","V2.1","V2.2","V2.3","V2.4")

annotation_col = data.frame(
  Treat = c(rep("Control",8),rep("Gilteritinib",8),rep("Venetoclax",8)),
  Knockout = c(rep("sc", 4), rep("CHD7",4),rep("sc", 4), rep("CHD7",4),rep("sc", 4), rep("CHD7",4)))

row.names(annotation_col) <- colnames(protein_data)           
# 重新排列数据框
pi3k_data <- pi3k_data[, new_order]

pi3k<- pheatmap(pi3k_data, 
         scale = "row",  # 按行标准化（基因的表达）
         cluster_rows = TRUE,  # 行聚类
         cluster_cols = F,  # 列聚类
         show_rownames = TRUE,  # 显示基因名
         show_colnames = TRUE,  # 显示样本名
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(50),  # 设置颜色
         main = "PI3K-Akt Signaling Pathway Protein Expression")
pi3k
ggsave(plot=pi3k,"03_result/DE/PI3K_genes.pdf",width = 8,height = 8)
row_order <- pi3k$tree_row$order  # 获取当前的行聚类顺序
reversed_row_order <- rev(row_order)  # 反转顺序
reversed_pi3k_data <- pi3k_data[rev(row_order), ]  # 反转行顺序
mapk_reversed <- pheatmap(reversed_pi3k_data, 
                          scale = "row",  
                          cluster_rows = FALSE,  # 禁用行聚类
                          cluster_cols = FALSE,  # 禁用列聚类
                          show_rownames = TRUE,  
                          show_colnames = TRUE,  
                          color = colorRampPalette(c("blue", "white", "red"))(50),  
                          main = "PI3K-Akt Signaling Pathway Protein Expression")
ggsave(plot = mapk_reversed, "03_result/DE/heatmap_reversed_pi3k.pdf", width = 8, height = 8)


mapk_data <- mapk_data[, new_order]
mapk <-pheatmap(mapk_data, 
                scale = "row",  # 按行标准化（基因的表达）
                cluster_rows = TRUE,  # 行聚类
                cluster_cols = F,  # 列聚类
                show_rownames = TRUE,  # 显示基因名
                show_colnames = TRUE,  # 显示样本名
                color = colorRampPalette(c("blue", "white", "red"))(50),  # 设置颜色
                main = "MAPK Signaling Pathway Protein Expression")
mapk
ggsave(plot=mapk,"03_result/DE/mapk_genes.pdf",width = 8,height = 8)

#row_order <- mapk$tree_row$order  # 获取当前的行聚类顺序
#reversed_row_order <- rev(row_order)  # 反转顺序
#reversed_mapk_data <- mapk_data[rev(row_order), ]  # 反转行顺序
#mapk_reversed <- pheatmap(reversed_mapk_data, 
                          scale = "row",  
                          cluster_rows = FALSE,  # 禁用行聚类
                          cluster_cols = FALSE,  # 禁用列聚类
                          show_rownames = TRUE,  
                          show_colnames = TRUE,  
                          color = colorRampPalette(c("blue", "white", "red"))(50),  
                          main = "MAPK Signaling Pathway Protein Expression")
ggsave(plot = mapk, "03_result/DE/heatmap_reversed_row_order.pdf", width = 8, height = 8)

rownames(pi3k_data) %in% targetgenes
com <- mapk_data[rownames(mapk_data) %in% targetgenes, ]
