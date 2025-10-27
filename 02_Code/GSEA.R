#  Rm workspace ----
rm(list=ls())
setwd("../AML_project/aml_analysis")

# 1. Library packages ----
#BiocManager::install("msigdbr")
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(curl)
library(msigdbr)
library(openxlsx)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
library(KEGG.db)

# 2. Set output path ----
dir <- "./03_Result/GSEA/OCI_M2_single_fill/High.vs.Ctrl/"
dir.create("./03_Result/GSEA/OCI_M2_single_fill/High.vs.Ctrl/C6/")

# 3. Preranked gene list ----
data_gene <- read.csv("./03_Result/DEP/OCI_AML2_single_fill/High_vs_Ctrl/result_DE.csv")
head(data_gene)
colnames(data_gene) 

# 提取基因名（symbol）
y <- data_gene$Genes
SYMBOL <- unlist(lapply(y,
                        function(y) strsplit(as.character(y),";")[[1]][1]))
id <- data.frame(SYMBOL)
diff_prote <- cbind(id,data_gene)

# 检查有无重复symbol
duplicates <- duplicated(SYMBOL)
a <- print(diff_prote[duplicates,])

#如果有则进行以下处理
# 只保留 SYMBOL 唯一且 AveExpr 最高的基因
diff_prote <- diff_prote[order(diff_prote$AveExpr, decreasing = TRUE), ]  # 按 AveExpr 降序排序
diff_prote <- diff_prote[!duplicated(diff_prote$SYMBOL), ]                # 去除重复 SYMBOL，保留第一个

# switch to ENTREZID names 
# Database set 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库

# gene ID转换 
ENTREZID_gene <- clusterProfiler::bitr(diff_prote$SYMBOL, 
                                       fromType = 'SYMBOL', 
                                       toType = 'ENTREZID', 
                                       OrgDb = GO_database)

# rank gene list 
data_all <- diff_prote %>% 
  inner_join(ENTREZID_gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$ENTREZID 

# 4. KEGG_GSEA ---- 
# msigdbr基因集中没有KEGG集，需要单独做
# database set 
KEGG_database="hsa"

# gsea analysis 
KEGG_GSEA <- gseKEGG(geneList, 
                     organism = KEGG_database, 
                     pvalueCutoff = 0.5,   # 实际上是adjp，为了显示更多p<0.05的通路，设置为0.5
                     eps = 0,
                     use_internal_data =  T,
                     nPermSimple = 10000) # 随机排列次数，提高结果准确性
KEGG_GSEA<- setReadable(KEGG_GSEA,        # 转换可读基因名
                        OrgDb=org.Hs.eg.db,
                        keyType = 'ENTREZID')

table(KEGG_GSEA@result$pvalue<0.05)      # 查看有多少个通路富集出来

# res output
KEGG_results <- as.data.frame(KEGG_GSEA)
write.xlsx(KEGG_results,file = paste0(dir,"KEGG/KEGG_results.xlsx"))
save(KEGG_GSEA, file = paste0(dir,"KEGG/KEGG_result.RData"))   # S4 res

# plot 
load("./03_Result/GSEA/MOLM13/VEN_VS_WT/kegg_result.RData")

pdf(paste0(dir,"gseaplot.pdf"),width = 11,height = 9)
#dotplot(KEGG_GSEA,showCategory = Inf,label_format = 100)
#ridgeplot(KEGG_GSEA,label_format = 100) # 1000*800
gseaplot2(KEGG_GSEA,c("hsa00565"),pvalue_table = T) 
dev.off()

# 5. MsigDBR_GSEA ---- 

## H ----
sigDBGR_H <- msigdbr(species = "Homo sapiens",
                     category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_h <- GSEA(geneList, 
                       TERM2GENE = sigDBGR_H,       # 预设基因集
                       minGSSize = 1,          # 最小基因集大小
                       maxGSSize = 500,        # 最大基因集大小
                       eps = 0,                # 防止p值过小
                       pvalueCutoff = 0.5,    # 实际上是adjp，为了显示更多p<0.05的通路，设置为0.5
                       pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_h <- setReadable(gsea_results_h, 
                              OrgDb = org.Hs.eg.db, 
                              keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_h@result$pvalue < 0.05)

# 结果导出
h_results <- as.data.frame(gsea_results_h)
write.xlsx(h_results,file = paste0(dir,"H/H_results.xlsx"))
save(gsea_results_h, file = paste0(dir,"H/H_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_h, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_h, label_format = 100)
#gseaplot2(gsea_results_h, 1, pvalue_table = TRUE)  # 这里的 `1` 是基因集的索引

## C1 ----
# 染色体位置
sigDBGR_C1 <- msigdbr(species = "Homo sapiens",
                      category = "C1") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C1 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C1,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 1,    # 实际上是adjp，为了显示更多p<0.05的通路，设置为0.5
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C1 <- setReadable(gsea_results_C1, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C1@result$pvalue < 0.05)

# 结果导出
C1_results <- as.data.frame(gsea_results_C1)
write.xlsx(C1_results,file = paste0(dir,"C1/C1_results.xlsx"))
save(gsea_results_C1, file = paste0(dir,"C1/C1_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_C1, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C1, label_format = 100)
#gseaplot2(gsea_results_C1, 1, pvalue_table = TRUE)

## C2 ----
# 来源于实验研究、文献和通路数据库（如 KEGG、Reactome、BioCarta）
colnames(msigdbr(species = "Homo sapiens", collection = "C2"))
sigDBGR_C2 <- msigdbr(species = "Homo sapiens",
                      collection = "C2") %>%
  dplyr::select(gs_name, ncbi_gene)
# GSEA分析
gsea_results_C2 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C2,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 0.5,    # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C2 <- setReadable(gsea_results_C2, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C2@result$pvalue < 0.05)

# 结果导出
C2_results <- as.data.frame(gsea_results_C2)
write.xlsx(C2_results,file = paste0(dir,"C2/C2_results.xlsx"))
save(gsea_results_C2, file = paste0(dir,"C2/C2_results.RData"))   # S4 res


# 可视化
#dotplot(gsea_results_C2, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C2, label_format = 100)
#gseaplot2(gsea_results_C2, 1, pvalue_table = TRUE)

## C3 ----
# motif富集分析
sigDBGR_C3 <- msigdbr(species = "Homo sapiens",
                      category = "C3") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C3 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C3,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 1,    # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C3 <- setReadable(gsea_results_C3, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C3@result$pvalue < 0.05)

# 结果导出
C3_results <- as.data.frame(gsea_results_C3)
write.xlsx(C3_results,file = paste0(dir,"C3/C3_results.xlsx"))
save(gsea_results_C3, file = paste0(dir,"C3/C3_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_C3, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C3, label_format = 100)
#gseaplot2(gsea_results_C3, 1, pvalue_table = TRUE)

## C4 ----
sigDBGR_C4 <- msigdbr(species = "Homo sapiens",
                      category = "C4") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C4 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C4,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 1,    # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C4 <- setReadable(gsea_results_C4, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C4@result$pvalue < 0.05)

# 结果导出
C4_results <- as.data.frame(gsea_results_C4)
write.xlsx(C4_results,file = paste0(dir,"C4/C4_results.xlsx"))
save(gsea_results_C4, file = paste0(dir,"C4/C4_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_C4, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C4, label_format = 100)
#gseaplot2(gsea_results_C4, 1, pvalue_table = TRUE)

## C5 ----
# GO 富集分析 ，分析时间较长
sigDBGR_C5 <- msigdbr(species = "Homo sapiens",
                      category = "C5") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C5 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C5,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 0.5,    # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C5 <- setReadable(gsea_results_C5, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C5@result$pvalue < 0.05)

# 结果导出
C5_results <- as.data.frame(gsea_results_C5)
write.xlsx(C5_results,file = paste0(dir,"C5/C5_results.xlsx"))
save(gsea_results_C5, file = paste0(dir,"C5/C5_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_C5, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C5, label_format = 100)
#gseaplot2(gsea_results_C5, 1, pvalue_table = TRUE)

## C6 ----
# 癌症相关表达谱
sigDBGR_C6 <- msigdbr(species = "Homo sapiens",
                      category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C6 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C6,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 1,      # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C6 <- setReadable(gsea_results_C6, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C6@result$pvalue < 0.05)

# 结果导出
C6_results <- as.data.frame(gsea_results_C6)
write.xlsx(C6_results,file = paste0(dir,"C6/C6_results.xlsx"))
save(gsea_results_C6, file = paste0(dir,"C6/C6_results.RData"))   # S4 res

# 可视化
#dotplot(gsea_results_C6, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C6, label_format = 100)
#gseaplot2(gsea_results_C6, 1, pvalue_table = TRUE)

## C7 ----
# 免疫相关表型
sigDBGR_C7 <- msigdbr(species = "Homo sapiens",
                      category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)
# GSEA分析
gsea_results_C7 <- GSEA(geneList, 
                        TERM2GENE = sigDBGR_C7,       # 预设基因集
                        minGSSize = 1,          # 最小基因集大小
                        maxGSSize = 500,        # 最大基因集大小
                        eps = 0,                # 防止p值过小
                        pvalueCutoff = 1,    # p 值阈值
                        pAdjustMethod = "BH")   # p 值调整方法
# 设置为可读格式 
gsea_results_C7 <- setReadable(gsea_results_C7, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENTREZID')

# 统计有多少个基因集富集
table(gsea_results_C7@result$pvalue < 0.05)

# 结果导出
C7_results <- as.data.frame(gsea_results_C7)
write.xlsx(C7_results,file = paste0(dir,"C7/C7_results.xlsx"))
save(gsea_results_C7, file = paste0(dir,"C7/C7_results.RData"))   # S4 res


# 可视化
#dotplot(gsea_results_C7, showCategory = Inf, label_format = 100)
#ridgeplot(gsea_results_C7, label_format = 100)
#gseaplot2(gsea_results_C7, 1, pvalue_table = TRUE)
