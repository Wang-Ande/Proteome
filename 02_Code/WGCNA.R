# WGCNA

# 加载包
if(!require(biomaRt ,quietly = TRUE)){BiocManager::install("biomaRt")}
library(biomaRt)
library(readr)
library(openxlsx)
library(readxl)
library(dplyr)
library(scales)
library(tidyr)
library(WGCNA)
library(forcats)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
# 0. 设置输出目录 ----
folder_path <- "./03_Result/WGCNA/ALL_gene"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

dir_WGCNA <- "./03_Result/WGCNA/ALL_gene/"

# 1.data input ----
# 差异基因取并集
DEG1 <- read.csv("./03_Result/DE/MOLM13/MOLM13_2W_vs_MOLM13_WT/result_DE.csv")
DEG2 <- read.csv("./03_Result/DE/MOLM13/MOLM13_6W_vs_MOLM13_WT/result_DE.csv")
DEG3 <- read.csv("./03_Result/DE/MV4_11/MV4_11_2W_vs_MV4_11_WT/result_DE.csv")
DEG4 <- read.csv("./03_Result/DE/MV4_11/MV4_11_6W_vs_MV4_11_WT/result_DE.csv")
DEG5 <- read.csv("./03_Result/DE/OCI_AML2/OCI_2W_vs_OCI_WT/result_DE.csv")
DEG6 <- read.csv("./03_Result/DE/OCI_AML2/OCI_6W_vs_OCI_WT/result_DE.csv")

DEG1 <- DEG1[DEG1$Sig != "stable", 3]
DEG2 <- DEG2[DEG2$Sig != "stable", 3]
DEG3 <- DEG3[DEG3$Sig != "stable", 3]
DEG4 <- DEG4[DEG4$Sig != "stable", 3]
DEG5 <- DEG5[DEG5$Sig != "stable", 3]
DEG6 <- DEG6[DEG6$Sig != "stable", 3]

# 取并集
union_genes <- Reduce(union, list(DEG1, DEG2, DEG3, 
                                  DEG4, DEG4, DEG6))

## 1.2 数据预处理 ----
# proteinExpr_input
Expr_raw <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv")
rownames(Expr_raw) <- Expr_raw$X
Expr_raw <- Expr_raw[,-1]
min(Expr_raw)
Expr_raw <- Expr_raw[rownames(Expr_raw)%in%union_genes,]

datExpr <- Expr_raw

## 1.3 标准化 log2+1 ----
datExpr <- log2(datExpr)

## 1.4 filter ---- 
# WGCNA可以不过滤，WGCNA包中有自己的过滤方案
# var计算每个基因方差，筛选基因变化较大的基因，此次选取前75%的基因
if(T){
vars_res <- apply(datExpr, 1, var)

# 计算百分位数截止值
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) # quantile生成分位数
per_res

upperGene <- datExpr[ which(vars_res > per_res[4]),]  # 选取方差位于前75%的基因
dim(upperGene)
datExpr <- data.matrix(upperGene)
}

# 检查是否有目的基因被过滤掉
targetgene_rows <- grep("BAK1", rownames(datExpr))
targetgene <- datExpr[targetgene_rows,]
# TP53\BCL2\MCL1\BAX\BAK1基因都存在

# 2 t转置 ----
# 为了执行 WGCNA，通常需要转置数据，使其符合 WGCNA 的标准格式
datExpr <- t(datExpr)                       # 让每一列代表一个基因

# 3. 检查数据好坏 ----
gsg <- goodSamplesGenes(datExpr, verbose = 3)

# 查看标记结果
gsg$allOK
#[1] TRUE

# 如果为false则运行下段
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", 
                     paste(names(datExpr)[!gsg$goodSamples], collapse = ", ")))
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

ncol(datExpr)

# 对样本聚类查看有无异常样本
sampleTree = hclust(dist(datExpr), method = "average")
png(paste0(dir_WGCNA, "sample_clustering.png"), 
    width = 1500, 
    height = 1000,
    res = 300)
par(cex = 0.6)
par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2)
#abline(h = 150, col = 'red')
dev.off()

# 剔除离群样本，无则跳过
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10) 
table(clust)               #显示剔除多少样本， 0代表剔除的个数，1代表保留的个数
datExpr = datExpr[clust == 1, ]

# or
datExpr <- datExpr[-grep("MV4_11_6W_1",rownames(datExpr)),] # delete MV4_6W_1
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 4. 选择软阈值 ----
# 定义1:30的power值，用于构建不同的基因共表达网络
powers <- c(c(1:10), 
            seq(from = 12, 
                to = 30,
                by = 2))

sft <-  pickSoftThreshold(datExpr,
                          powerVector = powers,
                          verbose = 5,
                          networkType = 'unsigned')
sft$powerEstimate
# 绘制拟合度图来选择合适的软阈值
pdf(file = paste0(dir_WGCNA, "sft_par.pdf"),width = 10,height = 6.5)
#sizeGrWindow(9, 5) 
par(mfrow = c(1,2))
#cex1 = 0.85
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     col="steelblue")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1],  sft$fitIndices[,5], 
     type="n", 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,  col="steelblue")
dev.off()

# 使用合适的软阈值（没有合适就选6）
power <- 16

# 5. 构建共表达矩阵 ----
net = blockwiseModules(datExpr,
                       power = power,           # 软阀值选为6
                       TOMType = "unsigned",    # 构建无尺度网络
                       minModuleSize = 30,      # 最小模块基因数为30
                       reassignThreshold = 0,   # 重分配阈值          
                       mergeCutHeight = 0.25,   # 模块合并阀值
                       numericLabels = F,       # 返回字符标签（如模块颜色名称）
                       pamRespectsDendro = FALSE,         
                       saveTOMs = FALSE,         # 不保存TOM矩阵
                       verbose = 3,
                       maxBlockSize = 20000)    # 可处理的最大模块基因数
# 显示所有模块个数和各个模块中基因数量
table(net$colors)

## 5.1 聚类分析 ----
# 使用层次聚类
geneTree = net$dendrograms[[1]] 
geneTree

png(paste0(dir_WGCNA, "Module_colors_clusterdendrogram.png"),
    width = 2000,
    height = 1500,
    res = 300,
    type = "cairo")
moduleColors = net$colors
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",             # 标题
                    dendroLabels = FALSE,        # 不显示基因标签
                    hang = 0.03,                 # 树枝悬挂高度（相对于根部）
                    addGuide = TRUE,             # 颜色引导条
                    guideHang = 0.05)            # 颜色引导条与树状图悬挂距离
dev.off()


# 6. 模块与表型的关系 ----

# 假设 `traitData` 是你样本的表型数据，可以是疾病分组、临床特征等
trait <- read.csv("./01_Data/trait.csv")
rownames(trait) <- trait$X
trait <- trait[,-1]
colnames(trait)


if(T){
  # 模块特征
  MEs0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  # 计算模块特征向量
  MEs = orderMEs(MEs0)                                       # 对模块特征向量排序
  
  # 计算模块特征向量与表型的相关系数矩阵
  moduleTraitCor <- cor(MEs,trait,use = "p")  
  
  # 计算相关系数矩阵的p值
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  
  
  # 构建绘图时用的文本矩阵
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  
  
  # 修改文本矩阵的维度，与相关系数矩阵相同
  dim(textMatrix)=dim(moduleTraitCor)  
}

### 6.2 相关性热图 ----
if(T){
  pdf(file = paste0(dir_WGCNA, "Module_Trait labeledHeatmap.pdf"),
      width = 12,height = 8)
  
  # mar（）分别代表图形边距的底部、左侧、顶部和右侧的边距
  par(mar = c(7, 7, 2, 2)) 
  labeledHeatmap(Matrix = moduleTraitCor,        # 绘制带标签的热图
                 xLabels = colnames(trait),     # x轴标签
                 yLabels = names(MEs),           # y轴标签
                 ySymbols = names(MEs),          # y轴符号
                 colorLabels = FALSE,            # 不显示颜色标签
                 colors = blueWhiteRed(50),      # 颜色范围
                 textMatrix = textMatrix,        # 显示文本矩阵
                 setStdMargins = FALSE,          # 不设置标准边距
                 cex.text = 0.5,                 # 文本大小
                 cex.lab.x = 0.7,                # X轴文本大小
                 zlim = c(-1,1),                 # 颜色映射范围
                 main = paste("Module-trait relationships"))  # 绘图标题
  dev.off()
}


df <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(df,file = paste0(dir_WGCNA,"module_gene.csv"))

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
#colnames(plot_data) <- paste0("Cluster", seq(ncol(plot_data)))
plot_data <- as.matrix(plot_data)
pdf(file =  "03_result/WGCNA/10.pathwayheatmap.pdf", width = 10, height = 10)
p<- pheatmap(plot_data,name =" -log(pvalue)",
             cluster_cols = T,
             cluster_rows = F,
             color = colorRampPalette(c("white", "blue"))(100))
             
### 6.3 基因与性状和重要模块的关系：基因重要性和模块成员 ----
IC50 = as.data.frame(trait$IC50)
names(Time) = "IC50"
modNames = substring(names(MEs), 3)                              # 提取模块名称
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")      # 修改列名
names(MMPvalue) = paste("p.MM", modNames, sep="")                # 修改列名

# 基因显著性
geneTraitSignificance = as.data.frame(cor(datExpr, IC50, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", names(IC50), sep="")
names(GSPvalue) = paste("p.GS.", names(IC50), sep="")

### 6.4 模内分析：鉴定具有高GS和MM的基因 ----
module = "blue"                      # 选择模块
column = match(module, modNames)     # 匹配模块名称
moduleGenes = moduleColors==module   # 提取模块内基因
table(moduleGenes)

pdf(paste0(dir_WGCNA,"Module_blue membership vs gene significance.pdf"),
    width = 6,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for IC50",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")
dev.off()

### 6.5 提取指定模块的基因做PPI网络分析 ----

adjacency = adjacency(datExpr, power =6)  # 计算邻接矩阵
TOM = TOMsimilarity(adjacency)             # 计算拓扑重叠矩阵（TOM）

module="royalblue"                               # 选择要导出的模块
probes = colnames(datExpr)                  # 获取基因名称
inModule = (moduleColors==module)           # 找到属于当前模块的基因
modProbes=probes[inModule]                  # 提取属于当前模块的基因名称
head(modProbes)                             # 显示基因名称前几行
modTOM = TOM[inModule,inModule]             # 提取属于当前模块的基因之间的TOM值
dimnames(modTOM)=list(modProbes,modProbes)  # 修改维度名称

# 绘制模块基因表达量箱线图
library(reshape2)
library(ggplot2)

datExpr_royalblue <- datExpr[,modProbes]                   #提取基因模块内基因 
datExpr_royalblue <- t(datExpr_royalblue)

# 添加分组信息
group <- traitData[,-3]
colnames(group)[2] <- "group"
colnames(group)[1] <- "id"
group$id <- rownames(group)
# 计算起始和结束位置
start = nchar(group$id) - 3  # 倒数第四的位置
end = nchar(group$id) - 2    # 倒数第二的位置
# 提取子字符串
result = substr(group$id, start, end)
print(result)
group$group <- result

table(group$group)
value_colour <- c("WT" = "#E64B35FF",# Experimental group
                  "2w" = "#4DBBD5FA",# other group1
                  "6w" = "#F2A200")# other group2
# 提取不同细胞系的时间分组
group_molm13 <- group[grepl("MOLM13",group$id),]
group_mv4_11 <- group[grepl("MV4_11",group$id),]
group_oci <- group[grepl("OCI_AML2",group$id),]

pdf(file = paste0(dir_WGCNA, "QC_boxplot_oci.pdf"),
    width = 6,
    height = 4)
QC_boxplot(2^datExpr_royalblue,data_group = group_oci,
           value_colour = value_colour,title = paste("Module",module))
dev.off()


# 这里只选了top100的基因
nTop=165                                       # 设置要选择的基因数目
IMConn = softConnectivity(datExpr[,modProbes])  # 计算当前模块中基因之间的相似性
top=(rank(-IMConn)<=nTop)                      # 找到相似性排名前nTop的基因
filterTOM=modTOM[top,top]                      # 提取相似性排名前nTop的基因之间的TOM值
# for visANT
vis = exportNetworkToVisANT(filterTOM,
                            file = paste("visANTinput-",module,".txt",sep = ""),
                            weighted = T,threshold = 0)  

# for cytoscape
cyt = exportNetworkToCytoscape(filterTOM,
                               edgeFile = paste("./03_result/WGCNA/Count/CytoscapeInput-edges-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               nodeFile = paste("./03_result/WGCNA/Count/CytoscapeInput-nodes-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               altNodeNames = 
                                 nodeNames = modProbes[top], 
                               nodeAttr = moduleColors[inModule][top])  
# 7. KEGG GO -------------------------------------------------------------------
## 7.1 Set output catagory----
#OCI_AML2/MV4_11/MOLM13
# 指定文件夹路径
dir.create("./03_Result/GO&KEGG/MOLM13/VEN_vs_WT/")
dir_WGCNA <- "./03_Result/WGCNA/"

## 7.2 WGCNA_res input ----
WGCNA_gene <- read.csv('./')
protein_anno <- read_xlsx("./01_Data/data_anno.xlsx")
protein_anno <- as.data.frame(protein_anno)
rownames(protein_anno) <- protein_anno$Protein.Group

WGCNA_gene <- merge(WGCNA_gene,protein_anno,by="row.names",all=FALSE)

# 转换基因名 
y <- WGCNA_gene$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
WGCNA_gene$gene <- gene

## 7.3 Module genes ----
Module_genes <- subset(WGCNA_gene, WGCNA_gene$module == "turquoise")

if(T){
  # 设置数据库 
  GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
  KEGG_database <- 'hsa'         # KEGG是hsa数据库
  
  # gene ID转换 
  gene <- bitr(Module_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
  
  ## 7.4 GO ----
  # GO富集分析
  kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                  OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                  keyType = "ENTREZID", # 设定读取的gene ID类型
                  ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)# 设定q值阈值
  
  ## 7.5 KEGG ----
  # KEGG富集分析
  kk <- enrichKEGG(gene = gene$ENTREZID,
                   keyType = "kegg",
                   organism = KEGG_database,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
  
  ## GO、KEGG结果整合 
  result <- list(enrichGO = kkd, enrichKEGG = kk)
  
  # 结果标记为下调 
  result_down <- result
  GO_res <- result_down$enrichGO
  KEGG_res <- result_down$enrichKEGG
  
  ## 7.6  res_output ----
  # 导出enrichGO 
  write.csv(GO_res@result, file = paste0(dir_WGCNA, "/GO_res.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_WGCNA, "/GO.pdf"), width = 6, height = 7)
  p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
  print(p1)
  dev.off()
  
  # 导出enrichKEGG
  write.csv(KEGG_res@result, file = paste0(dir_WGCNA, "/KEGG_res.csv"), quote = F, row.names = F)
  
  # dotplot
  pdf(file = paste0(dir_WGCNA, "/KEGG.pdf"), width = 6, height = 5)
  p2 <- dotplot(KEGG_gene,showCategory = 10)
  print(p2)
  dev.off()
}

kegg <- kk_up@result[kk_up@result$pvalue<0.05,]
p2 <- ggplot(kegg,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color= -log10(pvalue)))+
  theme_bw()+labs(y="",x="GeneRatio")+ 
  scale_color_gradient(low="blue",high="red")+
  #scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # 调整Y轴标签角度

ggsave(plot = p2,filename = "./03_Result/GO&KEGG/MOLM13/6W_vs_WT/downkegg.pdf")

## 7.7 统计通路的数量 ----
# GO 
table(GO_res@result$p.adjust<0.05)
# KEGG 
table(KEGG_res@result$p.adjust<0.05)
