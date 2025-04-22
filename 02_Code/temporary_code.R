# library(pathview)
# library(clusterprofiler)
kegg <- read.csv("./03_Result/GO&KEGG/OCI_AML2/Low_vs_Con/KEGG_down.csv")
gene_ids <- "4697/6392/9167/4706/27089/4724/4711/4709/4717/4725/4720/9550/4707/4710/4714/4701/4705/1345/1347/1350/4713/4729/7385/1349/515/4719/529/7386/4723/479/514/4700/10312/4704/9296/1355/535/64077/55967/4715/1353"
gene_list <- unlist(strsplit(gene_ids, "/"))
gene <- clusterProfiler::bitr(gene_list, fromType = 'ENTREZID', 
                              toType = 'SYMBOL', OrgDb = org.Hs.eg.db)


# logFC
DE_gene <- read.csv("./03_Result/DEP/OCI_AML2/Low_vs_Ctrl/result_DE.csv",row.names = 1)
# Gene_anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)

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
         pathway.id = "00190",       # 通路ID
         species = "hsa",            # 人类
         gene.idtype = "entrez",     # ID类型
         limit = list(gene=1, cpd=1),
         bins = list(gene=10, cpd=10),
         low = list(gene="#93A5CB", cpd="#66F1A9"),      # 深蓝（下调）
         mid = list(gene="#F7F7F7", cpd="#F0F0F0"),      # 白灰色（中性）
         high = list(gene="#C6133B", cpd="#FC4E2A")     # 深红（上调）
         )    


library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
expr <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv",row.names = 1)
anno <- read.xlsx("./01_Data/data_anno.xlsx",rowNames = TRUE)
expr <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)
targeted_gene <- expr[grep("TP53",expr$Genes),]
expr_targeted <- targeted_gene[2,grep("MV4_11",colnames(targeted_gene))]
expr_targeted <- t(log2(expr_targeted))
colnames(expr_targeted)[1] <- "expr"
expr_targeted <- as.data.frame(expr_targeted)

group <- read.xlsx("./01_Data/IC50_group.xlsx")
group_targeted <- group[grep("MV4_11",group$id),c(1,2)]
expr_targeted$group <- group_targeted$group

# 定义变量
gene <- "TP53"          # 目标基因
mycol <- c("#6388B4","#FFAE34","#EF6F6A")  # 配色

# 计算样本数（假设 group1 是分组列）
dt <- expr_targeted
num_normal <- sum(dt$group == "Ctrl")
num_low_resistance <- sum(dt$group == "Low")
num_high_resistance <- sum(dt$group == "High")

# 确保分组列是因子，并指定顺序
dt$group <- factor(dt$group, levels = c("Ctrl", "Low", "High"))
rownames(dt)
dt$pair <- c(1,1,1,2,3,2,3,2,3)
table(dt$pair)

dir <- "./03_Result/DEP/MV4_11/"

# 绘图
p1 <- ggplot(data = dt,
             aes(x = group, 
                 y = expr,  # 假设表达量列名为 expression
                 color = group, 
                 fill = group)) +
  geom_line(aes(group = pair), color = "gray50", linewidth = 0.5, alpha = 0.6) +  # 加入配对线
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +  # 隐藏箱线图异常值
  geom_jitter(width = 0.2, size = 2.8, shape = 21, 
              stroke = 0.8,aes(fill = group), 
              alpha = 0.7) +
  theme_classic() +  scale_color_manual(name = "Group", values = mycol) +
  scale_fill_manual(name = "Group" ,values = mycol) +
  labs(
    title = paste("Expression of", gene, "across three states"),
    x = "MV4_11 Group",
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
p1 <- p1 + stat_compare_means(comparisons = list(c("Ctrl", "Low"),
                                                 c("Low", "High"),
                                                 c("Ctrl", "High")),
                              method = "t.test",  # 非参数检验 wilcox.test, 参数检验 t.test 
                              label = "p.signif")  # 显示显著性标记（***, ​**, *）

ggsave(filename = paste0(dir,"/TP53_expr_boxplot.pdf"),
       plot = p1, device = cairo_pdf, 
       width = 6, height = 5)
  
