# 1. Packages ----
library(FactoMineR)
library(factoextra)
library(openxlsx)
library(readr)
library(readxl)
library(edgeR)
library(corrplot)

# 2. Data input ----
## 2.1 Group input ----
data_group <- read_excel("./01_Data/IC50_group.xlsx")
data_group <- as.data.frame(data_group)
data_group <- data_group[-23,]
data_group <- data_group[-grep("MV4_11_6W_1",data_group$id),-c(3,4)]
colnames(data_group)[2] <- "group"
table(data_group$group)
data_group <- data_group[grep("OCI", data_group$id),]
rownames(data_group) <- data_group$id


## 2.2 Expr input ----
exprSet_raw <- read.csv("./01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
exprSet_raw <- as.data.frame(exprSet_raw)
exprSet_raw <- exprSet_raw[,data_group$id]

## 2.3 PCA analysis ----
# log2 转换
expr_log <- log2(exprSet_raw + 1)

# Z-score 标准化（按列，即按基因）
expr_scaled <- t(scale(t(expr_log)))
# 这样得到的矩阵：每个基因均值=0，方差=1

# PCA 输入
expr_t <- t(expr_scaled)

# 确认分组顺序和样本对应
all(rownames(expr_t) == data_group$id)  # 应该返回 TRUE

# PCA
pca_res <- PCA(expr_t, graph = FALSE)

## 2.4 PCA plot ----
# 配色设置
value_colour <- c("High" = "#E64B35FF",
                  "Ctrl" = "#3C5488FF",
                  "Low" = "#F2A200")
# 控制图例顺序
data_group$group <- factor(data_group$group, levels = c("Ctrl", "Low", "High"))

dir <- "./03_Result/QC/OCI_AML2/"
pdf(file = paste0(dir,"QC_pca_1.pdf"),
    width = 5.5,
    height = 5)
fviz_pca_ind(pca_res,
             pointsize = 2,              # 点的大小
             mean.point = FALSE,           # 去除分组的中心点
             label = "none",               # 隐藏每个样本标签
             habillage = data_group$group, # 根据样本类型着色
             palette = value_colour,       # 三个组的颜色
             addEllipses = TRUE,           # 添加分组椭圆
             legend.title = "Groups",      # 图例标题
             ellipse.level = 0.95          # 椭圆置信区间（默认95%）
)+
  theme_classic()+
  # 修改图例
  theme(legend.title = element_text(size = 14, face = "plain"), # 图例标题大小
        legend.text  = element_text(size = 12),                # 图例文字大小
        legend.position = "right")+           
  guides(color = guide_legend(override.aes = list(size=5))) + # 放大图例点
  theme(legend.text = element_text(size = 12))
  
dev.off()

