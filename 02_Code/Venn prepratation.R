# 一共有三处需要更改的地方，
# 1、提取子集处的命名
# 2、最后输出文件的命名
# 3、输出路径dir

# 1. Set workspace -------------------------------------------------------------
setwd("../AML_project/aml_analysis")

# 2. Library packages ----------------------------------------------------------
library(dplyr)
library(openxlsx)

# 3. Creat output category -----------------------------------------------------
dir.create("03_result/Veen/Transcripton&Proteome/")
dir <- "./03_Result/Venn/GSEA/"

# 4. Data input ----------------------------------------------------------------
data1 <- read.xlsx("./03_Result/GSEA/MOLM13/2W_VS_WT/KEGG_results_selected.xlsx")
data2 <- read.xlsx("./03_Result/GSEA/MOLM13/6W_VS_WT/KEGG_results_selected.xlsx")
data3 <- read.xlsx("./03_Result/GSEA/MV4_11/2W_VS_WT/KEGG_results_selectd.xlsx")
data4 <- read.xlsx("./03_Result/GSEA/MV4_11/2W_VS_WT/KEGG_results_selectd.xlsx")
data5 <- read.xlsx("./03_Result/GSEA/OCI_M2/2W_VS_WT/KEGG_results_selected.xlsx")
data6 <- read.xlsx("./03_Result/GSEA/OCI_M2/6W_VS_WT/KEGG_results_selected.xlsx")


MOLM13_up <- data1[data1$NES > 0,c(2:3)]
MOLM13_down <- data1[data1$NES < 0,c(2:3)]
MOLM13_6w_up <- data2[data2$NES > 0,c(2:3)]
MOLM13_6w_down <- data2[data2$NES < 0,c(2:3)]

MV4_11_2w_up <- data3[data3$NES > 0,c(2:3)]
MV4_11_6w_up <- data4[data4$NES > 0,c(2:3)]
MV4_11_2w_down <- data3[data3$NES < 0,c(2:3)]
MV4_11_6w_down <- data4[data3$NES < 0,c(2:3)]

OCI_2w_up <- data5[data5$NES > 0,c(2:3)]
OCI_6w_up <- data6[data6$NES > 0,c(2:3)]
OCI_2w_down <- data5[data5$NES < 0,c(2:3)]
OCI_6w_down <- data6[data6$NES < 0,c(2:3)]

MOLM13_up_all <- bind_rows(MOLM13_up,MOLM13_6w_up)
MOLM13_down_all <- bind_rows(MOLM13_down,MOLM13_6w_down)
colnames(OCI_down_all)[2]<- "Group"
OCI_down_all[,2] <- "OCI_AML2"

MV4_11_up_all <- bind_rows(MV4_11_2w_up, MV4_11_6w_up)
MV4_11_down_all <- bind_rows(MV4_11_2w_down, MV4_11_6w_down)

OCI_up_all <- bind_rows(OCI_2w_up,OCI_6w_up)
OCI_down_all <- bind_rows(OCI_2w_down,OCI_6w_down)

All_up <- bind_rows(MOLM13_up_all, MV4_11_up_all, OCI_up_all)
All_down <- bind_rows(MOLM13_down_all, MV4_11_down_all, OCI_down_all)

write.xlsx(All_up, "./03_Result/Venn/GSEA/All_up.xlsx")
write.xlsx(All_down, "./03_Result/Venn/GSEA/All_down.xlsx")

if(T){
  # 5. Select gene names（symbol） -----------------------------------------------
  #Group1
  rownames(data1)<-data1$Row.names
  y <- rownames(data1)
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
  id<-data.frame(gene)
  data1<-cbind(id,data1)
  #Group2
  rownames(data2)<-data2$Row.names
  y <- rownames(data2)
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
  id<-data.frame(gene)
  data2<-cbind(id,data2)
  #Group3
  rownames(data3)<-data3$Row.names
  y <- rownames(data3)
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
  id<-data.frame(gene)
  data3<-cbind(id,data3)
  #Group4
  y <- data4$Genes
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
  id<-data.frame(gene)
  data4<-cbind(id,data4)
  #Group5
  y <- data5$Genes
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
  id<-data.frame(gene)
  data5<-cbind(id,data5)
  #Group6
  y <- data6$Genes
  gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
  id<-data.frame(gene)
  data6<-cbind(id,data6)
}

# 加Sig列
# logfc和pvalue需要更改！！！！！！！！！！！！！！！！！！！！！！！！！！！！
if(T){
  # group1
  logFC = 1
  P.Value = 0.05
  k1 <- (data1$P.Value < P.Value)&(data1$logFC < -logFC) 
  k2 <- (data1$P.Value < P.Value)&(data1$logFC > logFC)
  data1 <- mutate(data1,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
  # group2
  logFC = 1
  P.Value = 0.05
  k1 <- (data2$P.Value < P.Value)&(data2$logFC < -logFC) 
  k2 <- (data2$P.Value < P.Value)&(data2$logFC > logFC)
  data2 <- mutate(data2,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
  # group3
  logFC = 1
  P.Value = 0.05
  k1 <- (data3$P.Value < P.Value)&(data3$logFC < -logFC) 
  k2 <- (data3$P.Value < P.Value)&(data3$logFC > logFC)
  data3 <- mutate(data3,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
  # group4
  logFC = 0.585
  P.Value = 0.05
  k1 <- (data4$P.Value < P.Value)&(data4$logFC < -logFC) 
  k2 <- (data4$P.Value < P.Value)&(data4$logFC > logFC)
  data4 <- mutate(data4,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
  # group5
  logFC = 0.585
  P.Value = 0.05
  k1 <- (data5$P.Value < P.Value)&(data5$logFC < -logFC) 
  k2 <- (data5$P.Value < P.Value)&(data5$logFC > logFC)
  data5 <- mutate(data5,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
  # group6
  logFC = 0.585
  P.Value = 0.05
  k1 <- (data6$P.Value < P.Value)&(data6$logFC < -logFC) 
  k2 <- (data6$P.Value < P.Value)&(data6$logFC > logFC)
  data6 <- mutate(data6,Sig = ifelse(k1,"down",
                                     ifelse(k2,"up","stable")))
}
table(data1$Sig)
table(data4$Sig)

# 6. Extract subset ------------------------------------------------------------
# change名需要更改 ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
if(T){
  #Group1
  up_geneSet1 <- data1[data1$Sig == "up",c(1,ncol(data1))]
  down_geneSet1 <- data1[data1$Sig == "down",c(1,ncol(data1))]
  up_geneSet1$Sig <- "2W_UP"
  down_geneSet1$Sig <- "2W_DOWN"
  #Group2
  up_geneSet2 <- data2[data2$Sig == "up",c(1,ncol(data2))]
  down_geneSet2 <- data2[data2$Sig == "down",c(1,ncol(data2))]
  up_geneSet2$Sig <- "6W_UP"
  down_geneSet2$Sig <- "6W_DOWN"
  #Group3
  up_geneSet3 <- data3[data3$Sig == "up",c(1,ncol(data3))]
  down_geneSet3 <- data3[data3$Sig == "down",c(1,ncol(data3))]
  up_geneSet3$Sig <- "Transcription_VEN_vs_WT"
  down_geneSet3$Sig <- "Transcription_VEN_vs_WT" 
  #Group4
  up_geneSet4 <- data4[data4$Sig == "up",c(1,ncol(data4))]
  down_geneSet4 <- data4[data4$Sig == "down",c(1,ncol(data4))]
  up_geneSet4$Sig <- "2W_UP"
  down_geneSet4$Sig <- "2W_DOWN"
  #Group5
  up_geneSet5 <- data5[data5$Sig == "up",c(1,ncol(data5))]
  down_geneSet5 <- data5[data5$Sig == "down",c(1,ncol(data5))]
  up_geneSet5$Sig <- "6W_UP"
  down_geneSet5$Sig <- "6W_DOWN"
  #Group6
  up_geneSet6 <- data6[data6$Sig == "up",c(1,ncol(data6))]
  down_geneSet6 <- data6[data6$Sig == "down",c(1,ncol(data6))]
  up_geneSet6$Sig <- "Proteome_VEN_vs_WT"
  down_geneSet6$Sig <- "Proteome_VEN_vs_WT" 
}

# 7. Merge subset --------------------------------------------------------------
if(T){
  # 上调基因
  All_Up_geneSet <- bind_rows(up_geneSet4, up_geneSet5)
  All_Up_geneSet <- All_Up_geneSet %>% rename(Group=Sig)
  # 下调基因
  All_Down_geneSet <- bind_rows(down_geneSet4,down_geneSet5)
  All_Down_geneSet <- All_Down_geneSet %>% rename(Group=Sig)
  
  # 8. Result output -------------------------------------------------------------
  # 命名前方的前缀“1”表示差异基因为取cutoff=1的结果
  write.xlsx(All_Up_geneSet,file = paste0(dir,"All_UpGeneSet.xlsx"))
  write.xlsx(All_Down_geneSet,file = paste0(dir,"All_DownGeneSet.xlsx"))
}
table(All_Up_geneSet$Group)
table(All_Down_geneSet$Group)


