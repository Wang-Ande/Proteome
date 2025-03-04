run_DE <- function(data, data_group, data_anno=NULL,log2,
                   group_1,group_2,logfc_threshold, pvalue_threshold, 
                   qvalue_threshold = NULL,dir = getwd()) {
  library(limma)
  # check data
  # 设置输出目录，并生成目录
  # 检查行名是否一致
  if (!is.null(data_anno)){
    if (!all(rownames(data) == rownames(data_anno))) {
      stop("Error: The row names of 'data' and 'data_anno' do not match.")
    }
  }
  output_dir <- paste0(dir, group_1, "_vs_", group_2)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # 筛选对应分组
  limma_group <- subset(data_group,data_group$group %in% c(group_1,group_2))
  group <- factor(limma_group$group,levels = c(group_1,group_2))# "D","WT"

  limma_design <- model.matrix(~-1+group)

  rownames(limma_design) <- limma_group$id

  targeted_expr <- data[,colnames(data) %in% limma_group$id]
  
  if(log2){
    targeted_expr <- log2(targeted_expr)
  }

  contrast.matrix <- makeContrasts(contrasts = paste0("group",levels(group)[1],"-","group",levels(group)[2]),
                                   levels = limma_design)

  fit <- lmFit(targeted_expr,limma_design)

  fit1 <- contrasts.fit(fit, contrast.matrix)

  fit2 <- eBayes(fit1,0.01)

  tempOutput = topTable(fit2,
                        coef = paste0("group",levels(group)[1],"-","group",levels(group)[2]),
                        number = nrow(fit2),
                        lfc = log2(1),
                        adjust.method = "fdr")
  DE_res <- na.omit(tempOutput)
  
  # 加change列 ----
  # 标记上下调基因，可根据需求设定阈值
  logFC = logfc_threshold
  P.Value = pvalue_threshold
  k1 <- (DE_res$P.Value < P.Value) & (DE_res$logFC < -logFC)
  k2 <- (DE_res$P.Value < P.Value) & (DE_res$logFC > logFC)
  DE_res <- mutate(DE_res, 
                   Sig = ifelse(k1, "down", 
                                ifelse(k2, "up", "stable")))
  
  result_merge <- merge(2^targeted_expr,
                        DE_res,
                        by = "row.names")
  rownames(result_merge) <- result_merge$Row.names
  
  if (!is.null(data_anno)){
    result_merge <- merge(data_anno,
                          result_merge,
                          by = "row.names")
    }
  write.csv(result_merge,file = paste0(output_dir,"/result_DE.csv"))
  
  # 火山图 ----
  # change列因子化 ----
  DE_res$Sig <- factor(
    DE_res$Sig,
    levels = c("up", "down", "stable"))  # 强制按此顺序排列
  
  # 计算对称的横坐标范围（基于logFC绝对值的最大值）
  max_abs_logfc <- max(abs(DE_res$logFC), na.rm = TRUE)
  
  # 扩展5%的余量，避免点紧贴坐标轴边缘
  x_limit <- max_abs_logfc * 1.05
  
  # plot ----
  p <- ggplot(data = DE_res, 
              aes(x = logFC, 
                  y = -log10(P.Value))) +
    geom_point(alpha = 0.5, size = 3, 
               aes(color = Sig)) +
    ylab("-log10(Pvalue)")+
    # 按因子顺序指定颜色
    scale_color_manual(
      name = "Change",                              # 图例标题
      values = c(
        "up" = "#B30000",      # 红
        "stable" = "grey",     # 灰
        "down" = "#003366"),    # 蓝,
      labels = c(                                  # 显示上下调的个数
        "up" = paste0("Up ：", sum(DE_res$Sig == "up")),
        "stable" = paste0("Stable ：", sum(DE_res$Sig == "stable")),
        "down" = paste0("Down ：", sum(DE_res$Sig == "down"))))+
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
               col = "black", lwd = 0.8,alpha=0.4) +
    geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
               col = "black", lwd = 0.8,alpha=0.4) +
    labs(title = paste0(group_1,"-",group_2)) +
    xlim(-x_limit, x_limit)+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),  # 标题居中
          aspect.ratio = 1)  # 设置纵横比，调整为更高
  
  ggsave(filename = paste0(output_dir,"/volc.pdf"),
         plot = p, device = "pdf", 
         width = 6, height = 5)
  
  
  # 热图 
  return(list(result_merge = result_merge, DE_res = DE_res))
}
