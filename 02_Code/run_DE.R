#' Run Differential Expression Analysis using limma
#'
#' This function performs differential expression (DE) analysis between two groups
#' using the `limma` package. It supports both unpaired and paired analysis,
#' includes log transformation, volcano plot generation, and optional annotation merging.
#'
#' @param data A numeric expression matrix with rows as features and columns as samples.
#' @param data_group A data.frame with columns `id` (sample IDs) and `group` (group labels).
#' @param data_anno Optional. A data.frame containing annotation for each feature. Row names must match `data`.
#' @param log2 Logical. Whether to log2-transform the expression matrix.
#' @param group_1 Character. The first group to compare (e.g., "WT").
#' @param group_2 Character. The second group to compare (e.g., "KO").
#' @param logfc_threshold Numeric. Threshold for absolute log2 fold-change.
#' @param pvalue_threshold Numeric. P-value threshold for significance.
#' @param qvalue_threshold Optional. Not currently used. Leave as NULL.
#' @param dir Output directory to save results. Defaults to current working directory.
#' @param paired Logical. Whether the design is paired.
#' @param pair_col Character. Name of the column in `data_group` representing pairing ID (e.g., patient ID).
#'
#' @return A list containing:
#' \describe{
#'   \item{result_merge}{Merged expression matrix with DE results and annotation (if provided).}
#'   \item{DE_res}{A data.frame with DE results including logFC, P.Value, and significance labels.}
#' }
#' 
#' @import limma
#' @import ggplot2
#' @import dplyr
#' @importFrom utils write.csv
#' @importFrom stats model.matrix
#' 
#' @examples
#' run_DE(data = expr_matrix,
#'        data_group = group_info,
#'        group_1 = "Control",
#'        group_2 = "Treated",
#'        log2 = TRUE,
#'        logfc_threshold = 1,
#'        pvalue_threshold = 0.05,
#'        paired = TRUE,
#'        pair_col = "patient_id")
#'
#' @export

run_DE <- function(data, data_group, data_anno=NULL, log2,
                   group_1, group_2, logfc_threshold, pvalue_threshold, 
                   qvalue_threshold = NULL, paired = FALSE, pair_col = NULL, dir = getwd()) {
  library(limma)
  library(ggplot2)
  library(dplyr)
  
  if (!is.null(data_anno)) {
    if (!all(rownames(data) == rownames(data_anno))) {
      stop("Error: The row names of 'data' and 'data_anno' do not match.")
    }
  }
  
  output_dir <- paste0(dir, group_2, "_vs_", group_1)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  limma_group <- subset(data_group, data_group$group %in% c(group_1, group_2))
  group <- factor(limma_group$group, levels = c(group_1, group_2))
  
  targeted_expr <- data[, colnames(data) %in% limma_group$id]
  if(log2){
    targeted_expr <- log2(targeted_expr)
  }
  
  if (paired) {
    if (is.null(pair_col) || !pair_col %in% colnames(limma_group)) {
      stop("Please provide a valid 'pair_col' when 'paired = TRUE'")
    }
    pair_factor <- factor(limma_group[[pair_col]])
    design <- model.matrix(~pair_factor + group)
    fit <- lmFit(targeted_expr, design)
    fit2 <- eBayes(fit)
    coef_name <- paste0("group", group_2)
  } else {
    design <- model.matrix(~-1 + group)
    colnames(design) <- paste0("group", levels(group))
    rownames(design) <- limma_group$id
    contrast.matrix <- makeContrasts(contrasts = paste0("group", group_2, " - group", group_1),
                                     levels = design)
    fit <- lmFit(targeted_expr, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1, 0.01)
    coef_name <- paste0("group", group_2, " - group", group_1)
  }
  
  tempOutput <- topTable(fit2,
                         coef = coef_name,
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
  write.csv(result_merge, file = paste0(output_dir,"/result_DE.csv"))
  
  # 火山图 ----
  DE_res$Sig <- factor(DE_res$Sig, levels = c("up", "down", "stable"))
  max_abs_logfc <- max(abs(DE_res$logFC), na.rm = TRUE)
  x_limit <- max_abs_logfc * 1.05
  
  # plot ----
  p <- ggplot(data = DE_res, 
              aes(x = logFC, 
                  y = -log10(P.Value))) +
    geom_point(alpha = 0.5, size = 3, 
               aes(color = Sig)) +
    ylab("-log10(Pvalue)")+
    scale_color_manual(
      name = "Change",
      values = c("up" = "#B30000", "stable" = "grey", "down" = "#003366"),
      labels = c("up" = paste0("Up ：", sum(DE_res$Sig == "up")),
                 "stable" = paste0("Stable ：", sum(DE_res$Sig == "stable")),
                 "down" = paste0("Down ：", sum(DE_res$Sig == "down"))))+
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
               col = "black", lwd = 0.8, alpha = 0.4) +
    geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
               col = "black", lwd = 0.8, alpha = 0.4) +
    labs(title = paste0(group_2,"-",group_1)) +
    xlim(-x_limit, x_limit)+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  
  ggsave(filename = paste0(output_dir,"/volc.pdf"),
         plot = p, device = "pdf", 
         width = 6, height = 5)
  
  return(list(result_merge = result_merge, DE_res = DE_res))
}
