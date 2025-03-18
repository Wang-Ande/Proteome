extract_genes <- function(cm, cl, membership_threshold) {
  # 过滤出指定 cluster 的基因
  target_cl <- cm$wide.res %>%
    filter(cluster == cl) %>%
    pull(gene)
  
  # 如果需要筛选 membership，则进行筛选
  if (!is.null(membership_threshold)) {
    # 获取 target_cl 对应的 membership 值
    membership_values <- cm$wide.res %>%
      filter(cluster == cl) %>%
      pull(membership)
    
    # 选择 membership 大于阈值的基因
    target_cl <- target_cl[membership_values > membership_threshold]
  }
  
  return(target_cl)
}
