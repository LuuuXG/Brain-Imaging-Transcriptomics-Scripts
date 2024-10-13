# threshold: filter genes from pls_component

gsea_GO_from_pls <- function(pls_result_path, ont='ALL', threshold=1, output_dir, simplify=FALSE, simplify_cutoff=NULL){
  
  #---- 1. 导入数据，生成genelist
  genelist_input <- read.table(pls_result_path, header = T, sep = ',')
  # get the suffix of pls_result_path
  pls_result_path_suffix <- sub("^.*pls_component_([^.]*).*", "\\1", pls_result_path)
  
  genelist_input <- subset(genelist_input, `p.value` < threshold)
  
  genename <- as.character(genelist_input[,1])   #提取第一列基因名
  
  duplicate_genes <- genelist_input[duplicated(genelist_input$Gene), "Gene"]
  genelist_input <- genelist_input[!duplicated(genelist_input$Gene), ]
  
  gene_entrezid <- bitr(geneID = genename, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", # 转成ENTREZID
                        OrgDb = "org.Hs.eg.db"
  )
  
  # 去重复，前面是1:1 mapping的话应该没有影响
  non_duplicates_idx <- which(duplicated(gene_entrezid$SYMBOL) == FALSE)
  gene_entrezid <- gene_entrezid[non_duplicates_idx, ]
  
  colnames(gene_entrezid)[1]<-"Gene"
  
  gene_entrezid <- inner_join(gene_entrezid, genelist_input, by = "Gene")
  
  genelist <- gene_entrezid$Z.score
  names(genelist) = gene_entrezid$ENTREZID
  
  #---- 2. gseGO()
  Go_gseresult <- gseGO(genelist, 'org.Hs.eg.db', keyType = "ENTREZID", ont=ont, pvalueCutoff=0.05)
  # 细胞组成（cellular component，CC）
  # 生物过程（biological process，BP）
  # 分子功能（Molecular Function，MF）
  
  # 如果设置了simplify_cutoff但是没有设置simplify，则打印说明simplify_cutoff无效
  if (!simplify && !is.null(simplify_cutoff)){
    cat(warn("`simplify_cutoff` is set but `simplify` is not set. `simplify_cutoff` will be ignored.\n"))
  }
  
  if (simplify){
    Go_gseresult <- simplify_by_NES(Go_gseresult, cutoff=simplify_cutoff)
  }
  
  #保存富集分析结果
  go_results<-as.data.frame(Go_gseresult)
  
  if (simplify){
    gsea_go_output_path <- paste0(output_dir, "/Go_gseresult_", pls_result_path_suffix, "_ont=", ont, "_simplified.csv")
  } else {
    gsea_go_output_path <- paste0(output_dir, "/Go_gseresult_", pls_result_path_suffix, "_ont=", ont, ".csv")
  }

  write.csv (go_results, file = gsea_go_output_path)
  
  cat(blue$bold("GSEA_GO result has been saved to ", gsea_go_output_path, "\n"))
  
  return(Go_gseresult)
}