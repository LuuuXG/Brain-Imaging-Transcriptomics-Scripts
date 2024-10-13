# threshold: filter genes from pls_component

gsea_user_defined_from_pls <- function(pls_result_path, geneset, threshold=1, output_dir){
  
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
  
  #---- 2B. GSEA()
  gsea_res <- GSEA(genelist, 
                   TERM2GENE = geneset,
                   minGSSize = 5,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH"
  )
  
  result <- as.data.frame(gsea_res)
  
  gsea_msigdb_output_path <- paste0(output_dir, "/User_Defined_GeneSet_gseresult_", pls_result_path_suffix, ".csv")
  write.csv (result, file = gsea_msigdb_output_path)
  
  return(gsea_res)
}