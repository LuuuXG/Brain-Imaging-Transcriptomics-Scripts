# clusterProfiler::simplify
# choose the represent geneset by the absolyte value of NES

simplify_by_NES <- function(gsea_res, cutoff=0.7) {
  gsea_res@result$NES_original <- gsea_res@result$NES
  
  gsea_res@result$NES <- abs(gsea_res@result$NES)
  
  gsea_res_simplified <- clusterProfiler::simplify(gsea_res, 
                                                   by="NES", 
                                                   select_fun=max, 
                                                   cutoff=cutoff)
  
  gsea_res_simplified@result$NES <- gsea_res@result$NES_original[match(gsea_res_simplified@result$ID, gsea_res@result$ID)]
  
  gsea_res_simplified@result <- gsea_res_simplified@result[order(gsea_res_simplified@result$pvalue), ]
  
  gsea_res_simplified@result <- gsea_res_simplified@result %>%
    distinct(ONTOLOGY, core_enrichment, .keep_all = TRUE)
  
  cat(note("The GSEA results have been simplified by the ABSOLUTE NES value. The cutoff is set to", cutoff, ".\n"))
  
  return(gsea_res_simplified)
}
