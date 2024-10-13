# generate gene set data for GSEA analysis

fetch_cell_specific_data <- function(thresholds){
  load('../data/04_pSI_output/pSI_output.RData')
  
  final_gene_set_data <- data.frame(gs_name = character(), 
                                    ENTREZID = character(), 
                                    stringsAsFactors = FALSE)
  
  if (thresholds == 0.001){
    index <- 1
  } else if (thresholds == 0.01){
    index <- 2
  } else if (thresholds == 0.05){
    index <- 3
  } else if (thresholds == 0.0001){
    index <- 4
  }
  
  
  for (i in 1:length(gse_number)){
    path <- cell_specific_data_list[[i]][[2]][[index]]
    gse_id <- gse_number[i]
    
    pSI_data <- read.csv(path, stringsAsFactors = FALSE, header = F)
    
    gene_set_data <- pSI_data %>%
      pivot_longer(cols = -V1, names_to = NULL, values_to = "SYMBOL") %>%
      rename(gs_name = V1) %>%
      filter(!is.na(SYMBOL) & SYMBOL != "")
    
    gene_set_data <- gene_set_data %>%
      left_join(bitr(gene_set_data$SYMBOL, fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                by = c("SYMBOL" = "SYMBOL"))
    
    gene_set_data <- gene_set_data %>%
      filter(!is.na(ENTREZID)) %>%
      dplyr::select(gs_name, ENTREZID) %>%
      mutate(gs_name = paste0(gse_id, "_", gs_name))
    
    final_gene_set_data <- bind_rows(final_gene_set_data, gene_set_data)
  }
  
  return(final_gene_set_data)
}
