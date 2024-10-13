## load data for WGCNA analysis ##
# the main difference between load_data_raw.R and load_data.py is that the
# former does not perform z-score transform, and the later is designed for
# PLS+GSEA in python.

load_data_raw <- function(imaging_parc_path, gene_expression_path, 
                            nulls_path, roi_range=NA){
  data_list <- list()
  
  # load gene expression data
  gene_expression_data_raw <- read.csv(gene_expression_path, sep = ',')
  gene_expression_data_raw <- gene_expression_data_raw[, -1]
  rownames(gene_expression_data_raw) <- gene_expression_data_raw[, 1]
  gene_expression_data_raw <- gene_expression_data_raw[, -1]
  
  # load imaging data
  imaging_data_raw_temp <- read.csv(imaging_parc_path)
  imaging_data_raw <- imaging_data_raw_temp[, c(-1, -2)]
  imaging_data_raw <- as.data.frame(imaging_data_raw)
  rownames(imaging_data_raw) <- imaging_data_raw_temp[,2]
  colnames(imaging_data_raw) <- colnames(imaging_data_raw_temp)[3]
  
  # load spatial nulls data
  nulls_data_raw <- read.csv(nulls_path)
  rownames(nulls_data_raw) <- imaging_data_raw_temp[,2]
  
  if (is.na(roi_range[1]) == FALSE){
    imaging_data_raw <- imaging_data_raw[roi_range, , drop = FALSE]
    gene_expression_data_raw <- gene_expression_data_raw[roi_range,]
    nulls_data_raw <- nulls_data_raw[roi_range,]
  }
  
  data_list$imaging_data_raw <- imaging_data_raw
  data_list$gene_expression_data_raw <- gene_expression_data_raw
  data_list$nulls_data_raw <- nulls_data_raw
  
  return(data_list)
}