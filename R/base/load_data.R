# load data for analysis

load_data_R <- function(imaging_parc_path, gene_expression_path, nulls_path, roi_range=NULL) {
  imaging_data <- read.csv(imaging_parc_path)
  nulls_data <- read.csv(nulls_path)
  gene_exp_orig <- read.csv(gene_expression_path)
  
  if (!is.null(roi_range)) {
    imaging_data <- imaging_data[roi_range, ]
    nulls_data <- nulls_data[roi_range, ]
    gene_exp_orig <- gene_exp_orig[roi_range, ]
  }
  
  genelabels <- colnames(gene_exp_orig)[3:ncol(gene_exp_orig)]
  
  new_column <- imaging_data[, 3, drop = FALSE]
  colnames(new_column) <- colnames(imaging_data)[3]
  gene_exp_orig_for_analysis <- cbind(gene_exp_orig[, 1:2], new_column, gene_exp_orig[, 3:ncol(gene_exp_orig)])
  
  # z-score the data by column
  gene_exp_orig_for_analysis[, 3:ncol(gene_exp_orig_for_analysis)] <- scale(gene_exp_orig_for_analysis[, 3:ncol(gene_exp_orig_for_analysis)])
  nulls_data[, ] <- scale(nulls_data[, ])
  
  return(list(gene_exp_orig_for_analysis = gene_exp_orig_for_analysis, nulls_data = nulls_data, genelabels = genelabels))
}
