library(ggplot2)
library(data.table)
library(pheatmap)

#### Gene expression heatmap ####
source('./base/load_data.R')

load(file=file.path(output_root_folder, 'gene_exp+img_parc+nulls.RData'))

data_for_plot_list <- load_data_R(imaging_parc_path, gene_expression_path, nulls_path, roi_range=1:83)

data_for_plot <- data_for_plot_list$gene_exp_orig_for_analysis

data_for_plot <- data_for_plot[,c(-3, -1)]

rownames(data_for_plot) <- data_for_plot[,1]

data_for_plot <- data_for_plot[,-1]

heatmap_plot <- pheatmap(
  data_for_plot,          # 数据框
  cluster_rows = FALSE,   # 不聚类行
  cluster_cols = FALSE,   # 不聚类列
  scale = "none",         # 保持原始表达值
  color = colorRampPalette(c("blue", "white", "red"))(100), # 颜色梯度从蓝到红
  show_rownames = F,      # 去掉纵轴标签
  show_colnames = F,      # 去掉横轴标签
  legend = F
)

gene_expression_folder <- file.path(output_root_folder, 'expression')

#保存为jpeg
jpeg(file=file.path(gene_expression_folder, 'gene_expression_heatmap.jpg'), width=800, height=800)
print(heatmap_plot)
dev.off()
