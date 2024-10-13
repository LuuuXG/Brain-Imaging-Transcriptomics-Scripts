load(file.path(output_root_folder, 'gsea_cell_specific_results_data_pSI_threshold=0.001.RData'))
gsea_output_folder <- file.path(output_root_folder, 'gsea_output')

gsea_cell_specific_results_for_plot <- gsea_cell_specific_results_data[[1]]

result_df <- gsea_cell_specific_results_for_plot@result

# data preparation

cell_of_interest <- c('Astrocyte',
                      'Neuron',
                      'Oligodendrocyte',
                      'Endothelial',
                      'Microglia')

# 筛选result_df中Description列名称包含cell_of_interest中元素的行
result_df_filter <- result_df[grepl(paste(cell_of_interest, collapse='|'), 
                             result_df$Description),]

# Split the 'Description' column into 'GSE_ID' and 'Cell_Type'
result_df_filter$GSE_ID <- sapply(strsplit(result_df_filter$Description, "_"), `[`, 1)
result_df_filter$Cell_Type <- sapply(strsplit(result_df_filter$Description, "_"), `[`, 2)

# Select only the columns GSE_ID, Cell_Type, NES, and p.adjust
result_df_extracted <- result_df_filter[, c("GSE_ID", "Cell_Type", "NES", "p.adjust")]

# 创建一个占位数据行
placeholder_row <- data.frame(
  GSE_ID = "GSE73721",        # GSE编号
  Cell_Type = "Microglia",    # 细胞类型
  NES = NA,                   # 设定为 NA 或占位数值
  p.adjust = NA               # 设定为 NA 或占位数值
)

# 将占位数据行添加到现有数据框中
result_df_extracted <- rbind(result_df_extracted, placeholder_row)
rownames(result_df_extracted)[nrow(result_df_extracted)] <- "GSE73721_Microglia_0.001"

result_df_extracted$log_p_adjust <- -log10(result_df_extracted$p.adjust)

threshold1 <- -log10(0.05)
threshold2 <- -log10(0.01)

cell_specific_plot <- ggplot(result_df_extracted, aes(x = Cell_Type, y = GSE_ID, size = log_p_adjust, fill = NES)) +
  geom_point(
    aes(color = ifelse(log_p_adjust > threshold2, "#000000",
                       ifelse(log_p_adjust > threshold1, "#606060", "#C0C0C0")),
        stroke = ifelse(log_p_adjust > threshold2, 2,
                        ifelse(log_p_adjust > threshold1, 1.25, 0.1)),
        alpha = ifelse(log_p_adjust > threshold2, 1,
                       ifelse(log_p_adjust > threshold1, 0.8, 0.6))),
    shape = 21
  ) +
  scale_size_continuous(name = expression(-log[10] * "(" * italic(p) * "-adjust)"), range = c(1, 20), breaks = c(2, 4, 5)) +
  scale_fill_gradientn(name = "NES", colors = c('#33FFFF', '#66FFFF', '#99FFFF', '#CCFFFF', 'white', '#CCFFCC', '#99FF99', '#66FF66', '#33FF33')) +
  scale_color_identity() +  # 让自定义的颜色不显示图例
  scale_alpha_identity(guide = "none") +  # 关闭透明度图例
  #labs(x = "Cell Type") +
  theme_minimal() +
  scale_x_discrete(limits = c("Endothelial", "Neuron", "Oligodendrocyte", "Astrocyte", "Microglia")) +  # 设置横轴的顺序
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
    axis.text.y = element_text(size = 13, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    legend.position = "right",  # 控制图例整体位置
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm")
  ) +
  guides(
    size = guide_legend(order = 2),  # 设定 size 图例的顺序
    fill = guide_colorbar(order = 1)  # 设定 fill 图例的顺序
  )

jpeg(file.path(gsea_output_folder, 'GSEA_cell_specific_results.jpg'), width = 2600, height = 1600, res=300)
print(cell_specific_plot)
dev.off()

