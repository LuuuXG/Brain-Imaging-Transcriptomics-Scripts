load(file=file.path(output_root_folder, 'gene_exp+img_parc+nulls.RData'))

data_for_analysis_list <- load_data_R(imaging_parc_path, gene_expression_path, nulls_path, roi_range=1:83)

bootstrap_linear_regression <- function(target_gene_name, gene_exp_orig_for_analysis, nulls_data, n_perm=1000){
  # check if the target gene is in the gene expression data
  if (!(target_gene_name %in% colnames(gene_exp_orig_for_analysis))) {
    stop(paste("The target gene", target_gene_name, "is not in the gene expression data."))
  } else {
    cat(note("The target gene", target_gene_name, "is in the gene expression data."))
  }
  
  # Step 1: calculate the beta value of the target gene in the original data
  formula <- as.formula(paste(colnames(gene_exp_orig_for_analysis)[3], "~", target_gene_name))
  lm_orig <- lm(data = gene_exp_orig_for_analysis[, c(3, which(colnames(gene_exp_orig_for_analysis) == target_gene_name))], 
                formula = formula)
  beta_orig <- coef(lm_orig)[2]
  
  col_name1 <- colnames(gene_exp_orig_for_analysis)[3]
  col_name2 <- target_gene_name
  
  orig_data <- data.frame(
    imaging_value = gene_exp_orig_for_analysis[, 3],
    gene = gene_exp_orig_for_analysis[[target_gene_name]]
  )
  
  orig_data <- setNames(orig_data, c(col_name1, col_name2))
  
  # Step 2: calculate the beta value of the target gene in the nulls data
  beta_perm <- numeric(ncol(nulls_data))
  count <- 0
  
  for (i in 1:n_perm) {
    perm_data <- data.frame(
      imaging_value = nulls_data[, i],
      gene = gene_exp_orig_for_analysis[[target_gene_name]]
    )
    
    lm_perm <- lm(imaging_value ~ gene, data = perm_data)
    beta_perm[i] <- coef(lm_perm)[2]
    
    if ((beta_perm[i]*beta_orig > 0) && (abs(beta_perm[i]) > abs(beta_orig))){
      count <- count + 1
    }
  }
  
  # Step 3: calculate the z-statistic
  beta_sd_perm <- sd(beta_perm)
  z_stat <- beta_orig / beta_sd_perm
  
  # Step 4: calculate the p-value
  # use z-statistic to calculate p-value
  p_value_z <- 2 * (1 - pnorm(abs(z_stat)))
  
  # use permutation to calculate p-value
  p_value_perm <- count / n_perm
  
  return(list(orig_data = orig_data, target_gene_name = target_gene_name, beta_orig = beta_orig, beta_sd_perm = beta_sd_perm, z_stat = z_stat, p_value_z = p_value_z, p_value_perm = p_value_perm))
}

target_gene_name <- "WNT7A"

regression_result_list <- bootstrap_linear_regression(target_gene_name, data_for_analysis_list$gene_exp_orig_for_analysis, data_for_analysis_list$nulls_data, n_perm=1000)

data_for_plot <- regression_result_list$orig_data

y_col <- colnames(gene_exp_orig_for_analysis)[3]
x_col <- target_gene_name

x_max <- max(data_for_plot[[x_col]], na.rm = TRUE)
y_max <- max(data_for_plot[[y_col]], na.rm = TRUE)

x_min <- min(data_for_plot[[x_col]], na.rm = TRUE)
y_min <- min(data_for_plot[[y_col]], na.rm = TRUE)

x_distance <- x_max - x_min
y_distance <- y_max - y_min

p_value_color <- ifelse(regression_result_list$p_value_z < 0.05, "#FF3333", "black")

ggplot(data_for_plot, aes(x = !!sym(x_col), y = !!sym(y_col))) +
  geom_point(size = 1.75) +
  geom_smooth(method = "lm", se = T, color = 'skyblue', alpha = 0.25, linewidth = 0.5) +
  labs(x = x_col, y = y_col) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.75),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  # 添加 beta 值
  annotate(
    "text", x = x_max * 0.98, y = y_max - y_distance * 0.02, 
    label = paste("β =", round(regression_result_list$beta_orig, 3)),
    hjust = 1, vjust = 1, size = 6,
    color = '#606060', fontface = 'bold'
  ) +
  # 添加 p 值
  annotate(
    "text", x = x_max * 0.98, y = y_max - y_distance * 0.08, 
    label = paste0("p = ", round(regression_result_list$p_value_z, 3)),
    hjust = 1, vjust = 1, size = 6,
    color = p_value_color, fontface = 'bold'
  )
