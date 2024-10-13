library(ggplot2)
library(cowplot)

load(file=file.path(output_root_folder, 'boot_pls+boot_genes.RData'))

pls_output_folder <- file.path(output_root_folder, 'pls_output')

pls_variance_df <- read.csv(boot_pls_result)

# 取前n个成分
n_components <- 10
pls_variance_df <- pls_variance_df[1:n_components,]

############################################################################
################################### PLOT ###################################
############################################################################

#### variance only ####
plot_variance <- ggplot(pls_variance_df, aes(x = Component)) + 
  geom_bar(aes(y = Explained.Variance), stat = "identity", fill = "skyblue", width = 0.7) + 
  geom_line(aes(y = Cumulative.R2 * max(Explained.Variance), group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = Cumulative.R2 * max(Explained.Variance)), color = "red", size = 2) +
  scale_y_continuous(
    name = "Explained Variance",
    sec.axis = sec_axis(~./max(pls_variance_df$Explained.Variance), name = "Cumulative Variance"),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = pls_variance_df$Component,
    expand = c(0.02, 0)  # 确保与热图对齐
  ) +
  labs(x = "Component", title = "PLS Explained Variance and Cumulative Variance") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # 外边距
  )

#### variance + correlation ####
plot_variance_and_correlation <- ggplot(pls_variance_df, aes(x = Component)) + 
  geom_bar(aes(y = Explained.Variance, fill = Score.Imaging.Correlation), stat = "identity", width = 0.7) + 
  geom_line(aes(y = Cumulative.R2 * max(Explained.Variance), group = 1), color = "#666666", linewidth = 1) +
  geom_point(aes(y = Cumulative.R2 * max(Explained.Variance)), color = "#333333", size = 2) +
  geom_text(aes(y = Explained.Variance+0.008, 
                label = paste0("bold(italic(R)^2 == ", round(Explained.Variance, 2), ")")),
            vjust = -0.5, color = "black", size = 4, parse = TRUE) +
  geom_text(aes(y = Explained.Variance-0.001, 
                label = paste0("bold(italic(corr) == ", round(Score.Imaging.Correlation, 2), ")")),
            vjust = -0.5, color = "black", size = 4, parse = TRUE) +
  scale_fill_gradient2(low = "#6600FF", mid = "white", high = "#FF0066", midpoint = 0, name = "Correlation") +  # 设置颜色渐变
  scale_y_continuous(
    name = "Explained Variance",
    sec.axis = sec_axis(~./max(pls_variance_df$Explained.Variance), name = "Cumulative Variance"),
    expand = c(0, 0, 0.1, 0)
  ) +
  #ylim(0, max(pls_variance_df$Explained.Variance) * 1.05) +
  scale_x_continuous(
    breaks = pls_variance_df$Component,
    expand = c(0.025, 0, 0.025, 0)  # 确保与热图对齐
  ) +
  labs(x = "PLS Component") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 20),
    axis.ticks.y = element_line(linewidth = 0.75),
    plot.title = element_text(color = "black", hjust = 0.5, size = 20),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

print(plot_variance_and_correlation)

############################################################################
################################### SAVE ###################################
############################################################################

jpeg(file=file.path(pls_output_folder, 
                    paste0('plot_pls_variance_n=', n_components, '.jpeg')), 
     width=3600, height=2000, res=300)
print(plot_variance_and_correlation)
dev.off()
