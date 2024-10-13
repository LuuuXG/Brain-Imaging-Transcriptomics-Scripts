# plot parcellated imaging data
# available atlases: DK

#install.packages('ggseg')
library(ggseg)
library(ggplot2)
library(dplyr)
library(tidyr)

# # Enable this universe
# options(repos = c(
#   ggseg = 'https://ggseg.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
# 
# # Install some packages
# install.packages('ggseg3d')
library(ggseg3d)

#library(plotly)

library(htmlwidgets)

############################################################################
############################# Data preparation #############################
############################################################################

load(file=file.path(output_root_folder, 'gene_exp+img_parc+nulls.RData'))

imaging_parc_folder <- file.path(output_root_folder, 'imaging_parcellation')

imaging_parc_data <- read.csv(imaging_parc_path)

color_max <- max(imaging_parc_data[3])
color_min <- min(imaging_parc_data[3])

measures_name <- colnames(imaging_parc_data)[3]
colnames(imaging_parc_data)[3] <- "value"

# palette_magma <- c("#FCFDBF"=color_min, 
#                    "#FBD581"=color_min+1*(color_max-color_min)/9, 
#                    "#F89E4F"=color_min+2*(color_max-color_min)/9, 
#                    "#EC663C"=color_min+3*(color_max-color_min)/9, 
#                    "#CC375E"=color_min+4*(color_max-color_min)/9, 
#                    "#A0206F"=color_min+5*(color_max-color_min)/9, 
#                    "#6D1972"=color_min+6*(color_max-color_min)/9, 
#                    "#41166C"=color_min+7*(color_max-color_min)/9, 
#                    "#1A0F53"=color_min+8*(color_max-color_min)/9, 
#                    "#000004"=color_max)

palette_coolwarm <- c("#99CCFF" = color_min,
                      "#3399FF" = color_min*0.75,
                      "#0066CC" = color_min*0.5,
                      "#003366" = color_min*0.25,
                      "#001933" = color_min*0.125,
                      "black" = 0,
                      "#330000" = color_max*0.125,
                      "#660000" = color_max*0.25,
                      "#CC0000" = color_max*0.5,
                      "#FF3333" = color_max*0.75,
                      "#FF9999" = color_max)

#### cortical ####
imaging_parc_cortical_data <- imaging_parc_data[c(1:34, 42:75),]

imaging_parc_cortical_data$label <- ifelse(grepl("_L$", imaging_parc_cortical_data$label),
                                           paste0("lh_", sub("_L$", "", imaging_parc_cortical_data$label)),
                                           imaging_parc_cortical_data$label)

imaging_parc_cortical_data$label <- ifelse(grepl("_R$", imaging_parc_cortical_data$label),
                                           paste0("rh_", sub("_R$", "", imaging_parc_cortical_data$label)),
                                           imaging_parc_cortical_data$label)

#### subcortical ####
imaging_parc_subcortical_data <- imaging_parc_data[c(35:38, 40:41, 76:79, 81:82, 83),] # 83: brainstem

label_rename <- c('Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', 'Left-Hippocampus', 'Left-Amygdala',
                  'Right-Thalamus-Proper', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala',
                  'Brain-Stem')

imaging_parc_subcortical_data$label <- label_rename

############################################################################
################################### PLOT ###################################
############################################################################

#### 1. cortical parcellation (2D) ####
plot_cortical_2d <- ggplot(imaging_parc_cortical_data) +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = value)) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_void()

jpeg(file = file.path(imaging_parc_folder, "cortical_parcellation_2D.jpg"), 
     res = 300, width = 3000, height = 3000)
print(plot_cortical_2d)
dev.off()

#### 2. subcortical parcellation (2D) ####
plot_subcortical_2d <- ggplot(imaging_parc_subcortical_data) +
  geom_brain(atlas = aseg, 
             aes(fill = value)) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_void()

jpeg(file = file.path(imaging_parc_folder, "subcortical_parcellation_2D.jpg"), 
     res = 300, width = 3000, height = 3000)
print(plot_subcortical_2d)
dev.off()

#### 3. cortical parcellation (3D) ####

# right hemisphere
someData_rh = dk_3d %>% 
  filter(surf == "inflated" & hemi == "right") %>% 
  unnest(ggseg_3d) %>% 
  ungroup() %>% 
  dplyr::select(c(region, label)) %>% 
  na.omit()

someData_rh <- merge(someData_rh, imaging_parc_cortical_data, by = "label")

plot_cortical_3d_rh_lateralview <- ggseg3d(.data = someData_rh, 
                   atlas = dk_3d,
                   colour = "value", text = "value",
                   hemisphere = "right",
                   palette = palette_coolwarm,
                   options.legend = list(title=list(text=measures_name))) %>% 
  pan_camera("right lateral") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

plot_cortical_3d_rh_medialview <- ggseg3d(.data = someData_rh, 
                   atlas = dk_3d,
                   colour = "value", text = "value",
                   hemisphere = "right",
                   palette = palette_coolwarm,
                   options.legend = list(title=list(text=measures_name))) %>% 
  pan_camera("right medial") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

htmlwidgets::saveWidget(plot_cortical_3d_rh_lateralview, file = file.path(imaging_parc_folder, "cortical_parcellation_rh_3D_lateral.html"))
htmlwidgets::saveWidget(plot_cortical_3d_rh_medialview, file = file.path(imaging_parc_folder, "cortical_parcellation_rh_3D_medial.html"))

# left hemisphere
someData_lh = dk_3d %>% 
  filter(surf == "inflated" & hemi == "left") %>% 
  unnest(ggseg_3d) %>% 
  ungroup() %>% 
  dplyr::select(c(region, label)) %>% 
  na.omit()

someData_lh <- merge(someData_lh, imaging_parc_cortical_data, by = "label")

plot_cortical_3d_lh_lateralview <- ggseg3d(.data = someData_lh, 
                atlas = dk_3d,
                colour = "value", text = "value",
                hemisphere = "left",
                palette = palette_coolwarm,
                options.legend = list(title=list(text=measures_name))) %>% 
  pan_camera("left lateral") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

plot_cortical_3d_lh_medialview <- ggseg3d(.data = someData_lh, 
                atlas = dk_3d,
                colour = "value", text = "value",
                hemisphere = "left",
                palette = palette_coolwarm,
                options.legend = list(title=list(text=measures_name))) %>% 
  pan_camera("left medial") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

htmlwidgets::saveWidget(plot_cortical_3d_lh_lateralview, file = file.path(imaging_parc_folder, "cortical_parcellation_lh_3D_lateral.html"))
htmlwidgets::saveWidget(plot_cortical_3d_lh_medialview, file = file.path(imaging_parc_folder, "cortical_parcellation_lh_3D_medial.html"))

# subcortical parcellation (3D)
someData_aseg = aseg_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  dplyr::select(label)

someData_aseg <- merge(someData_aseg, imaging_parc_subcortical_data, by = "label")

plot_subcortical_3d_leftview <- ggseg3d(.data = someData_aseg, atlas = aseg_3d, 
        colour = "value", text = "value", 
        na.alpha= .5,
        palette = palette_coolwarm,
        options.legend = list(title=list(text=measures_name))) %>% 
  add_glassbrain() %>%
  pan_camera("left lateral") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

plot_subcortical_3d_rightview <- ggseg3d(.data = someData_aseg, atlas = aseg_3d, 
        colour = "value", text = "value", 
        na.alpha= .5,
        palette = palette_coolwarm,
        options.legend = list(title=list(text=measures_name))) %>% 
  add_glassbrain() %>%
  pan_camera("right lateral") %>%
  plotly::layout(
    scene = list(
      xaxis = list(visible = FALSE),  
      yaxis = list(visible = FALSE),  
      zaxis = list(visible = FALSE)   
    )
  )

htmlwidgets::saveWidget(plot_subcortical_3d_leftview, file = file.path(imaging_parc_folder, "subcortical_parcellation_3D_leftview.html"))
htmlwidgets::saveWidget(plot_subcortical_3d_rightview, file = file.path(imaging_parc_folder, "subcortical_parcellation_3D_rightview.html"))

#### plot the colorbar if needed ####
color_data <- data.frame(
  value = c(color_min, color_min*0.75, color_min*0.5, color_min*0.25, color_min*0.125, 0, color_max*0.125, color_max*0.25, color_max*0.5, color_max*0.75, color_max),
  color = c("#99CCFF", "#3399FF", "#0066CC", "#003366", "#001933", "black", "#330000", "#660000", "#CC0000", "#FF3333", "#FF9999")
)

gradient_data <- data.frame(
  x = 1,  # 固定 x 位置
  y = seq(color_min, color_max, length.out = 1000),  # 连续的 y 值
  value = seq(color_min, color_max, length.out = 1000)  # 连续的 value 值
)

# 使用 ggplot2 绘制自定义颜色条
colorbar_plot <- ggplot(gradient_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = color_data$color, 
    values = scales::rescale(color_data$value, to = c(0, 1)), # 将颜色值范围缩放到 0-1
    guide = "colorbar",
    na.value = "grey50"
  ) +
  theme_void() +  # 移除背景
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 1.5,  # 控制图例的宽度
      barheight = 20   # 控制图例的高度
    )
  ) +
  labs(fill = "Value")

# 显示颜色条
print(colorbar_plot)
