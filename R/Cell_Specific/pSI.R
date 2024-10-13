# This script generates cell-specific gene sets from the pSI output files. The 
# ouput files are saved in the `../data/04_pSI_output` directory.

# If you want to replicate this script, you need to install the pSI package.
# Please download the source code from 
# https://sites.wustl.edu/doughertylab/psi_package-page/ and install the package
# from the downloaded .gz file.

# For example, if you have downloaded the pSI_1.1.tar_.gz and 
# pSI.data_1.0.tar_.gz in the current directory, you can install the package by 
# running the following:

# install.packages('./pSI_1.1.tar_.gz', repos = NULL, type="source")
# install.packages('./pSI.data_1.0.tar_.gz', repos = NULL, type="source")

library(pSI)
library(pSI.data)

wd <- getwd()

geo_rawdata_path <- '../data/03_GEO_Rawdata'

gse_number <- c('GSE52564', 'GSE67835', 'GSE73721')

cell_type <- list(
  GSE52564 = c('Astrocyte', 'Neuron', 'OPC', 'NFO', 'Oligodendrocyte', 
               'Microglia', 'Endothelial'),
  GSE67835 = c('Oligodendrocyte', 'Astrocyte', 'OPC', 'Microglia', 'Neuron', 
               'Endothelial'),
  GSE73721 = c('Astrocyte', 'Neuron', 'Oligodendrocyte', 'Myeloid', 
               'Endothelial')
)

cell_averaged_data_path <- file.path(geo_rawdata_path, 
                                     gse_number, 
                                     paste0(gse_number, 
                                            '_fpkm_cell_type_averaged.csv'))

cell_specific_data_list <- list()

pSI_output_path <- '../data/04_pSI_output'

for (i in 1:length(gse_number)){
  cell_averaged_data <- read.csv(cell_averaged_data_path[i])
  
  # remove rows with duplicated gene names
  cell_averaged_data <- 
    cell_averaged_data[!duplicated(cell_averaged_data$gene),]
  
  rownames(cell_averaged_data) <- cell_averaged_data[,1]
  cell_averaged_data <- cell_averaged_data[,-1]
  
  colnames(cell_averaged_data) <- cell_type[[gse_number[i]]]
  
  pSI.output <- specificity.index(pSI.in=cell_averaged_data, e_min=0.3)
  pSI.out.list <- pSI.list(pSIs=pSI.output, write.csv=T)
  file.rename(file.path(wd, 'pSi_0.001.csv'), 
              file.path(pSI_output_path, 
                        paste0(gse_number[i], '_pSi_0.001.csv')))
  file.rename(file.path(wd, 'pSi_0.01.csv'), 
              file.path(pSI_output_path, 
                        paste0(gse_number[i], '_pSi_0.01.csv')))
  file.rename(file.path(wd, 'pSi_0.05.csv'), 
              file.path(pSI_output_path, 
                        paste0(gse_number[i], '_pSi_0.05.csv')))
  file.rename(file.path(wd, 'pSi_1e-04.csv'), 
              file.path(pSI_output_path, 
                        paste0(gse_number[i], '_pSi_0.0001.csv')))
  
  cell_specific_data_list[[gse_number[i]]] <- list(
    cell_averaged_data_path = cell_averaged_data_path[i],
    pSI_output_path = list(
      pSi_0.001 = file.path(pSI_output_path, paste0(gse_number[i], 
                                                    '_pSi_0.001.csv')),
      pSi_0.01 = file.path(pSI_output_path, paste0(gse_number[i], 
                                                   '_pSi_0.01.csv')),
      pSi_0.05 = file.path(pSI_output_path, paste0(gse_number[i], 
                                                   '_pSi_0.05.csv')),
      pSi_0.0001 = file.path(pSI_output_path, paste0(gse_number[i], 
                                                     '_pSi_0.0001.csv'))
    )
  )
}

save(gse_number, cell_specific_data_list, 
     file = '../data/04_pSI_output/pSI_output.RData')