# Written by Youjie Wang (wangyoujie2002@163.com)
# Department of Neurology, West China Hospital, Sichuan University, China
# Last update: 2024-09-26

############################################################################
# make sure the working directory is where this file (`main.R`) is located #
############################################################################

# If directly open the `main.R` file, the working directory is where this file 
# is located. If not, please set the working directory to where this file is 
# located.

setwd('E:/Codes/Imaging-Transcriptomics/R') # [user defined]
getwd()
rm(list = ls())

############################################################################
################### 0. Load basic packages and settings ####################
############################################################################

source('./base/load_packages.R')
# 'crayon' settings
source('./base/crayon_color_settings.R')
source('./base/create_folder.R')

# set python path
#use_python("D:/Anaconda/python.exe", required = TRUE) # [user defined]
use_python('E:/Anaconda3/python.exe', required = TRUE) # [user defined]

# set output root path
output_root_folder <- 
  'F:/QSM_process/YC+YT_QSM_BIDS/derivatives/IMT/VBA_v240929_TFCE' # [user defined]
create_folder(output_root_folder)

############################################################################
##################### 1. Generate gene expression data ##################### 
############################################################################

#### 1.1. Define a brain atlas ####
source_python('../python/fetch_atlas.py')

result_temp <- fetch_atlas_data('dk') # [user defined]

atlas_1mm <- result_temp[[1]]
atlas_2mm <- result_temp[[2]]
atlas_info <- result_temp[[3]]

#### 1.2. Generate gene expression data (use `abagen`) ####
gene_expression_folder <- file.path(output_root_folder, 'expression')
create_folder(gene_expression_folder)
source_python('../python/gene_expression.py')

atlas_suffix <- 'DK_1mm' # [user defined]
missing <- 'interpolate' # [user defined]
result_temp <- generate_gene_expression(atlas_1mm, 
                                        atlas_info, 
                                        suffix=atlas_suffix, 
                                        output_dir=gene_expression_folder, 
                                        missing=missing)

gene_expression_path <- result_temp[[1]]
abagen_report <- result_temp[[2]]

############################################################################
################ 2. Generate imaging data and spatial nulls ################
############################################################################

source_python('../python/fetch_spatial_nulls.py')
#### 2.1. Parcellate imaging data ####
imaging_parc_folder <- file.path(output_root_folder, 'imaging_parcellation')
create_folder(imaging_parc_folder)

raw_imaging_folder <- file.path(output_root_folder, 'raw_imaging')
create_folder(raw_imaging_folder)

imaging_file_name <- 'TFCE_0001.nii' # [user defined]: please put the imaging
                                     # file in the `raw_imaging` folder
imaging_path <- file.path(raw_imaging_folder, imaging_file_name)
parcellation_mni <- atlas_2mm

imaging_suffix <- 'QSM_TFCE_PT-HC' # [user defined]: used in output file names
measures_name <- 'TFCE' # [user defined]: used in the column name of the output
                        # pacellated imaging file (.csv)
imaging_parc_path <- parcellate_imaging(imaging_path, 
                                        parcellation_mni, 
                                        atlas_info=atlas_info, 
                                        suffix=imaging_suffix, 
                                        output_dir=imaging_parc_folder,
                                        measures_name=measures_name)

#### 2.2. Generate spatial nulls ####
nulls_folder <- file.path(output_root_folder, 'spatial_nulls')
create_folder(nulls_folder)

n_perm_nulls <- 1000 # [user defined]
nulls_path = generate_nulls(imaging_parc_path, parcellation_mni, 
                            n_perm=n_perm_nulls, suffix=imaging_suffix,
                            output_dir=nulls_folder,
                            measures_name=measures_name) 

#### 2.3. Save RData for future use ####
save(imaging_parc_path, gene_expression_path, nulls_path, 
     file=file.path(output_root_folder, 'gene_exp+img_parc+nulls.RData'))

############################################################################
##################### 3. Link genes with neuroimaging ######################
############################################################################

#### 3.0(A). Load data path (From .RData) ####
load(file=file.path(output_root_folder, 'gene_exp+img_parc+nulls.RData'))

#### 3.0(B). Load data path (From python process) ####
# todo

#### 3.1. Load data (imaging, gene expression, and spatial nulls) ####
source_python('../python/load_data.py')
result_temp <- load_data(imaging_parc_path, 
                         gene_expression_path, 
                         nulls_path, 
                         roi_range=1:83)

imaging_data <- result_temp[[1]]
nulls_data <- result_temp[[2]]
gene_exp_orig <- result_temp[[3]]
gene_labels <- result_temp[[4]]

#### 3.2(A) PLS + GSEA ####
### 3.2(A).1. PLS Regression ###
source_python('../python/bootstrap.py')
pls_output_folder <- file.path(output_root_folder, 'pls_output')
create_folder(pls_output_folder)

n_components <- 10 # [user defined]
n_iter <- 1000 # [user defined]: no more than spatial nulls

## pls components
boot_pls_result <- boot_pls(imaging_data, gene_exp_orig, nulls_data, 
                                      n_components, n_iter, pls_output_folder, 
                                      suffix=imaging_suffix)

## genes importance
# You can check the .csv file generated in the `pls_output_folder` to see the
# explained variance of each component. This can help to determine the number
# of components to use in the following gene importance analysis, which can 
# save the processing time.
n_components <- 3 # [user defined]: reset n_components
boot_genes_result <- boot_genes(imaging_data, gene_exp_orig, gene_labels, 
                                nulls_data, n_components, n_iter, 
                                pls_output_folder, suffix=imaging_suffix)

save(boot_pls_result, boot_genes_result, 
     file=file.path(output_root_folder, 'boot_pls+boot_genes.RData'))

### 3.2(A).2 GSEA ###
gc()
load(file=file.path(output_root_folder, 'boot_pls+boot_genes.RData'))

gsea_output_folder <- file.path(output_root_folder, 'gsea_output')
create_folder(gsea_output_folder)

source('./GSEA/GSEA_GO.R')
source('./GSEA/GSEA_MSigDB.R')
source('./GSEA/GSEA_user_defined_GeneSet.R')
source('./GSEA/simplify_by_NES.R')

## GSEA GO
gsea_GO_results_data <- list()
ont <- 'ALL' # [user defined]: 'BP', 'CC', 'MF', 'ALL'

for (i in 1:length(boot_genes_result)){
  gsea_GO_results_data[[i]] <- gsea_GO_from_pls(boot_genes_result[[i]],
                                                ont=ont, threshold=1,
                                                output_dir=gsea_output_folder,
                                                simplify=T,
                                                simplify_cutoff=0.6)
}

save(gsea_GO_results_data, 
     file=file.path(output_root_folder, 
                    paste0('gsea_GO_results_data_ont=', ont, '.RData')))

## GSEA MsigDB
gsea_MsigDB_results_data <- list()
category <- 'C5' # [user defined]
subcategory <- NULL # [user defined]

for (i in 1:length(boot_genes_result)){
  gsea_MsigDB_results_data[[i]] <- gsea_MSigDB_from_pls(
    boot_genes_result[[i]], 
    category=category,
    subcategory=subcategory,
    output_dir=gsea_output_folder)
}

save(gsea_MsigDB_results_data, 
     file=file.path(output_root_folder, 
                    paste0('gsea_MsigDB_results_data_category=', 
                           category, '_subcategory=', subcategory, '.RData')))

## GSEA user defined gene set (cell-specific data for example)
source('./Cell_Specific/pSI_geneset_threshold.R')

pSI_threshold <- 0.001 # [user defined]
cell_specific_geneset <- fetch_cell_specific_data(pSI_threshold)

gsea_cell_specific_results_data <- list()

for (i in 1:length(boot_genes_result)){
  gsea_cell_specific_results_data[[i]] <- gsea_user_defined_from_pls(
    boot_genes_result[[i]],
    cell_specific_geneset, threshold=1,
    output_dir=gsea_output_folder)
}

save(gsea_cell_specific_results_data, 
     file=file.path(output_root_folder, 
                    paste0('gsea_cell_specific_results_data_pSI_threshold=', 
                           pSI_threshold, '.RData')))

#### 3.2(B) WGCNA ####
WGCNA_output_folder <- file.path(output_root_folder, 'WGCNA_output')
create_folder(WGCNA_output_folder)

source('./WGCNA/load_data_WGCNA.R')
result_temp <- load_data_raw(imaging_parc_path, gene_expression_path,
                             nulls_path, 1:83)

imaging_data_raw <- result_temp[[1]]
gene_expression_data_raw <- result_temp[[2]]
nulls_data_raw <- result_temp[[3]]

# todo
