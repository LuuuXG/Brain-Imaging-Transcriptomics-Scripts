#%% set the root output directory
import os
output_root_folder = r'..\data\output_test1'

#%% choose an atlas
from fetch_atlas import fetch_atlas_data

atlas_1mm, atlas_2mm, atlas_info = fetch_atlas_data('dk')

atlas_suffix = 'DK_1mm'

#%% generate gene expression data
from gene_expression import generate_gene_expression

gene_expression_folder = os.path.join(output_root_folder, 'expression')
os.makedirs(gene_expression_folder, exist_ok=True)

gene_expression_path, abagen_report = generate_gene_expression(atlas_1mm, atlas_info, suffix=atlas_suffix, missing='interpolate', output_dir=gene_expression_folder)

#%% create permuted imaging data
from fetch_spatial_nulls import generate_nulls, parcellate_imaging

#imaging：影像数据
#parcellation_mni：选择的图谱（parcellation/atlas），请务必确保parcellation_mni是2mm分辨率的！
imaging_parc_folder = os.path.join(output_root_folder, 'imaging_parcellation')
os.makedirs(imaging_parc_folder, exist_ok=True)

raw_imaging_folder = os.path.join(output_root_folder, 'raw_imaging')
os.makedirs(raw_imaging_folder, exist_ok=True)

imaging_path = os.path.join(raw_imaging_folder, 'mean_QSM_PT.nii.gz')
parcellation_mni = atlas_2mm

imaging_suffix = 'QSM_Chi_mean_PT'
measures_name = 'Chi'

imaging_parc_path = parcellate_imaging(imaging_path, parcellation_mni, atlas_info=atlas_info, suffix=imaging_suffix, output_dir=imaging_parc_folder, measures_name=measures_name)

nulls_folder = os.path.join(output_root_folder, 'spatial_nulls')
os.makedirs(nulls_folder, exist_ok=True)

nulls_path = generate_nulls(imaging_parc_path, parcellation_mni, 1000, suffix=imaging_suffix, output_dir=nulls_folder, measures_name=measures_name)

#%%
from load_data import load_data

imaging_data, nulls_data, gene_exp_orig, gene_labels = load_data(imaging_parc_path, gene_expression_path, nulls_path, 1:83)

#%% component bootstrap
from bootstrap import boot_pls, boot_genes
#from pyls import pls_regression
pls_output_folder = os.path.join(output_root_folder, 'pls_output')
os.makedirs(pls_output_folder, exist_ok=True)

n_components = 3
n_iter = 1000 # 小于等于前面设置的nulls的数目

boot_pls_result = boot_pls(imaging_data, gene_exp_orig, nulls_data, n_components, n_iter, pls_output_folder, suffix=imaging_suffix)

boot_genes_result = boot_genes(imaging_data, gene_exp_orig, gene_labels, nulls_data, n_components, n_iter, pls_output_folder, suffix=imaging_suffix)

