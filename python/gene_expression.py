#%% generate gene expression data
import abagen
import pandas as pd

def generate_gene_expression(atlas_for_AHBA_image, atlas_for_AHBA_info, suffix, output_dir, missing=None):
    files = abagen.fetch_microarray(donors='all', verbose=0)

    # get dense gene expression data
    expression, report = abagen.get_expression_data(atlas_for_AHBA_image, atlas_for_AHBA_info, missing=missing, return_report=True)

    # 输出路径
    gene_expression_path = output_dir + f'/gene_expression_data_{suffix}.csv'
    expression.to_csv(gene_expression_path, index=True)

    # 添加Region列（和imt保持一致）
    gene_expression_data = pd.read_csv(gene_expression_path)
    atlas_data = pd.read_csv(atlas_for_AHBA_info)

    # 根据'label'和'id'合并数据
    gene_expression_data_merged = gene_expression_data.merge(atlas_data[['id', 'label']], left_on='label', right_on='id')
    gene_expression_data_merged.rename(columns={'label_x': 'label', 'label_y': 'Region'}, inplace=True)
    gene_expression_data_merged.drop(columns=['id'], inplace=True)
    column_order = ['label', 'Region'] + [col for col in gene_expression_data_merged.columns if col not in ['label', 'Region']]
    gene_expression_data_merged = gene_expression_data_merged[column_order]

    # 保存处理后的数据到新的CSV文件
    gene_expression_data_merged.to_csv(gene_expression_path, index=False)
    print(f"Processed data saved to {gene_expression_path}")

    return gene_expression_path, report
