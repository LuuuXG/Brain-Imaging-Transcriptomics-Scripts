import pandas as pd
from scipy.stats import zscore

def load_data(imaging_parc_path, gene_expression_path, nulls_path, roi_range=None):
    imaging_data_df = pd.read_csv(imaging_parc_path)
    gene_expression_data_df = pd.read_csv(gene_expression_path)
    nulls_df = pd.read_csv(nulls_path)

    # 提取数值数据
    imaging_data_orig = imaging_data_df.iloc[:, 2:].values
    gene_exp = gene_expression_data_df.iloc[:, 2:].values
    gene_labels = gene_expression_data_df.columns[2:]
    nulls = nulls_df.values

    if roi_range is not None:
        adjusted_range = [r - 1 for r in roi_range]  # 将用户输入的范围减去1来调整为正确的索引
        imaging_data_orig = imaging_data_orig[adjusted_range, :]
        gene_exp = gene_exp[adjusted_range, :]
        nulls = nulls[adjusted_range, :]

    # 标准化数据
    imaging_data = zscore(imaging_data_orig, axis=0, ddof=1)
    nulls_data = zscore(nulls, axis=0, ddof=1)
    gene_exp = zscore(gene_exp, ddof=1)
    gene_exp_orig = gene_exp # because of this bug (https://github.com/netneurolab/pypyls/issues/63)

    print(f"Successfully loaded data from {imaging_parc_path}, {gene_expression_path}, and {nulls_path}")

    return imaging_data, nulls_data, gene_exp_orig, gene_labels