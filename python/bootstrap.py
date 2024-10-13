import pandas as pd
import numpy as np
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
from pyls import pls_regression

def correlate(c1, c2):
    """Return the MATLAB style correlation between two vectors."""
    return np.corrcoef(np.hstack((c1, c2)), rowvar=False)[:-1, -1]

def boot_pls(imaging_data, gene_exp_orig, nulls_data, n_components, n_iter, output_dir, suffix):
    n_components = int(n_components)
    n_iter = int(n_iter)

    # 初始化存储结果的数组
    cumulative_r2 = np.zeros(n_components)
    single_explained_var = np.zeros(n_components)
    single_p_val = np.zeros(n_components)
    cumulative_p_val = np.zeros(n_components)

    # 计算原始数据的累计解释方差
    imaging_data_copy = np.array(imaging_data, copy=True)
    gene_exp_copy = np.array(gene_exp_orig, copy=True)

    _res = pls_regression(gene_exp_copy, imaging_data_copy.reshape(imaging_data.shape[0], 1),
                          n_components=n_components, n_perm=0, n_boot=0)
    original_exp_var = _res.get("varexp")
    r1 = correlate(_res.get("x_scores"), imaging_data_copy.reshape(imaging_data.shape[0], 1))

    # 计算累计解释方差
    original_cumulative_exp_var = np.array([np.sum(original_exp_var[:i + 1]) for i in range(n_components)])

    # 置换测试
    permuted_cumulative_exp_var = np.zeros((n_iter, n_components))
    permuted_explained_var = np.zeros((n_iter, n_components))

    for i in range(n_iter):
        if (i+1) % 100 == 0:
            print(f"Bootstrap iteration: {i + 1}/{n_iter}")

        y_data = nulls_data[:, i].reshape(imaging_data.shape[0], 1)
        gene_exp_copy = np.array(gene_exp_orig, copy=True)
        y_data_copy = np.array(y_data, copy=True)
        _res = pls_regression(gene_exp_copy, y_data_copy, n_components=n_components, n_perm=0, n_boot=0)
        permuted_exp_var = _res.get("varexp")
        permuted_cumulative_exp_var[i, :] = [np.sum(permuted_exp_var[:j + 1]) for j in range(n_components)]
        permuted_explained_var[i, :] = permuted_exp_var[:n_components]

    # 计算p值
    for component in range(n_components):
        cumulative_r2[component] = original_cumulative_exp_var[component]
        cumulative_p_val[component] = np.sum(
            permuted_cumulative_exp_var[:, component] >= original_cumulative_exp_var[component]) / n_iter
        single_explained_var[component] = original_exp_var[component]
        single_p_val[component] = np.sum(
            permuted_explained_var[:, component] >= original_exp_var[component]) / n_iter

    # 保存结果
    results_df = pd.DataFrame({
        'Component': np.arange(1, n_components + 1),
        'Cumulative R2': cumulative_r2,
        'Explained Variance': single_explained_var,
        'P-value (Cumulative)': cumulative_p_val,
        'P-value (Single)': single_p_val,
        'Score-Imaging Correlation': r1
    })

    # 指定输出文件路径
    output_file = f"{output_dir}/pls_bootstrap_results_{suffix}.csv"
    results_df.to_csv(output_file, index=False)

    print(f"Results saved to {output_file}")

    return output_file

def boot_genes(imaging_data, gene_exp_orig, gene_labels, nulls_data, n_components, n_iter, output_dir, suffix):
    n_components = int(n_components)
    n_iter = int(n_iter)

    # PLS回归
    gene_exp = gene_exp_orig
    n_genes = gene_exp.shape[1]

    orig_data = {
        'weights': np.zeros((n_components, n_genes)),
        'genes': np.zeros((n_components, n_genes), dtype=object),
        'index': np.zeros((n_components, n_genes), dtype=np.int32),
        'zscored': np.zeros((n_components, n_genes))
    }

    # 初始化 boot 数据存储字典
    boot_data = {
        'weights': np.zeros((n_components, n_genes, n_iter)),
        'std': np.zeros((n_components, n_genes)),
        'z_score': np.zeros((n_components, n_genes)),
        'pval': np.zeros((n_components, n_genes)),
        'pval_corr': np.zeros((n_components, n_genes)),
        'genes': np.zeros((n_components, n_genes), dtype=object)
    }

    imaging_data_copy = np.array(imaging_data, copy=True)
    gene_exp_copy = np.array(gene_exp_orig, copy=True)
    res = pls_regression(gene_exp_copy, imaging_data_copy.reshape(imaging_data.shape[0], 1), n_components=n_components, n_boot=0,
                         n_perm=0)

    permuted_imaging = nulls_data

    n_roi = imaging_data.shape[0]

    #r1 = correlate(res.get("x_scores"), imaging_data_orig.reshape(n_roi, 1))
    r1 = correlate(res.get("x_scores"), imaging_data.reshape(n_roi, 1))

    #print('correlation with original imaging data: R^2 = {0}'.format(r0))
    print('correlation of each component with z-scored imaging data: R^2 = {0}'.format(r1))

    weights = res.get("x_weights")
    scores = res.get("x_scores")

    # 处理权重方向
    for i in range(r1.size):
        if r1[i] < 0:
            weights[:, i] *= -1
            scores[:, i] *= -1

    # 排序并存储原始结果
    for idx in range(n_components):
        sorted_indexes = np.argsort(weights[:, idx], kind="mergesort")[::-1]
        orig_data['index'][idx, :] = sorted_indexes
        orig_data['genes'][idx, :] = gene_labels[sorted_indexes]
        orig_data['weights'][idx, :] = weights[sorted_indexes, idx]
        orig_data['zscored'][idx, :] = zscore(orig_data['weights'][idx, :], ddof=1)

    permuted_imaging_copy = np.array(permuted_imaging, copy=True)
    # Bootstrap
    for iter in range(n_iter):
        if (iter+1) % 100 == 0:
            print(f"Bootstrap iteration: {iter + 1}/{n_iter}")

        gene_exp_copy = np.array(gene_exp_orig, copy=True)
        res_i = pls_regression(gene_exp_copy, permuted_imaging_copy[:, iter].reshape(imaging_data.shape[0], 1), n_components=n_components, n_boot=0, n_perm=0)
        weights_i = res_i.get("x_weights")

        for comp in range(n_components):
            temp_weights = weights_i[:, comp]
            new_weights = temp_weights[orig_data['index'][comp, :]]

            corr = correlate(
                orig_data['weights'][comp, :].reshape(orig_data['weights'][comp, :].size, 1),
                new_weights.reshape(new_weights.shape[0], 1)
            )

            if corr < 0:
                new_weights *= -1
            boot_data['weights'][comp, :, iter] = new_weights

    # 计算统计量
    for component in range(n_components):
        boot_data['std'][component, :] = boot_data['weights'][component, :, :].std(axis=1, ddof=1)
        z = orig_data['weights'][component, :] / boot_data['std'][component, :]
        boot_data['z_score'][component, :] = np.sort(z, axis=0, kind="mergesort")[::-1]
        index = np.argsort(z, axis=0, kind="mergesort")[::-1]
        boot_data['genes'][component, :] = orig_data['genes'][component, index]
        p_val = norm.sf(abs(boot_data['z_score'][component, :]))
        boot_data['pval'][component, :] = p_val

        _, p_corr, _, _ = multipletests(p_val, method='fdr_bh', is_sorted=False)

        boot_data['pval_corr'][component, :] = p_corr

    # 保存结果
    result_paths = []  # 初始化存储结果路径的列表
    for component in range(n_components):
        result_df = pd.DataFrame({
            'Gene': boot_data['genes'][component, :],
            'Z-score': boot_data['z_score'][component, :],
            'p-value': boot_data['pval'][component, :],
            'p-value (FDR)': boot_data['pval_corr'][component, :]
        })
        result_file_path = f'{output_dir}/pls_component_{component + 1}_{suffix}.csv'
        result_df.to_csv(result_file_path, index=False)
        result_paths.append(result_file_path)

        print(f"Results saved to {result_file_path}")

    return result_paths