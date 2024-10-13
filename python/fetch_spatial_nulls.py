#%% create permuted imaging data
from neuromaps.images import load_nifti
from neuromaps.parcellate import Parcellater
import pandas as pd

def parcellate_imaging(imaging_path, parcellation_mni, atlas_info, suffix, output_dir, measures_name='Measures'):
    # 加载NIfTI文件
    imaging = load_nifti(imaging_path)
    parcellater_mni = Parcellater(parcellation_mni, 'mni152')
    imaging_parc = parcellater_mni.fit_transform(imaging, space='mni152', ignore_background_data=True, background_value=0)

    # 加载atlas_info文件，并获取前两列
    atlas_info_df = pd.read_csv(atlas_info)
    atlas_info_columns = atlas_info_df.iloc[:, :2]

    # 将成像数据进行分区并转换为DataFrame
    imaging_parc_reshaped = imaging_parc.reshape(-1, 1)
    imaging_parc_df = pd.DataFrame(imaging_parc_reshaped, columns=[measures_name])

    # 合并atlas_info的前两列和imaging_parc_df
    merged_df = pd.concat([atlas_info_columns, imaging_parc_df], axis=1)

    # 保存结果到CSV文件
    imaging_parc_path = output_dir + f'/imaging_data_{suffix}.csv'
    merged_df.to_csv(imaging_parc_path, index=False)

    print(f"Processed data saved to {imaging_parc_path}")

    return imaging_parc_path

def generate_nulls(imaging_parc_path, parcellation_mni, n_perm, suffix, output_dir, measures_name='Measures'):
    #imaging：影像数据
    #parcellation_mni：选择的图谱（parcellation/atlas），请务必确保parcellation_mni是2mm分辨率的！

    imaging_parc = pd.read_csv(imaging_parc_path)[measures_name].values.reshape(-1, 1)

    n_perm = int(n_perm) # avoid error using `reticulate` in R

    # 这一步大约数十分钟，和n_perm关系不大，主要就是parcellation的分辨率
    from neuromaps import nulls
    nulls = nulls.burt2020(imaging_parc, atlas='MNI152', density='2mm',n_perm=n_perm, seed=1234, parcellation=parcellation_mni)

    # save nulls to .csv file
    nulls_df = pd.DataFrame(nulls)
    nulls_path = output_dir + f'/nulls_{suffix}.csv'
    nulls_df.to_csv(nulls_path, index=False)

    print(f"Processed data saved to {nulls_path}")

    return nulls_path
