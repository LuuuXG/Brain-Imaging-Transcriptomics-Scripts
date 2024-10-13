from nilearn import datasets
import abagen
from abagen import images
import nibabel as nib
import numpy as np
from neuromaps import transforms

# 读取 NIfTI 文件
img = nib.load(r'D:\WYJ\Codes\Imaging-Transcriptomics\data\01_Atlas\AAL_v4\ROI_MNI_V4_1mm_MNI152.nii')  # 替换为你的 NIfTI 文件路径
data = img.get_fdata()

# 获取所有唯一值
unique_values = np.unique(data)

# 打印所有唯一值
print("Unique values in the atlas:")
print(unique_values)

if __name__ == '__main__':
    # AAL from nilearn
    #aal = datasets.fetch_atlas_aal()
    #aal_vol = aal['maps']
    # todo: info

    ###### AAL_V4 from nilearn downloads ######
    aal_v4_vol_2mm = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4.nii' # 2mm resolution
    aal_v4_vol_info_2mm = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4.csv'
    #atlas_check = images.check_atlas(aal_v4_vol, aal_v4_vol_info)

    # 得到1mm的AAL模板
    aal_2mm = nib.load(aal_v4_vol_2mm)
    aal_1mm = transforms.mni152_to_mni152(aal_2mm, '1mm', method='nearest')
    nib.save(aal_1mm, r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4_1mm_MNI152.nii')
    aal_v4_vol_1mm = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4_1mm_MNI152.nii'

    ###### DK from https://github.com/alegiac95/Imaging-transcriptomics ######
    # abagen得到的DK模板和imt自带的DK模板不同，但如果用DK模板，获得基因表达矩阵的结果是一样的
    DK = abagen.fetch_desikan_killiany()
    DK_vol_1mm = DK['image']
    DK_vol_info_1mm = DK['info']

    # 得到2mm的DK模板
    #DK_1mm = nib.load(r'..\data\atlas\desikan_killiany\atlas-desikankilliany_1mm_MNI152.nii')
    #DK_2mm = transforms.mni152_to_mni152(DK_1mm, '2mm')
    #nib.save(DK_2mm, r'..\data\atlas\desikan_killiany\atlas-desikankilliany_2mm_MNI152.nii')
    DK_vol_2mm = r'..\data\atlas\desikan_killiany\atlas-desikankilliany_2mm_MNI152.nii'

    # todo: harvard_oxford
    #harvard_oxford_vol = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr50-2mm')

    ###### Allen atlas ######
    import nibabel as nib
    import pandas as pd
    import numpy as np

    # 读取CSV文件
    csv_path = r'E:\Atlas_Templetes_Masks\Atlas\Allen\label.csv'
    df = pd.read_csv(csv_path)

    # 创建映射字典
    value_map = dict(zip(df['id'], df['graph_order']))

    # 读取NIfTI文件
    nii_path = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full.nii'  # 替换为实际路径
    img = nib.load(nii_path)
    data = img.get_fdata()

    # 创建新的数据数组，用于存储重新赋值的区域
    new_data = np.copy(data)

    # 遍历映射字典，修改数据数组中的值
    for old_value, new_value in value_map.items():
        #print(old_value, new_value)
        new_data[data == old_value] = new_value

    coordinates = np.where(new_data == 141)
    new_data[194, 174, 15]
    data[194, 174, 15]

    # 创建新的NIfTI图像
    new_img = nib.Nifti1Image(new_data, img.affine, img.header)

    # 保存修改后的NIfTI文件
    new_nii_path = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified.nii'  # 替换为实际路径
    nib.save(new_img, new_nii_path)

    print("NIfTI文件重新赋值并保存成功。")

    ## 区分左右半球
    # 读取NIfTI文件
    input_filepath = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified.nii'
    output_filepath = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified_LR.nii'

    # 加载NIfTI文件
    img = nib.load(input_filepath)
    data = img.get_fdata()

    # 计算数据的中线
    midline = data.shape[0] // 2 - 1

    midline_data = data[midline, :, :]

    data[midline, :, :] = 0

    # 找出不为0的值及其位置
    non_zero_indices = np.argwhere(midline_data != 0)
    non_zero_values = np.unique(midline_data[midline_data != 0])

    # 获取所有非零标签的最大值
    max_label = int(np.round(np.max(data[data > 0])))

    # 创建一个新的数组来存储整数结果
    new_data = np.zeros_like(data, dtype=int)

    # 复制左半球的数据
    new_data[:midline, :, :] = np.round(data[:midline, :, :]).astype(int)

    # 修改右半球的非零标签值
    new_data[midline + 1:, :, :] = np.where(data[midline + 1:, :, :] > 0,
                                            np.round(data[midline + 1:, :, :]).astype(int) + max_label, 0)

    # 创建新的NIfTI图像
    new_img = nib.Nifti1Image(new_data, img.affine, img.header)

    # 保存修改后的文件
    nib.save(new_img, output_filepath)

    ## create 2mm atlas
    Allen_05mm = nib.load(r'E:\Codes\Imaging-Transcriptomics\data\atlas\allen_human_reference_atlas_subcortical\allen_human_reference_atlas_subcortical.nii')
    Allen_2mm = transforms.mni152_to_mni152(Allen_05mm, '1mm')
    nib.save(Allen_2mm, r'E:\Codes\Imaging-Transcriptomics\data\atlas\allen_human_reference_atlas_subcortical\allen_human_reference_atlas_subcortical_1mm.nii')

    ## subcortical
    # 读取 atlas 数据
    atlas_path = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified_LR.nii'
    atlas_img = nib.load(atlas_path)
    atlas_data = atlas_img.get_fdata()

    # 读取 CSV 文件数据
    csv_data = pd.read_csv(r'E:\Atlas_Templetes_Masks\Atlas\Allen\allen_human_reference_atlas.csv')

    # 获取需要转换的 id 列和 transform 列
    id_values = csv_data['id']
    transform_values = csv_data['transform']

    # 创建一个新的数组用于存储转换后的数据，并确保数据类型为整数
    new_atlas_data = np.zeros_like(atlas_data, dtype=int)

    # 遍历 id 和 transform 列
    for i, transform in zip(id_values, transform_values):
        new_atlas_data[atlas_data == i] = int(transform)

    # 保存新的 atlas 文件
    new_atlas_img = nib.Nifti1Image(new_atlas_data, atlas_img.affine, atlas_img.header)
    new_atlas_path = r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified_LR_sub.nii'
    nib.save(new_atlas_img, new_atlas_path)

    # check nii value
    img = nib.load(r'E:\Atlas_Templetes_Masks\Atlas\Allen\annotation_full_modified_LR_sub.nii')
    data = img.get_fdata()

    np.any(data == 1)
    value = np.unique(data)

    ###### harvard_oxford ######
    harvard_oxford_vol = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-1mm')

    ###### subcortical ######
    subcortical_1mm = r'..\data\atlas\subcortical_nuclei\subcortical_nuclei.nii'
    subcortical_2mm = transforms.mni152_to_mni152(subcortical_1mm, '2mm')
    nib.save(subcortical_2mm, r'..\data\atlas\subcortical_nuclei\subcortical_nuclei_2mm.nii')