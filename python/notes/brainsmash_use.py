import nibabel as nib
import numpy as np

# 加载 NIfTI 文件
t1_image = nib.load(r'F:\Codes\AHBA\data\T1_biascorr.nii')
brain_mask = nib.load(r'F:\Codes\AHBA\data\T1_biascorr_brain_mask.nii')

# 获取数据
t1_data = t1_image.get_fdata()
mask_data = brain_mask.get_fdata()

# 找出掩膜内的体素
indices = np.where(mask_data == 1)  # 假设掩膜中的脑区域标记为 1
coords = np.array(indices).T  # 转换坐标格式
t1_values = t1_data[indices]  # 提取T1图像中掩膜区域的值

# 保存坐标和值到文本文件
np.savetxt(r'F:\Codes\AHBA\data\voxel_coordinates.txt', coords, fmt='%d')
np.savetxt(r'F:\Codes\AHBA\data\brain_map.txt', t1_values, fmt='%f')

print('Files have been saved: voxel_coordinates.txt and brain_map_values.txt')

#%%
from brainsmash.workbench.geo import volume

coord_file = r'F:\Codes\AHBA\data\voxel_coordinates.txt'
output_dir = r'F:\Codes\AHBA\data'

filenames = volume(coord_file, output_dir)