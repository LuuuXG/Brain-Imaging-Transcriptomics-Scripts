import nibabel as nib
import pandas as pd
import os

# 路径配置
nii_path = r'E:\Neuroimage\QSM_test\aparc+aseg.mni152.nii.gz'
lut_path = r'E:\Neuroimage\QSM_test\FreeSurferColorLUT.txt'
output_excel_path = r'E:\Neuroimage\QSM_test\brain_regions.xlsx'

# 读取NIfTI文件，提取所有唯一标签值
nii_data = nib.load(nii_path).get_fdata()
unique_labels = sorted(list(set(nii_data.flat)))

# 读取颜色查找表，建立标签值到脑区名称的映射
label_to_region = {}
with open(lut_path, 'r') as lut_file:
    for line in lut_file:
        if not line.startswith('#') and line.strip():
            parts = line.split()
            label = int(parts[0])
            region_name = ' '.join(parts[1:-1])
            label_to_region[label] = region_name

# 提取脑区名称
regions = [label_to_region.get(label, 'Unknown') for label in unique_labels]

# 保存到Excel
df = pd.DataFrame({'Label': unique_labels, 'Region Name': regions})
df.to_excel(output_excel_path, index=False)

print(f'Done. The brain regions and their labels have been saved to {output_excel_path}')
