#%% create permuted imaging data
from neuromaps.images import load_nifti, load_gifti, relabel_gifti
from neuromaps.parcellate import Parcellater
from neuromaps import nulls, datasets, resampling, transforms
from abagen import fetch_desikan_killiany
from neuromaps.datasets import fetch_annotation

# for test
satterthwaite = datasets.fetch_annotation(source='satterthwaite2014')
neurosynth = datasets.fetch_annotation(source='neurosynth')
abagen = datasets.fetch_annotation(source='abagen')

# imaging for analysis
#imaging_path = r'E:\Codes\Imaging-Transcriptomics\data\HC45_mni_qsm_smooth_3mm.nii'
imaging_path = r'E:\Codes\Imaging-Transcriptomics\data\atlas-desikankilliany_1mm_MNI152.nii'
imaging = load_nifti(imaging_path)

# Test: DK atlas
# MNI-DK-1: from Imaging Transcriptomics Scripts
parcellation_mni = r'E:\Codes\Imaging-Transcriptomics\data\atlas-desikankilliany_1mm_MNI152.nii'

parcellater_mni = Parcellater(parcellation_mni, 'mni152')
imaging_parc = parcellater_mni.fit_transform(imaging, space='mni152')

#nulls = nulls.burt2020(imaging_parc, atlas='MNI152', density='2mm', n_perm=10, seed=1234, parcellation=parcellation_mni)

#print(nulls.shape)

# MNI-DK-2: from abagen
DK_mni = fetch_desikan_killiany()['image']
DK_mni_parcellater = Parcellater(DK_mni, 'mni152')
imaging_mni_parc_2 = DK_mni_parcellater.fit_transform(imaging, space='mni152')

# surface-DK-1: transform imaging to fsaverage
imaging_surface = transforms.mni152_to_fsaverage(imaging, '10k')

imaging_surface_lh, imaging_surface_rh = imaging_surface

# fetch DK surface
DK_surface = fetch_desikan_killiany(surface=True)

DK_surface_lh = load_gifti(DK_surface['image'][0])
DK_surface_rh = load_gifti(DK_surface['image'][1])

parcellation = relabel_gifti((DK_surface_lh, DK_surface_rh), background=['Medial_wall'])

DK = Parcellater(parcellation, 'fsaverage').fit()

imaging_surface_parc = DK.transform(imaging_surface, 'fsaverage')

# surface-DK-2: resample imaging
imaging_surface_2, DK_surface = resampling.resample_images(imaging_path, DK_surface['image'], 'MNI152', 'fsaverage')

imaging_surface_parc_2 = DK.transform(imaging_surface_2, 'fsaverage')

# surface-DK-3: try different space
imaging_surface_3 = transforms.mni152_to_fslr(imaging, '32k')

DK_surface = fetch_desikan_killiany(surface=True)
DK_surface_3 = transforms.fsaverage_to_fslr(DK_surface['image'], '32k')

abagen = fetch_annotation(source='abagen')
fslr = transforms.fsaverage_to_fslr(abagen, '32k')

#%% compare different parcellation
import numpy as np
import matplotlib.pyplot as plt

def plot_data(*args, labels=None, title='Comparison of Segmentation Methods', xlabel='Region Index', ylabel='Values'):
    # 计算数据序列的数量
    num_datasets = len(args)
    if labels is None:
        labels = [f'Dataset {i + 1}' for i in range(num_datasets)]

    # 确保每个数据都被扁平化并截取到前34个值
    data = [np.array(dataset).flatten()[:34] for dataset in args]

    # 位置数组，根据数据集数量动态调整位置和宽度
    positions = np.arange(1, 35)  # 1到34
    width = 0.8 / num_datasets  # 动态调整宽度以适应所有柱子

    # 创建图形和轴
    fig, ax = plt.subplots(figsize=(14, 8))

    # 绘制柱状图
    for i, dataset in enumerate(data):
        ax.bar(positions + (i - num_datasets / 2) * width, dataset, width=width, label=labels[i], align='center')

    # 添加标题和标签
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(positions)
    ax.legend()

    # 显示图形
    plt.show()


# 使用示例
plot_data(imaging_parc, imaging_mni_parc_2, imaging_surface_parc, imaging_surface_parc_2,
          labels=['imaging_parc', 'imaging_mni_parc_2', 'imaging_surface_parc', 'imaging_surface_parc_2'])
