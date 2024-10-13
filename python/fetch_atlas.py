#%% fetch atlas (2mm)
# 只有2mm数据才能计算出模拟数据，1mm目前无法算出
def fetch_atlas_data(atlas='dk'):
    if atlas.lower() == 'aal':
        # AAL atlas
        atlas_1mm = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4_1mm_MNI152.nii'  # 1mm resolution
        atlas_2mm = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4.nii'
        atlas_info = r'..\data\01_Atlas\AAL_v4\ROI_MNI_V4.csv'

    elif atlas.lower() == 'dk':
        # DK atlas
        atlas_1mm = r'..\data\01_Atlas\desikan_killiany\atlas-desikankilliany_1mm_MNI152.nii'
        atlas_2mm = r'..\data\01_Atlas\desikan_killiany\atlas-desikankilliany_2mm_MNI152.nii'
        atlas_info = r'..\data\01_Atlas\desikan_killiany\atlas-desikankilliany.csv'

    elif atlas.lower() == 'allen_sub':
        # Allen Human Reference Atlas
        atlas_1mm = r'..\data\01_Atlas\allen_human_reference_atlas_subcortical\allen_human_reference_atlas_subcortical_1mm.nii'
        atlas_2mm = r'..\data\01_Atlas\allen_human_reference_atlas_subcortical\allen_human_reference_atlas_subcortical_2mm.nii'
        atlas_info = r'..\data\01_Atlas\allen_human_reference_atlas_subcortical\allen_human_reference_atlas_subcortical.csv'

    elif atlas.lower() == 'subcortical_nuclei':
        # Subcortical nuclei atlas
        atlas_1mm = r'..\data\01_Atlas\subcortical_nuclei\subcortical_nuclei.nii'
        atlas_2mm = r'..\data\01_Atlas\subcortical_nuclei\subcortical_nuclei_2mm.nii'
        atlas_info = r'..\data\01_Atlas\subcortical_nuclei\subcortical_nuclei.csv'

    else:
        raise ValueError("Atlas must be either 'aal' or 'dk' or 'allen_sub' or 'subcortical_nuclei'!")

    return atlas_1mm, atlas_2mm, atlas_info

