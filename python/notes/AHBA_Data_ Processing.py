# use abagen to process AHBA data

#%%
import abagen

## 1. The Allen Human Brain Atlas dataset ##
# download and fetching the AHBA data to F:\Public_dataset\AHBA
files = abagen.fetch_microarray(data_dir=r'F:\Public_dataset\AHBA', donors='all', verbose=1)

# loading the AHBA data
# this can check the donors and the samples
print(files.keys())
print(sorted(files['9861']))
data = files['9861']
annotation = abagen.io.read_annotation(data['annotation'])
print(annotation)

## 2. Defining a parcellation ##
# parcellation means 'atlas'
# in abagen:
#   (1) a NIFTI image in MNI space
#   (2) a tuple of GIFTI images in fsaverage space (and with fsaverage5 resolution!)

# for demonstration purposes: Desikan-Killiany atlas
atlas = abagen.fetch_desikan_killiany()
print(atlas['image'])
print(atlas['info'])

# or the surface version
atlas = abagen.fetch_desikan_killiany(surface=True)
print(atlas['image'])

# or individualized atlas (in donor-native space)
atlas = abagen.fetch_desikan_killiany(native=True)
print(atlas['image'].keys())
print(atlas['image']['9861'])

# and the surface verison of individualized atlas
atlas = abagen.fetch_desikan_killiany(native=True, surface=True)
print(atlas['image'].keys())
print(atlas['image']['9861'])

# User defined atlas
atlas = {
    'image': r'F:\Codes\AHBA\data\AAL3_ROI_MNI_V7_1mm.nii',
    'info': r'F:\Codes\AHBA\data\AAL3_ROI_MNI_V7_1mm.csv'
}
print(atlas['image'])
print(atlas['info'])
#atlas = abagen.images.check_atlas(atlas['image'], atlas['info']);

# non-standard atlas
# check atlas (when define custom atlas?)
# atlas = abagen.images.check_atlas(atlas['image'], atlas['info'])

## 3. Parcellating expression data ##
expression = abagen.get_expression_data(atlas['image'], atlas['info'], data_dir=r'F:\Public_dataset\AHBA')
print(expression)

# if using the default command, there will be two missing rows (right frontal pole (label 72) and right temporal pole (label 73))
print(expression.loc[[72, 73]])

# 2 ways to get dense expression data
expression = abagen.get_expression_data(atlas['image'], atlas['info'], data_dir=r'F:\Public_dataset\AHBA',
                                        missing='centroids')

expression = abagen.get_expression_data(atlas['image'], atlas['info'], data_dir=r'F:\Public_dataset\AHBA',
                                        missing='interpolate')

# save the expression data to csv file
expression.to_csv(r'F:\Public_dataset\AHBA\expression_dense_AAL.csv')