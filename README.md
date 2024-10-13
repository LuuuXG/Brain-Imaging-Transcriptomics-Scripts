**The README is still being built**

# Brain Imaging Transcriptomics Scripts

This is a script for **human** imaging transcriptomics analysis using Allen Human Brain Atlas (AHBA). The script is used in `R` and some `python` scripts are called via the `reticulate` package.

The codes were tested on Windows 10 system with `R` 4.4.1 and `python` 3.9.7.

## Download

The preprocessed data (brain atlas, RNA-seq data etc.) can be downloaded from [here](https://drive.google.com/file/d/1HmgRQwb_udEsU32j-6DthM5YJzrouFTp/view?usp=sharing). Users should put the `data` folder in the root path of this script.

`python` packages requirement:
- [pyls](https://pyls.readthedocs.io/en/latest/) (should be installed from github)
- [neuromaps](https://netneurolab.github.io/neuromaps/installation.html) (should be installed from github; [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) is needed)
- [brainsmash](https://brainsmash.readthedocs.io/en/latest/installation.html) ([Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) is needed)
- [abagen](https://abagen.readthedocs.io/en/stable/installation.html)

`R` packages needed can be found in `./R/base/load_packages.R` for main analysis and in corresponding scripts in `./R/Plots` for visualization.

## Usage

What you need is to prepare a single 3D brain imaging file in Nifti format (`.nii` or `.nii.gz`). It can either be a statistic map or other imaging phenotypes. You can click and open the `main.R` file in the `./R` folder to perform an imaging transcriptomics analysis starting with the acquisition of the regional gene expression matrix.

## Possible Issues

1. `reticulate` crash or similar probelms: [py_run_file_impl() crashing since .v1.27](https://github.com/rstudio/reticulate/issues/1422): Reinstall `R` might help.
2. `Error in serialize(data, node$con) : error writing to connection`: This may be related to insufficient memory, especially for GSEA analysis of multiple PLS components.

## Similar Works

1. [Tight fitting genes: finding relations between statistical maps and gene expression patterns](https://www.semanticscholar.org/paper/Tight-fitting-genes%3A-finding-relations-between-maps-Gorgolewski-Fox/ae3a76f517bc921664d717a5032b3dc75fe9d7ae?utm_source=direct_link&sort=relevance&page=5) (repo: https://github.com/chrisgorgo/alleninf)
2. [Integrating neuroimaging and gene expression data using the imaging transcriptomics toolbox](https://linkinghub.elsevier.com/retrieve/pii/S2666166722001952) (repo: https://github.com/alegiac95/Imaging-transcriptomics)

**Highlight:** Compared with these toolboxes, our processing workflow allows the use of any user-provided brain atlas/parcellations in the MNI space, and offer a variaty of statistic methods to link spatial gene expression with brain imaging such as partial least squares regression (PLSR), linear regression, and weighted gene coexpression network analysis (WGCNA). At the same time, we also used a spatial null model to increase the statistical power.