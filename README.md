**The README is still being built**

# Brain Imaging Transcriptomics Scripts

This is a script for imaging transcriptomics analysis using Allen Human Brain Atlas (AHBA). The script is used in `R` and some `python` scripts are called via the `reticulate` package.

## Download

The preprocessed data (brain atlas, RNA-seq data etc.) can be downloaded from [here](https://drive.google.com/file/d/1HmgRQwb_udEsU32j-6DthM5YJzrouFTp/view?usp=sharing). Users should put the `data` folder in the root path of this script.

python packages requirement:
- pyls
- neuromaps
- brainsmash
- nilearn
- abagen
- nibabel

## Usage

You can click and open the `main.R` file in the `./R` folder to perform an imaging transcriptomics analysis starting with the acquisition of gene expression matrix.

## Possible Issues

1. `reticulate` crash or similar probelms: [py_run_file_impl() crashing since .v1.27](https://github.com/rstudio/reticulate/issues/1422)
Reinstall `R` might help.
2. `Error in serialize(data, node$con) : error writing to connection`
This may be related to insufficient memory, especially for GSEA analysis of multiple PLS components.

## Similar Works

1. [Tight fitting genes: finding relations between statistical maps and gene expression patterns](https://www.semanticscholar.org/paper/Tight-fitting-genes%3A-finding-relations-between-maps-Gorgolewski-Fox/ae3a76f517bc921664d717a5032b3dc75fe9d7ae?utm_source=direct_link&sort=relevance&page=5) (repo: https://github.com/chrisgorgo/alleninf)
2. [Integrating neuroimaging and gene expression data using the imaging transcriptomics toolbox](https://linkinghub.elsevier.com/retrieve/pii/S2666166722001952) (repo: https://github.com/alegiac95/Imaging-transcriptomics)