# Load packages needed for the analysis.

# It is recommended to install the packages manually as running this script
# might not work as expected.

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
  
# If packages are not installed, install them
options(repos = BiocManager::repositories())

if (!requireNamespace("HDO.db", quietly = TRUE)){
  BiocManager::install("HDO.db", force = TRUE)
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)){
  devtools::install_github("YuLab-SMU/clusterProfiler")
}

if (!requireNamespace("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt", force = TRUE)
}

if (!requireNamespace("ReactomePA", quietly = TRUE)){
  BiocManager::install("ReactomePA", force = TRUE)
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db", force = TRUE)
}

if (!requireNamespace("mixOmics", quietly = TRUE)){
  BiocManager::install("mixOmics", force = TRUE)
}

if (!requireNamespace("impute", quietly = TRUE)){
  BiocManager::install("impute", force = TRUE)
}

if (!requireNamespace("preprocessCore", quietly = TRUE)){
  BiocManager::install("preprocessCore")
}

if (!requireNamespace("simplifyEnrichment", quietly = TRUE)){
  BiocManager::install("simplifyEnrichment")
}

options(repos = c(CRAN = "https://cran.rstudio.com/"))

if (!requireNamespace("shadowtext", quietly = TRUE)){
  install.packages("shadowtext")
}

if (!requireNamespace("msigdbr", quietly = TRUE)){
  install.packages("msigdbr")
}

if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse", force = TRUE)
}

if (!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")
}

if (!requireNamespace("enrichplot", quietly = TRUE)){
  install.packages("enrichplot")
}

if (!requireNamespace("crayon", quietly = TRUE)){
  install.packages("crayon")
}

if (!requireNamespace("reticulate", quietly = TRUE)){
  install.packages("reticulate")
}

if (!requireNamespace("WGCNA", quietly = TRUE)){
  install.packages("WGCNA", force = TRUE)
}

if (!requireNamespace("ggridges", quietly = TRUE)){
  install.packages("ggridges")
}

# load packages
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(msigdbr)
library(crayon)
library(reticulate)
library(WGCNA)
library(simplifyEnrichment)
