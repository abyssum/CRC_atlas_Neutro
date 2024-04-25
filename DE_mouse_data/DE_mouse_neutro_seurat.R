# load("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
# filename <- file.choose("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
# mouse.seurat <- readRDS("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
library(Seurat)
library(sceasy)
library(reticulate)

# The following needs a python virtual env containing anndata package
# You can create this by running $ ~/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user anndata
# and then sourcing the virtual env with reticulate

h5ad_file <- "/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.h5ad"
mouse.seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")




