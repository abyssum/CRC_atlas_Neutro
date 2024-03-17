#! /bin/bash/

Rscript ./run_DESeq2.R \
--count_mat=./tables/tumor_blood_Neutrophil_count_mat.csv \
--colData=./tables/tumor_blood_Neutrophil_colData.csv \
--rowData=./tables/tumor_blood_Neutrophil_rowData.csv \
--prefix=MUI \
--outdir=./tables/ \
--sample_col=sample \
--cond_col=sample_type \
--covariate_formula='patient_id +' \
--c1=tumor \
--c2=blood \
--cpus=2
