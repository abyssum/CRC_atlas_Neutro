# load("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
# filename <- file.choose("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
# mouse.seurat <- readRDS("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")
library(Seurat)
library(ggplot2)
library(sceasy)
library(reticulate)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library("stringr")

# The following needs a python virtual env containing anndata package
# You can create this by running $ ~/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user anndata
# and then sourcing the virtual env with reticulate
h5ad_file <- "/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.h5ad"
mouse.seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

# Plot - sanity check
DimPlot(mouse.seurat, group.by = c("tissue","condition"))

# Normalize data
mouse.seurat <- NormalizeData(mouse.seurat)

# DE analysis - Blood
mouse.seurat2 <- mouse.seurat

Idents(mouse.seurat2) <- mouse.seurat2$condition

mouse.neutro.Blood_de <- FindMarkers(mouse.seurat2, ident.1 = "blood_tumor", ident.2 = "blood_wt", verbose = FALSE)
head(mouse.neutro.Blood_de, n = 10)

mouse.neutro.Blood_de_MAST <- FindMarkers(mouse.seurat2, ident.1 = "blood_tumor", ident.2 = "blood_wt", test.use = "MAST")
head(mouse.neutro.Blood_de_MAST, n = 10)


# Volcano plot
EnhancedVolcano(mouse.neutro.Blood_de_MAST,
                lab = rownames(mouse.neutro.Blood_de_MAST),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Blood tumor vs. Blood healthy - Zurich mice dataset',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

# KEGG pathway
mouse.neutro.Blood_de_MAST_sig <- filter(mouse.neutro.Blood_de_MAST, p_val<=0.05)

geneList_blood <- c(mouse.neutro.Blood_de_MAST_sig$avg_log2FC)
genes_blood <- c(row.names(mouse.neutro.Blood_de_MAST_sig))
geneListNames_blood <- mapIds(org.Mm.eg.db, genes_blood, 'ENTREZID', 'SYMBOL')
names(geneList_blood) <- geneListNames_blood
geneList_blood <- sort(geneList_blood, decreasing = TRUE)

geneList_blood <- na.omit(geneList_blood)

kegg.res_blood <- gseKEGG(geneList_blood,
                    organism = "mmu",
                    keyType = "kegg",
                    # exponent = 1,
                    # minGSSize = 10,
                    # maxGSSize = 500,
                    # eps = 1e-10,
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    verbose = TRUE,
                    use_internal_data = FALSE,
                    seed = FALSE,
                    by = "fgsea"
)

kegg.res_blood <- setReadable(kegg.res_blood, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

dotplot(kegg.res_blood, showCategory=30)

kegg.res_blood@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", 
                                          kegg.res_blood@result$Description, fixed = T)

dotplot(kegg.res_blood, showCategory=30, split=".sign") +
  facet_grid(.~.sign) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))


# DE analysis - BM
mouse.neutro.BM_de <- FindMarkers(mouse.seurat2, ident.1 = "BM_tumor", ident.2 = "BM_wt", verbose = FALSE)
head(mouse.neutro.BM_de, n = 10)
  
mouse.neutro.BM_de_MAST <- FindMarkers(mouse.seurat2, ident.1 = "BM_tumor", ident.2 = "BM_wt", test.use = "MAST")
head(mouse.neutro.BM_de_MAST, n = 10)
  
  
# Volcano plots
EnhancedVolcano(mouse.neutro.BM_de_MAST,
                  lab = rownames(mouse.neutro.BM_de_MAST),
                  x = 'avg_log2FC',
                  y = 'p_val',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  title = 'Bone marrow tumor vs. Bone marrow healthy - Zurich mice dataset',
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  pointSize = 4.0,
                  labSize = 6.0,
                  colAlpha = 1,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
  
# KEGG pathway
mouse.neutro.BM_de_MAST_sig <- filter(mouse.neutro.BM_de_MAST, p_val<=0.05)
  
geneList_bm <- c(mouse.neutro.BM_de_MAST_sig$avg_log2FC)
genes_bm <- c(row.names(mouse.neutro.BM_de_MAST_sig))
geneListNames_bm <- mapIds(org.Mm.eg.db, genes_bm, 'ENTREZID', 'SYMBOL')
names(geneList_bm) <- geneListNames_bm
geneList_bm <- sort(geneList_bm, decreasing = TRUE)
  
geneList_bm <- na.omit(geneList_bm)
  
kegg.res_bm <- gseKEGG(geneList_bm,
                      organism = "mmu",
                      keyType = "kegg",
                      # exponent = 1,
                      # minGSSize = 10,
                      # maxGSSize = 500,
                      # eps = 1e-10,
                      pvalueCutoff = 0.5,
                      pAdjustMethod = "BH",
                      verbose = TRUE,
                      use_internal_data = FALSE,
                      seed = FALSE,
                      by = "fgsea"
  )
  
kegg.res_bm <- setReadable(kegg.res_bm, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

kegg.res_bm@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                    replacement = "", kegg.res_bm@result$Description, fixed = T)

dotplot(kegg.res_bm, showCategory=30)

plot <- c("Cell cycle", "DNA replication",
          "p53 signaling pathway", "Systemic lupus erythematosus", 
          "Nucleotide excision repair", "Mismatch repair", 
          "Base excision repair", "Motor proteins",
          "Fanconi anemia pathway", "Homologous recombination",
          "Platinum drug resistance", "Nucleotide metabolism",
          "Viral carcinogenesis", "Neutrophil extracellular trap formation",
          "Cellular senescence", "Pyrimidine metabolism", 
          "SNARE interactions in vesicular transport",
          "ATP-dependent chromatin remodeling", "Ribosome", 
          "Arginine and proline metabolism", "Glutathione metabolism", "Purine metabolism",
          "FoxO signaling pathway","Colorectal cancer",
          "Polycomb repressive complex", "RNA degradation",
          "Glycine, serine and threonine metabolism", "Nucleocytoplasmic transport",
          "Intestinal immune network for IgA production", "Cushing syndrome", 
          "Apoptosis - multiple species", "Citrate cycle (TCA cycle)",
          "Non-homologous end-joining", "Gap junction",
          "Cysteine and methionine metabolism")

dotplot(kegg.res_bm, showCategory=plot, split=".sign") +
    facet_grid(.~.sign) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=80))

# Do not use - initial
plot <- c("Cell cycle", "DNA replication", 
          "p53 signaling pathway", "Base excision repair", 
          "Mismatch repair", "Nucleotide excision repair", 
          "Systemic lupus erythematosus", 
          "SNARE interactions in vesicular transport", 
          "Platinum drug resistance", "Nucleotide metabolism", 
          "Motor proteins", "Homologous recombination", 
          "Fanconi anemia pathway", "Small cell lung cancer", 
          "Neutrophil extracellular trap formation", 
          "Human T-cell leukemia virus 1 infection", 
          "Drug metabolism - other enzymes", 
          "Viral carcinogenesis", "Glutathione metabolism", 
          "Pyrimidine metabolism", 
          "Endocrine and other factor-regulated calcium reabsorption")

plot <- c(
  "Cell cycle", "DNA replication", 
  "p53 signaling pathway", "Nucleotide excision repair", 
  "Mismatch repair", "Base excision repair", 
  "Systemic lupus erythematosus", "Motor proteins", 
  "Fanconi anemia pathway", "Platinum drug resistance", 
  "Homologous recombination", "Nucleotide metabolism", 
  "Human T-cell leukemia virus 1 infection", "Viral carcinogenesis", 
  "Neutrophil extracellular trap formation", 
  "Ribosome biogenesis in eukaryotes", "Cellular senescence", 
  "Pyrimidine metabolism", "Small cell lung cancer", 
  "SNARE interactions in vesicular transport", 
  "Drug metabolism - other enzymes", 
  "Ribosome", "Autophagy - other", 
  "ATP-dependent chromatin remodeling", "Mineral absorption", 
  "Arginine and proline metabolism", "Purine metabolism", 
  "Glutathione metabolism", "FoxO signaling pathway", 
  "Antifolate resistance", "Colorectal cancer", 
  "Polycomb repressive complex", "Cocaine addiction", 
  "Pathways in cancer", "B cell receptor signaling pathway", 
  "RNA degradation", "Nucleocytoplasmic transport", 
  "Endocrine and other factor-regulated calcium reabsorption", 
  "Cushing syndrome", "Rheumatoid arthritis", 
  "Carbohydrate digestion and absorption", 
  "Apoptosis - multiple species", "Citrate cycle (TCA cycle)", 
  "Intestinal immune network for IgA production", 
  "Non-homologous end-joining"
)