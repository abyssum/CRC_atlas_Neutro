library(Seurat)
library(ggplot2)
library(sceasy)
library(reticulate)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)

# Do not run
# Returns error reading from connection - possible file corruption (talk with Kris)
# mouse.seurat <- readRDS("/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.rds")

# Use the provided H5AD file.
# The following needs a python virtual env containing anndata package
# You can create this by running $ ~/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user anndata
# and then sourcing the virtual env with reticulate
h5ad_file <- "/data/projects/2022/CRCA/data/own_datasets/arnold_lab_mouse/tumor_bm_blood_neutrophils.h5ad"
mouse.seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

# To be used later for BM clusters comparison
mouse.seurat3 <- mouse.seurat

# Plot - sanity check
DimPlot(mouse.seurat, group.by = c("tissue","condition"))

# DE analysis - Blood
mouse.seurat2 <- mouse.seurat

# Normalize data
mouse.seurat2 <- NormalizeData(mouse.seurat2)

# all.genes <- rownames(mouse.seurat2)
# mouse.seurat2 <- ScaleData(mouse.seurat2, features = all.genes)

Idents(mouse.seurat2) <- mouse.seurat2$condition

# Do not run -> used for testing/comparing results
# mouse.neutro.Blood_de <- FindMarkers(mouse.seurat2, ident.1 = "blood_tumor", ident.2 = "blood_wt", verbose = FALSE)
# head(mouse.neutro.Blood_de, n = 10)

# Perform DE with MAST
mouse.neutro.Blood_de_MAST <- FindMarkers(mouse.seurat2, ident.1 = "blood_tumor", ident.2 = "blood_wt", test.use = "MAST")
# Check top 10
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

# Export filtered DEGs to CSV
write.csv(mouse.neutro.Blood_de_MAST_sig, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/DE_mouse_blood.csv")

# Create the named genes list
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

# Set genes to symbols (instead of ENTREZID)
kegg.res_blood <- setReadable(kegg.res_blood, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

# Plot the top 30 terms
dotplot(kegg.res_blood, showCategory=30)

# Remove the " - Mus musculus (house mouse)" from results
kegg.res_blood@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", 
                                          kegg.res_blood@result$Description, fixed = T)

# Plot the top up-/down-regulated terms
dotplot(kegg.res_blood, showCategory=30, split=".sign") +
  facet_grid(.~.sign) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

# Esport the results of the GSEA to CSV
write.csv(kegg.res_blood@result, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/GSEA_mouse_blood.csv")

# DE analysis - BM
# Do not run -> used for testing/comparing results
# mouse.neutro.BM_de <- FindMarkers(mouse.seurat2, ident.1 = "BM_tumor", ident.2 = "BM_wt", verbose = FALSE)
# head(mouse.neutro.BM_de, n = 10)

# Perform DE with MAST
mouse.neutro.BM_de_MAST <- FindMarkers(mouse.seurat2, ident.1 = "BM_tumor", ident.2 = "BM_wt", test.use = "MAST")
# Check top 10
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
# Filter for p-val
mouse.neutro.BM_de_MAST_sig <- filter(mouse.neutro.BM_de_MAST, p_val<=0.05)

# Export filtered DEGs to CSV
write.csv(mouse.neutro.BM_de_MAST_sig, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/DE_mouse_bm.csv")

# Create the named genes list
geneList_bm <- c(mouse.neutro.BM_de_MAST_sig$avg_log2FC)
genes_bm <- c(row.names(mouse.neutro.BM_de_MAST_sig))
geneListNames_bm <- mapIds(org.Mm.eg.db, genes_bm, 'ENTREZID', 'SYMBOL')
names(geneList_bm) <- geneListNames_bm
geneList_bm <- sort(geneList_bm, decreasing = TRUE)
  
geneList_bm <- na.omit(geneList_bm)

# Run GSEA with KEGG terms
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

# Set genes to symbols (instead of ENTREZID)
kegg.res_bm <- setReadable(kegg.res_bm, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

# Remove the " - Mus musculus (house mouse)" from GSEA results
kegg.res_bm@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                    replacement = "", kegg.res_bm@result$Description, fixed = T)

# Plot the top 30 terms
dotplot(kegg.res_bm, showCategory=30)

# Chose specific terms to be ploted - requires to check the results
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

# Plot the chosen KEGG terms
dotplot(kegg.res_bm, showCategory=plot, split=".sign") +
    facet_grid(.~.sign) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=80))

# Export the results of the GSEA to CSV
write.csv(kegg.res_bm@result, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/GSEA_mouse_bm.csv")


# DE between BM 10 vs. 6,8,9 clusters
# subset - only the BM clusters
mouse.seurat.bm <- subset(x = mouse.seurat3, subset = tissue == "BM")

# Further subset for the specific BM subclusters
mouse.seurat.bm2 <- subset(x = mouse.seurat.bm, subset = seurat_clusters == c(10,6,8,9))

# Normalize data - Note: the initial seurat object used for subseting is not normalised
mouse.seurat.bm2 <- NormalizeData(mouse.seurat.bm2)

# Set the idents to the seurat clusters - for the DE analysis
Idents(mouse.seurat.bm2) <- mouse.seurat.bm2$seurat_clusters

# Do not run - just for sanity check
# DimPlot(mouse.seurat.bm2, group.by = c("seurat_clusters"), label = T)

# Perform the DE analysis using MAST
mouse.neutro.BM_sub_MAST <- FindMarkers(mouse.seurat.bm2, ident.1 = '10', test.use = "MAST")

# Since there are several genes with reported 0 p-val -> makes it difficult to plot the (log10) y axis (will result in -Inf)
# Need to replace the 0 with p-vals lower than the lowest reported p-val from the rest of the genes.
# This is the lowest (not 0) reported p-val 
new_pval <- 3.952525e-323
# this loop will iterate through the indexes of all rows with p_val < 1.976263e-323 and replace them with lower p-vals
# this way we avoid duplicate p-values while ensuring that the genes remain the top DE genes
for (i in c(which(mouse.neutro.BM_sub_MAST$p_val < 1.976263e-323))){
  mouse.neutro.BM_sub_MAST[i,]$p_val <- new_pval
  # print(new_pval)
  new_pval <- new_pval + 1.976263e-323
  print(new_pval)
}

# Volcano plot
EnhancedVolcano(mouse.neutro.BM_sub_MAST,
                lab = rownames(mouse.neutro.BM_sub_MAST),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                # ylab = bquote(~-Log[10] ~ italic(P)),
                # ylim = c(0, max(-log10(mouse.neutro.BM_sub_sig$p_val), na.rm = TRUE) + 5),
                title = 'BM cluster 10 vs. BM clusters 6,8,9 - Zurich mice dataset',
                pCutoff = 10e-4,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = F,
                widthConnectors = 0.75)

# KEGG pathway
# Filter for p-val adjusted (due to the sheer amount of DE genes with p-val < 0.05)
mouse.neutro.BM_sub_sig <- filter(mouse.neutro.BM_sub_MAST, p_val_adj<=0.05)

# Export filtered DEGs to CSV
write.csv(mouse.neutro.BM_sub_sig, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/DE_mouse_bm_sub.csv")

# Create the named genes list
geneList_bm_sub <- c(mouse.neutro.BM_sub_sig$avg_log2FC)
genes_bm_sub <- c(row.names(mouse.neutro.BM_sub_sig))
geneListNames_bm_sub <- mapIds(org.Mm.eg.db, genes_bm_sub, 'ENTREZID', 'SYMBOL')
names(geneList_bm_sub) <- geneListNames_bm_sub
geneList_bm_sub <- sort(geneList_bm_sub, decreasing = TRUE)

geneList_bm_sub <- na.omit(geneList_bm_sub)

# Perform the GSEA analysis using KEGG terms
kegg.res_bm_sub <- gseKEGG(geneList_bm_sub,
                           organism = "mmu",
                           keyType = "kegg",
                           # exponent = 1,
                           # minGSSize = 10,
                           # maxGSSize = 500,
                           # eps = 1e-10,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = TRUE,
                           use_internal_data = FALSE,
                           seed = FALSE,
                           by = "fgsea"
)

# Set genes to symbols (instead of ENTREZID)
kegg.res_bm_sub <- setReadable(kegg.res_bm_sub, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

# Remove the " - Mus musculus (house mouse)" from GSEA results
kegg.res_bm_sub@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                           replacement = "", kegg.res_bm_sub@result$Description, fixed = T)

# Plot the top 30 KEGG terms
dotplot(kegg.res_bm_sub, showCategory=30)

# Plot the top 30 up-/done-regulated KEGG terms
dotplot(kegg.res_bm_sub, showCategory=30, split=".sign") +
  facet_grid(.~.sign) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))

# Export the results of the GSEA to CSV
write.csv(kegg.res_bm_sub@result, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/GSEA_mouse_bm_sub.csv")

# DE between BM 10 vs. rest clusters
# subset - only the BM clusters
mouse.seurat.bm_all <- subset(x = mouse.seurat3, subset = tissue == "BM")

# Normalize data - Note: the initial seurat object used for subseting is not normalised
mouse.seurat.bm_all <- NormalizeData(mouse.seurat.bm_all)

# Set the idents to the seurat clusters - for the DE analysis
Idents(mouse.seurat.bm_all) <- mouse.seurat.bm_all$seurat_clusters

# Do not run - just for sanity check
# DimPlot(mouse.seurat.bm2, group.by = c("seurat_clusters"), label = T)

# Perform the DE analysis using MAST
mouse.neutro.BM_sub_all_MAST <- FindMarkers(mouse.seurat.bm_all, ident.1 = '10', test.use = "MAST")

# Since there are several genes with reported 0 p-val -> makes it difficult to plot the (log10) y axis (will result in -Inf)
# Need to replace the 0 with p-vals lower than the lowest reported p-val from the rest of the genes.
# This is the lowest (not 0) reported p-val 
new_pval2 <- 3.952525e-323
# this loop will iterate through the indexes of all rows with p_val < 1.976263e-323 and replace them with lower p-vals
# this way we avoid duplicate p-values while ensuring that the genes remain the top DE genes
for (i in c(which(mouse.neutro.BM_sub_all_MAST$p_val < 1.976263e-323))){
  mouse.neutro.BM_sub_all_MAST[i,]$p_val <- new_pval2
  # print(new_pval)
  new_pval2 <- new_pval2 + 1.976263e-323
  print(new_pval2)
}

# Volcano plot
EnhancedVolcano(mouse.neutro.BM_sub_all_MAST,
                lab = rownames(mouse.neutro.BM_sub_all_MAST),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                # ylab = bquote(~-Log[10] ~ italic(P)),
                # ylim = c(0, max(-log10(mouse.neutro.BM_sub_sig$p_val), na.rm = TRUE) + 5),
                title = 'BM cluster 10 vs. rest - Zurich mice dataset',
                pCutoff = 10e-4,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = F,
                widthConnectors = 0.75)

# KEGG pathway
# Filter for p-val adjusted (due to the sheer amount of DE genes with p-val < 0.05)
mouse.neutro.BM_sub_all_sig <- filter(mouse.neutro.BM_sub_all_MAST, p_val_adj<=0.05)

# Export filtered DEGs to CSV
write.csv(mouse.neutro.BM_sub_all_sig, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/DE_mouse_bm_sub_all.csv")

# Create the named genes list
geneList_bm_sub_all <- c(mouse.neutro.BM_sub_all_sig$avg_log2FC)
genes_bm_sub_all <- c(row.names(mouse.neutro.BM_sub_all_sig))
geneListNames_bm_sub_all <- mapIds(org.Mm.eg.db, genes_bm_sub_all, 'ENTREZID', 'SYMBOL')
names(geneList_bm_sub_all) <- geneListNames_bm_sub_all
geneList_bm_sub_all <- sort(geneList_bm_sub_all, decreasing = TRUE)

geneList_bm_sub_all <- na.omit(geneList_bm_sub_all)

# Perform the GSEA analysis using KEGG terms
kegg.res_bm_sub_all <- gseKEGG(geneList_bm_sub_all,
                           organism = "mmu",
                           keyType = "kegg",
                           # exponent = 1,
                           # minGSSize = 10,
                           # maxGSSize = 500,
                           # eps = 1e-10,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = TRUE,
                           use_internal_data = FALSE,
                           seed = FALSE,
                           by = "fgsea"
)

# Set genes to symbols (instead of ENTREZID)
kegg.res_bm_sub_all <- setReadable(kegg.res_bm_sub_all, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

# Remove the " - Mus musculus (house mouse)" from GSEA results
kegg.res_bm_sub_all@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", 
                                           replacement = "", kegg.res_bm_sub_all@result$Description, fixed = T)

# Plot the top 30 KEGG terms
dotplot(kegg.res_bm_sub_all, showCategory=30)

# Plot the top 30 up-/done-regulated KEGG terms
dotplot(kegg.res_bm_sub_all, showCategory=30, split=".sign") +
  facet_grid(.~.sign) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))

# Export the results of the GSEA to CSV
write.csv(kegg.res_bm_sub_all@result, "/home/fotakis/myScratch/CRC_atlas_Neutro/DE_analysis/tables/GSEA_mouse_bm_sub_all.csv")


DimPlot(mouse.seurat.bm_all, group.by = c("tissue","condition"))
mouse.seurat.bm_all

# All markers
library(dplyr)
mouse.seurat4 <- mouse.seurat

mouse.seurat4 <- NormalizeData(mouse.seurat4)

Idents(mouse.seurat4) <- mouse.seurat4$seurat_clusters

mouse.markers <- FindAllMarkers(mouse.seurat4, only.pos = TRUE)
mouse.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

mouse.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
all.genes <- rownames(mouse.seurat4)
mouse.seurat4 <- ScaleData(mouse.seurat4, features = all.genes)
DoHeatmap(mouse.seurat4, features = top10$gene) + NoLegend()
