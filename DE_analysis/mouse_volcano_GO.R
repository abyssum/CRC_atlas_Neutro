library(EnhancedVolcano)

data_dir <- "~/myScratch/CRC_atlas_Neutro/DE_analysis/tables/"

sample_res <- "bm_wt_vs_bm_tumor.tsv"

res <- read.table(paste0(data_dir, sample_res),
                  sep = "\t",
                  header = T,
                  row.names = 1)
# rownames(res) <- res$gene_id

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Bone marrow healthy vs. Bone marrow tumor - Zurich mice dataset',
                #pCutoff = 0.5,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

