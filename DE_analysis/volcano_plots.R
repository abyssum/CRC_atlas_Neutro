library(EnhancedVolcano)

data_dir <- "/home/fotakis/myScratch/CRC_atlas_Neutro_back/DE_analysis/tables/"

sample_res <- "MUI-tumor_vs_blood-DESeq2_result.tsv"

res <- read.table(paste0(data_dir, sample_res),
                sep = "\t",
                header = T,
                row.names = 1)
# rownames(res) <- res$gene_id

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Tumor vs. Whole Blood - MUI dataset',
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

