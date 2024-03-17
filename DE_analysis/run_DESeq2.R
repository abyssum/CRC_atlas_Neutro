#!/usr/bin/env Rscript
'
Usage:
DESeq2_DGEA.R --count_mat=<count_mat> --colData=<colData> --rowData=<rowData> --prefix=<prefix> --outdir=<outdir> --sample_col=<sample_col> --cond_col=<cond_col> --covariate_formula=<covariate_formula> --c1=<c1> --c2=<c2> --cpus=<cpus> [options]

Mandatory arguments:
  --count_mat=<count_mat>
  --colData=<colData>
  --rowData=<rowData>
  --prefix=<prefix>
  --outdir=<outdir>
  --sample_col=<sample_col>
  --cond_col=<cond_col>
  --covariate_formula=<covariate_formula>
  --c1=<c1>
  --c2=<c2>
  --cpus=<cpus>
' -> doc

# load required packages
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

suppressPackageStartupMessages({
library(BiocParallel)
library(conflicted)
library(readr)
library(tibble)
library(dplyr)
library(stringr)
library(forcats)
library(DESeq2)
library(IHW)
library(limma)
})

# Load parameters
prefix <- arguments$prefix
count_mat <- read_csv(arguments$count_mat)
colData <- read_csv(arguments$colData)
rowData <- read_csv(arguments$rowData)
outdir <- arguments$outdir
sample_col <- arguments$sample_col
cond_col <- arguments$cond_col
covariate_formula <- arguments$covariate_formula
c1 <- arguments$c1
c2 <- arguments$c2
n_cpus <- arguments$cpus


# design_formula <- as.formula(paste0("~", cond_col))

design_formula <- as.formula(paste0("~", covariate_formula, " ", cond_col))


register(MulticoreParam(workers = n_cpus))


dds <- DESeqDataSetFromMatrix(
  countData = count_mat |> column_to_rownames(var = "gene_id") |> ceiling(),
  colData = colData |> column_to_rownames(var = sample_col),
  rowData = rowData,
  design = design_formula
)
# define reference level (not really necessary when uisng contrasts)
dds[[cond_col]] <- relevel(dds[[cond_col]], ref = c2)

## keep only genes where we have >= 10 reads per samplecondition in at least 2 samples
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 2
dds <- dds[keep, ]
#rownames(dds) <- rowData(dds)$GeneSymbol

# save normalized filtered count file
norm_mat <- counts(dds, normalized = TRUE) |> as_tibble(rownames = "gene_id")
write_tsv(norm_mat, file.path(outdir, paste0(prefix, "_NormalizedCounts.tsv")))

# save normalized batch corrected filtered count file
vst <- vst(dds, blind = FALSE)
batch <- gsub("\\+", "", covariate_formula) |> str_squish()
assay(vst) <- limma::removeBatchEffect(x = assay(vst), batch = vst[[batch]])
write_tsv(assay(vst) |> as_tibble(rownames = "gene_id"), file.path(outdir, paste0(prefix, "_vst_batch_corrected_NormalizedCounts.tsv")))

# run DESeq
dds <- DESeq(dds, parallel = (n_cpus > 1))

# set names of contrasts
contrasts <- list(c(cond_col, c1, c2))
names(contrasts) <- sprintf("%s_vs_%s", c1, c2)

## IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- lapply(names(contrasts), function(name) {
  contrast <- contrasts[[name]]
  results(dds, filterFun = ihw, contrast = contrast) |>
    as_tibble(rownames = "gene_id") |>
    mutate(comparison = name) |>
    arrange(pvalue)
}) |> bind_rows()

write_tsv(resIHW, file.path(outdir, paste0(prefix, "-", names(contrasts), "-DESeq2_result.tsv")))