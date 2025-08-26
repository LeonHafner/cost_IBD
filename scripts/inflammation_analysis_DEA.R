#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(anndataR)
  library(argparser)
})

set.seed(0)

# args
p <- arg_parser("DESeq2 DEA per cell type × tissue (inflamed vs Non-inflamed)")
p <- add_argument(p, "--input", help = "input pseudobulk .h5ad")
p <- add_argument(p, "--outdir", help = "output root dir", default = "results/region_inflammation_analysis")
argv <- parse_args(p)

input  <- argv$input
outdir <- argv$outdir

# load counts and make integers
sce <- read_h5ad(input, to = "SingleCellExperiment")
assays(sce)$counts <- round(assays(sce)$X)

# factors (expect: dataset, sample, condition, annotation.coarse, tissue, disease.status)
cd <- as.data.frame(colData(sce))
# disease status levels: Non-inflamed as reference, inflamed as treatment
cd$disease.status <- factor(cd$disease.status, levels = c("Non-inflamed", "inflamed"))
colData(sce) <- S4Vectors::DataFrame(cd)

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# iterate per cell type × tissue, fit once, then extract the contrast
celltypes <- unique(colData(sce)$annotation.coarse)
tissues   <- unique(colData(sce)$tissue)

for (ct in celltypes) {
  for (tx in tissues) {
    print(paste("Processing cell type", ct, "and tissue", tx))
    # subset only by cell type and tissue; keep both disease.status groups
    sce_ct_tx <- sce[, colData(sce)$annotation.coarse == ct & colData(sce)$tissue == tx, drop = FALSE]

    if (ncol(sce_ct_tx) == 0) next

    dds <- DESeqDataSetFromMatrix(
      countData = assays(sce_ct_tx)$counts,
      colData   = as.data.frame(colData(sce_ct_tx)),
      design    = ~ dataset + disease.status
    )
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)

    # single comparison: inflamed vs Non-inflamed
    cond1 <- "Inflamed"; cond2 <- "Non-inflamed"
    res <- results(dds, contrast = c("disease.status", cond1, cond2))

    # make a clean table with gene first and a final 'significant' column (padj < 0.05, NA-propagating)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
    res_df$significant <- ifelse(is.na(res_df$padj), NA, res_df$padj < 0.05)

    # write to results/region_inflammation_analysis/<CellType>_<cond1>_vs_<cond2>/results.tsv
    # to avoid collisions across tissues, put each tissue in its own subfolder
    dir <- file.path(outdir, sanitize(tx), paste0(sanitize(ct), "_", cond1, "_vs_", cond2))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    write.table(res_df,
      file = file.path(dir, "results.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}
