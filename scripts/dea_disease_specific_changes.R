#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(anndataR)
  library(argparser)
})

set.seed(0)

# args
p <- arg_parser("DESeq2 DEA per cell type")
p <- add_argument(p, "--input", help = "input pseudobulk.h5ad")
p <- add_argument(p, "--outdir", help = "output root dir", default = "../results/disease_specific_changes/dea")
argv <- parse_args(p)

input  <- argv$input
outdir <- argv$outdir

# io
sce <- read_h5ad(input, to = "SingleCellExperiment")
assays(sce)$counts <- round(assays(sce)$X)  # ensure integers for DESeq2

# simple helpers
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# expect these columns to exist in colData(sce):
# dataset, sample, condition, annotation:coarse
cd <- as.data.frame(colData(sce))
cd$condition <- factor(cd$condition, levels = c("HC", "UC", "CD"))
colData(sce) <- S4Vectors::DataFrame(cd)

comparisons <- list(
  c("UC", "HC"),
  c("CD", "HC"),
  c("UC", "CD")
)

celltypes <- unique(colData(sce)$annotation.coarse)

print(paste("Found cell types:", paste(celltypes, collapse = ", ")))

for (ct in celltypes) {
  sce_ct <- sce[, colData(sce)$annotation.coarse == ct, drop = FALSE]

  for (cmp in comparisons) {
    print(paste("Processing", ct, paste(cmp, collapse = " vs ")))
    cond1 <- cmp[[1]]
    cond2 <- cmp[[2]]

    keep <- sce_ct$condition %in% c(cond1, cond2)
    sub  <- sce_ct[, keep, drop = FALSE]
    sub$condition <- factor(sub$condition, levels = c(cond2, cond1))  # base = cond2

    dds <- DESeqDataSetFromMatrix(
      countData = assays(sub)$counts,
      colData   = as.data.frame(colData(sub)),
      design    = ~ dataset + condition
    )
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)

    dir <- file.path(outdir, paste0(sanitize(ct), "_", cond1, "vs", cond2))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    write.table(res, file = file.path(dir, "results.tsv"),
                sep = "\t", quote = FALSE, row.names = TRUE)
  }
}
