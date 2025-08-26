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
p <- add_argument(p, "--outdir", help = "output root dir", default = "../../results/disease_specific_changes/dea")
argv <- parse_args(p)

input  <- argv$input
outdir <- argv$outdir

# load counts and make integers
sce <- read_h5ad(input, to = "SingleCellExperiment")
assays(sce)$counts <- round(assays(sce)$X)

# factors
cd <- as.data.frame(colData(sce))
cd$condition <- factor(cd$condition, levels = c("HC", "UC", "CD"))
colData(sce) <- S4Vectors::DataFrame(cd)

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
comparisons <- list(c("CD","UC"), c("HC","CD"), c("HC","UC"))
celltypes <- unique(colData(sce)$annotation.coarse)

for (ct in celltypes) {
  # subset only by cell type; keep all conditions
  sce_ct <- sce[, colData(sce)$annotation.coarse == ct, drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = assays(sce_ct)$counts,
    colData   = as.data.frame(colData(sce_ct)),
    design    = ~ dataset + condition
  )
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds)

  for (cmp in comparisons) {
    cond1 <- cmp[[1]]; cond2 <- cmp[[2]]
    res <- results(dds, contrast = c("condition", cond1, cond2))

    # convert to data.frame and add gene column
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

    # add significance column
    res_df$significant <- ifelse(is.na(res_df$padj), NA, res_df$padj < 0.05)

    dir <- file.path(outdir, paste0(sanitize(ct), "_", cond1, "vs", cond2))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    write.table(res_df,
      file = file.path(dir, "results.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}
