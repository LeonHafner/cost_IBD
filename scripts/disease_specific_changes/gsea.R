library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(ggplot2)
library(tools)

dea_dir  <- "../../results/disease_specific_changes/dea"
gsea_dir <- "../../results/disease_specific_changes/gsea"

# Find nested results.tsv files
files <- list.files(dea_dir, pattern = "^results\\.tsv$", full.names = TRUE, recursive = TRUE)

# Hallmark gene sets (new arg name)
msigdbr_df <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene  <- msigdbr_df[, c("gs_name", "gene_symbol")]
colnames(term2gene) <- c("term", "gene")

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

for (file in files) {
  cat("Processing:", file, "\n")

  # use the parent directory name, which encodes <celltype>_<cond1>vs<cond2>
  parent_dir <- basename(dirname(file))
  out_dir <- file.path(gsea_dir, parent_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  res <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)

  # drop NA stats
  res <- res[!is.na(res$stat), ]

  # build ranked list; de-dup gene names if needed
  gl <- res$stat
  names(gl) <- res$gene
  keep <- !is.na(gl)
  gl <- gl[keep]

  # if duplicate gene names exist, keep max |stat|
  if (any(duplicated(names(gl)))) {
    o <- order(abs(gl), decreasing = TRUE)
    gl <- gl[o][!duplicated(names(gl[o]))]
  }

  gl <- sort(gl, decreasing = TRUE)

  gseaRes <- tryCatch(
    GSEA(gl, TERM2GENE = term2gene, pvalueCutoff = 0.05, verbose = FALSE, eps = 0),
    error = function(e) { cat("GSEA error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(gseaRes)) next

  # save table
  write.csv(as.data.frame(gseaRes), file = file.path(out_dir, "gsea_results.csv"), row.names = FALSE)

  # save one plot per term
  gsea_df <- as.data.frame(gseaRes)
  if (nrow(gsea_df) == 0) {
    cat("No significant gene sets for", parent_dir, "\n")
    next
  }
  for (i in seq_len(nrow(gsea_df))) {
    term <- gsea_df$ID[i]
    p <- gseaplot2(gseaRes, geneSetID = term, title = term)
    ggsave(file.path(out_dir, paste0(sanitize(term), ".png")), plot = p, width = 8, height = 6)
  }
}
