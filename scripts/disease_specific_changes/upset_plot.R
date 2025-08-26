library(ComplexUpset)
library(ggplot2)
library(dplyr)

conditions <- c("CDvsUC", "HCvsCD", "HCvsUC")

dirs <- list.dirs(path = "../../results/disease_specific_changes/dea", recursive = FALSE, full.names = TRUE)

gene_sets <- list()

for (dir in dirs) {
  dir_name <- basename(dir)
  parts <- strsplit(dir_name, "_")[[1]]
  if (length(parts) < 2) next
  
  cond <- tail(parts, 1)
  cell_type <- paste(head(parts, length(parts) - 1), collapse = "_")
  if (!(cond %in% conditions)) next
  
  file_path <- file.path(dir, "results.tsv")
  if (!file.exists(file_path)) next
  
  res <- read.delim(file_path, header = TRUE, stringsAsFactor = FALSE)
  
  if (!("padj" %in% colnames(res)) || !("log2FoldChange" %in% colnames(res))) {
    warning(paste("File", file_path, "missing required columns. Skipping."))
    next
  }
  
  res_filt <- res %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
  
  if ("gene" %in% colnames(res_filt)) {
    genes <- unique(res_filt$gene)
  } else {
    genes <- rownames(res_filt)
  }
  
  if (!cond %in% names(gene_sets)) {
    gene_sets[[cond]] <- list()
  }
  
  gene_sets[[cond]][[cell_type]] <- genes
}

create_membership_df <- function(gene_list) {
  all_genes <- unique(unlist(gene_list))
  df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  for (ct in names(gene_list)) {
    df[[ct]] <- df$gene %in% gene_list[[ct]]
  }
  return(df)
}

for (cond in conditions) {
  if (!cond %in% names(gene_sets)) {
    message(paste("No data for condition", cond, "- skipping plot."))
    next
  }
  
  membership_df <- create_membership_df(gene_sets[[cond]])
  
  p <- upset(membership_df,
             intersect = setdiff(names(membership_df), "gene"),
             set_sizes = upset_set_size() + theme(axis.text = element_text(size=15)),
             base_annotations = list('Intersection size' = intersection_size(
               text = list(
                 size=5
               )
             ) +
               theme(axis.title = element_text(size=15)))
             ) +
    ggtitle(paste("UpSet Plot for", cond, "\n(overlap of significantly expressed genes p < 0.05 and |l2FC| > 1)")) +
    theme(text = element_text(size = 15),
          title = element_text(size=12),
          legend.text = element_text(size = 12),     # Legend labels
          legend.title = element_text(size = 14),
          axis.text = element_text(size=15))
  ggsave(filename = paste0("../../results/disease_specific_changes/dea/upset_", cond, ".png"), plot = p)
}
