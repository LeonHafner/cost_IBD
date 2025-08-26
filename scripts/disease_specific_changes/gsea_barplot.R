library(data.table)
library(ggplot2)
library(cowplot)


files <- list.files(path = "../../results/disease_specific_changes/gsea", pattern = ".csv", full.names = TRUE, recursive = TRUE)

process_file <- function(file) {
  df <- fread(file)
  df <- df[, c("Description", "enrichmentScore", "p.adjust")]
  df$nlogpadj <- -log10(df$p.adjust)
  df$Description <- gsub("^HALLMARK_", "", df$Description)
  
  p <- ggplot(df, aes(x = reorder(Description, -enrichmentScore), y = enrichmentScore, fill=nlogpadj)) +
        geom_bar(stat="identity") +
        theme_cowplot() +
        labs(x = "", y = "Enrichment Score", fill = "-log10(padj)") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(10, 10, 10, 100)
        )
  
  output_dir <- file.path("../../results/disease_specific_changes/gsea", tail(unlist(strsplit(dirname(file), split = "/")), 1))
  
  ggsave(filename = file.path(output_dir, "gsea_barplot.png"), plot = p)
  
}

lapply(files, process_file)
