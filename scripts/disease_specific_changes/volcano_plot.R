library(ggplot2)
library(data.table)
library(ggrepel)
library(plotly)
library(data.table)
library(fs)
library(tools)


files <- list.files(path = "../../results/disease_specific_changes/dea", pattern = ".tsv", full.names = TRUE, recursive = TRUE)

process_file <- function(file) {
  df <- fread(file)
  df <- df[!is.na(padj)]
  df <- df[order(padj)]
  
  df$significant <- df$padj < 0.05 & abs(df$log2FoldChange) > 1
  
  p1 <- ggplot(df, aes(text = paste("Gene:", gene), x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), color="gray", linetype=2) +
    geom_vline(xintercept = -1, color="gray", linetype=2) +
    geom_vline(xintercept = 1, color="gray", linetype=2) +
    geom_text_repel(data = subset(df, significant), aes(label = gene), size = 5, max.overlaps = 100, color="#00BFC4") +
    ggtitle(strsplit(file, "/")[[1]][2]) +
    theme_bw() +
    theme(axis.title = element_text(size=15),
          legend.title = element_text(size=15),
          legend.text = element_text(size=12))
  
  p2 <- ggplotly(p1)
  
  output_dir <- file.path("..", "..", "results", "disease_specific_changes", "dea", tail(unlist(strsplit(dirname(file), split = "/")), 1))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  ggsave(filename = file.path(output_dir, "volcano_static.png"), plot = p1, width = 8, height = 6)
  htmlwidgets::saveWidget(as_widget(p2), file.path(output_dir, "volcano_interactive.html"))
}


lapply(files, process_file)

