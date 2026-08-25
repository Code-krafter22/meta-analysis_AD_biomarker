suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
input_dir <- file.path(project_root, "results", "figures", "intermediate", "figure1")
output_dir <- file.path(project_root, "results", "figures", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- file.path(input_dir, c(
  "meta_effects_37_genes.csv",
  "stouffer_37_genes.csv",
  "metavolcano_37_genes.csv"
))
if (!all(file.exists(required))) {
  stop("Run 00_prepare_37_gene_volcano_inputs.R first.")
}

make_panel <- function(data, x_col, p_col, fdr_col, x_label, panel_label) {
  data$x <- suppressWarnings(as.numeric(data[[x_col]]))
  data$p <- pmax(suppressWarnings(as.numeric(data[[p_col]])), .Machine$double.xmin)
  data$fdr <- suppressWarnings(as.numeric(data[[fdr_col]]) )
  data$neg_log10_p <- -log10(data$p)
  data$direction <- ifelse(data$x > 0, "Up", "Down")

  ggplot(data, aes(x = x, y = neg_log10_p, color = direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, color = "grey35") +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_text_repel(aes(label = SYMBOL), size = 3.2, max.overlaps = Inf,
                    box.padding = 0.35, point.padding = 0.2,
                    min.segment.length = 0, seed = 42, show.legend = FALSE) +
    scale_color_manual(values = c(Down = "#3C78A8", Up = "#D64F43")) +
    labs(title = panel_label, x = x_label, y = expression(-log[10](P)), color = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "right")
}

meta_effects <- read.csv(required[1], check.names = FALSE)
stouffer <- read.csv(required[2], check.names = FALSE)
metavolcano <- read.csv(required[3], check.names = FALSE)

p1 <- make_panel(meta_effects, "meta_logFC", "meta_p", "meta_FDR",
                 expression(Meta~log[2]~fold~change), "i")
p2 <- make_panel(stouffer, "Z", "P", "FDR", "Stouffer Z", "ii")
p3 <- make_panel(metavolcano, "metafc", "metap", "FDR",
                 "MetaVolcanoR meta fold change", "iii")

combined <- p1 | p2 | p3
ggsave(file.path(output_dir, "Figure1_volcano_37_genes.png"), combined,
       width = 15, height = 6, dpi = 600, bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "Figure1_volcano_37_genes.pdf"), combined,
       width = 15, height = 6, bg = "white", limitsize = FALSE)

capture.output(sessionInfo(), file = file.path(output_dir, "Figure1_volcano_sessionInfo.txt"))
cat("Saved Figure 1 volcano panels to:", output_dir, "\n")
