suppressPackageStartupMessages({ library(here); library(png); library(grid); library(gridExtra) })
intermediate <- here("results", "figures", "intermediate", "figure1")
main_dir <- here("results", "figures", "main"); dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
heatmap_file <- file.path(intermediate, "Figure1_heatmap_37genes.png")
volcano_file <- file.path(main_dir, "Figure1_volcano_37_genes.png")
stopifnot(file.exists(heatmap_file), file.exists(volcano_file))
heatmap <- rasterGrob(readPNG(heatmap_file), interpolate = TRUE)
volcano <- rasterGrob(readPNG(volcano_file), interpolate = TRUE)
figure <- arrangeGrob(
  textGrob("(A)", x = 0, just = "left", gp = gpar(fontsize = 16, fontface = "bold")),
  heatmap,
  textGrob("(B)", x = 0, just = "left", gp = gpar(fontsize = 16, fontface = "bold")),
  volcano,
  ncol = 1, heights = unit.c(unit(.3, "in"), unit(6.4, "in"), unit(.3, "in"), unit(4.2, "in"))
)
png(file.path(main_dir, "Figure1_heatmap_and_volcano.png"), width = 4200, height = 3400, res = 300)
grid.draw(figure); dev.off()
pdf(file.path(main_dir, "Figure1_heatmap_and_volcano.pdf"), width = 14, height = 11.3, useDingbats = FALSE)
grid.draw(figure); dev.off()
message("Saved reproducible Figure 1 (heatmap + volcano).")
