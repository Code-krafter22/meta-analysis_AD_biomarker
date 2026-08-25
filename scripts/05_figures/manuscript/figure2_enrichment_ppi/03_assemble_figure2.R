suppressPackageStartupMessages({ library(here); library(png); library(grid); library(gridExtra) })
int <- here("results", "figures", "intermediate", "figure2")
main <- here("results", "figures", "main"); dir.create(main, recursive = TRUE, showWarnings = FALSE)
enrichment <- rasterGrob(readPNG(file.path(int, "Figure2_enrichment.png")), interpolate = TRUE)
ppi <- rasterGrob(readPNG(file.path(int, "Figure2_STRING_PPI.png")), interpolate = TRUE)
figure <- arrangeGrob(
  arrangeGrob(textGrob("(A)", x = 0, just = "left", gp = gpar(fontsize = 16, fontface = "bold")), enrichment,
              ncol = 1, heights = c(.4, 12)),
  arrangeGrob(textGrob("(B)", x = 0, just = "left", gp = gpar(fontsize = 16, fontface = "bold")), ppi,
              ncol = 1, heights = c(.4, 12)),
  ncol = 2
)
png(file.path(main, "Figure2_enrichment_and_STRING_PPI.png"), width = 4800, height = 2700, res = 300)
grid.draw(figure); dev.off()
pdf(file.path(main, "Figure2_enrichment_and_STRING_PPI.pdf"), width = 16, height = 9, useDingbats = FALSE)
grid.draw(figure); dev.off()
message("Saved reproducible Figure 2 (enrichment + STRING PPI).")
