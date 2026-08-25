# HISTORICAL PROVENANCE COPY — DO NOT RUN. The reproducible Figure 1 excludes
# the Venn diagram and uses run_figure1.R.
library(VennDiagram)
library(grid)

gene_sets <- list(
  "i" = c("KANK2","CHML","PRELP","ADAMTS2","CCDC102A","HSPB1","AEBP1",
          "RAB3C","GPR22","KLF15","ELOVL4","ITFG1","NAP1L5","PIK3R5",
          "GFAP","RPH3A","TCEA3","HPRT1","NEUROD6","MRGPRF","NRN1",
          "MAS1","FOXJ1","NUPR1","ADAM33","SERPINI1","CLDN9","COL25A1",
          "HSPB7","SCG2","PRX","STAT4","TMPRSS5","RAB3B","NCALD",
          "RGS4","OPN3","GAD2","TRIM36","GAD1"),
  "ii" = c("CHML","KANK2","PRELP","HSPB1","PIK3R5","AEBP1","GFAP",
           "ITFG1","TCEA3","CCDC102A","ELOVL4","KLF15","NAP1L5",
           "NEUROD6","HPRT1","NRN1","MAS1","RPH3A","HSPB7","NUPR1",
           "ADAM33","COL25A1","RAB3C","SERPINI1","STAT4","SCG2",
           "TMPRSS5","RGS4","PRX","RAB3B","NCALD","OPN3","CLDN9",
           "TRIM36","GAD2","MRGPRF","GAD1"),
  "iii" = c("CHML","FOXJ1","HSPB1","MAS1","KANK2","ADAMTS2","PRELP","CCDC102A",
            "NEUROD6","AEBP1","RPH3A","GPR22","RAB3C","NAP1L5","GFAP","MRGPRF",
            "PIK3R5","ELOVL4","KLF15","NRN1","RGS4","TCEA3","COL25A1","HPRT1",
            "ADAM33","NUPR1","ITFG1","CLDN9","STAT4","HSPB7","SERPINI1","PRX",
            "SCG2","TMPRSS5","RAB3B","GAD2","NCALD","TNFRSF10D","OPN3","GAD1",
            "TRIM36","ARL4D","ARHGAP42","LINC01561","PLPP4","VWA7","DNAJC5G")
)

library(extrafont)

# Set a universal font for all plots
universal_font <- "Arial"

# Apply globally for ggplot2
theme_set(theme_bw(base_size = 18, base_family = universal_font))

col1 <- "#4E79A7"
col2 <- "#E15759"
col3 <- "#76B7B2"

# Calculate the triple intersection count explicitly
triple_count <- length(Reduce(intersect, gene_sets))

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  fill     = c(col1, col2, col3),
  alpha    = 0.45,
  col      = c("black", "black", "black"),
  lwd      = 0.2,
  print.mode   = "raw",
  cex          = 1,
  fontfamily   = "Arial",          # <-- numbers inside circles
  cat.cex         = 1,
  cat.fontfamily  = "Arial",
  cat.col         = c("black", "black", "black"),
  cat.default.pos = "outer",
  cat.dist        = c(0.12, 0.12, 0.12),
  cat.pos         = c(-20, 20, 180),
  scaled  = FALSE,
  euler.d = FALSE,
  margin  = 0.02
)

# Find all text grobs and blank everything EXCEPT the triple intersection value
for (i in seq_along(venn.plot)) {
  grob <- venn.plot[[i]]
  if (inherits(grob, "text") && !is.null(grob$label)) {
    # Keep only the grob whose label matches the triple intersection count
    # AND is not a category label (i, ii, iii)
    if (grob$label != as.character(triple_count) ||
        grob$label %in% c("i", "ii", "iii")) {
      if (!grob$label %in% c("i", "ii", "iii")) {
        venn.plot[[i]]$label <- ""
      }
    }
  }
}

jpeg("C:/Users/kashvichirag/Box/Dr. ZHANG/alzheimer_meta/ad_meta/FIGURES/37-genes/VennDiagram_MetaAnalysis.jpeg",
     width = 1000, height = 700, res = 600, bg = "white")
grid.newpage()
grid.draw(venn.plot)
dev.off()
