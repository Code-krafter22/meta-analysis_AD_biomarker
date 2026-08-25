source(file.path("scripts", "04_dataset_downstream", "GSE278723", "00_config.R"))

has_venn <- requireNamespace("VennDiagram", quietly = TRUE)

up_sets <- lapply(up_files, read_gene_list)
down_sets <- lapply(down_files, read_gene_list)
common_up <- Reduce(intersect, up_sets)
common_down <- Reduce(intersect, down_sets)

writeLines(common_up, file.path(overlap_dir, "common_up_all_5_regions.txt"))
writeLines(common_down, file.path(overlap_dir, "common_down_all_5_regions.txt"))

counts <- data.frame(
  region = regions,
  upregulated = lengths(up_sets),
  downregulated = lengths(down_sets)
)
write.csv(counts, file.path(overlap_dir, "regional_DEG_counts.csv"), row.names = FALSE)

pairwise_counts <- function(sets) {
  pairs <- utils::combn(names(sets), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(x) {
    data.frame(region_1 = x[1], region_2 = x[2],
               intersection_n = length(intersect(sets[[x[1]]], sets[[x[2]]])))
  }))
}
write.csv(pairwise_counts(up_sets), file.path(overlap_dir, "pairwise_up_counts.csv"), row.names = FALSE)
write.csv(pairwise_counts(down_sets), file.path(overlap_dir, "pairwise_down_counts.csv"), row.names = FALSE)

draw_venn <- function(sets, common, filename, title) {
  grDevices::pdf(filename, width = 8, height = 8)
  grid::grid.newpage()
  plot <- VennDiagram::venn.diagram(
    x = sets, filename = NULL,
    fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"),
    alpha = 0.45, cex = rep(0, 31), cat.cex = 1.1, margin = 0.08
  )
  grid::grid.draw(plot)
  grid::grid.text(length(common), x = 0.5, y = 0.5,
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  grid::grid.text(title, x = 0.5, y = 0.97,
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  grDevices::dev.off()
}

if (has_venn) {
  draw_venn(up_sets, common_up, file.path(overlap_dir, "venn_up_5_regions.pdf"),
            "GSE278723: upregulated genes")
  draw_venn(down_sets, common_down, file.path(overlap_dir, "venn_down_5_regions.pdf"),
            "GSE278723: downregulated genes")
} else {
  warning("VennDiagram is not installed; overlap tables were created but Venn PDFs were skipped. ",
          "Install it with install.packages('VennDiagram') and rerun this script.")
}

message("Common UP across all five regions: ", length(common_up))
message("Common DOWN across all five regions: ", length(common_down))
save_session_info()
