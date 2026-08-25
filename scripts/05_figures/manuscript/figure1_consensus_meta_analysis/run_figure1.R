script_dir <- file.path(here::here(), "scripts", "05_figures", "manuscript", "figure1_consensus_meta_analysis")
for (script in c("00_prepare_37_gene_volcano_inputs.R", "02_reproducible_heatmap.R",
                 "03_meta_volcanoes_from_csv.R", "04_assemble_reproducible_figure1.R")) {
  message("\nRunning Figure 1 step: ", script)
  sys.source(file.path(script_dir, script), envir = globalenv())
}
