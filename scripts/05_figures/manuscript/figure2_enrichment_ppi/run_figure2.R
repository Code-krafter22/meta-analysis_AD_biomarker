script_dir <- file.path(here::here(), "scripts", "05_figures", "manuscript", "figure2_enrichment_ppi")
for (script in c("00_prepare_signature.R", "01_reproducible_enrichment.R",
                 "02_reproducible_ppi_network.R", "03_assemble_figure2.R")) {
  message("\nRunning Figure 2 step: ", script)
  sys.source(file.path(script_dir, script), envir = globalenv())
}
