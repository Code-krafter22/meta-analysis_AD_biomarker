root <- file.path(here::here(), "scripts", "05_figures", "manuscript")
scripts <- c(
  file.path("shared_braak", "00_prepare_braak_analysis.R"),
  file.path("figure3_braak_overall", "01_overall_waterfall.R"),
  file.path("figure3_braak_overall", "02_reproducible_forest_plot.R"),
  file.path("figure3_braak_overall", "03_assemble_reproducible_figure3.R"),
  file.path("figure4_braak_regional", "01_regional_waterfall.R")
)
for (script in scripts) {
  message("Running: ", script)
  sys.source(file.path(root, script), envir = globalenv())
}
message("Figures 3 and 4 completed from one shared Braak preparation.")
