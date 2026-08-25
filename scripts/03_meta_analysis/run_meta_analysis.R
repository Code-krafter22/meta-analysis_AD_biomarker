# INTERNAL SHARED ENGINE. Do not run this file directly.
# Public entry points are 00_run_all_meta_analyses.R and the two numbered
# branch runners. Statistical stages are intentionally shared.
if (!exists("analysis_mode", inherits = FALSE)) {
  stop("Set analysis_mode before sourcing run_meta_analysis.R")
}

script_dir <- file.path(here::here(), "scripts", "03_meta_analysis")
stage_files <- c(
  "00_config_and_inputs.R",
  file.path("stages", "01_harmonize_inputs.R"),
  file.path("stages", "02_collapse_regions_and_qc.R"),
  file.path("stages", "03_stouffer_meta_analysis.R"),
  file.path("stages", "04_effect_size_meta_analysis.R"),
  file.path("stages", "05_metavolcano_and_diagnostics.R")
)

for (stage in stage_files) {
  message("\n========== Running ", stage, " ==========")
  sys.source(file.path(script_dir, stage), envir = globalenv())
}

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
message("\nCompleted ", analysis_mode, " meta-analysis.")
