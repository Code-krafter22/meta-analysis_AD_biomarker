required_packages <- c(
  "dplyr", "readr", "stringr", "purrr", "tidyr", "tibble", "ggplot2",
  "AnnotationDbi", "org.Hs.eg.db", "metafor", "MetaVolcanoR", "here",
  "pheatmap", "reshape2", "plotly", "htmlwidgets"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       ". Install them before running the analysis.")
}

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr); library(tidyr)
  library(AnnotationDbi); library(org.Hs.eg.db); library(metafor)
  library(MetaVolcanoR); library(ggplot2); library(here)
})

if (!exists("analysis_mode", inherits = FALSE)) {
  stop("Set analysis_mode to 'filtered_unadjusted' or 'unfiltered_all_genes'.")
}
analysis_mode <- match.arg(analysis_mode,
                           c("filtered_unadjusted", "unfiltered_all_genes",
                             "logfc_threshold"))

# Screening knobs. Only the fold-change sensitivity branch sets them; for the
# two prespecified branches they are NULL so stage 01 applies no extra screen
# (this also clears them if a previous branch ran in the same session).
if (analysis_mode == "logfc_threshold") {
  if (!exists("screen_lfc_cut", inherits = FALSE) || !is.finite(screen_lfc_cut)) {
    stop("Set screen_lfc_cut (|log2FC| gate) before running the logfc_threshold branch.")
  }
  if (!exists("screen_p_cut", inherits = FALSE)) screen_p_cut <- 0.05
} else {
  screen_lfc_cut <- NULL
  screen_p_cut   <- NULL
}

project_root <- here::here()
deg_root <- file.path(project_root, "results", "differential_expression")
output_tag <- if (analysis_mode == "logfc_threshold") {
  paste0("logfc_gt_", formatC(screen_lfc_cut, format = "f", digits = 2))
} else {
  analysis_mode
}
default_output <- file.path(project_root, "results", "meta_analysis", output_tag)
output_dir <- Sys.getenv("META_ANALYSIS_OUTPUT_DIR", unset = default_output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Per-dataset input table. Overridable so the FDR-screened branch can read the
# adjusted DEG tables (adj.P.Val <= 0.05 & |log2FC| > 0.58) without a new mode.
default_input <- if (analysis_mode == "filtered_unadjusted") {
  "significant_DEGs_unadjusted.csv"
} else {
  "full_DE_results.csv"
}
input_file <- Sys.getenv("META_ANALYSIS_INPUT_FILE", unset = default_input)

files <- tibble::tribble(
  ~acc,             ~region,                  ~path,
  "GSE278723_CA1", "hippocampus_CA1",        file.path(deg_root, "GSE278723", "CA1", input_file),
  "GSE278723_CA2", "hippocampus_CA2",        file.path(deg_root, "GSE278723", "CA2", input_file),
  "GSE278723_CA3", "hippocampus_CA3",        file.path(deg_root, "GSE278723", "CA3", input_file),
  "GSE278723_CA4", "hippocampus_CA4",        file.path(deg_root, "GSE278723", "CA4", input_file),
  "GSE278723_NX",  "entorhinal_cortex",      file.path(deg_root, "GSE278723", "NX", input_file),
  "GSE67333",       "hippocampus",            file.path(deg_root, "GSE67333", input_file),
  "GSE95587",       "fusiform_gyrus",         file.path(deg_root, "GSE95587", input_file),
  "GSE159699",      "lateral_temporal_lobe",  file.path(deg_root, "GSE159699", input_file),
  "GSE203206",      "occipital_lobe",         file.path(deg_root, "GSE203206", input_file)
)

missing_inputs <- files$path[!file.exists(files$path)]
if (length(missing_inputs)) {
  stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}

resolve_base <- function(rel_paths) {
  if (!length(rel_paths)) return(character())
  vapply(rel_paths, function(p) {
    if (file.exists(p)) return(normalizePath(p, mustWork = TRUE))
    candidate <- file.path(project_root, p)
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
    stop("Input file not found: ", p)
  }, character(1))
}

message("Meta-analysis mode: ", analysis_mode)
message("Input type: ", input_file)
message("Output directory: ", output_dir)
