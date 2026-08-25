# ============================================================
# META-ANALYSIS: RUN BOTH PRESPECIFIED ANALYSIS BRANCHES
# ============================================================
# Run from the project root:
# source("scripts/03_meta_analysis/00_run_all_meta_analyses.R")
#
# Branch 1 (primary/historical): filtered DEGs using unadjusted P < 0.05
# Branch 2 (reviewer analysis):  all tested genes, with no DEG prefilter
#
# Both branches use the same implementation in stages/. They differ only in
# input selection and output directory, which are set by analysis_mode.

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}

runner_dir <- file.path(here::here(), "scripts", "03_meta_analysis")

message("\n============================================================")
message("BRANCH 1 OF 2: FILTERED DEGs / UNADJUSTED P-VALUE")
message("============================================================")
sys.source(file.path(runner_dir, "01_meta_analysis_filtered_unadjusted.R"), envir = globalenv())

message("\n============================================================")
message("BRANCH 2 OF 2: UNFILTERED / ALL TESTED GENES")
message("============================================================")
sys.source(file.path(runner_dir, "02_meta_analysis_unfiltered.R"), envir = globalenv())

message("\nBoth prespecified meta-analysis branches completed.")
