# REVIEWER-REQUESTED SENSITIVITY BRANCH
# Inputs: all tested genes, with no P-value or logFC DEG prefilter.
# Purpose: evaluates robustness without threshold-based input selection.
analysis_mode <- "unfiltered_all_genes"
source(file.path(here::here(), "scripts", "03_meta_analysis", "run_meta_analysis.R"))
