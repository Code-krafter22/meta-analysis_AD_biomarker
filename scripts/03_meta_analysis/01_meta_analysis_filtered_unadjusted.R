# PRIMARY/HISTORICAL BRANCH
# Inputs: dataset DEG tables filtered with unadjusted P.Value < 0.05.
# Purpose: reproduces the analysis from which the 37-gene signature arose.
analysis_mode <- "filtered_unadjusted"
source(file.path(here::here(), "scripts", "03_meta_analysis", "run_meta_analysis.R"))
