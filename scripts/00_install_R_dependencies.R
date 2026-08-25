# Install packages required by the current reproducible R workflows.
cran <- c(
  "broom", "dplyr", "ggplot2", "ggrepel", "ggraph", "gridExtra", "here",
  "igraph", "matrixStats", "metafor", "patchwork", "pheatmap", "png",
  "purrr", "readr", "readxl", "reshape2", "stringr", "tibble", "tidyr",
  "tidygraph", "VennDiagram"
)
bioc <- c(
  "AnnotationDbi", "clusterProfiler", "edgeR", "enrichplot", "GEOquery",
  "limma", "org.Hs.eg.db", "sva"
)

missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

message("R dependency check complete.")
