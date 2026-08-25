project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "scripts", "04_dataset_downstream"))) {
  stop("Run from the project root (the directory containing scripts/ and results/).")
}

regions <- c("CA1", "CA2", "CA3", "CA4", "NX")
deg_root <- file.path(project_root, "results", "differential_expression", "GSE278723")
out_root <- file.path(project_root, "results", "dataset_downstream", "GSE278723")

overlap_dir <- file.path(out_root, "01_regional_overlap")
heatmap_dir <- file.path(out_root, "02_common_gene_heatmap")
enrichment_dir <- file.path(out_root, "03_common_gene_enrichment")
invisible(lapply(c(overlap_dir, heatmap_dir, enrichment_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

up_files <- setNames(file.path(deg_root, regions, "upregulated_genes.txt"), regions)
down_files <- setNames(file.path(deg_root, regions, "downregulated_genes.txt"), regions)
logcpm_files <- setNames(file.path(deg_root, regions,
                                   paste0(regions, "_filtered_logcpm.csv")), regions)

required_inputs <- c(up_files, down_files, logcpm_files)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing GSE278723 DEG input(s):\n", paste(missing_inputs, collapse = "\n"))
}

read_gene_list <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  unique(x[nzchar(x)])
}

save_session_info <- function() {
  capture.output(sessionInfo(), file = file.path(out_root, "sessionInfo.txt"))
}
