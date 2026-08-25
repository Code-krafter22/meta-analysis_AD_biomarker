project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(project_root, "results", "meta_analysis"))) {
  stop("Run from the project root containing results/meta_analysis.")
}

meta_dir <- file.path(project_root, "results", "meta_analysis", "filtered_unadjusted")
out_dir <- file.path(project_root, "results", "figures", "intermediate", "figure1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

signature_37 <- c(
  "CHML", "KANK2", "PRELP", "HSPB1", "PIK3R5", "AEBP1", "GFAP",
  "ITFG1", "TCEA3", "CCDC102A", "ELOVL4", "KLF15", "NAP1L5",
  "NEUROD6", "HPRT1", "NRN1", "MAS1", "RPH3A", "HSPB7", "NUPR1",
  "ADAM33", "COL25A1", "RAB3C", "SERPINI1", "STAT4", "SCG2",
  "TMPRSS5", "RGS4", "PRX", "RAB3B", "NCALD", "OPN3", "CLDN9",
  "TRIM36", "GAD2", "MRGPRF", "GAD1"
)
stopifnot(length(signature_37) == 37L, !anyDuplicated(signature_37))

inputs <- c(
  meta_effects = file.path(meta_dir, "meta_effects_meta_results.csv"),
  stouffer = file.path(meta_dir, "stouffer_meta_results.csv"),
  metavolcano = file.path(meta_dir, "metavolcano", "updated_metavolcano_combining_results.csv")
)
missing_files <- inputs[!file.exists(inputs)]
if (length(missing_files)) stop("Missing meta-analysis file(s):\n", paste(missing_files, collapse = "\n"))

tables <- lapply(inputs, read.csv, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- list(
  meta_effects = c("SYMBOL", "meta_logFC", "meta_p", "meta_FDR"),
  stouffer = c("SYMBOL", "Z", "P", "FDR"),
  metavolcano = c("SYMBOL", "metafc", "metap", "FDR")
)
for (name in names(tables)) {
  absent <- setdiff(required_columns[[name]], names(tables[[name]]))
  if (length(absent)) stop(name, " is missing column(s): ", paste(absent, collapse = ", "))
}

subset_signature <- function(x, method) {
  result <- x[x$SYMBOL %in% signature_37, , drop = FALSE]
  result <- result[match(signature_37[signature_37 %in% result$SYMBOL], result$SYMBOL), , drop = FALSE]
  missing_genes <- setdiff(signature_37, result$SYMBOL)
  if (length(missing_genes)) warning(method, " missing: ", paste(missing_genes, collapse = ", "))
  result
}

subsets <- Map(subset_signature, tables, names(tables))
for (name in names(subsets)) {
  write.csv(subsets[[name]], file.path(out_dir, paste0(name, "_37_genes.csv")), row.names = FALSE)
}

availability <- data.frame(
  SYMBOL = signature_37,
  meta_effects = signature_37 %in% tables$meta_effects$SYMBOL,
  stouffer = signature_37 %in% tables$stouffer$SYMBOL,
  metavolcano = signature_37 %in% tables$metavolcano$SYMBOL
)
write.csv(availability, file.path(out_dir, "signature_37_method_availability.csv"), row.names = FALSE)
writeLines(signature_37, file.path(out_dir, "signature_37_genes.txt"))

cat("Prepared Figure 1 volcano inputs in:", out_dir, "\n")
print(colSums(availability[-1]))
