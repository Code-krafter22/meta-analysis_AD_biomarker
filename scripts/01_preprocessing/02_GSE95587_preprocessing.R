# Recreate the GSE95587 processed count table from GEO downloads.

library(here)
library(readxl)

dataset_id <- "GSE95587"
raw_dir <- here::here("data", "raw", dataset_id)
metadata_file <- here::here("data", "metadata", dataset_id, "metadata.xlsx")
output_dir <- here::here("data", "processed", dataset_id)
output_file <- file.path(output_dir, "Edited_raw_counts.csv")

counts_candidates <- file.path(
  raw_dir,
  c(
    "GSE95587_raw_counts_GRCh38.p13_NCBI.tsv",
    "GSE95587_raw_counts.tsv"
  )
)
counts_file <- counts_candidates[file.exists(counts_candidates)][1]
annotation_file <- file.path(raw_dir, "Human.GRCh38.p13.annot.tsv")

if (is.na(counts_file)) {
  stop(
    "GSE95587 count table not found. Expected one of:\n",
    paste0("- ", counts_candidates, collapse = "\n"),
    "\nSee data/raw/GSE95587/README.md."
  )
}

required_files <- c(annotation_file, metadata_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Required input file(s) not found:\n", paste0("- ", missing_files, collapse = "\n"))
}

counts_raw <- read.delim(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation <- read.delim(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- readxl::read_excel(metadata_file)

if (!identical(colnames(counts_raw)[1], "GeneID")) {
  stop("The first count-table column must be named GeneID.")
}
if (!all(c("GeneID", "Symbol") %in% colnames(annotation))) {
  stop("The annotation table must contain GeneID and Symbol columns.")
}
if (!all(c("GSM", "IDENTIFIER") %in% colnames(metadata))) {
  stop("The metadata workbook must contain GSM and IDENTIFIER columns.")
}

symbol_index <- match(counts_raw$GeneID, annotation$GeneID)
if (anyNA(symbol_index)) {
  stop(sum(is.na(symbol_index)), " GeneID value(s) lack an annotation match.")
}

sample_columns <- colnames(counts_raw)[-1]
sample_index <- match(sample_columns, metadata$GSM)
if (anyNA(sample_index)) {
  stop("Metadata are missing for: ", paste(sample_columns[is.na(sample_index)], collapse = ", "))
}

analysis_names <- as.character(metadata$IDENTIFIER[sample_index])
if (anyNA(analysis_names) || anyDuplicated(analysis_names)) {
  stop("Analysis identifiers must be complete and unique.")
}

processed <- data.frame(
  Symbol = annotation$Symbol[symbol_index],
  counts_raw[, -1, drop = FALSE],
  check.names = FALSE
)
colnames(processed)[-1] <- analysis_names

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(processed, output_file, row.names = FALSE)
message("Saved ", output_file, " (", nrow(processed), " genes; ", length(analysis_names), " samples)")
