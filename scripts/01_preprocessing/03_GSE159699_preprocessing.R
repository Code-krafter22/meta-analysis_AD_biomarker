# Recreate the GSE159699 processed count table from the GEO supplementary table.

library(here)
library(readxl)

dataset_id <- "GSE159699"
counts_file <- here::here(
  "data",
  "raw",
  dataset_id,
  "GSE159699_summary_count.star.txt"
)
metadata_file <- here::here("data", "metadata", dataset_id, "metadata.xlsx")
output_dir <- here::here("data", "processed", dataset_id)
output_file <- file.path(output_dir, "EDITED_COUNTS.csv")

required_files <- c(counts_file, metadata_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste0("- ", missing_files, collapse = "\n"),
    "\nDownload the GEO file described in data/raw/GSE159699/README.md."
  )
}

counts_raw <- read.delim(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- readxl::read_excel(metadata_file)

if (!identical(colnames(counts_raw)[1], "refGene")) {
  stop("The first count-table column must be named refGene.")
}
if (!all(c("Sample ID", "identifier") %in% colnames(metadata))) {
  stop("The metadata workbook must contain Sample ID and identifier columns.")
}

sample_columns <- colnames(counts_raw)[-1]
sample_ids <- sub("-.*$", "", sample_columns)
metadata_ids <- as.character(metadata$`Sample ID`)
sample_index <- match(sample_ids, metadata_ids)

if (anyNA(sample_index)) {
  stop("Metadata are missing for source sample(s): ", paste(sample_columns[is.na(sample_index)], collapse = ", "))
}

analysis_names <- as.character(metadata$identifier[sample_index])
if (anyNA(analysis_names) || anyDuplicated(analysis_names)) {
  stop("Analysis identifiers must be complete and unique.")
}

processed <- counts_raw
colnames(processed)[-1] <- analysis_names

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(processed, output_file, row.names = FALSE)
message("Saved ", output_file, " (", nrow(processed), " genes; ", length(analysis_names), " samples)")
