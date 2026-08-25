# Preprocess GSE278723 counts for region-specific differential expression

library(here)

dataset_id <- "GSE278723"
input_file <- here::here("data", "raw", dataset_id, "GSE278723_Counts.txt")
output_dir <- here::here("data", "processed", dataset_id)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(
    "Raw count file not found: ", input_file,
    "\nDownload the GEO count file and place it at this location."
  )
}

counts_raw <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Column selection retained exactly from the original preprocessing script.
selected_columns <- c(
  1,
  12:25,
  31:40,
  51:65,
  71:85,
  91:95,
  101:105,
  131:164,
  170:174,
  185:199,
  205:209,
  220:229,
  240:244
)

if (max(selected_columns) > ncol(counts_raw)) {
  stop(
    "The raw table contains ", ncol(counts_raw),
    " columns, but the retained selection requires column ",
    max(selected_columns), "."
  )
}

selected <- counts_raw[, selected_columns, drop = FALSE]
colnames(selected) <- sub("^X", "", colnames(selected))
colnames(selected)[1] <- "Gene"

if (anyDuplicated(colnames(selected))) {
  stop("Duplicate column names were found after cleaning sample names.")
}

regions <- c("CA1", "CA2", "CA3", "CA4", "NX")

for (region in regions) {
  sample_columns <- grep(
    paste0("\\.", region, "$"),
    colnames(selected),
    value = TRUE
  )

  if (length(sample_columns) == 0L) {
    stop("No sample columns were found for region ", region, ".")
  }

  region_counts <- selected[, c("Gene", sample_columns), drop = FALSE]
  output_file <- file.path(output_dir, paste0("output_", region, ".csv"))

  write.csv(region_counts, output_file, row.names = FALSE)
  message("Saved ", output_file, " (", length(sample_columns), " samples)")
}
