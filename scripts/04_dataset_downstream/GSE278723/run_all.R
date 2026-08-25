scripts <- file.path("scripts", "04_dataset_downstream", "GSE278723", c(
  "01_regional_overlap.R",
  "02_common_gene_heatmap.R",
  "03_common_gene_enrichment.R"
))
for (script in scripts) {
  message("Running ", script)
  source(script, local = new.env(parent = globalenv()))
}
message("GSE278723 downstream workflow complete.")
