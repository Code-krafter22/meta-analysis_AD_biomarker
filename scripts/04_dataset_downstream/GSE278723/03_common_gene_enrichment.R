source(file.path("scripts", "04_dataset_downstream", "GSE278723", "00_config.R"))

needed <- c("clusterProfiler", "org.Hs.eg.db")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required Bioconductor package(s): ", paste(missing_packages, collapse = ", "))
}

run_enrichment <- function(gene_file, direction) {
  genes <- read_gene_list(gene_file)
  mapping <- clusterProfiler::bitr(
    genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  mapping <- mapping[!duplicated(mapping$ENTREZID), , drop = FALSE]
  write.csv(mapping, file.path(enrichment_dir, paste0(direction, "_symbol_to_entrez.csv")),
            row.names = FALSE)

  if (nrow(mapping) == 0) {
    warning("No Entrez IDs mapped for ", direction)
    return(invisible(NULL))
  }

  for (ontology in c("BP", "CC", "MF")) {
    result <- clusterProfiler::enrichGO(
      gene = mapping$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID", ont = ontology, pAdjustMethod = "BH",
      pvalueCutoff = 0.05, qvalueCutoff = 0.20, readable = TRUE
    )
    write.csv(as.data.frame(result),
              file.path(enrichment_dir, paste0(direction, "_GO_", ontology, ".csv")),
              row.names = FALSE)
  }

  kegg <- clusterProfiler::enrichKEGG(
    gene = mapping$ENTREZID, organism = "hsa", pvalueCutoff = 0.05,
    pAdjustMethod = "BH", qvalueCutoff = 0.20
  )
  write.csv(as.data.frame(kegg), file.path(enrichment_dir, paste0(direction, "_KEGG.csv")),
            row.names = FALSE)
}

up_file <- file.path(overlap_dir, "common_up_all_5_regions.txt")
down_file <- file.path(overlap_dir, "common_down_all_5_regions.txt")
if (!all(file.exists(c(up_file, down_file)))) stop("Run 01_regional_overlap.R first.")

run_enrichment(up_file, "common_up")
run_enrichment(down_file, "common_down")
save_session_info()
