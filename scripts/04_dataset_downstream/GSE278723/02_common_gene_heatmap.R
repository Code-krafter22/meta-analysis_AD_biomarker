source(file.path("scripts", "04_dataset_downstream", "GSE278723", "00_config.R"))

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Install pheatmap first: install.packages('pheatmap')")
}

common_up_file <- file.path(overlap_dir, "common_up_all_5_regions.txt")
common_down_file <- file.path(overlap_dir, "common_down_all_5_regions.txt")
if (!all(file.exists(c(common_up_file, common_down_file)))) {
  stop("Run 01_regional_overlap.R first.")
}
common_genes <- unique(c(read_gene_list(common_up_file), read_gene_list(common_down_file)))
if (!length(common_genes)) stop("No genes were common to all five regions.")

read_logcpm <- function(path) {
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  colnames(x) <- sub("\\..*$", "", colnames(x))
  x
}

mats <- lapply(logcpm_files, read_logcpm)
shared_genes <- Reduce(intersect, c(list(common_genes), lapply(mats, rownames)))
shared_subjects <- Reduce(intersect, lapply(mats, colnames))
if (!length(shared_genes) || !length(shared_subjects)) {
  stop("No shared genes or subjects remain across the five regional matrices.")
}

mats <- lapply(mats, function(x) x[shared_genes, shared_subjects, drop = FALSE])
mean_logcpm <- Reduce(`+`, mats) / length(mats)
write.csv(mean_logcpm, file.path(heatmap_dir, "common_genes_mean_logCPM_across_regions.csv"))

z <- t(scale(t(mean_logcpm)))
z[!is.finite(z)] <- 0
top_n <- min(50L, nrow(z))
variable_genes <- names(sort(apply(z, 1, var), decreasing = TRUE))[seq_len(top_n)]

subject_number <- suppressWarnings(as.integer(shared_subjects))
condition <- ifelse(subject_number <= 14, "Control", "AD")
annotation <- data.frame(Condition = factor(condition, levels = c("Control", "AD")))
rownames(annotation) <- shared_subjects
ord <- order(annotation$Condition, subject_number)

grDevices::pdf(file.path(heatmap_dir, "common_genes_top50_heatmap.pdf"), width = 10, height = 8)
pheatmap::pheatmap(
  z[variable_genes, ord, drop = FALSE],
  annotation_col = annotation[ord, , drop = FALSE],
  cluster_cols = FALSE, cluster_rows = TRUE,
  show_colnames = FALSE, border_color = NA,
  color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)
grDevices::dev.off()
save_session_info()
