suppressPackageStartupMessages({
  library(here); library(readxl); library(pheatmap); library(sva)
})
input_dir <- here("data", "figure_inputs", "figure1")
out_dir <- here("results", "figures", "intermediate", "figure1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
expression_file <- file.path(input_dir, "heatmap_data_5dataset.xlsx")
phenotype_file <- file.path(input_dir, "pheno_batch_updated.csv")
stopifnot(file.exists(expression_file), file.exists(phenotype_file))

df <- as.data.frame(read_xlsx(expression_file, sheet = 1))
df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]
genes <- as.character(df[[1]])
stopifnot(nrow(df) == 37L, !anyNA(genes), !anyDuplicated(genes))
mat <- as.matrix(df[, -1, drop = FALSE]); mode(mat) <- "numeric"; rownames(mat) <- genes

# The last 26 source columns are GSE278723: 12 controls then 14 AD samples.
stopifnot(ncol(mat) >= 26L)
idx <- (ncol(mat) - 25L):ncol(mat)
colnames(mat)[idx] <- c(paste0("GSE278723_CTRL", 1:12), paste0("GSE278723_AD", 1:14))

pheno <- read.csv(phenotype_file, stringsAsFactors = FALSE)
stopifnot(all(c("SampleID", "diag", "batch") %in% names(pheno)),
          !anyDuplicated(pheno$SampleID))
if (!setequal(colnames(mat), pheno$SampleID)) stop("Expression and phenotype samples differ.")
mat <- mat[, pheno$SampleID, drop = FALSE]
if (anyNA(mat)) stop("Expression matrix contains missing values.")

design <- model.matrix(~ diag, data = pheno)
mat_combat <- ComBat(dat = mat, batch = pheno$batch, mod = design,
                     par.prior = TRUE, prior.plots = FALSE)
mat_z <- t(scale(t(mat_combat))); mat_z[!is.finite(mat_z)] <- 0
condition <- factor(pheno$diag, levels = c("AD", "CTRL"), labels = c("AD", "Control"))
ord <- order(condition, pheno$batch, pheno$SampleID)
mat_z <- mat_z[, ord, drop = FALSE]
ann <- data.frame(Condition = condition[ord], row.names = colnames(mat_z))
gap <- sum(ann$Condition == "AD")

write.csv(mat_combat, file.path(out_dir, "Figure1_37gene_ComBat_log2CPM.csv"))
write.csv(mat_z, file.path(out_dir, "Figure1_37gene_ComBat_row_zscore.csv"))
cols <- colorRampPalette(c("deepskyblue2", "#F7F7F7", "deeppink2"))(201)
breaks <- seq(-3, 3, length.out = 202)
mat_plot <- mat_z
mat_plot[mat_plot > 3] <- 3
mat_plot[mat_plot < -3] <- -3

draw <- function(path, type) {
  if (type == "png") png(path, 3600, 1900, res = 300)
  if (type == "pdf") pdf(path, 12, 6.4, useDingbats = FALSE)
  pheatmap(mat_plot, color = cols, breaks = breaks, annotation_col = ann,
    annotation_colors = list(Condition = c(AD = "aquamarine2", Control = "darkorange2")),
    show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
    gaps_col = gap, border_color = NA, fontsize_row = 7, fontface_row = "italic",
    main = "37 consensus genes across five datasets")
  dev.off()
}
draw(file.path(out_dir, "Figure1_heatmap_37genes.png"), "png")
draw(file.path(out_dir, "Figure1_heatmap_37genes.pdf"), "pdf")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Figure1_heatmap_sessionInfo.txt"))
message("Saved heatmap for 37 genes and ", ncol(mat_z), " samples.")
