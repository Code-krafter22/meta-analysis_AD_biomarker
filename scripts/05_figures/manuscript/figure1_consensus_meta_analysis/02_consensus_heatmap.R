# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use 02_reproducible_heatmap.R.
library(readxl)
library(pheatmap)
library(sva)
library(RColorBrewer)

# ============================================================
# 1) Read expression data (log2-CPM, genes x samples)
# ============================================================
infile <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\heatmap_data_5dataset.xlsx"
df <- as.data.frame(read_xlsx(infile, sheet = 1))
df <- df[rowSums(is.na(df)) < ncol(df), , drop = FALSE]

genes <- as.character(df[[1]])
genes[is.na(genes) | genes == ""] <- paste0("GENE_", which(is.na(genes) | genes == ""))

mat <- as.matrix(df[, -1, drop = FALSE])
mode(mat) <- "numeric"
rownames(mat) <- make.unique(genes)
colnames(mat) <- make.unique(colnames(mat))

cnames <- colnames(mat)

ctrl_xlsx_names <- paste0("CTRL...", 203:214)  # 12 controls
ad_xlsx_names   <- paste0("AD...", 215:228)    # 14 AD

ctrl_new_names <- paste0("GSE278723_CTRL", 1:12)
ad_new_names   <- paste0("GSE278723_AD", 1:14)

for (i in seq_along(ctrl_xlsx_names)) {
  idx <- which(cnames == ctrl_xlsx_names[i])
  if (length(idx) == 1) cnames[idx] <- ctrl_new_names[i]
}
for (i in seq_along(ad_xlsx_names)) {
  idx <- which(cnames == ad_xlsx_names[i])
  if (length(idx) == 1) cnames[idx] <- ad_new_names[i]
}

colnames(mat) <- make.unique(cnames)

cat("Expression matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")


# ============================================================
# 2) Read updated phenotype/batch file (with GSE278723 fixed)
# ============================================================
pheno <- read.csv("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\Data\\expression_data\\pheno_batch_updated.csv", stringsAsFactors = FALSE)

# ============================================================
# 3) Align samples between expression and pheno
# ============================================================
common <- intersect(colnames(mat), pheno$SampleID)
mat   <- mat[, common, drop = FALSE]
pheno <- pheno[match(common, pheno$SampleID), ]
rownames(pheno) <- pheno$SampleID

cat("Aligned samples:", ncol(mat), "| Genes:", nrow(mat),
    "| Batches:", length(unique(pheno$batch)), "\n")
cat("Samples per batch:\n")
print(table(pheno$batch))

# ============================================================
# 4) ComBat batch correction
# ============================================================
mod <- model.matrix(~ diag, data = pheno)
mat_combat <- ComBat(dat = mat, batch = pheno$batch, mod = mod)

# ============================================================
# 5) Condition assignment
# ============================================================
ann <- data.frame(
  Condition = factor(
    ifelse(pheno$diag == "AD", "AD", "Control"),
    levels = c("AD", "Control")
  )
)
rownames(ann) <- pheno$SampleID

# ============================================================
# 6) Order columns: AD together, then Control
# ============================================================
ord <- order(ann$Condition)
mat_combat <- mat_combat[, ord, drop = FALSE]
ann <- ann[colnames(mat_combat), , drop = FALSE]

n_ad <- sum(ann$Condition == "AD")
gaps <- if (n_ad > 0 && n_ad < ncol(mat_combat)) n_ad else NULL

# ============================================================
# 7) Row Z-score scaling (for visualization only)
# ============================================================
mat_z <- t(scale(t(mat_combat)))
mat_z[is.na(mat_z)] <- 0
write.csv(mat_z, "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\Data\\expression_data\\combat_zscore_5_datasets.csv")
# ============================================================
# 8) Color settings
# ============================================================
hm_cols <- colorRampPalette(c("deepskyblue2", "#F7F7F7", "deeppink2"))(201)

lim <- 3
mat_z[mat_z >  lim] <-  lim
mat_z[mat_z < -lim] <- -lim
hm_breaks <- seq(-lim, lim, length.out = length(hm_cols) + 1)

ann_colors <- list(
  Condition = c(AD = "aquamarine2", Control = "darkorange2")
)
library(extrafont)

# Set a universal font for all plots
universal_font <- "Arial"

# Apply globally for ggplot2
theme_set(theme_bw(base_size = 18, base_family = universal_font))
# ============================================================
# 9) Heatmap (no dendrograms)
# ============================================================
jpeg("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\heatmap_combat_zscore_5datasets.jpeg",
    width = 3400, height = 1800, res = 600, bg = "white")
pheatmap(mat_z,
         color      = hm_cols,
         breaks     = hm_breaks,
         annotation_col    = ann,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         cluster_cols  = FALSE,
         cluster_rows  = FALSE,
         gaps_col      = gaps,
         border_color  = NA,
         fontsize_row  = 6,
         fontface_row = "bold",
         fontsize      = 6,
         fontfamily    = universal_font
         )
dev.off()

cat("Done! Saved: heatmap_combat_zscore_5datasets.png\n")

library(grid)

gene_order <- rownames(mat_z)
n_genes    <- length(gene_order)

# ---- TUNABLES: nudge until the wedge meets the heatmap edge -----------------
BODY_R <- 0.545   # right edge of the heatmap BODY (not the gtable)
BODY_T <- 0.905   # top edge of body (just under the Condition annotation bar)
BODY_B <- 0.150   # bottom edge of body
COL_X  <- 0.735   # left edge of the gene column
TXT_T  <- 0.965; TXT_B <- 0.135
GUIDES <- TRUE    # set FALSE once aligned

outfile <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\heatmap_combat_zscore_5datasets_magnified.jpeg"

ph <- pheatmap(mat_z,
               color = hm_cols, breaks = hm_breaks,
               annotation_col = ann, annotation_colors = ann_colors,
               show_colnames = FALSE, show_rownames = FALSE,
               cluster_cols = FALSE, cluster_rows = FALSE,
               gaps_col = gaps, border_color = NA,
               legend = FALSE, annotation_legend = FALSE,
               fontsize = 6, fontfamily = universal_font,
               silent = TRUE)

jpeg(outfile, width = 4800, height = 3000, res = 600, bg = "white")
grid.newpage()

pushViewport(viewport(x = 0, width = 0.60, just = "left"))
grid.draw(ph$gtable)
popViewport()

# --- single magnifier wedge --------------------------------------------------
grid.polygon(x = c(BODY_R, BODY_R, COL_X - 0.012, COL_X - 0.012),
             y = c(BODY_T, BODY_B, TXT_B - 0.01, TXT_T + 0.01),
             gp = gpar(fill = rgb(0.55, 0.55, 0.55, 0.16), col = NA))
grid.lines(x = c(BODY_R, COL_X - 0.012), y = c(BODY_T, TXT_T + 0.01),
           arrow = arrow(ends = "last", length = unit(0.5, "mm")),
           gp = gpar(col = "grey45", lwd = 0.7))
grid.lines(x = c(BODY_R, COL_X - 0.012), y = c(BODY_B, TXT_B - 0.01),
           arrow = arrow(ends = "last", length = unit(0.5, "mm")),
           gp = gpar(col = "grey45", lwd = 0.7))

# --- gene labels (all 37, current order, not bold) ---------------------------
ys <- seq(TXT_T, TXT_B, length.out = n_genes)
grid.text(gene_order, x = COL_X, y = ys, just = "left",
          gp = gpar(fontsize = 8, fontface = "italic",
                    fontfamily = universal_font))

# --- embedded key ------------------------------------------------------------
grid.raster(matrix(hm_cols, nrow = 1), x = 0.06, y = 0.075, just = "left",
            width = 0.20, height = 0.024, interpolate = TRUE)
grid.text(c("-3", "0", "+3"), x = c(0.06, 0.16, 0.26), y = 0.042,
          gp = gpar(fontsize = 6, fontfamily = universal_font))
grid.text("Row z-score", x = 0.06, y = 0.108, just = "left",
          gp = gpar(fontsize = 6, fontfamily = universal_font))
grid.rect(x = c(0.34, 0.44), y = 0.075, width = 0.016, height = 0.024,
          gp = gpar(fill = c("aquamarine2", "darkorange2"), col = NA))
grid.text(c("AD", "Control"), x = c(0.360, 0.460), y = 0.075, just = "left",
          gp = gpar(fontsize = 6, fontfamily = universal_font))

if (GUIDES) {
  grid.lines(x = c(BODY_R, BODY_R), y = c(0, 1), gp = gpar(col = "red", lty = 3, lwd = 0.5))
  grid.lines(x = c(0, 1), y = c(BODY_T, BODY_T), gp = gpar(col = "red", lty = 3, lwd = 0.5))
  grid.lines(x = c(0, 1), y = c(BODY_B, BODY_B), gp = gpar(col = "red", lty = 3, lwd = 0.5))
}
dev.off()
