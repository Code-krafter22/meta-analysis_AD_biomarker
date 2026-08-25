# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use
# 04_assemble_reproducible_figure1.R through run_figure1.R.
library(magick)

# ----------------------------
# 1) Load images
# ----------------------------
venn    <- image_read("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\VennDiagram_MetaAnalysis.jpeg")
heatmap <- image_read("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\heatmap_combat_zscore_5datasets_magnified.jpeg")
volcano <- image_read("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\volcano_outputs\\volcano_all_three.jpeg")

# ----------------------------
# 2) Trim whitespace
# ----------------------------
venn    <- image_trim(venn)
heatmap <- image_trim(heatmap)
volcano <- image_trim(volcano)

# ----------------------------
# 3) Resize — Venn smaller, Heatmap larger, Volcano wider
# ----------------------------
venn    <- image_resize(venn,    "1000x1000")
heatmap <- image_resize(heatmap, "3800x3800")
volcano <- image_resize(volcano, "4000x1300")   # wider and taller

# ----------------------------
# 4) Add minimal padding
# ----------------------------
add_padding <- function(img, pad = 30) {
  info  <- image_info(img)
  new_w <- info$width + 2 * pad
  new_h <- info$height + 2 * pad
  image_extent(img, paste0(new_w, "x", new_h), gravity = "center", color = "white")
}

venn    <- add_padding(venn, 20)
heatmap <- add_padding(heatmap, 20)
volcano <- add_padding(volcano, 5)    # minimal padding on volcano

# ----------------------------
# 5) Equalize heights for top row
# ----------------------------
top_h <- max(image_info(venn)$height, image_info(heatmap)$height)

venn    <- image_extent(venn,    paste0(image_info(venn)$width, "x", top_h),
                        gravity = "center", color = "white")
heatmap <- image_extent(heatmap, paste0(image_info(heatmap)$width, "x", top_h),
                        gravity = "center", color = "white")

# ----------------------------
# 6) Build top row with smaller gap
# ----------------------------
gap_px <- 60       # reduced from 120
h_gap  <- image_blank(gap_px, top_h, color = "white")
row1   <- image_append(c(venn, h_gap, heatmap), stack = FALSE)

row1_w <- image_info(row1)$width

# ----------------------------
# 7) Stretch volcano to full top-row width (no side white space)
# ----------------------------
volcano <- image_resize(volcano, paste0(row1_w, "x"))   # scale width to match, height auto
# Crop height if it got too tall, or extent if slightly short
volcano <- image_extent(volcano, paste0(row1_w, "x", image_info(volcano)$height),
                        gravity = "center", color = "white")

# ----------------------------
# 8) Create label strips (thinner)
# ----------------------------
label_h <- 100     # reduced from 120
font_sz <- 85

venn_w <- image_info(venn)$width

label_top <- image_blank(row1_w, label_h, color = "white")
label_top <- image_annotate(label_top, text = "(A)", size = font_sz, font = "Arial",
                            weight = 700, color = "black",
                            location = "+10+5")

x_b <- venn_w + gap_px + 10
label_top <- image_annotate(label_top, text = "(B)", size = font_sz, font = "Arial",
                            weight = 700, color = "black",
                            location = paste0("+", round(x_b), "+5"))

label_bot <- image_blank(row1_w, label_h, color = "white")
label_bot <- image_annotate(label_bot, text = "(C)", size = font_sz, font = "Arial",
                            weight = 700, color = "black",
                            location = "+10+5")

# ----------------------------
# 9) Stack with minimal vertical gap
# ----------------------------
v_gap <- image_blank(row1_w, 30, color = "white")   # reduced from 80

figure <- image_append(c(label_top, row1, v_gap,
                         label_bot, volcano),
                       stack = TRUE)

figure <- image_background(figure, color = "white", flatten = TRUE)

# ----------------------------
# 10) Save
# ----------------------------

figure <- image_resize(figure, "4800x5000!")
output_path <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Figure_37_genes_new1.pdf"
image_write(
  figure,
  path    = output_path,
  format  = "pdf",
  density = "900x900"
)
cat("Saved:", output_path, "\n")

library(magick)

figure <- image_read("heatmap_combat_zscore_5datasets_magnified.jpeg")
figure <- image_resize(figure, "4800x")   # width-only: aspect preserved

base <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Figure_37_genes_new1"

image_write(figure, path = paste0(base, ".pdf"),
            format = "pdf", density = "900x900")
image_write(figure, path = paste0(base, ".png"),
            format = "png", quality = 600, density = "900x900")

cat("Saved:", base, ".pdf / .jpeg\n")
