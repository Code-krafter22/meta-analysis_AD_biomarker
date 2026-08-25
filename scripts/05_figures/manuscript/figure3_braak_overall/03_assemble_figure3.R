# HISTORICAL PROVENANCE COPY — DO NOT RUN.
# Use 03_assemble_reproducible_figure3.R through run_figure3.R.
library(magick)

# Read the images
img_a <- image_read("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Final_Waterfall_Braak_AllGenes_FC.jpeg")
img_b <- image_read("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Forest_Plot_OR_Braak_SideBySide.jpeg")

# Scale both to the same height
target_height <- max(image_info(img_a)$height, image_info(img_b)$height)
img_a <- image_scale(img_a, paste0("x", target_height))
img_b <- image_scale(img_b, paste0("x", target_height))

# Create white label strips above each image
label_height <- 150
label_a <- image_blank(width = image_info(img_a)$width, height = label_height, color = "white")
label_a <- image_annotate(label_a, "(A)", size = 150, weight = 700,
                          color = "black", location = "+30+10")

label_b <- image_blank(width = image_info(img_b)$width, height = label_height, color = "white")
label_b <- image_annotate(label_b, "(B)", size = 150, weight = 700,
                          color = "black", location = "+30+10")

# Stack label on top of each panel
panel_a <- image_append(c(label_a, img_a), stack = TRUE)
panel_b <- image_append(c(label_b, img_b), stack = TRUE)

# Add a small white spacer between panels
spacer <- image_blank(width = 40, height = image_info(panel_a)$height, color = "white")

# Merge horizontally
combined <- image_append(c(panel_a, spacer, panel_b), stack = FALSE)

# Save outputs
image_write(combined, path = "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Figure_Combined_AB.pdf",
            format = "pdf", density = "900x900")

image_write(combined, path = "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\Figure_Combined_AB.tiff",
            format = "tiff", density = "900x900")

cat("Done! Saved as Figure_Combined_AB.tiff and Figure_Combined_AB.png\n")
