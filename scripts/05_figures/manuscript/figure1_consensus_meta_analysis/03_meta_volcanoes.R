# =========================
# Volcano plots (3 sheets) + label ALL genes
# Input: volcanodata.xlsx
# Output: volcano_outputs/*.png
# =========================

# install.packages(c("readxl", "ggplot2", "ggrepel", "patchwork"))
# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use 03_meta_volcanoes_from_csv.R.
library(readxl)
library(ggplot2)
library(ggrepel)
library(patchwork)

infile <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\volcanodata.xlsx"          # <-- change path if needed
outdir <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\37-genes\\volcano_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

library(extrafont)

# Set a universal font for all plots
universal_font <- "Arial"

# Apply globally for ggplot2
theme_set(theme_bw(base_size = 18, base_family = universal_font))

make_volcano <- function(df, gene_col, x_col, p_col, title,
                         p_cutoff = 0.05, p_floor = 1e-300,
                         show_legend = FALSE) {
  
  d <- df[, c(gene_col, x_col, p_col)]
  colnames(d) <- c("GENE", "X", "P")
  
  d$GENE <- as.character(d$GENE)
  d$X <- suppressWarnings(as.numeric(d$X))
  d$P <- suppressWarnings(as.numeric(d$P))
  d <- d[!is.na(d$GENE) & !is.na(d$X) & !is.na(d$P), ]
  
  d$P[d$P < p_floor] <- p_floor
  d$neglog10P <- -log10(d$P)
  
  # Color groups
  d$status <- "NS"
  d$status[d$P < p_cutoff & d$X > 0] <- "Up"
  d$status[d$P < p_cutoff & d$X < 0] <- "Down"
  d$status <- factor(d$status, levels = c("Down", "NS", "Up"))
  
  p <- ggplot(d, aes(X, neglog10P)) +
    geom_point(aes(color = status), size = 6, alpha = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.7, linetype = "solid") +
    geom_hline(yintercept = -log10(p_cutoff), linewidth = 0.7, linetype = "dashed") +
    
    geom_text_repel(
      aes(label = GENE, color = status),
      size = 8.5,                  # ← increased from 7.0
      max.overlaps = Inf,
      box.padding = 1.0,           # ← slightly more padding to reduce overlap
      point.padding = 0.4,
      force = 5,
      force_pull = 0.5,
      direction = "both",
      min.segment.length = 0,
      segment.alpha = 0.4,
      max.iter = 2e6,
      max.time = 15,
      seed = 123,
      bg.color = "white",
      bg.r = 0.15,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = c(
      "Down" = "#1f77b4",
      "NS"   = "grey70",
      "Up"   = "#d62728"
    )) +
    labs(
      title = title,
      x = x_col,
      y = paste0("-log10(", p_col, ")"),
      color = ""
    ) +
    
    theme_classic(base_size = 18, base_family = universal_font) +
    theme(
      plot.title   = element_text(size = 48, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 42, face = "bold"),
      axis.title.y = element_text(size = 42, face = "bold"),
      axis.text.x  = element_text(size = 36, face = "bold"),
      axis.text.y  = element_text(size = 36, face = "bold"),
      axis.ticks   = element_line(linewidth = 0.9),
      
      # Legend: only show on the designated panel
      legend.position = if (show_legend) "right" else "none",
      legend.text  = element_text(size = 42, face = "bold"),
      legend.title = element_text(size = 42, face = "bold"),
      legend.key.size  = unit(1.5, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      
      plot.margin = margin(10, 15, 10, 10)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}


# -------------------------
# Sheet 1: meta_effects
# -------------------------
df1 <- read_xlsx(infile, sheet = "meta_effects")
p1 <- make_volcano(df1,
             gene_col = "SYMBOL",
             x_col    = "meta_logFC",
             p_col    = "meta_p",
             title    = "i",
             show_legend = FALSE)

# -------------------------
# Sheet 2: stouffer
# -------------------------
df2 <- read_xlsx(infile, sheet = "stouffer")
p2 <- make_volcano(df2,
             gene_col = "SYMBOL",
             x_col    = "Z",
             p_col    = "P",
             title    = "ii",
             show_legend = FALSE)

# -------------------------
# Sheet 3: metavolcanoR  ← legend only on this (rightmost) panel
# -------------------------
df3 <- read_xlsx(infile, sheet = "metavolcanoR")
p3 <- make_volcano(df3,
             gene_col = "SYMBOL",
             x_col    = "metafc",
             p_col    = "metap",
             title    = "iii",
             show_legend = TRUE)     # ← single legend here

# -------------------------
# Save individual plots
# -------------------------
ggsave(file.path(outdir, "volcano_meta_effects.jpeg"),
       plot = p1, width = 18, height = 10, dpi = 600, bg = "white")
ggsave(file.path(outdir, "volcano_stouffer.jpeg"),
       plot = p2, width = 18, height = 10, dpi = 600, bg = "white")
ggsave(file.path(outdir, "volcano_metavolcanoR.jpeg"),
       plot = p3, width = 18, height = 10, dpi = 600, bg = "white")

# -------------------------
# Combined plot: single legend from p3 on the right
# -------------------------
combined <- (p1 | p2 | p3) +
  plot_layout(widths = c(1, 1, 1.15))   # slightly wider panel iii to accommodate legend

combined_outfile <- file.path(outdir, "volcano_all_three.jpeg")
ggsave(combined_outfile, plot = combined, width = 30, height = 18, dpi = 600, bg = "white", limitsize = FALSE)

cat("Done. Saved plots in:", normalizePath(outdir), "\n")
