# HISTORICAL PROVENANCE COPY — DO NOT RUN.
# Use 02_reproducible_forest_plot.R through run_figure3.R.
#================================================================
# FOREST PLOT: ODDS RATIOS FROM LOGISTIC REGRESSION
# 
# Purpose: For each of the 37 consensus genes, fit a logistic
#          regression predicting Braak High(>3) vs Low(<=3),
#          adjusting for REGION, AGE, SEX. Extract odds ratios
#          and 95% CIs, then display as a forest plot.
#
# Pipeline: Load data → Standardize labels → ComBat batch 
#           correction → Z-score normalization → Per-gene 
#           logistic regression → FDR correction → Forest plot
#
# Input:  BRAAK_DATA_4studies.xlsx (log2-CPM, 4 studies)
# Output: Forest plot (PNG/JPEG) + OR table (CSV)
#================================================================

# ============================================================
# 0) SETUP
# ============================================================
library(readxl)
library(tidyverse)
library(sva)
library(broom)
library(ggplot2)

# File paths — UPDATE THESE TO YOUR LOCAL PATHS
infile  <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\Data\\stats\\BRAAK_DATA_4studies.xlsx"
outdir  <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) LOAD DATA
# ============================================================
df <- read_excel(infile)
gene_cols <- setdiff(colnames(df),
                     c("SampleID", "SEX", "AGE", "BRAAK1",
                       "PATHDX", "REGION", "BATCH"))
cat("Number of genes:", length(gene_cols), "\n")

# ============================================================
# 2) STANDARDIZE LABELS
# ============================================================
df$SEX    <- recode(toupper(df$SEX),    "MALE" = "M", "FEMALE" = "F")
df$PATHDX <- recode(toupper(df$PATHDX), "ALZHEIMERS" = "AD", "HEALTHY" = "CON")

# ============================================================
# 3) COMBAT BATCH CORRECTION → Z-SCORE
# ============================================================
expr_mat <- as.matrix(df[, gene_cols])

# Z-score first, then ComBat
expr_z <- scale(expr_mat)
rownames(expr_z) <- df$SampleID

mod <- model.matrix(~ PATHDX + AGE + SEX, data = df)
expr_combat_z <- t(ComBat(dat = t(expr_z), batch = df$BATCH, mod = mod))
rownames(expr_combat_z) <- df$SampleID

# ============================================================
# 4) CREATE BRAAK BINARY OUTCOME
# ============================================================
df <- df %>%
  mutate(Braak_bin = factor(ifelse(BRAAK1 <= 3, "Low", "High"),
                            levels = c("Low", "High")))

cat("\nBraak binary distribution:\n")
print(table(df$Braak_bin))

# ============================================================
# 5) DEFINE CONSENSUS GENE LIST
# ============================================================
genes <- c("CHML","KANK2","PRELP","HSPB1","PIK3R5","AEBP1","GFAP",
           "ITFG1","TCEA3","CCDC102A","ELOVL4","KLF15","NAP1L5",
           "NEUROD6","HPRT1","NRN1","MAS1","RPH3A","HSPB7","NUPR1",
           "ADAM33","COL25A1","RAB3C","SERPINI1","STAT4","SCG2",
           "TMPRSS5","RGS4","PRX","RAB3B","NCALD","OPN3","CLDN9",
           "TRIM36","GAD2","MRGPRF","GAD1")

genes_use <- intersect(genes, colnames(expr_combat_z))
cat("Using", length(genes_use), "of", length(genes), "genes\n")

# ============================================================
# 6) LOGISTIC REGRESSION: PER-GENE ODDS RATIOS
#    - Outcome: Braak_bin (High vs Low)
#    - Predictor: each gene (ComBat + z-scored)
#    - Covariates: REGION, AGE, SEX
# ============================================================
or_results <- map_dfr(genes_use, function(g) {
  
  tmp <- df %>%
    mutate(gene_expr = expr_combat_z[SampleID, g])
  
  fit <- glm(Braak_bin ~ gene_expr + REGION + AGE + SEX,
             data   = tmp,
             family = binomial(link = "logit"))
  
  # Extract the gene term
  tidy_fit <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  gene_row <- tidy_fit %>% filter(term == "gene_expr")
  
  tibble(
    Gene     = g,
    OR       = gene_row$estimate,
    CI_lower = gene_row$conf.low,
    CI_upper = gene_row$conf.high,
    p_value  = gene_row$p.value,
    z_stat   = gene_row$statistic
  )
})

# FDR correction
or_results <- or_results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

# Significance labels
or_results <- or_results %>%
  mutate(sig_label = case_when(
    FDR <= 0.001 ~ "***",
    FDR <= 0.01  ~ "**",
    FDR <= 0.05  ~ "*",
    TRUE         ~ ""
  ))

cat("\n=== LOGISTIC REGRESSION RESULTS (Top 10) ===\n")
print(head(or_results, 10))
cat("\nGenes with FDR <= 0.05:", sum(or_results$FDR <= 0.05),
    "out of", nrow(or_results), "\n")

# Save table
write.csv(or_results,
          file.path(outdir, "Forest_Plot_OR_Braak_LogisticRegression.csv"),
          row.names = FALSE)

# ============================================================
# 7) FOREST PLOT
# ============================================================

library(patchwork)   # for side-by-side layout
 
# --- Classify and prepare ---
plot_df <- or_results %>%
  mutate(
    significant = FDR <= 0.05,
    direction   = ifelse(OR >= 1, "Up (Risk)", "Down (Protective)")
  )
 
# Split into two data frames
df_up   <- plot_df %>% filter(direction == "Up (Risk)")   %>% mutate(Gene = fct_reorder(Gene, OR))
df_down <- plot_df %>% filter(direction == "Down (Protective)") %>% mutate(Gene = fct_reorder(Gene, OR))
 
# --- Per-panel x-axis ranges (no wasted space) ---
# UP panel: data is all >= 1, show a bit left of 1 and stretch right
x_lo_up <- max(0, min(df_up$CI_lower, na.rm = TRUE) * 0.85)
x_hi_up <- 10
 
# DOWN panel: data is all < 1, stretch left for stars and just past 1 on right
x_lo_down <- max(0, min(df_down$CI_lower, na.rm = TRUE) * 0.70)
x_hi_down <- 1.2  # slightly past the dashed line — no empty 2-8 range
 
# --- Shared theme ---
forest_theme <- theme_minimal(base_size = 18) +
  theme(
    plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.y      = element_text(size = 16, face = "bold"),
    axis.text.x      = element_text(size = 15),
    axis.title.x     = element_text(size = 18, face = "bold"),
    axis.title.y     = element_blank(),
    legend.position  = "none",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )
 
# --- LEFT PANEL: Up-regulated (OR >= 1, risk) ---
p_up <- ggplot(df_up, aes(x = OR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.3, linewidth = 0.8, color = "tomato3") +
  geom_point(aes(shape = significant), size = 4, color = "tomato3") +
  geom_text(aes(x = CI_upper, label = sig_label),
            hjust = -0.3, size = 5.5, fontface = "bold") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16)) +
  coord_cartesian(xlim = c(x_lo_up, x_hi_up)) +
  labs(x = "Odds Ratio (95% CI)") +
  forest_theme
 
ggsave(file.path(outdir, "Forest_Plot_OR_Braak_up.jpeg"),
       p_up, width = 7.5, height = 6, dpi = 600, bg = "white")

# --- RIGHT PANEL: Down-regulated (OR < 1, protective) ---
p_down <- ggplot(df_down, aes(x = OR, y = Gene)) +
  geom_vline(xintercept = 1.2, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.3, linewidth = 0.8, color = "steelblue3") +
  geom_point(aes(shape = significant), size = 4, color = "steelblue3") +
  geom_text(aes(x = CI_upper, label = sig_label),
            hjust = -0.3, size = 5.5, fontface = "bold") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16)) +
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
  coord_cartesian(xlim = c(x_lo_down, x_hi_down)) +
  labs(x = "Odds Ratio (95% CI)") +
  forest_theme
 
# --- Combine side by side ---
p_combined <- p_up + p_down +
  plot_annotation(
    #title    = "Forest Plot: Gene-level Odds Ratios for Braak High(>3) vs Low(\u22643)",
    #subtitle = "Logistic regression adjusted for Region, Age, Sex  |  * FDR<0.05  ** FDR<0.01  *** FDR<0.001",
    theme = theme(
      plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "grey30")
    )
  )
 
# Save — compact size per user preference
ggsave(file.path(outdir, "Forest_Plot_OR_Braak_SideBySide.pdf"),
       p_combined, width = 7.5, height = 6, dpi = 600, bg = "white")

cat("\n✓ Forest plot saved.\n")
cat("✓ OR table saved.\n")
cat("Done!\n")
