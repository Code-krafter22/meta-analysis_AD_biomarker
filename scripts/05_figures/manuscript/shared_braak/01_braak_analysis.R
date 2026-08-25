# HISTORICAL PROVENANCE COPY — DO NOT RUN.
# Use shared_braak/00_prepare_braak_analysis.R and the current Figure 3/4
# runners. This file retains the original absolute paths and duplicate blocks.
#================================================================
# BRAAK STAGE ASSOCIATION WITH GENE EXPRESSION
# 
# Purpose: Test whether the 37-gene consensus signature is
#          associated with Braak neuropathological staging
#          across 4 combined transcriptomic datasets.
#
# Pipeline: Load data → Standardize labels → ComBat batch 
#           correction → Z-score normalization → ANCOVA 
#           (Braak Low vs High) → Per-region ANCOVA → Plots
#
# ANOVA & Boxplots: ComBat + z-score normalized
# Waterfall (Fold Change): ComBat log2-CPM only
# Waterfall asterisks: From z-score ANCOVA
#
# Input:  BRAAK_DATA_4studies.xlsx (log2-CPM, 4 studies)
# Output: Tables (CSV) + Figures (PNG)
#================================================================

# ============================================================
# 0) SETUP
# ============================================================
library(readxl)
library(tidyverse)

library(broom)
library(sva)
library(ggplot2)
library(RColorBrewer)

# File paths
infile <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\Data\\stats\\BRAAK_DATA_4studies.xlsx"
outdir <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) LOAD DATA
# ============================================================
df <- read_excel(infile)
glimpse(df)

gene_cols <- setdiff(colnames(df), c("SampleID", "SEX", "AGE", "BRAAK1", "PATHDX", "REGION", "BATCH"))
cat("Number of genes:", length(gene_cols), "\n")

# ============================================================
# 2) STANDARDIZE LABELS
# ============================================================
df$SEX <- toupper(df$SEX)
df$SEX <- recode(df$SEX, "MALE" = "M", "FEMALE" = "F")

df$PATHDX <- toupper(df$PATHDX)
df$PATHDX <- recode(df$PATHDX, "ALZHEIMERS" = "AD", "HEALTHY" = "CON")

cat("\n--- Label verification ---\n")
cat("SEX levels:\n"); print(table(df$SEX))
cat("PATHDX levels:\n"); print(table(df$PATHDX))
cat("BATCH levels:\n"); print(table(df$BATCH))

# ============================================================
# 3) PCA BEFORE COMBAT
# ============================================================
expr_mat <- as.matrix(df[, gene_cols])

pca_before <- prcomp(expr_mat, scale. = TRUE)
pca_before_df <- as.data.frame(pca_before$x) %>%
  bind_cols(df[, c("BATCH", "PATHDX", "REGION", "BRAAK1")])

p_pca_before <- ggplot(pca_before_df, aes(PC1, PC2, color = BATCH, shape = PATHDX)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Before Batch Correction") +
  theme_minimal()

# ============================================================
# 4) Z-SCORE NORMALIZATION (before ComBat)
#    - Standardize each gene to mean=0, sd=1 FIRST
# ============================================================
expr_z <- scale(expr_mat)
rownames(expr_z) <- df$SampleID

# ============================================================
# 5) COMBAT BATCH CORRECTION (on z-scored data)
#    - Preserves biological signal (PATHDX, AGE, SEX)
# ============================================================
mod <- model.matrix(~ PATHDX + AGE + SEX, data = df)
expr_combat_z <- ComBat(dat = t(expr_z), batch = df$BATCH, mod = mod)
expr_combat_z <- t(expr_combat_z)
rownames(expr_combat_z) <- df$SampleID

# ComBat on original log2-CPM (for fold change calculations)
expr_combat <- ComBat(dat = t(expr_mat), batch = df$BATCH, mod = mod)
expr_combat <- t(expr_combat)
rownames(expr_combat) <- df$SampleID

# ============================================================
# 6) PCA AFTER COMBAT
# ============================================================
pca_after <- prcomp(expr_combat, scale. = TRUE)
pca_after_df <- as.data.frame(pca_after$x) %>%
  bind_cols(df[, c("BATCH", "PATHDX", "REGION")])

p_pca_after <- ggplot(pca_after_df, aes(PC1, PC2, color = BATCH, shape = PATHDX)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA After Batch Correction (ComBat)") +
  theme_minimal()

# ============================================================
# 7) CREATE BRAAK GROUPING VARIABLES
# ============================================================
df <- df %>%
  mutate(Braak_group = factor(BRAAK1, levels = c("0", "1", "2", "3", "4", "5", "6")))

df <- df %>%
  mutate(Braak_bin = ifelse(BRAAK1 <= 3, "Low(<=3)", "High(>3)")) %>%
  mutate(Braak_bin = factor(Braak_bin, levels = c("Low(<=3)", "High(>3)")))

cat("\nBraak group distribution:\n"); print(table(df$Braak_group))
cat("\nBraak binary distribution:\n"); print(table(df$Braak_bin))

# ============================================================
# 8) DEFINE CONSENSUS GENE LIST
# ============================================================
genes <- c("CHML","KANK2","PRELP","HSPB1","PIK3R5","AEBP1","GFAP","ITFG1","TCEA3",
           "CCDC102A","ELOVL4","KLF15","NAP1L5","NEUROD6","HPRT1","NRN1","MAS1",
           "RPH3A","HSPB7","NUPR1","ADAM33","COL25A1","RAB3C","SERPINI1","STAT4",
           "SCG2","TMPRSS5","RGS4","PRX","RAB3B","NCALD","OPN3","CLDN9","TRIM36",
           "GAD2","MRGPRF","GAD1")

genes_use <- intersect(genes, colnames(expr_combat))
missing_genes <- setdiff(genes, genes_use)
if (length(missing_genes) > 0) {
  message("Missing genes (not in data): ", paste(missing_genes, collapse = ", "))
}
cat("Using", length(genes_use), "of", length(genes), "genes\n")

# ============================================================
# 9) ANCOVA: BRAAK LOW vs HIGH (z-score)
#    - Adjusts for REGION, AGE, SEX
#    - PATHDX excluded (confounded with Braak)
# ============================================================
aov_results <- map_dfr(colnames(expr_combat_z), function(g) {
  m <- aov(expr_combat_z[df$SampleID, g] ~ Braak_bin + REGION + AGE + SEX, data = df)
  tidy(m) %>%
    filter(term == "Braak_bin") %>%
    mutate(Gene = g)
}) %>%
  mutate(FDR = p.adjust(p.value, "fdr")) %>%
  arrange(FDR)

cat("\n=== ANCOVA RESULTS (Top 10) ===\n")
print(head(aov_results, 10))
cat("Genes with FDR <= 0.05:", sum(aov_results$FDR <= 0.05), "out of", nrow(aov_results), "\n")

write.csv(aov_results, file.path(outdir, "Final_ANCOVA_Braak_results_new.csv"), row.names = FALSE)

# ============================================================
# 10) PER-REGION ANCOVA (z-score)
#     - Adjusts for AGE, SEX
#     - FDR corrected within each region
# ============================================================
region_models <- df %>%
  group_split(REGION) %>%
  map_dfr(function(sub_df) {
    reg <- unique(sub_df$REGION)
    map_dfr(genes_use, function(g) {
      y <- expr_combat_z[sub_df$SampleID, g]
      m <- aov(y ~ Braak_bin + AGE + SEX, data = sub_df)
      tidy(m) %>%
        filter(term == "Braak_bin") %>%
        transmute(Gene = g, Region = reg, statistic, p.value)
    })
  }) %>%
  group_by(Region) %>%
  mutate(FDR = p.adjust(p.value, "fdr")) %>%
  ungroup() %>%
  arrange(Region, FDR)

write.csv(region_models, file.path(outdir, "Final_ANCOVA_Braak_per_region_new.csv"), row.names = FALSE)

# ============================================================
# 10b) SPEARMAN CORRELATION: GENE vs BRAAK STAGE (z-score)
#      - Tests monotonic association with numeric Braak (0-6)
#      - No covariates (non-parametric)
# ============================================================
corr_results <- map_dfr(genes_use, function(gene) {
  test <- cor.test(expr_combat_z[df$SampleID, gene], df$BRAAK1, method = "spearman")
  tibble(Gene = gene, rho = test$estimate, p_value = test$p.value)
})

corr_results <- corr_results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

cat("\n=== SPEARMAN CORRELATION (Top 10) ===\n")
print(head(corr_results, 10))
cat("Genes with FDR <= 0.05:", sum(corr_results$FDR <= 0.05), "out of", nrow(corr_results), "\n")

write.csv(corr_results, file.path(outdir, "Final_Spearman_Braak_zscore_new.csv"), row.names = FALSE)

# Bar plot of Spearman rho
corr_plot_df <- corr_results %>%
  mutate(Gene = fct_reorder(Gene, rho),
         direction = ifelse(rho > 0, "Positive", "Negative"),
         sig_label = case_when(FDR <= 0.001 ~ "***",
                               FDR <= 0.01 ~ "**",
                               FDR <= 0.05 ~ "*",
                               TRUE ~ ""))

p_corr <- ggplot(corr_plot_df, aes(x = Gene, y = rho, fill = direction)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = sig_label,
                hjust = ifelse(rho >= 0, -0.3, 1.3)), size = 2) +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "aquamarine3", "Negative" = "lightsalmon3"),
                    labels = c("Positive" = "Positive correlation", "Negative" = "Negative correlation")) +
  labs(x = "", y = "Spearman rho",
       fill = "") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10),
      axis.title.x = element_text(size = 6, face = "bold"),
      axis.text.y = element_text(size = 5, face = "bold"),
      axis.text.x = element_text(size = 8),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "top",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, -5, 0),
      plot.margin = margin(2, 2, 2, 2))

ggsave(file.path(outdir, "Final_Spearman_Braak_zscore.jpeg"), p_corr,
       width = 3.5, height = 3, dpi = 600, bg = "white")
# ============================================================
# 11) PREPARE PLOTTING DATA
# ============================================================
plot_df <- df %>%
  select(SampleID, BRAAK1, Braak_bin, REGION, AGE, SEX) %>%
  mutate(SampleID = as.character(SampleID))

# Z-score data for boxplots
expr_long_genes_z <- as.data.frame(expr_combat_z[df$SampleID, genes_use, drop = FALSE]) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Gene", values_to = "Expression") %>%
  left_join(plot_df, by = "SampleID")

# ComBat log2-CPM data for waterfall fold change
expr_long_genes <- as.data.frame(expr_combat[df$SampleID, genes_use, drop = FALSE]) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Gene", values_to = "Expression") %>%
  left_join(plot_df, by = "SampleID")

# ============================================================
# 12) BOXPLOT: TOP 20 GENES (z-score)
# ============================================================
top20_genes <- aov_results %>%
  arrange(FDR) %>%
  slice_head(n = 20) %>%
  pull(Gene)

expr_long_top20_z <- expr_long_genes_z %>%
  filter(Gene %in% top20_genes) %>%
  mutate(Gene = factor(Gene, levels = top20_genes))

p_top20 <- ggplot(expr_long_top20_z, aes(x = Braak_bin, y = Expression, fill = Braak_bin)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.7) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Low(<=3)" = "yellowgreen", "High(>3)" = "tomato4")) +
  labs(x = "Braak group", y = "ComBat-corrected expression (z-score)") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "Final_Top20_Braak_Low_vs_High.jpeg"), p_top20,
       width = 14, height = 10, dpi = 600, bg = "white")

# ============================================================
# 13) WATERFALL PLOT: FOLD CHANGE (ComBat log2-CPM)
#     - Bars: fold change from ComBat log2-CPM
#     - Asterisks: from per-region ANCOVA (z-score)
# ============================================================

# Gene order based on overall log2FC across all regions
gene_order <- expr_long_genes %>%
  group_by(Gene, Braak_bin) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Braak_bin, values_from = mean_expr) %>%
  mutate(log2FC = `High(>3)` - `Low(<=3)`) %>%
  arrange(desc(log2FC)) %>%
  pull(Gene)

# Plot 1: Waterfall without significance markers
waterfall_df <- expr_long_genes %>%
  group_by(Gene, Region = REGION, Braak_bin) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Braak_bin, values_from = mean_expr) %>%
  mutate(log2FC = `High(>3)` - `Low(<=3)`,
         FC = ifelse(log2FC >= 0, 2^log2FC, -(2^abs(log2FC))),
         Gene = factor(Gene, levels = gene_order),
         direction = ifelse(log2FC > 0, "Up", "Down"))

gene_order <- expr_long_genes %>%
  group_by(Gene, Braak_bin) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Braak_bin, values_from = mean_expr) %>%
  mutate(log2FC = `High(>3)` - `Low(<=3)`,
         direction = ifelse(log2FC > 0, "Up", "Down")) %>%
  left_join(region_models %>%
              group_by(Gene) %>%
              summarise(min_FDR = min(FDR), .groups = "drop"),
            by = "Gene") %>%
  mutate(sig_tier = case_when(min_FDR <= 0.001 ~ 1,   # ***
                              min_FDR <= 0.01  ~ 2,   # **
                              min_FDR <= 0.05  ~ 3,   # *
                              TRUE ~ 4)) %>%           # ns
  arrange(direction, sig_tier, desc(abs(log2FC))) %>%
  pull(Gene)


# region -specifc analysis figure
# ============================================================
# 13) WATERFALL PLOT: FOLD CHANGE (ComBat log2-CPM)
# ============================================================

# Compute per-region fold changes
waterfall_df <- expr_long_genes %>%
  group_by(Gene, Region = REGION, Braak_bin) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Braak_bin, values_from = mean_expr) %>%
  mutate(log2FC = `High(>3)` - `Low(<=3)`,
         FC = ifelse(log2FC >= 0, 2^log2FC, -(2^abs(log2FC))),
         direction = ifelse(log2FC > 0, "Up", "Down"))

# Add significance from per-region ANCOVA
waterfall_sig_df <- waterfall_df %>%
  left_join(region_models %>% select(Gene, Region, FDR), by = c("Gene", "Region")) %>%
  mutate(sig_label = case_when(FDR <= 0.001 ~ "***",
                               FDR <= 0.01  ~ "**",
                               FDR <= 0.05  ~ "*",
                               TRUE ~ ""),
         sig_tier = case_when(FDR <= 0.001 ~ 1,
                              FDR <= 0.01  ~ 2,
                              FDR <= 0.05  ~ 3,
                              TRUE ~ 4),
         label_vjust = ifelse(FC >= 0, -0.5, 1.5))

# Gene order: direction → sig_tier → FC magnitude
# Use overall direction and best (min) sig_tier across regions per gene
gene_summary <- waterfall_sig_df %>%
  group_by(Gene) %>%
  summarise(overall_log2FC = mean(log2FC),
            best_sig = min(sig_tier),
            .groups = "drop") %>%
  mutate(direction = ifelse(overall_log2FC > 0, "Up", "Down"))

# Up genes: most significant LAST (toward center)
up_genes <- gene_summary %>%
  filter(direction == "Up") %>%
  arrange(desc(best_sig), abs(overall_log2FC)) %>%
  pull(Gene)

# Down genes: most significant FIRST (toward center)
down_genes <- gene_summary %>%
  filter(direction == "Down") %>%
  arrange(desc(best_sig), abs(overall_log2FC)) %>%
  pull(Gene)

gene_order <- c(up_genes, down_genes)
# Apply gene order
waterfall_sig_df <- waterfall_sig_df %>%
  mutate(Gene = factor(Gene, levels = gene_order))

# Plot 2: Waterfall with ANCOVA significance markers
p_waterfall2 <- ggplot(waterfall_sig_df, aes(x = Gene, y = FC, fill = direction)) +
  geom_col(width = 0.3) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = sig_label, vjust = label_vjust), size = 3.5) +
  facet_wrap(~ Region, scales = "free_y", ncol = 2) +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  scale_fill_manual(values = c("Up" = "skyblue", "Down" = "yellowgreen"),
                    labels = c("Up" = "Upregulated", "Down" = "Downregulated")) +
  labs(x = "", y = "Fold Change", fill = "") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "top",
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.2, "cm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

ggsave(file.path(outdir, "Final_Waterfall_Braak_by_region_FC_with_sig.pdf"), p_waterfall2,
       width =12.5, height = 6, dpi = 900, bg = "white")


# Add significance markers from per-region ANCOVA (z-score)
waterfall_sig_df <- waterfall_df %>%
  left_join(region_models %>% select(Gene, Region, FDR), by = c("Gene", "Region")) %>%
  mutate(sig_label = case_when(FDR <= 0.001 ~ "***",
                               FDR <= 0.01  ~ "**",
                               FDR <= 0.05  ~ "*",
                               TRUE ~ ""),
         label_vjust = ifelse(FC >= 0, -0.5, 1.5))%>%
  arrange(direction, desc(sig_tier), desc(abs(log2FC))) %>%
  mutate(Gene = factor(Gene, levels = unique(Gene)))


# p_waterfall1 <- ggplot(waterfall_df, aes(x = Gene, y = FC, fill = direction)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  facet_wrap(~ Region, scales = "free_y", ncol = 2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("Up" = "tomato4", "Down" = "yellowgreen"),
                    labels = c("Up" = "Upregulated", "Down" = "Downregulated")) +
  labs(title = "Fold change: Braak High(>3) vs Low(≤3) by Region",
       x = "", y = "Fold Change (High / Low)",
       fill = "") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "top")

# ggsave(file.path(outdir, "Final_Waterfall_Braak_by_region_FC.png"), p_waterfall1,
       # width = 16, height = 10, dpi = 600, bg = "white")

# Plot 2: Waterfall with ANCOVA significance markers
p_waterfall2 <- ggplot(waterfall_sig_df, aes(x = Gene, y = FC, fill = direction)) +
  geom_col(width = 0.3) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = sig_label, vjust = label_vjust), size = 2.5) +
  facet_wrap(~ Region, scales = "free_y", ncol = 2) +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  scale_fill_manual(values = c("Up" = "skyblue", "Down" = "yellowgreen"),
                    labels = c("Up" = "Upregulated", "Down" = "Downregulated")) +
  labs(x = "", y = "Fold Change", fill = "") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 9),
        legend.position = "top",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.2, "cm"),
        panel.spacing = unit(0.2, "cm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

ggsave(file.path(outdir, "Final_Waterfall_Braak_by_region_FC_with_sig.jpg"), p_waterfall2,
       width = 7.5, height = 4.5, dpi = 600, bg = "white")

# ============================================================
# 14) WATERFALL PLOT: ALL GENES – OVERALL FOLD CHANGE
#     - Bars: fold change from ComBat log2-CPM (High vs Low)
#     - Asterisks: from overall ANCOVA (z-score)
# ============================================================

# Compute overall fold change per gene (across all regions)
overall_fc_df <- expr_long_genes %>%
  group_by(Gene, Braak_bin) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Braak_bin, values_from = mean_expr) %>%
  mutate(log2FC = `High(>3)` - `Low(<=3)`,
         FC = ifelse(log2FC >= 0, 2^log2FC, -(2^abs(log2FC))),
         direction = ifelse(log2FC > 0, "Up", "Down"))

# Merge ANCOVA significance from overall results (section 9)
overall_fc_df <- overall_fc_df %>%
  left_join(aov_results %>% select(Gene, FDR), by = "Gene") %>%
  mutate(sig_label = case_when(FDR <= 0.001 ~ "***",
                               FDR <= 0.01  ~ "**",
                               FDR <= 0.05  ~ "*",
                               TRUE ~ ""),
         Gene = fct_reorder(Gene, log2FC),
         label_vjust = ifelse(FC >= 0, -0.5, 1.5))

p_waterfall_all <- ggplot(overall_fc_df, aes(x = Gene, y = FC, fill = direction)) +
  geom_col(width = 0.85) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = sig_label, vjust = label_vjust), size = 4) +
  scale_fill_manual(values = c("Up" = "tomato4", "Down" = "yellowgreen"),
                    labels = c("Up" = "Upregulated", "Down" = "Downregulated")) +
  labs(x = "", y = "Fold Change",
       fill = "") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 10, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "top")

ggsave(file.path(outdir, "Final_Waterfall_Braak_AllGenes_FC.jpg"), p_waterfall_all,
       width = 8, height=5.5, dpi = 600, bg = "white")
