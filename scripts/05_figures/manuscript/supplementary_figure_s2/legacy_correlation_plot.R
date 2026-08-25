#================================================================
# CORRELATION OF BRAAK AND GENE EXPRESSION
#================================================================
library(pacman)
pacman::p_load(tidyverse, broom, readxl, sva)

# Load the data 
df <- read_xlsx("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\Data\\stats\\BRAAK_DATA_4studies.xlsx")

# Define gene columns
gene_cols <- setdiff(colnames(df), c("SampleID", "BRAAK1", "PATHDX", "AGE", "SEX", "REGION", "BATCH"))

# ---- ComBat batch correction ----
expr_mat <- as.matrix(df[, gene_cols])
rownames(expr_mat) <- df$SampleID

# Standardize PATHDX
df$PATHDX <- toupper(df$PATHDX)
df$PATHDX <- recode(df$PATHDX, "ALZHEIMERS" = "AD", "HEALTHY" = "CON")

# Standardize SEX
df$SEX <- toupper(df$SEX)
df$SEX <- recode(df$SEX, "FEMALE" = "F", "MALE" = "M")

# Verify
print(table(df$BATCH, df$PATHDX))
print(table(df$BATCH, df$SEX))

# Now run ComBat
mod <- model.matrix(~ PATHDX + AGE + SEX, data = df)
expr_combat <- t(ComBat(dat = t(expr_mat), batch = df$BATCH, mod = mod))

# Replace raw expression with ComBat-corrected values
df[, gene_cols] <- as.data.frame(expr_combat)

# ---- Spearman correlation for each gene vs BRAAK1 ----
corr_results <- map_dfr(gene_cols, function(gene) {
  test <- cor.test(df[[gene]], df$BRAAK1, method = "spearman")
  tibble(Gene = gene, rho = test$estimate, p_value = test$p.value)
})

corr_results <- corr_results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

head(corr_results)

# ---- Bar plot ----
library(ggplot2)

top_corr <- corr_results %>%
  arrange(desc(rho)) %>%
  slice_head(n = 37) %>%
  bind_rows(corr_results %>% arrange(rho) %>% slice_head(n = 37))

pos_col <- "aquamarine3"
neg_col <- "lightsalmon3"

p <- ggplot(top_corr, aes(x = reorder(Gene, rho), y = rho, fill = rho > 0)) +
  geom_col(width = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = neg_col, "TRUE" = pos_col)) +
  labs(title = "Top Genes Correlated with Braak Stage (ComBat-corrected)",
       x = "Gene", y = "Spearman rho") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )

outdir <- "C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(outdir, "Braak_top_corr_barplot_combat.png"),
  plot = p,
  width = 8, height = 8, units = "in",
  dpi = 600, bg = "white"
)


