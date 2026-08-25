# Give all subregions from the same dataset a shared study_id
# e.g., GSE278723_CA1 -> GSE278723; others stay as-is
all_deg <- all_deg %>%
  dplyr::mutate(study_id = sub("_.*$", "", acc))

## ---- GSE278723: between-subregion heterogeneity + design-effect SE ----
gse278723_regions <- all_deg %>%
  dplyr::filter(study_id == "GSE278723",
                is.finite(logFC), is.finite(SE), SE > 0)

het_one_gene <- function(d) {
  m <- nrow(d)
  if (m < 2) return(tibble::tibble())
  wi <- 1 / (d$SE^2)
  yi <- d$logFC
  beta_pooled <- sum(wi * yi) / sum(wi)
  se_pooled   <- sqrt(1 / sum(wi))
  Q  <- sum(wi * (yi - beta_pooled)^2)
  df <- m - 1
  tibble::tibble(
    SYMBOL = d$SYMBOL[1],
    m_regions = m,
    beta_pooled = beta_pooled,
    se_pooled = se_pooled,
    Q = Q, Q_df = df,
    Q_p = pchisq(Q, df, lower.tail = FALSE),
    I2 = max(0, (Q - df) / Q) * 100,
    se_rho0.0 = se_pooled,
    se_rho0.3 = se_pooled * sqrt(1 + (m - 1) * 0.3),
    se_rho0.5 = se_pooled * sqrt(1 + (m - 1) * 0.5),
    se_rho0.7 = se_pooled * sqrt(1 + (m - 1) * 0.7)
  )
}

gse278723_het <- gse278723_regions %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::group_modify(~ het_one_gene(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Q_FDR = p.adjust(Q_p, "BH"))

message("GSE278723 genes with all 5 regions: ",
        sum(gse278723_het$m_regions == 5))
print(summary(gse278723_het$I2))
message("Fraction homogeneous (Q_FDR > 0.05): ",
        round(mean(gse278723_het$Q_FDR > 0.05, na.rm = TRUE), 3))
readr::write_csv(gse278723_het, file.path(output_dir, "GSE278723_subregion_heterogeneity.csv"))


# Collapse multiple regions within a dataset to ONE study-level row per gene
collapse_study_gene <- function(df_one_gene_one_study){
  # Expect columns: SYMBOL, logFC, SE, P.Value
  df <- df_one_gene_one_study %>% dplyr::filter(is.finite(logFC))
  # Use fixed-effect pooling if any SE are available
  if (any(is.finite(df$SE))) {
    df2 <- df %>% dplyr::filter(is.finite(SE))
    wi  <- 1/(df2$SE^2)
    yi  <- df2$logFC
    meta_logFC <- sum(wi*yi) / sum(wi)
    meta_var   <- 1 / sum(wi)
    meta_se    <- sqrt(meta_var)
    z          <- meta_logFC / meta_se
    p          <- 2*pnorm(abs(z), lower.tail = FALSE)
    tibble::tibble(SYMBOL = df2$SYMBOL[1],
                   study_id = df2$study_id[1],
                   logFC = meta_logFC, SE = meta_se, P.Value = p)
  } else {
    # Fallback if SE missing: direction-aware Stouffer within this dataset
    d <- df %>% dplyr::filter(is.finite(P.Value))
    if (nrow(d) == 0) return(tibble::tibble())
    sgn <- ifelse(d$logFC > 0, 1L, ifelse(d$logFC < 0, -1L, 0L))
    maj <- ifelse(mean(sgn==1, na.rm=TRUE) >= mean(sgn==-1, na.rm=TRUE), 1L, -1L)
    keep <- (sgn == maj) & (sgn != 0L); if (!any(keep)) keep <- sgn != 0L
    pvec <- pmin(pmax(d$P.Value[keep], .Machine$double.eps), 1)
    zmag <- qnorm(pvec/2, lower.tail = FALSE)
    Z    <- sum(zmag) / sqrt(length(zmag))
    p    <- 2*pnorm(abs(Z), lower.tail = FALSE)
    # Use robust central tendency for effect; no SE available
    lfc  <- stats::median(d$logFC[keep], na.rm = TRUE) * maj
    tibble::tibble(SYMBOL = d$SYMBOL[1],
                   study_id = d$study_id[1],
                   logFC = lfc, SE = NA_real_, P.Value = p)
  }
}

# Build the collapsed, study-level table
by_study_deg <- all_deg %>%
  dplyr::group_by(study_id, SYMBOL) %>%
  dplyr::group_modify(~ collapse_study_gene(.x)) %>%
  dplyr::ungroup()

gse278723_collapsed_genes <- by_study_deg %>%
  filter(study_id == "GSE278723")
# Save the full results table (gene, logFC, SE, P.Value) to a CSV file
write_csv(gse278723_collapsed_genes, file.path(output_dir, "GSE278723_collapsed_results.csv"))
write_csv(all_deg, file.path(output_dir, "ALL_DATASET_COMBINED.csv"))

library(ggplot2)

## Overall logFC distribution (all regions, all studies)
ggplot(all_deg, aes(x = logFC)) +
  geom_histogram(bins = 60) +
  labs(title = "Distribution of logFC across all datasets",
       x = "log2 Fold Change", y = "Count")

## Per-study density of logFC (after collapsing regions to study_id)
ggplot(by_study_deg, aes(x = logFC)) +
  geom_density() +
  facet_wrap(~ study_id, scales = "free_y") +
  labs(title = "logFC density by study",
       x = "log2 Fold Change", y = "Density")
ggplot(all_deg, aes(x = -log10(P.Value))) +
  geom_histogram(bins = 60) +
  labs(title = "P-value distribution across all datasets",
       x = expression(-log[10](p)), y = "Count")

ggplot(by_study_deg, aes(x = study_id, y = logFC)) +
  geom_boxplot(outlier.size = 0.3) +
  coord_flip() +
  labs(title = "logFC distribution by study",
       x = "Study", y = "log2 Fold Change")
ggplot(all_deg, aes(x = acc, y = logFC)) +
  geom_boxplot(outlier.size = 0.3) +
  coord_flip() +
  labs(title = "logFC distribution by region/study (pre-collapse)",
       x = "Accession / Region", y = "log2 Fold Change")


# (Optional) sanity check
message("Rows before collapsing: ", nrow(all_deg),
        " | after: ", nrow(by_study_deg),
        " | studies: ", dplyr::n_distinct(by_study_deg$study_id))
# --- Sanity Check: Which studies are missing Standard Error (SE)? ---
se_availability_summary <- all_deg |>
  # Create a column that is TRUE if SE is present, FALSE if it's missing (NA)
  mutate(has_SE = !is.na(SE)) |>
  # Group by the original study/region accession
  group_by(acc, region) |>
  # For each study, calculate the percentage of its genes that have an SE
  summarise(
    total_genes = n(),
    percent_with_SE = mean(has_SE) * 100
  ) |>
  # Arrange the results to see the ones with the least data first
  arrange(percent_with_SE)

# Print the summary table to your console
print(se_availability_summary)

# Save this summary to a file for your records
readr::write_csv(se_availability_summary, file.path(output_dir, "se_availability_by_study.csv"))
