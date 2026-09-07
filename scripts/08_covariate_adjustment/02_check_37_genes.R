# Reviewer 2, Major Comment 2 - part 2.
#
# Input:  the DE tables written by 01_covariate_adjusted_DE.R
#         (Results_covariate_adjusted/de_tables/*_{unadjusted,adjusted}.csv)
# Output: the 37-gene comparison tables and a plain-text findings summary.
#
# This script re-reads files only. It never refits a model, so it can be rerun
# cheaply after editing thresholds or the gene list.

suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db); library(here)
})

outdir <- here::here("Results_covariate_adjusted")
de_dir <- file.path(outdir, "de_tables")
stopifnot(dir.exists(de_dir))

P_CUT <- 0.05
# The primary pipeline calls a DEG on |logFC| > 0.58 AND unadjusted P < 0.05.
# Both criteria are reported: Sig_* is the P-only flag the reviewer asked for,
# DEG_* additionally applies the fold-change cut so the counts are comparable
# to the published per-dataset DEG tables.
LFC_CUT <- 0.58

# ------------------------------------------------------- the 37 signature ----
# Typos in earlier source files: CCD102A -> CCDC102A, HSBP7 -> HSPB7.
genes37 <- c(
  "ADAM33", "AEBP1", "CLDN9", "CCDC102A", "HSPB1", "HSPB7", "GFAP", "KANK2",
  "KLF15", "MRGPRF", "NUPR1", "PIK3R5", "PRELP", "PRX", "TCEA3", "TMPRSS5",
  "CHML", "COL25A1", "ELOVL4", "GAD1", "GAD2", "HPRT1", "ITFG1", "MAS1",
  "NAP1L5", "NCALD", "NEUROD6", "NRN1", "OPN3", "RAB3B", "RAB3C", "RGS4",
  "RPH3A", "SCG2", "SERPINI1", "STAT4", "TRIM36")
stopifnot(length(genes37) == 37L, !any(duplicated(genes37)))

# ------------------------------------------------------------- unit index ----
units <- data.frame(
  Dataset = c("GSE67333", "GSE95587", "GSE159699", "GSE203206",
              rep("GSE278723", 5)),
  Region  = c(NA, NA, NA, NA, "CA1", "CA2", "CA3", "CA4", "NX"),
  stringsAsFactors = FALSE)
units$Tag <- ifelse(is.na(units$Region), units$Dataset,
                    paste0(units$Dataset, "_", units$Region))

model_summary <- read.csv(file.path(outdir, "model_summary.csv"),
                          stringsAsFactors = FALSE)

# ---------------------------------------------------------------- helpers ----

# Same harmonisation as scripts/03_meta_analysis/stages/01_harmonize_inputs.R:
# strip Ensembl version suffix, map to SYMBOL, uppercase, one row per symbol
# (lowest P, ties broken by largest |logFC|).
looks_like_ensembl <- function(x)
  mean(grepl("^ENSG\\d{11}$", sub("\\.\\d+$", "", x)), na.rm = TRUE) > 0.5

harmonize <- function(df) {
  g <- trimws(as.character(df$Gene))
  if (looks_like_ensembl(g)) {
    clean <- sub("\\.\\d+$", "", g)
    conv <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = unique(clean), keytype = "ENSEMBL",
      columns = c("ENSEMBL", "SYMBOL")))
    conv <- conv[!duplicated(conv$ENSEMBL), ]
    df$SYMBOL <- conv$SYMBOL[match(clean, conv$ENSEMBL)]
  } else {
    df$SYMBOL <- g
  }
  df$SYMBOL <- toupper(df$SYMBOL)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[order(df$P.Value, -abs(df$logFC)), ]
  df[!duplicated(df$SYMBOL), ]
}

read_fit <- function(tag, which) {
  f <- file.path(de_dir, paste0(tag, "_", which, ".csv"))
  if (!file.exists(f)) return(NULL)
  harmonize(read.csv(f, stringsAsFactors = FALSE))
}

# Inverse-variance fixed-effect pooling; identical to collapse_study_gene() in
# scripts/03_meta_analysis/stages/02_collapse_regions_and_qc.R.
pool_fe <- function(logFC, SE) {
  ok <- is.finite(logFC) & is.finite(SE) & SE > 0
  if (!any(ok)) return(c(logFC = NA, SE = NA, P.Value = NA))
  w <- 1 / SE[ok]^2
  b <- sum(w * logFC[ok]) / sum(w)
  s <- sqrt(1 / sum(w))
  c(logFC = b, SE = s, P.Value = 2 * pnorm(abs(b / s), lower.tail = FALSE))
}

# --------------------------------------------------- per-gene comparison -----
rows <- list()
for (i in seq_len(nrow(units))) {
  tag <- units$Tag[i]
  un <- read_fit(tag, "unadjusted")
  if (is.null(un)) { warning("missing unadjusted table for ", tag); next }
  ad <- read_fit(tag, "adjusted")   # NULL when no adjusted model was fitted
  ms <- model_summary[model_summary$Dataset == units$Dataset[i] &
                        (is.na(units$Region[i]) | model_summary$Region %in% units$Region[i]), ][1, ]

  iu <- match(genes37, un$SYMBOL)
  ia <- if (is.null(ad)) rep(NA_integer_, 37) else match(genes37, ad$SYMBOL)

  lfc_u <- un$logFC[iu];  p_u <- un$P.Value[iu];  q_u <- un$adj.P.Val[iu]
  lfc_a <- if (is.null(ad)) NA_real_ else ad$logFC[ia]
  p_a   <- if (is.null(ad)) NA_real_ else ad$P.Value[ia]
  se_a  <- if (is.null(ad)) NA_real_ else ad$SE[ia]
  q_a   <- if (is.null(ad)) NA_real_ else ad$adj.P.Val[ia]

  rows[[tag]] <- data.frame(
    Gene = genes37,
    Dataset = units$Dataset[i],
    Region = units$Region[i],
    Level = if (is.na(units$Region[i])) "study" else "subregion",
    Covariates_in_model = ms$Covariates_in_model,
    logFC_unadjusted = lfc_u,
    P_unadjusted = p_u,
    FDR_unadjusted = q_u,
    SE_unadjusted = un$SE[iu],
    logFC_adjusted = lfc_a,
    P_adjusted = p_a,
    FDR_adjusted = q_a,
    SE_adjusted = se_a,
    Direction_concordant = sign(lfc_u) == sign(lfc_a),
    Sig_unadjusted = p_u < P_CUT,
    Sig_adjusted = p_a < P_CUT,
    DEG_unadjusted = p_u < P_CUT & abs(lfc_u) > LFC_CUT,
    DEG_adjusted = p_a < P_CUT & abs(lfc_a) > LFC_CUT,
    FDRsig_unadjusted = q_u < P_CUT,
    FDRsig_adjusted = q_a < P_CUT,
    FDR_DEG_unadjusted = q_u < P_CUT & abs(lfc_u) > LFC_CUT,
    FDR_DEG_adjusted = q_a < P_CUT & abs(lfc_a) > LFC_CUT,
    logFC_change = lfc_a - lfc_u,
    residual_df = ms$residual_df_adjusted,
    stringsAsFactors = FALSE)
}
per_gene <- do.call(rbind, rows)

# GSE278723: collapse the five subregions to one study-level row per gene.
# Pooling is done genome-wide over every gene tested in all five subregions so
# that the pooled FDR is a genuine transcriptome-wide correction, not a BH
# adjustment over just the 37 signature genes. The 37 are subset afterwards.
pool_all <- function(which) {
  fits <- lapply(c("CA1", "CA2", "CA3", "CA4", "NX"),
                 function(r) read_fit(paste0("GSE278723_", r), which))
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  # Union, not intersection: a gene expressed in four of five subregions is
  # pooled over those four rather than dropped. Weights are 1/SE^2, so absent
  # subregions simply contribute nothing.
  all_sym <- sort(unique(unlist(lapply(fits, function(f) f$SYMBOL))))
  lfc <- vapply(fits, function(f) f$logFC[match(all_sym, f$SYMBOL)], numeric(length(all_sym)))
  se  <- vapply(fits, function(f) f$SE[match(all_sym, f$SYMBOL)],    numeric(length(all_sym)))
  w   <- ifelse(is.finite(se) & se > 0, 1 / se^2, 0)
  lfc[!is.finite(lfc)] <- 0
  sw  <- rowSums(w)
  b   <- ifelse(sw > 0, rowSums(w * lfc) / sw, NA_real_)
  sp  <- ifelse(sw > 0, sqrt(1 / sw), NA_real_)
  pv  <- 2 * pnorm(abs(b / sp), lower.tail = FALSE)
  data.frame(SYMBOL = all_sym, logFC = b, SE = sp, P.Value = pv,
             adj.P.Val = p.adjust(pv, method = "BH"),
             n_subregions = rowSums(w > 0), stringsAsFactors = FALSE)
}
pu <- pool_all("unadjusted")
pa <- pool_all("adjusted")
sub <- per_gene[per_gene$Dataset == "GSE278723", ]
iu <- match(genes37, pu$SYMBOL); ia <- match(genes37, pa$SYMBOL)
coll <- data.frame(
  Gene = genes37, Dataset = "GSE278723", Region = NA_character_,
  Level = "study (inverse-variance FE pool of 5 subregions)",
  Covariates_in_model = sub$Covariates_in_model[1],
  logFC_unadjusted = pu$logFC[iu],
  P_unadjusted = pu$P.Value[iu],
  FDR_unadjusted = pu$adj.P.Val[iu],
  SE_unadjusted = pu$SE[iu],
  logFC_adjusted = pa$logFC[ia],
  P_adjusted = pa$P.Value[ia],
  FDR_adjusted = pa$adj.P.Val[ia],
  SE_adjusted = pa$SE[ia],
  Direction_concordant = sign(pu$logFC[iu]) == sign(pa$logFC[ia]),
  Sig_unadjusted = pu$P.Value[iu] < P_CUT,
  Sig_adjusted = pa$P.Value[ia] < P_CUT,
  DEG_unadjusted = pu$P.Value[iu] < P_CUT & abs(pu$logFC[iu]) > LFC_CUT,
  DEG_adjusted = pa$P.Value[ia] < P_CUT & abs(pa$logFC[ia]) > LFC_CUT,
  FDRsig_unadjusted = pu$adj.P.Val[iu] < P_CUT,
  FDRsig_adjusted = pa$adj.P.Val[ia] < P_CUT,
  FDR_DEG_unadjusted = pu$adj.P.Val[iu] < P_CUT & abs(pu$logFC[iu]) > LFC_CUT,
  FDR_DEG_adjusted = pa$adj.P.Val[ia] < P_CUT & abs(pa$logFC[ia]) > LFC_CUT,
  logFC_change = pa$logFC[ia] - pu$logFC[iu],
  residual_df = NA_integer_,
  stringsAsFactors = FALSE)
per_gene <- rbind(per_gene, coll)
per_gene <- per_gene[order(per_gene$Gene, per_gene$Dataset, per_gene$Level), ]
write.csv(per_gene, file.path(outdir, "table1_per_gene_comparison.csv"),
          row.names = FALSE)

# ------------------------------------------------------ per-dataset summary --
avail <- read.csv(file.path(outdir, "covariate_availability_long.csv"),
                  stringsAsFactors = FALSE)

summ <- do.call(rbind, lapply(unique(units$Dataset), function(ds) {
  ms <- model_summary[model_summary$Dataset == ds, ]
  d  <- per_gene[per_gene$Dataset == ds & per_gene$Level != "subregion", ]
  a  <- avail[avail$Dataset == ds, ]
  excl <- a[a$Status != "Included", ]
  # No adjusted fit for this dataset -> the adjusted columns are undefined, not 0.
  fitted <- any(!is.na(d$logFC_adjusted))
  na_if <- function(x) if (fitted) x else NA_integer_
  data.frame(
    Dataset = ds,
    n_samples = paste(unique(ms$n_samples), collapse = "/"),
    n_AD = paste(unique(ms$n_AD), collapse = "/"),
    n_control = paste(unique(ms$n_control), collapse = "/"),
    Covariates_included = paste(a$Covariate[a$Status == "Included"], collapse = ", "),
    Covariates_excluded = paste(sprintf("%s (%s)", excl$Covariate, excl$Reason),
                                collapse = "; "),
    residual_df_unadjusted = paste(unique(ms$residual_df_unadjusted), collapse = "/"),
    residual_df_adjusted = paste(unique(ms$residual_df_adjusted), collapse = "/"),
    n_of_37_present = sum(!is.na(d$logFC_unadjusted)),
    n_directionally_concordant = na_if(sum(d$Direction_concordant, na.rm = TRUE)),
    n_sig_unadjusted = sum(d$Sig_unadjusted, na.rm = TRUE),
    n_sig_adjusted = na_if(sum(d$Sig_adjusted, na.rm = TRUE)),
    n_retaining_significance = na_if(sum(d$Sig_unadjusted & d$Sig_adjusted, na.rm = TRUE)),
    n_DEG_unadjusted = sum(d$DEG_unadjusted, na.rm = TRUE),
    n_DEG_adjusted = na_if(sum(d$DEG_adjusted, na.rm = TRUE)),
    n_retaining_DEG = na_if(sum(d$DEG_unadjusted & d$DEG_adjusted, na.rm = TRUE)),
    n_FDRsig_unadjusted = sum(d$FDRsig_unadjusted, na.rm = TRUE),
    n_FDRsig_adjusted = na_if(sum(d$FDRsig_adjusted, na.rm = TRUE)),
    n_retaining_FDRsig = na_if(sum(d$FDRsig_unadjusted & d$FDRsig_adjusted, na.rm = TRUE)),
    n_FDR_DEG_unadjusted = sum(d$FDR_DEG_unadjusted, na.rm = TRUE),
    n_FDR_DEG_adjusted = na_if(sum(d$FDR_DEG_adjusted, na.rm = TRUE)),
    n_retaining_FDR_DEG = na_if(sum(d$FDR_DEG_unadjusted & d$FDR_DEG_adjusted, na.rm = TRUE)),
    median_abs_logFC_unadjusted = round(median(abs(d$logFC_unadjusted), na.rm = TRUE), 4),
    median_abs_logFC_adjusted = round(median(abs(d$logFC_adjusted), na.rm = TRUE), 4),
    median_abs_logFC_change = round(median(abs(d$logFC_change), na.rm = TRUE), 4),
    stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(outdir, "table2_per_dataset_summary.csv"),
          row.names = FALSE)

# ------------------------------------------------------------ findings -------
# Study-level rows only - subregion rows would count GSE278723 five times over.
adj <- per_gene[!is.na(per_gene$logFC_adjusted) & per_gene$Level != "subregion", ]
flipped <- adj[!adj$Direction_concordant, c("Gene", "Dataset",
                                            "logFC_unadjusted", "logFC_adjusted",
                                            "P_unadjusted", "P_adjusted",
                                            "Sig_unadjusted", "Sig_adjusted")]
# A sign flip only matters if at least one of the two estimates was significant.
flipped$Material <- flipped$Sig_unadjusted | flipped$Sig_adjusted
lost <- adj[adj$Sig_unadjusted & !adj$Sig_adjusted, ]
lost_ct <- sort(table(lost$Gene), decreasing = TRUE)
lost_multi <- names(lost_ct)[lost_ct > 1]

out <- c(
  "Covariate-adjusted sensitivity analysis - findings summary",
  format(Sys.Date()), strrep("=", 70), "",
  sprintf("37-gene signature checked in %d datasets. Significance = P < %.2f.",
          length(unique(units$Dataset)), P_CUT),
  "Fit A (~ diagnosis) reproduced the deposited per-dataset DE tables exactly",
  "(see fitA_reproducibility_check.csv).", "",
  "Per dataset:", "")
for (i in seq_len(nrow(summ))) {
  s <- summ[i, ]
  out <- c(out,
    sprintf("* %s (n=%s; %s AD / %s control)", s$Dataset, s$n_samples, s$n_AD, s$n_control),
    sprintf("    included:  %s", ifelse(nzchar(s$Covariates_included), s$Covariates_included, "none")),
    sprintf("    excluded:  %s", s$Covariates_excluded),
    sprintf("    residual df: %s unadjusted / %s adjusted",
            s$residual_df_unadjusted, s$residual_df_adjusted),
    sprintf("    37 genes present %d; directionally concordant %d",
            s$n_of_37_present, s$n_directionally_concordant),
    sprintf("    P < 0.05 only:            %d -> %d (retained %d)",
            s$n_sig_unadjusted, s$n_sig_adjusted, s$n_retaining_significance),
    sprintf("    P < 0.05 AND |logFC|>0.58: %d -> %d (retained %d)",
            s$n_DEG_unadjusted, s$n_DEG_adjusted, s$n_retaining_DEG),
    sprintf("    FDR < 0.05 only:          %d -> %d (retained %d)",
            s$n_FDRsig_unadjusted, s$n_FDRsig_adjusted, s$n_retaining_FDRsig),
    sprintf("    FDR < 0.05 AND |logFC|>0.58: %d -> %d (retained %d)",
            s$n_FDR_DEG_unadjusted, s$n_FDR_DEG_adjusted, s$n_retaining_FDR_DEG),
    sprintf("    median |logFC| %s -> %s; median |dlogFC| %s",
            format(s$median_abs_logFC_unadjusted),
            format(s$median_abs_logFC_adjusted),
            format(s$median_abs_logFC_change)),
    "")
}
out <- c(out, strrep("-", 70), "",
  if (nrow(flipped) == 0) "No gene changed sign in any dataset after adjustment." else
    c("Genes changing sign after adjustment (name these explicitly).",
      "Material = TRUE means the gene was significant in at least one of the two fits;",
      "flips among non-significant estimates are noise-level sign changes.",
      capture.output(print(flipped, row.names = FALSE))),
  "",
  if (!length(lost_multi)) "No gene lost significance in more than one dataset." else
    c("Genes losing significance in more than one dataset (name these explicitly):",
      paste0("  ", lost_multi, " (", lost_ct[lost_multi], " datasets)")),
  "",
  "GSE278723: the GEO record deposits only tissue, subfield and diagnostic",
  "group. Sex, age and PMI were taken from Supplementary Table 1 of the source",
  "paper (Li C et al. 2024, Aging Dis) and joined by donor ID; RIN, sequencing",
  "batch and APOE are recorded nowhere. Both fits are run per subregion and",
  "pooled to one study-level estimate by inverse-variance fixed-effect",
  "weighting, matching the primary pipeline.",
  "",
  "GSE67333: PMI is likewise absent from GEO and was recovered from the source",
  "paper's supplement (Magistri et al. 2015, J Alzheimers Dis 48:647) by",
  "matching donors on age, gender, APOE and Braak stage. With n = 8 the adjusted",
  "model leaves 3 residual df, so its estimates are power-limited rather than",
  "informative about confounding.",
  "",
  "Donor grouping for GSE278723 follows Supplementary Table 1 (NC = 01-13,",
  "n = 13; AD = 14-28, n = 15), matching the corrected primary pipeline.")
writeLines(out, file.path(outdir, "findings_summary.txt"))
cat(paste(out, collapse = "\n"), "\n")

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_02.txt"))
