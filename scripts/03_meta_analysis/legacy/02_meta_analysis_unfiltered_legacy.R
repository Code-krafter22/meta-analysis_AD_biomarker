suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr); library(tidyr)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(metafor)
  library(MetaVolcanoR)
  library(ggplot2)
  library(here)
})

project_root <- here::here()
deg_root <- file.path(project_root, "results", "differential_expression")
output_dir <- file.path(project_root, "results", "meta_analysis", "unfiltered_all_genes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Tools you may need
#install.packages(c("remotes", "metafor"))
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 2) Bioconductor core deps often needed by MetaVolcanoR
#BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))

# 3) Install MetaVolcanoR from its source repo

# 4) Load
library(MetaVolcanoR)


files <- tibble::tribble(
  ~acc,              ~region,                 ~path,
  # GSE278723 subregions (hippocampal???entorhinal system)
  "GSE278723_CA1",   "hippocampus_CA1",       file.path(deg_root, "GSE278723", "CA1", "full_DE_results.csv"),
  "GSE278723_CA2",   "hippocampus_CA2",       file.path(deg_root, "GSE278723", "CA2", "full_DE_results.csv"),
  "GSE278723_CA3",   "hippocampus_CA3",       file.path(deg_root, "GSE278723", "CA3", "full_DE_results.csv"),
  "GSE278723_CA4",   "hippocampus_CA4",       file.path(deg_root, "GSE278723", "CA4", "full_DE_results.csv"),
  "GSE278723_NX",    "entorhinal_cortex",     file.path(deg_root, "GSE278723", "NX", "full_DE_results.csv"),
  
  # Other GEOs
  "GSE67333",        "hippocampus",           file.path(deg_root, "GSE67333", "full_DE_results.csv"),
  "GSE95587",        "fusiform_gyrus",        file.path(deg_root, "GSE95587", "full_DE_results.csv"),
  "GSE159699",       "lateral_temporal_lobe", file.path(deg_root, "GSE159699", "full_DE_results.csv"),
  "GSE203206",       "occipital_lobe",        file.path(deg_root, "GSE203206", "full_DE_results.csv")
)


resolve_base <- function(rel_paths) {
  cands <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  # pick the first directory where ALL files exist; else first where ANY exist; else wd
  for (d in cands) if (all(file.exists(file.path(d, rel_paths)))) return(d)
  for (d in cands) if (any(file.exists(file.path(d, rel_paths)))) return(d)
  getwd()
}
missing <- files$path[!file.exists(files$path)]
if (length(missing)) {
  stop(
    "These paths were not found:\n - ", paste(missing, collapse = "\n - "),
    "\nProject root: ", project_root
  )
}

## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
guess_gene_col <- function(df){
  nms <- names(df)
  cand <- nms[grepl("^(gene|genes|symbol|hgnc|gene_name)$", tolower(nms))]
  if (length(cand) >= 1) return(cand[1])
  cand2 <- nms[grepl("ensembl|gene_id|ensembl_gene", tolower(nms))]
  if (length(cand2) >= 1) return(cand2[1])
  # avoid unnamed index-like first column if it???s just 1:n
  if (ncol(df) && nms[1] %in% c("...1","X1","X","Unnamed: 0")) {
    v <- df[[1]]
    if (is.numeric(v) && all(v == seq_len(nrow(df)))) {
      if (ncol(df) >= 2) return(nms[2])
    }
  }
  nms[1]
}

looks_like_ensembl <- function(x){
  x <- as.character(x)
  x0 <- stringr::str_replace(x, "\\.\\d+$", "")
  mean(grepl("^ENSG\\d{11}$", x0), na.rm = TRUE) > 0.5
}

ensembl_to_symbol <- function(ensembl_vec){
  ens_clean <- stringr::str_replace(ensembl_vec, "\\.\\d+$", "")
  conv <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ens_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL","SYMBOL")
  )
  tibble::tibble(ENSEMBL_clean = ens_clean) |>
    dplyr::left_join(conv, by = c("ENSEMBL_clean" = "ENSEMBL"))
}

standardize_stat_cols <- function(df){
  pick_first <- function(cols){
    nm <- intersect(cols, names(df))[1]
    if (length(nm) == 0 || is.na(nm)) return(rep(NA_real_, nrow(df)))
    df[[nm]]
  }
  logFC   <- pick_first(c("logFC","log2FoldChange","LFC","log2fc"))
  P.Value <- pick_first(c("P.Value","pvalue","p_val","P","adj.P.Val","padj"))
  t       <- suppressWarnings(pick_first(c("t","stat","T")))
  lfcSE   <- suppressWarnings(pick_first(c("lfcSE","SE","se")))
  
  SE <- dplyr::case_when(
    !is.na(lfcSE)             ~ as.numeric(lfcSE),
    !is.na(logFC) & !is.na(t) ~ abs(as.numeric(logFC) / as.numeric(t)),
    TRUE                      ~ NA_real_
  )
  
  dplyr::mutate(df,
                logFC = logFC, P.Value = P.Value, t = t, lfcSE = lfcSE, SE = SE
  )
}

read_one_deg <- function(path, acc, region){
  ext <- tools::file_ext(path)
  df_raw <- if (tolower(ext) %in% c("tsv","txt")) {
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE, guess_max = 1e5))
  } else {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE, guess_max = 1e5))
  }
  
  # Drop a pure row-index column if present (e.g., ...1 / X1 / Unnamed: 0)
  if (ncol(df_raw) > 0 && names(df_raw)[1] %in% c("...1","X1","X","Unnamed: 0")) {
    v <- df_raw[[1]]
    if (is.numeric(v) && all(v == seq_len(nrow(df_raw)))) {
      df_raw <- dplyr::select(df_raw, -1)
    }
  }
  
  gene_col <- guess_gene_col(df_raw)
  df_raw[[gene_col]] <- trimws(as.character(df_raw[[gene_col]]))
  
  if (looks_like_ensembl(df_raw[[gene_col]])) {
    mapdf <- ensembl_to_symbol(df_raw[[gene_col]])
    df <- df_raw |>
      dplyr::mutate(ENSEMBL_clean = stringr::str_replace(.data[[gene_col]], "\\.\\d+$", "")) |>
      dplyr::left_join(mapdf, by = "ENSEMBL_clean")
  } else {
    df <- df_raw |> dplyr::mutate(SYMBOL = !!rlang::sym(gene_col))
  }
  
  df <- df |>
    dplyr::mutate(SYMBOL = toupper(SYMBOL)) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    standardize_stat_cols()
  
  # One row per SYMBOL per study (lowest p; tie ??? largest |logFC|)
  df |>
    dplyr::mutate(
      .p_aux = suppressWarnings(as.numeric(P.Value)),
      .p_aux = ifelse(is.na(.p_aux), 1, .p_aux),
      .absfc = suppressWarnings(abs(as.numeric(logFC)))
    ) |>
    dplyr::group_by(SYMBOL) |>
    dplyr::arrange(.p_aux, dplyr::desc(.absfc), .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(SYMBOL, logFC, SE, P.Value) |>
    dplyr::mutate(
      acc = acc,
      region = region,
      logFC   = suppressWarnings(as.numeric(logFC)),
      SE      = suppressWarnings(as.numeric(SE)),
      P.Value = suppressWarnings(as.numeric(P.Value))
    ) |>
    dplyr::filter(!is.na(logFC), !is.na(P.Value))
}

## ------------------------------------------------------------------
## Run
## ------------------------------------------------------------------
dir.create(output_dir, showWarnings = FALSE)

reader_safe <- purrr::possibly(read_one_deg, otherwise = tibble::tibble(), quiet = FALSE)

deg_std_list <- purrr::pmap(files, function(acc, region, path){
  message("Processing ", acc, " (", region, ") ...")
  reader_safe(path, acc, region)
}) |> purrr::set_names(files$acc)

all_deg <- dplyr::bind_rows(deg_std_list)

readr::write_csv(all_deg, file.path(output_dir, "all_studies_harmonized.csv"))
message(
  "Harmonized rows: ", nrow(all_deg),
  " | Genes: ", dplyr::n_distinct(all_deg$SYMBOL),
  " | Studies: ", dplyr::n_distinct(all_deg$acc)
)

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


# =======================
# Direction-aware Stouffer meta (weighted)
# =======================

dir.create(output_dir, showWarnings = FALSE)

# (1) Optional: set study weights (default = 1)
#     If you know per-study sample sizes, put them here and weights will be sqrt(n_total).
#     Otherwise everything is weighted equally.
study_weights <- tibble::tibble(
  study_id = unique(by_study_deg$study_id),
  n_total  = NA_real_   # fill with AD+Ctrl per dataset if you have it
) %>%
  dplyr::mutate(w = ifelse(is.finite(n_total), sqrt(n_total), 1))


# (2) Helper to convert two-sided p to signed Z
to_signed_z <- function(p, sign) {
  # clamp p into (0,1] to avoid Inf
  p <- pmin(pmax(p, .Machine$double.eps), 1)
  zmag <- qnorm(p/2, lower.tail = FALSE)   # >= 0
  sign * zmag
}

stouffer_one_gene <- function(d) {
  d <- d %>%
    dplyr::filter(is.finite(P.Value), P.Value > 0, P.Value <= 1, is.finite(logFC)) %>%
    dplyr::mutate(sign = dplyr::case_when(logFC > 0 ~ 1L, logFC < 0 ~ -1L, TRUE ~ 0L))
  
  if (nrow(d) == 0) return(tibble::tibble(K=0L, K_used=0L, Z=NA_real_, P=NA_real_, Dir=NA_character_))
  
  n_up <- sum(d$sign ==  1L); n_down <- sum(d$sign == -1L)
  majority <- ifelse(n_up >= n_down, 1L, -1L)
  keep <- d$sign == majority
  if (!any(keep, na.rm = TRUE)) keep <- d$sign != 0L
  d2 <- d[keep & d$sign != 0L, , drop = FALSE]
  
  if (nrow(d2) == 0) return(tibble::tibble(K=nrow(d), K_used=0L, Z=NA_real_, P=NA_real_, Dir=ifelse(majority>0,"up","down")))
  
  d2 <- d2 %>%
    dplyr::left_join(study_weights %>% dplyr::select(study_id, w), by = "study_id") %>%
    dplyr::mutate(w = dplyr::if_else(is.finite(w), w, 1))
  
  zi <- to_signed_z(d2$P.Value, d2$sign)
  wi <- d2$w
  Z  <- sum(wi * zi) / sqrt(sum(wi^2))
  P  <- 2 * pnorm(abs(Z), lower.tail = FALSE)
  
  tibble::tibble(K=nrow(d), K_used=nrow(d2), Z=Z, P=P, Dir=ifelse(Z>=0,"up","down"))
}


# (3) Run across all genes
stouffer_res <- by_study_deg %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::group_modify(~ stouffer_one_gene(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(FDR = p.adjust(P, "BH")) %>%
  dplyr::arrange(FDR, P)

# (4) Save + convenient subsets
readr::write_csv(stouffer_res, file.path(output_dir, "stouffer_meta_results.csv"))

stouffer_sig  <- stouffer_res %>% dplyr::filter(FDR <= 0.05)
stouffer_topN <- stouffer_res %>% dplyr::slice_head(n = 200)

 readr::write_csv(stouffer_sig,  file.path(output_dir, "stouffer_meta_significant_fdr05.csv"))
 readr::write_csv(stouffer_topN, file.path(output_dir, "stouffer_meta_top200.csv"))


#4a) Keep genes seen in ???3 studies (tuneable)
# keep genes seen in ???3 datasets (now study-level, not region-level)
keep_genes <- by_study_deg %>%
  dplyr::distinct(SYMBOL, study_id) %>%
  dplyr::count(SYMBOL, name = "k") %>%
  dplyr::filter(k >= 3)

meta_in <- by_study_deg %>% dplyr::semi_join(keep_genes, by = "SYMBOL")
meta_ready <- meta_in %>% dplyr::filter(!is.na(SE)) %>% dplyr::mutate(vi = SE^2)

# ... then your meta_one_gene() and meta_eff code exactly as before,
# except use 'study_id' instead of 'acc' for labels if you like:
# meta_ready should have: SYMBOL, study_id, logFC, vi ( = SE^2 )
stopifnot(all(c("SYMBOL","study_id","logFC","vi") %in% names(meta_ready)))

meta_one_gene <- function(d){
  # require at least 2 studies after filtering
  d <- d %>% dplyr::filter(is.finite(logFC), is.finite(vi))
  if (nrow(d) < 2) return(tibble::tibble())
  # guard against zero/negative variances
  d$vi <- pmax(d$vi, .Machine$double.eps)
  
  fit <- tryCatch(
    metafor::rma.uni(yi = d$logFC, vi = d$vi, method = "DL", slab = d$study_id),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble::tibble())
  
  tibble::tibble(
    SYMBOL     = d$SYMBOL[1],
    k          = nrow(d),
    meta_logFC = as.numeric(fit$b),
    meta_se    = fit$se,
    meta_z     = fit$zval,
    meta_p     = fit$pval,
    ci_lb      = fit$ci.lb,
    ci_ub      = fit$ci.ub,
    tau2       = fit$tau2,
    I2         = fit$I2
  )
}

meta_eff <- meta_ready %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::group_modify(~ meta_one_gene(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(meta_FDR = p.adjust(meta_p, "BH")) %>%
  dplyr::arrange(meta_FDR, meta_p)

 readr::write_csv(meta_eff, file.path(output_dir, "meta_effects_meta_results.csv"))
library(tidyr)
library(dplyr)
library(pheatmap)
library(tibble)

## choose top 50 significant genes ??? adjust N as you like
top_genes <- meta_eff %>%
  arrange(meta_FDR, meta_p) %>%
  filter(meta_FDR <= 0.05) %>%
  slice_head(n = 50) %>%
  pull(SYMBOL)

## wide matrix: rows = genes, cols = studies, values = logFC
mat_logfc <- by_study_deg %>%
  filter(SYMBOL %in% top_genes) %>%
  dplyr::select(SYMBOL, study_id, logFC) %>%
  pivot_wider(names_from = study_id, values_from = logFC) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

## optional: remove genes with too many NAs
keep_rows <- rowSums(is.na(mat_logfc)) < ncol(mat_logfc)  # seen in at least 1 study
mat_logfc <- mat_logfc[keep_rows, , drop = FALSE]

# Any non-finite values at all?
any(!is.finite(mat_logfc))

# Where are they?
which(!is.finite(mat_logfc), arr.ind = TRUE)
## heatmap (scaled per gene)
mat_logfc_imp <- mat_logfc
mat_logfc_imp[!is.finite(mat_logfc_imp)] <- 0  # replace NA/NaN/Inf with 0

# drop zero-variance rows after imputation
var0_rows <- apply(mat_logfc_imp, 1, function(x) sd(x, na.rm = TRUE) == 0)
mat_logfc_imp <- mat_logfc_imp[!var0_rows, ]

pheatmap(
  mat_logfc_imp,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top meta-DEGs (logFC across studies)"
)

#4c) Region moderation (optional)

#meta_one_gene_mod <- function(df) {
# required cols
# if (!all(c("yi","vi","region") %in% names(df)))
#  return(tibble::tibble(Qm_p = NA_real_, k = nrow(df), tau2 = NA_real_, method = NA, note = "missing yi/vi/region"))

# clean
#df <- df |>
# dplyr::mutate(yi = as.numeric(yi),
#              vi = as.numeric(vi),
#             region = as.factor(region)) |>
#dplyr::filter(is.finite(yi), is.finite(vi))
#if (nrow(df) < 2)
# return(tibble::tibble(Qm_p = NA_real_, k = nrow(df), tau2 = NA_real_, method = NA, note = "k<2 after cleaning"))

# strictly positive variances
#df$vi[df$vi <= 0] <- .Machine$double.eps

#  k <- nrow(df); L <- nlevels(df$region)
#method <- if (k <= L) "FE" else "REML"   # <- key change

#  fit <- try(metafor::rma(yi = yi, vi = vi, mods = ~ region, method = method, data = df), silent = TRUE)
#if (inherits(fit, "try-error"))
# return(tibble::tibble(Qm_p = NA_real_, k = k, tau2 = NA_real_, method = method,
#                     note = paste0("rma error: ", conditionMessage(fit))))

#  tibble::tibble(
#  Qm_p = unname(fit$QMp),                 # omnibus moderator p-value
#    k    = k,
#  tau2 = if (method == "FE") NA_real_ else unname(fit$tau2),
#  method = method,
#  note = NA_character_
# )
#}

#res <- meta_ready |>
# dplyr::group_by(SYMBOL) |>
#  dplyr::group_modify(~ meta_one_gene_mod(.x)) |>
# dplyr::ungroup()

#res |> dplyr::count(method, note, sort = TRUE) |> print(n = 20)
#region_mod <- res |> dplyr::mutate(Qm_FDR = p.adjust(Qm_p, method = "BH"))
#readr::write_csv(region_mod, file.path(output_dir, "meta_region_moderator.csv"))

#4d) LOO sensitivity
loo_one_gene <- function(d){
  d <- d %>% dplyr::filter(is.finite(SE))  # need SE for metafor
  if (nrow(d) < 3) return(tibble::tibble())  # need ???3 studies
  
  res <- purrr::map_df(unique(d$study_id), function(s){
    dd <- d %>% dplyr::filter(study_id != s)
    fit <- tryCatch(metafor::rma.uni(yi = dd$logFC, vi = dd$SE^2, method = "DL"),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble::tibble(omit = s, p = NA_real_, est = NA_real_))
    tibble::tibble(omit = s, p = fit$pval, est = as.numeric(fit$b))
  })
  
  tibble::tibble(
    SYMBOL    = d$SYMBOL[1],
    max_loo_p = suppressWarnings(max(res$p, na.rm = TRUE)),
    min_loo_p = suppressWarnings(min(res$p, na.rm = TRUE))
  )
}

meta_loo <- meta_ready %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::group_modify(~ loo_one_gene(.x)) %>%
  dplyr::ungroup()

 readr::write_csv(meta_loo, file.path(output_dir, "meta_leave_one_out.csv"))

#4e) Call consensus meta-DEGs
meta_sig <- meta_eff |> filter(meta_FDR <= 0.05, abs(meta_logFC) >= 0.3)
meta_up   <- meta_sig |> filter(meta_logFC > 0) |> arrange(meta_FDR)
meta_down <- meta_sig |> filter(meta_logFC < 0) |> arrange(meta_FDR)
 readr::write_lines(meta_up$SYMBOL,   file.path(output_dir, "consensus_up.txt"))
 readr::write_lines(meta_down$SYMBOL, file.path(output_dir, "consensus_down.txt"))

#5) Rank-based meta (uses everyone, even without SE)
# rank_lists <- by_study_deg %>%
#   dplyr::mutate(P.Value = dplyr::coalesce(P.Value, 1)) %>%
#   dplyr::group_by(study_id) %>%
#   dplyr::arrange(P.Value, .by_group = TRUE) %>%
#   dplyr::mutate(rank = dplyr::row_number()) %>%
#   dplyr::ungroup()
# 
# lists_for_rra <- rank_lists %>%
#   dplyr::group_by(study_id) %>%
#   dplyr::summarise(genes = list(SYMBOL), .groups = "drop") %>%
#   dplyr::pull(genes)
# 
# library(rlang)
# rra <- RobustRankAggreg::aggregateRanks(glist = lists_for_rra, method = "RRA") %>%
#   as_tibble() %>%
#   { if (!is.null(rownames(.)) && !any(c("Name","Item") %in% names(.)))
#     tibble::rownames_to_column(., var = "Name") else . } %>%
#   { names(.) <- tolower(names(.)); . } %>%
#   {
#     nm <- names(.)
#     id_col <- intersect(c("symbol","gene","name","item"), nm)[1]
#     p_col  <- intersect(c("rra_p","pvalue","p.value","p_value","p","score"), nm)[1]
#     if (is.na(id_col) || is.na(p_col)) {
#       stop(paste0("Couldn't find ID or p-value columns. Available columns: ",
#                   paste(nm, collapse = ", ")))
#     }
#     dplyr::rename(., SYMBOL = dplyr::all_of(id_col),
#                   RRA_P  = dplyr::all_of(p_col))
#   } %>%
#   dplyr::mutate(RRA_FDR = p.adjust(RRA_P, method = "BH")) %>%
#   dplyr::arrange(RRA_P)
# rra
# 
# readr::write_csv(rra, file.path(output_dir, "rra_results.csv"))
# 
# rra <- rra %>%
#   mutate(RRA_rank = row_number()) %>%
#   relocate(RRA_rank, .before = SYMBOL)
# 
# top50 <- rra %>% slice_head(n = 50)
# sig   <- rra %>% filter(RRA_FDR <= 0.05)
# readr::write_csv(sig, "rra_significant_fdr05.csv")

# library(ggplot2)
# ggplot(rra, aes(x = -log10(RRA_P))) +
#   geom_histogram(bins = 50) +
#   labs(x = expression(-log[10]~RRA[p]), y = "Count",
#        title = "RRA p-value distribution")


# MetaVolcanoR (if most studies have logFC + p)
mv_list <- split(by_study_deg, by_study_deg$study_id) %>%
  purrr::map(~ .x %>%
               dplyr::select(SYMBOL, logFC, P.Value) %>%
               dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
               dplyr::rename(p = P.Value, effectsize = logFC) %>%
               dplyr::mutate(p = as.numeric(p), effectsize = as.numeric(effectsize)) %>%
               dplyr::filter(is.finite(p), is.finite(effectsize))
  )
library(purrr)
library(dplyr)
# Patch ggplot2 class to restore compatibility
library(ggplot2)
# stack to long, count valid p's per gene
mv_long <- purrr::imap_dfr(mv_list, ~ dplyr::mutate(.x, study = .y))
k_tbl <- mv_long %>%
  filter(is.finite(p)) %>%
  count(SYMBOL, name = "k_valid")

keep_syms <- k_tbl %>% filter(k_valid >= 3) %>% pull(SYMBOL)

# filter each study to those genes
mv_list2 <- purrr::imap(mv_list, ~ dplyr::filter(.x, SYMBOL %in% keep_syms))

#remotes::install_version("ggplot2", version = "3.5.2")
dir.create(file.path(output_dir, "metavolcano"), recursive = TRUE, showWarnings = FALSE)
source(file.path(project_root, "scripts", "03_meta_analysis", "metavolcano_compat.R"))

mv <- combining_mv_compatible(
  diffexp       = mv_list2,
  pcriteria     = "p",
  foldchangecol = "effectsize",
  genenamecol   = "SYMBOL",
  metafc        = "Mean",
  metathr       = 0.05,
  jobname       = "AD_brain_meta",
  outputfolder  = file.path(output_dir, "metavolcano"),
  draw          = "HTML"
)



# results + how many studies contributed
meta_tbl <- mv$metaresult %>%
  left_join(k_tbl, by = "SYMBOL") %>%
  mutate(FDR = p.adjust(metap, method = "BH"))
 readr::write_csv(meta_tbl, file.path(output_dir, "metavolcano/updated_metavolcano_combining_results.csv"))
#6) Quick sanity checks you???ll want
# How many per study?
table(all_deg$acc)

# Which studies lacked SE (will be omitted from effect-size meta)?
all_deg |> mutate(hasSE = !is.na(SE)) |>
  group_by(acc) |> summarise(with_SE = mean(hasSE), n = n())

# Direction agreement for meta-significant genes
consistency <- meta_in |>
  semi_join(meta_sig, by="SYMBOL") |>
  mutate(sign = sign(logFC)) |>
  group_by(SYMBOL) |>
  summarise(k = n(), agree_frac = max(mean(sign==1), mean(sign==-1)), .groups="drop") |>
  arrange(desc(agree_frac))
readr::write_csv(consistency, file.path(output_dir, "direction_consistency.csv"))
