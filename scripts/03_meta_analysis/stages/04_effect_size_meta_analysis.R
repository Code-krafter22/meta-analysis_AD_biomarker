#4a) Keep genes seen in ≥3 studies (tuneable)
# keep genes seen in ≥3 datasets (now study-level, not region-level)
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

## choose top 50 significant genes – adjust N as you like
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
  main = "Top meta-DEGs (logFC across studies)",
  filename = file.path(output_dir, "top50_meta_DEGs_heatmap.pdf"),
  width = 8,
  height = 10
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
  if (nrow(d) < 3) return(tibble::tibble())  # need ≥3 studies
  
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
