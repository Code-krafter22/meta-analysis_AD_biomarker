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
