## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
guess_gene_col <- function(df){
  nms <- names(df)
  cand <- nms[grepl("^(gene|genes|symbol|hgnc|gene_name)$", tolower(nms))]
  if (length(cand) >= 1) return(cand[1])
  cand2 <- nms[grepl("ensembl|gene_id|ensembl_gene", tolower(nms))]
  if (length(cand2) >= 1) return(cand2[1])
  # avoid unnamed index-like first column if it’s just 1:n
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
  
  ## Optional per-dataset screening (fold-change sensitivity branch only).
  ## Applied HERE — on raw rows, per study arm, before the one-row-per-SYMBOL
  ## collapse and before GSE278723 subregion pooling. That ordering matters: the
  ## primary branch reads tables that were already screened on disk, so its
  ## collapse picks the best row among *passing* rows. Screening after the
  ## collapse instead would pick a representative row that may then be
  ## discarded, and does not reproduce the primary run. NULL for the two
  ## prespecified branches, which are unaffected.
  if (exists("screen_lfc_cut", inherits = TRUE) && !is.null(screen_lfc_cut)) {
    df <- df |>
      dplyr::filter(is.finite(suppressWarnings(as.numeric(P.Value))),
                    is.finite(suppressWarnings(as.numeric(logFC))),
                    suppressWarnings(as.numeric(P.Value)) <= screen_p_cut,
                    abs(suppressWarnings(as.numeric(logFC))) > screen_lfc_cut)
  }
  
  # One row per SYMBOL per study (lowest p; tie → largest |logFC|)
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

if (exists("screen_lfc_cut", inherits = TRUE) && !is.null(screen_lfc_cut)) {
  message("Per-dataset screen applied: P <= ", screen_p_cut,
          " & |logFC| > ", screen_lfc_cut)
  print(dplyr::count(all_deg, acc, name = "eligible_genes"))
}

readr::write_csv(all_deg, file.path(output_dir, "all_studies_harmonized.csv"))
message(
  "Harmonized rows: ", nrow(all_deg),
  " | Genes: ", dplyr::n_distinct(all_deg$SYMBOL),
  " | Studies: ", dplyr::n_distinct(all_deg$acc)
)
