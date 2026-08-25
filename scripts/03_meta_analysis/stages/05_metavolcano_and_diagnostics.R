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
#6) Quick sanity checks you’ll want
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
