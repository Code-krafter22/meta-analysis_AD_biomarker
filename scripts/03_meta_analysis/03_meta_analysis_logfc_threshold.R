# =============================================================================
# SENSITIVITY BRANCH 3 — FOLD-CHANGE THRESHOLD
#
# Reviewer-requested sensitivity analysis that HOLDS the significance criterion
# at the primary setting (unadjusted P <= 0.05) and varies ONLY the fold-change
# threshold. No combined FDR + logFC run is produced: the point is to move one
# knob at a time.
#
#   Primary       : unadj. P <= 0.05 & |log2FC| > 0.58
#   Sensitivity 3a: unadj. P <= 0.05 & |log2FC| > 1.00   (2-fold)
#   Sensitivity 3b: unadj. P <= 0.05 & |log2FC| > 0.25   (~1.19-fold)
#
# INPUTS. Both runs are driven from results/differential_expression/**/
# full_DE_results.csv — the unfiltered per-dataset tables retained for the
# existing no-screening sensitivity branch. The pre-filtered
# significant_DEGs_unadjusted.csv tables are NOT used: they were written to disk
# already screened at P <= 0.05 & |logFC| > 0.58 (verified: min |logFC| = 0.580
# in all nine study arms), so the 0.25 variant is unrecoverable from them.
# Driving both thresholds off the unfiltered tables avoids the subset/re-run
# asymmetry entirely.
#
# PIPELINE. Runs the same shared engine as the primary analysis
# (run_meta_analysis.R -> stages 01-05). Nothing downstream is reimplemented:
# per-dataset screening, GSE278723 inverse-variance subregion collapsing,
# direction-aware Stouffer, DerSimonian-Laird random effects, MetaVolcanoR,
# and leave-one-dataset-out all come from the shared stages.
#
# MetaVolcanoR / ggplot2 4.0.x. The engine's stage 05 sources the existing
# compatibility shim scripts/03_meta_analysis/metavolcano_compat.R, which is
# reused unchanged. NOTE: that shim sidesteps the S4 validity failure by calling
# MetaVolcanoR:::plot_mv directly rather than by unlocking ".__C__MetaVolcano"
# and redefining its slots as "ANY"; no unlockBinding patch exists in this repo.
#
# SEED. No stage is stochastic (DL and Stouffer are closed-form, LOO is
# exhaustive). A seed is set anyway so the run is bit-reproducible if a future
# stage introduces sampling.
#
# Package versions are not hard-coded here; each run writes its own
# sessionInfo.txt via the engine, and this script also writes a combined one.
# =============================================================================

set.seed(20260903)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(here)
})

project_root <- here::here()
script_dir   <- file.path(project_root, "scripts", "03_meta_analysis")
engine       <- file.path(script_dir, "run_meta_analysis.R")


# ----------------------------------------------------------- the 37 signature
# Typos carried by some downstream source files: CCD102A -> CCDC102A,
# HSBP7 -> HSPB7. The per-dataset DE tables themselves use the correct symbols,
# but the map is applied to every list read here so matching cannot silently
# miss them. COL25A1 is present: one validation script elsewhere carries only 36.
symbol_fix <- c(CCD102A = "CCDC102A", HSBP7 = "HSPB7")
harmonise <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(x %in% names(symbol_fix), symbol_fix[x], x)
}

genes37 <- harmonise(readr::read_lines(
  file.path(project_root, "results", "figures", "intermediate", "figure1",
            "signature_37_genes.txt")))
genes37 <- genes37[nzchar(genes37)]
stopifnot(length(genes37) == 37L, !anyDuplicated(genes37), "COL25A1" %in% genes37)

primary_direction <- readr::read_csv(
  file.path(project_root, "results", "figures", "intermediate", "figure2",
            "consensus_37_direction.csv"), show_col_types = FALSE) |>
  dplyr::mutate(SYMBOL = harmonise(SYMBOL),
                primary_dir = tolower(Direction)) |>
  dplyr::select(SYMBOL, primary_dir)
stopifnot(setequal(primary_direction$SYMBOL, genes37))

# ------------------------------------------------- consensus calling criteria
# Read off the primary analysis, not assumed. Applied to
# results/meta_analysis/filtered_unadjusted this returns 38 genes: all 37 of the
# published signature plus COL8A2 (see the audit written by this script).
#
#   1. quantified and eligible in >= 4 of the 5 datasets (GSE278723 counts once)
#   2. DL random effects: meta_FDR <= 0.05 AND |meta_logFC| >= 0.3
#      (the gate coded in stages/04_effect_size_meta_analysis.R)
#   3. Stouffer FDR <= 0.05
#   4. MetaVolcanoR combining-method FDR <= 0.05
#   5. I2 < 40%
#   6. max raw leave-one-dataset-out P < 0.05
#
# NOTE on the k gate. The shared engine pools any gene seen in >= 3 studies;
# the >= 4 requirement is a consensus-calling gate, not a pooling gate. Both are
# preserved exactly as the primary analysis has them and both are reported.
K_MIN_POOL      <- 3L
K_MIN_CONSENSUS <- 4L
FDR_CUT         <- 0.05
ABS_LFC_CUT     <- 0.30
I2_CUT          <- 40
LOO_P_CUT       <- 0.05

collect_run <- function(dir) {
  eff <- readr::read_csv(file.path(dir, "meta_effects_meta_results.csv"),
                         show_col_types = FALSE) |>
    dplyr::mutate(SYMBOL = harmonise(SYMBOL))
  sto <- readr::read_csv(file.path(dir, "stouffer_meta_results.csv"),
                         show_col_types = FALSE) |>
    dplyr::mutate(SYMBOL = harmonise(SYMBOL)) |>
    dplyr::select(SYMBOL, stouffer_Z = Z, stouffer_P = P, stouffer_FDR = FDR,
                  stouffer_Dir = Dir)
  mvc <- readr::read_csv(file.path(dir, "metavolcano",
                                   "updated_metavolcano_combining_results.csv"),
                         show_col_types = FALSE) |>
    dplyr::mutate(SYMBOL = harmonise(SYMBOL)) |>
    dplyr::select(SYMBOL, mv_P = metap, mv_FDR = FDR)
  loo <- readr::read_csv(file.path(dir, "meta_leave_one_out.csv"),
                         show_col_types = FALSE) |>
    dplyr::mutate(SYMBOL = harmonise(SYMBOL)) |>
    dplyr::select(SYMBOL, max_loo_p, min_loo_p)

  # k = number of distinct STUDIES contributing (GSE278723 subregions collapse
  # to one study before this table is written).
  combined <- readr::read_csv(file.path(dir, "ALL_DATASET_COMBINED.csv"),
                              show_col_types = FALSE) |>
    dplyr::mutate(SYMBOL = harmonise(SYMBOL),
                  study_id = sub("_.*$", "", acc))
  k_tbl <- combined |>
    dplyr::distinct(SYMBOL, study_id) |>
    dplyr::count(SYMBOL, name = "k_studies")

  stats <- eff |>
    dplyr::left_join(sto, by = "SYMBOL") |>
    dplyr::left_join(mvc, by = "SYMBOL") |>
    dplyr::left_join(loo, by = "SYMBOL") |>
    dplyr::left_join(k_tbl, by = "SYMBOL") |>
    dplyr::mutate(
      direction   = ifelse(meta_logFC > 0, "up", "down"),
      ci_95       = sprintf("[%.3f, %.3f]", ci_lb, ci_ub),
      pass_k      = k_studies >= K_MIN_POOL,
      pass_k4     = k_studies >= K_MIN_CONSENSUS,
      pass_dl     = meta_FDR <= FDR_CUT & abs(meta_logFC) >= ABS_LFC_CUT,
      pass_stou   = !is.na(stouffer_FDR) & stouffer_FDR <= FDR_CUT,
      pass_mv     = !is.na(mv_FDR) & mv_FDR <= FDR_CUT,
      pass_i2     = I2 < I2_CUT,
      pass_loo    = !is.na(max_loo_p) & max_loo_p < LOO_P_CUT,
      # I2 is NOT part of the consensus call: Panel A reports it as a separate
      # column, and the two genuinely differ (see the no-screening row, where
      # 37 are retained but only 22 have I2 < 40%).
      consensus   = pass_k & pass_dl & pass_stou & pass_mv & pass_loo
    )

  list(stats = stats, combined = combined)
}

# k per study for genes that did NOT reach the meta-analysis at all, so dropouts
# can be attributed to the screening step rather than to pooling.
eligibility_of <- function(combined, symbols) {
  combined |>
    dplyr::filter(SYMBOL %in% symbols) |>
    dplyr::distinct(SYMBOL, study_id) |>
    dplyr::count(SYMBOL, name = "k_eligible") |>
    dplyr::right_join(tibble::tibble(SYMBOL = symbols), by = "SYMBOL") |>
    dplyr::mutate(k_eligible = tidyr::replace_na(k_eligible, 0L))
}

# ------------------------------------------------------------------ run both
# 0.58 is re-run as a harness control: it must reproduce the primary analysis
# exactly. 1.00 and 0.25 are the reviewer-requested thresholds.
# 0.00 is a P-only reference run (unadjusted P <= 0.05, no magnitude gate). It
# is not a reported specification: it exists so a dropped gene can be attributed
# exactly to the magnitude gate rather than to the significance gate. A gene with
# k_P >= 4 but k < 4 at the threshold was removed BY the fold-change cut; a gene
# with k_P < 4 was never eligible in 4 cohorts on significance grounds at any
# threshold.
thresholds <- c(0.58, 1.00, 0.25, 0.00)
run_dirs <- setNames(
  file.path(project_root, "results", "meta_analysis",
            paste0("logfc_gt_", formatC(thresholds, format = "f", digits = 2))),
  formatC(thresholds, format = "f", digits = 2))

for (cut in thresholds) {
  message("\n#################### |log2FC| > ", cut, " ####################")
  analysis_mode  <- "logfc_threshold"
  screen_lfc_cut <- cut
  screen_p_cut   <- 0.05
  Sys.setenv(META_ANALYSIS_OUTPUT_DIR =
               run_dirs[[formatC(cut, format = "f", digits = 2)]])
  source(engine)
  Sys.unsetenv("META_ANALYSIS_OUTPUT_DIR")
}

# ------------------------------------------------------- assemble the outputs
out_root <- file.path(project_root, "results", "meta_analysis",
                      "logfc_threshold_sensitivity")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

primary_dir <- file.path(project_root, "results", "meta_analysis",
                         "filtered_unadjusted")
ponly_dir   <- run_dirs[["0.00"]]

all_dirs <- list(
  "primary"    = primary_dir,
  "s1"         = file.path(project_root, "results", "meta_analysis",
                           "unfiltered_all_genes"),
  "s2"         = file.path(project_root, "results", "meta_analysis",
                           "fdr_screened"),
  "s3a"        = run_dirs[["1.00"]],
  "s3b"        = run_dirs[["0.25"]],
  "control058" = run_dirs[["0.58"]],
  "ponly"      = ponly_dir
)
# Drop any prespecified branch that has not been generated.
all_dirs <- all_dirs[vapply(all_dirs, function(d)
  file.exists(file.path(d, "meta_effects_meta_results.csv")), logical(1))]
spec_meta <- list(
  primary    = list(name = "Primary analysis",
                    crit = "Unadjusted P <= 0.05 and |log2FC| > 0.58",
                    where = "Main text"),
  control058 = list(name = "Harness control (re-run of primary at 0.58)",
                    crit = "Unadjusted P <= 0.05 and |log2FC| > 0.58",
                    where = "Not reported (control)"),
  s1         = list(name = "Sensitivity analysis 1 (no screening)",
                    crit = "No P-value or fold-change screening",
                    where = "Panel C"),
  s2         = list(name = "Sensitivity analysis 2 (FDR)",
                    crit = "BH FDR <= 0.05 and |log2FC| > 0.58",
                    where = "Panel C"),
  s3a        = list(name = "Sensitivity analysis 3a",
                    crit = "Unadjusted P <= 0.05 and |log2FC| > 1.0",
                    where = "Panel C"),
  s3b        = list(name = "Sensitivity analysis 3b",
                    crit = "Unadjusted P <= 0.05 and |log2FC| > 0.25",
                    where = "Panel C"),
  ponly      = list(name = "P-only reference (no magnitude gate)",
                    crit = "Unadjusted P <= 0.05, no |log2FC| criterion",
                    where = "Not reported (attribution reference)")
)

panelA <- list(); panelC <- list(); dropout <- list()
retained_sets <- list(); consensus_sets <- list()

# Datasets in which each gene passes the significance gate alone.
k_pvalue_only <- readr::read_csv(file.path(ponly_dir, "ALL_DATASET_COMBINED.csv"),
                                 show_col_types = FALSE) |>
  dplyr::mutate(SYMBOL = harmonise(SYMBOL), study_id = sub("_.*$", "", acc)) |>
  dplyr::distinct(SYMBOL, study_id) |>
  dplyr::count(SYMBOL, name = "k_p_only")

for (nm in names(all_dirs)) {
  r <- collect_run(all_dirs[[nm]])
  st <- r$stats; cb <- r$combined

  cons <- st |> dplyr::filter(consensus)
  retained <- intersect(genes37, cons$SYMBOL)

  crit      <- spec_meta[[nm]]$crit
  spec_name <- spec_meta[[nm]]$name
  retained_sets[[nm]] <- retained
  consensus_sets[[nm]] <- sort(cons$SYMBOL)

  panelA[[nm]] <- tibble::tibble(
    Specification = spec_name,
    `Per-study eligibility criteria` = crit,
    `Genes entering meta-analysis` = sum(st$k_studies >= K_MIN_POOL, na.rm = TRUE),
    `Genes eligible in >= 4 of 5 datasets` = sum(st$pass_k4, na.rm = TRUE),
    `Consensus genes retained (of 37)` = length(retained),
    `Retained with I2 < 40%` = sum(st$SYMBOL %in% retained & st$pass_i2),
    `Reported in` = spec_meta[[nm]]$where
  )

  # Panel C: per-gene results for the retained consensus genes.
  panelC[[nm]] <- st |>
    dplyr::filter(SYMBOL %in% retained) |>
    dplyr::left_join(primary_direction, by = "SYMBOL") |>
    dplyr::transmute(
      Specification = spec_name,
      Gene = SYMBOL,
      `Datasets contributing (k)` = k_studies,
      `Pooled log2FC` = meta_logFC,
      SE = meta_se,
      `95% CI` = ci_95,
      Z = meta_z,
      `Pooled P (raw)` = meta_p,
      `Pooled P (FDR)` = meta_FDR,
      Direction = direction,
      `Direction matches primary` = direction == primary_dir,
      `I2 (%)` = I2,
      `Max raw LOO P` = max_loo_p
    ) |>
    dplyr::arrange(`Pooled P (FDR)`)

  # Full consensus audit (all genes passing, not only the 37) — this is where
  # the primary-run COL8A2 discrepancy is visible.
  readr::write_csv(cons, file.path(out_root, paste0(
    "consensus_audit_", gsub("[^0-9a-z]+", "_", tolower(nm)), ".csv")))

  # Dropout attribution for the 37, staged. The P-only reference describes the
  # significance/magnitude split for the fold-change specifications; it is not
  # a meaningful contrast for the no-screening or FDR branches.
  if (nm %in% c("s1", "s2")) next

  elig <- eligibility_of(cb, genes37) |>
    dplyr::left_join(k_pvalue_only, by = "SYMBOL") |>
    dplyr::mutate(k_p_only = tidyr::replace_na(k_p_only, 0L))
  dropout[[nm]] <- elig |>
    dplyr::left_join(st, by = "SYMBOL") |>
    dplyr::transmute(
      Specification = spec_name,
      Gene = SYMBOL,
      `Datasets eligible after screening` = k_eligible,
      Retained = SYMBOL %in% retained,
      `Datasets passing P alone` = k_p_only,
      `Datasets lost to magnitude gate` = k_p_only - k_eligible,
      Stage = dplyr::case_when(
        SYMBOL %in% retained ~ "retained",
        k_eligible < K_MIN_POOL & k_p_only >= K_MIN_CONSENSUS ~
          "dropped at per-dataset screening - failed the MAGNITUDE gate",
        k_eligible < K_MIN_POOL ~
          "dropped at per-dataset screening - not P-significant in >= 4 cohorts",
        is.na(meta_p) ~ "dropped at pooling (no poolable estimate)",
        TRUE ~ "dropped at pooling (failed a significance/heterogeneity gate)"),
      `Failed gate` = dplyr::case_when(
        SYMBOL %in% retained ~ NA_character_,
        k_eligible < K_MIN_POOL ~ NA_character_,
        TRUE ~ trimws(paste(
          ifelse(!pass_dl %in% TRUE, "DL", ""),
          ifelse(!pass_stou %in% TRUE, "Stouffer", ""),
          ifelse(!pass_mv %in% TRUE, "MetaVolcano", ""),
          ifelse(!pass_i2 %in% TRUE, "I2>=40", ""),
          ifelse(!pass_loo %in% TRUE, "LOO", "")))),
      `Pooled log2FC` = meta_logFC,
      `Pooled P (FDR)` = meta_FDR,
      `I2 (%)` = I2,
      `Max raw LOO P` = max_loo_p
    ) |>
    dplyr::arrange(Retained, `Datasets eligible after screening`, Gene)
}

# Harness control: re-running the engine at the primary 0.58 gate off the
# unfiltered tables must reproduce the primary analysis exactly.
if (!identical(consensus_sets$control058, consensus_sets$primary)) {
  stop("Harness control failed: 0.58 re-run does not reproduce the ",
       "primary consensus set.")
}
message("Harness control PASSED: 0.58 re-run reproduces the primary analysis (",
        length(consensus_sets$primary), " consensus genes).")

panelA_tbl <- dplyr::bind_rows(panelA)
readr::write_csv(panelA_tbl, file.path(out_root, "TableS6c_PanelA.csv"))
readr::write_csv(dplyr::bind_rows(panelC), file.path(out_root, "TableS6c_PanelC.csv"))
readr::write_csv(dplyr::bind_rows(dropout), file.path(out_root, "TableS6c_dropout_stage.csv"))

writeLines(capture.output(sessionInfo()),
           file.path(out_root, "sessionInfo.txt"))

message("\n===================== Supplementary Table S6c, Panel A =====================")
print(as.data.frame(panelA_tbl), row.names = FALSE)
message("\nWrote: ", out_root)
