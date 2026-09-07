# Reviewer 2, Major Comment 2: covariate-adjusted sensitivity analysis.
#
# For each discovery dataset (and each GSE278723 subregion) fit two limma-voom
# models on the SAME samples and the SAME filtered/TMM-normalised count matrix:
#
#   Fit A (unadjusted): ~ Condition            <- reproduces the primary pipeline
#   Fit B (adjusted):   ~ Condition + <covariates recorded for that dataset>
#
# Nothing in results/ is touched. Everything is written under
# Results_covariate_adjusted/. Run 02_check_37_genes.R afterwards.

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(readxl); library(here); library(GEOquery)
})

outdir  <- here::here("Results_covariate_adjusted")
de_dir  <- file.path(outdir, "de_tables")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

COMPLETENESS_MIN <- 0.90   # covariate must be recorded for >= 90% of samples

# ---------------------------------------------------------------- helpers ----

# Covariate availability log; one row per dataset x covariate.
avail <- list()
log_cov <- function(dataset, covariate, status, reason = "") {
  avail[[length(avail) + 1L]] <<- data.frame(
    Dataset = dataset, Covariate = covariate, Status = status, Reason = reason
  )
}

# Decide whether a covariate column is usable. Returns TRUE/FALSE and logs why.
usable <- function(dataset, name, x, not_deposited = FALSE) {
  if (not_deposited || is.null(x)) {
    log_cov(dataset, name, "Not deposited", "not present in GEO record")
    return(FALSE)
  }
  complete <- mean(!is.na(x))
  if (complete < COMPLETENESS_MIN) {
    log_cov(dataset, name, "Recorded but not used",
            sprintf("only %.0f%% complete (< %.0f%% threshold)",
                    100 * complete, 100 * COMPLETENESS_MIN))
    return(FALSE)
  }
  if (length(unique(x[!is.na(x)])) < 2) {
    log_cov(dataset, name, "Recorded but not used", "no variation")
    return(FALSE)
  }
  log_cov(dataset, name, "Included",
          sprintf("%.0f%% complete", 100 * complete))
  TRUE
}

# Some curated metadata workbooks under data/metadata/ hold only a subset of the
# fields GEO actually deposits (GSE95587 omits PMI and RIN, for example), so
# characteristics are pulled from the series matrix and cached alongside them.
geo_characteristics <- function(acc) {
  cache <- here::here("data", "metadata", acc,
                      paste0(acc, "_geo_characteristics.csv"))
  if (!file.exists(cache)) {
    message("  fetching GEO characteristics for ", acc)
    gse <- GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)
    p <- Biobase::pData(gse[[1]])
    p <- p[, c("geo_accession", "title", grep(":ch1$", names(p), value = TRUE)),
           drop = FALSE]
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    write.csv(p, cache, row.names = FALSE)
  }
  read.csv(cache, check.names = FALSE, stringsAsFactors = FALSE)
}

# GSE67333 PMI lives only in the source paper's supplement (Magistri et al.
# 2015, J Alzheimers Dis 48:647, Supplementary Table 1) - it is not in the GEO
# record. Donors are matched to the eight sequenced samples on the four fields
# GEO does deposit: age, gender, APOE genotype and Braak stage, restricted to
# donors with hippocampal RNA available. Six of eight match uniquely; LOAD1 and
# LOAD4 are indistinguishable on all four fields (both female, 82, 3/3, Braak
# VI) so that pair is assigned in file order. Both are AD and their PMI differs
# by 0.67 h, so the ambiguity cannot affect the diagnosis coefficient.
gse67333_pmi <- function() {
  cache <- here::here("data", "metadata", "GSE67333", "GSE67333_donor_PMI.csv")
  if (!file.exists(cache)) {
    supp <- Sys.getenv("GSE67333_SUPP",
                       "~/Downloads/jad-48-3-jad150398-s001.xlsx")
    supp <- path.expand(supp)
    if (!file.exists(supp)) {
      warning("GSE67333 supplement not found at ", supp,
              " - PMI will be unavailable. Set GSE67333_SUPP.")
      return(NULL)
    }
    s <- as.data.frame(read_excel(supp, sheet = 1, col_names = FALSE))
    names(s) <- as.character(unlist(s[1, ])); s <- s[-1, ]
    s <- s[grepl("Hipp", s$`RNA tissue available`), ]
    key <- function(age, sex, apoe, braak)
      paste(trimws(age), tolower(substr(trimws(sex), 1, 1)),
            tolower(trimws(apoe)), toupper(trimws(braak)), sep = "|")
    s$key <- key(s$Expired_age, s$Gender, s$ApoE, s$`Braak Stage`)

    geo <- geo_characteristics("GSE67333")
    geo$key <- key(geo$`age (yrs):ch1`, geo$`gender:ch1`,
                   geo$`apoe genotype:ch1`, geo$`braak stage:ch1`)
    # Consume supplement rows in order so the LOAD1/LOAD4 tie resolves stably.
    used <- rep(FALSE, nrow(s))
    pick <- vapply(geo$key, function(k) {
      i <- which(s$key == k & !used)
      if (!length(i)) return(NA_integer_)
      used[i[1]] <<- TRUE
      i[1]
    }, integer(1))
    if (anyNA(pick))
      stop("GSE67333: could not match ", sum(is.na(pick)),
           " GEO samples to the supplement")
    out <- data.frame(Sample = geo$title, CaseID = s$CaseID[pick],
                      PMI = suppressWarnings(as.numeric(s$PMI[pick])),
                      stringsAsFactors = FALSE)
    write.csv(out, cache, row.names = FALSE)
  }
  read.csv(cache, stringsAsFactors = FALSE)
}

read_counts_csv <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(raw)[1] <- "Gene"
  colnames(raw) <- sub("^X", "", colnames(raw))
  raw$Gene <- make.unique(as.character(raw$Gene))
  raw[is.na(raw)] <- 0
  m <- as.matrix(sapply(raw[, -1, drop = FALSE], as.numeric))
  rownames(m) <- raw$Gene
  m
}

# Filtering + TMM exactly as in scripts/02_differential_expression/*.Rmd.
build_dge <- function(counts_mat, metadata) {
  dge <- DGEList(counts = counts_mat, samples = metadata,
                 group = metadata$Condition)
  dge <- calcNormFactors(dge, method = "TMM")
  min_n <- min(table(metadata$Condition))
  keep  <- rowSums(cpm(dge, log = TRUE, prior.count = 1) > 1) >= min_n
  dge   <- dge[keep, , keep.lib.sizes = FALSE]
  calcNormFactors(dge, method = "TMM")
}

fit_model <- function(dge, metadata, formula_rhs) {
  design <- model.matrix(as.formula(paste("~", formula_rhs)), data = metadata)
  rownames(design) <- colnames(dge)
  if (qr(design)$rank < ncol(design))
    return(list(estimable = FALSE, design = design))
  if (nrow(design) - ncol(design) < 1L)
    return(list(estimable = FALSE, design = design))
  v   <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "ConditionAD", number = Inf, sort.by = "none")
  tt$SE <- sqrt(fit$s2.post) * fit$stdev.unscaled[rownames(tt), "ConditionAD"]
  tt$Gene <- rownames(tt)
  list(estimable = TRUE, table = tt,
       residual_df = nrow(design) - ncol(design), design = design)
}

# Analyse one (dataset, region) unit: writes both DE tables, returns a summary row.
run_unit <- function(dataset, region, counts_mat, metadata, covars) {
  tag <- if (is.na(region)) dataset else paste0(dataset, "_", region)
  message("\n=== ", tag, " ===")

  # Fit A: all samples with a Condition label.
  metadata$Condition <- factor(metadata$Condition, levels = c("Control", "AD"))
  metadata <- metadata[!is.na(metadata$Condition), , drop = FALSE]
  counts_mat <- counts_mat[, metadata$Sample, drop = FALSE]

  dge_a <- build_dge(counts_mat, metadata)
  fa <- fit_model(dge_a, metadata, "Condition")
  stopifnot(fa$estimable)
  write.csv(fa$table[, c("Gene", "logFC", "AveExpr", "t", "P.Value",
                         "adj.P.Val", "SE")],
            file.path(de_dir, paste0(tag, "_unadjusted.csv")), row.names = FALSE)

  base <- data.frame(
    Dataset = dataset, Region = region,
    n_samples = nrow(metadata),
    n_AD = sum(metadata$Condition == "AD"),
    n_control = sum(metadata$Condition == "Control"),
    residual_df_unadjusted = fa$residual_df,
    stringsAsFactors = FALSE
  )

  if (!length(covars)) {
    message("  no covariates available - adjusted model not fitted")
    return(cbind(base, Covariates_in_model = "none (none deposited)",
                 residual_df_adjusted = NA_integer_,
                 n_samples_adjusted = NA_integer_,
                 Adjusted_model_status = "not fitted: no covariates deposited"))
  }

  # Fit B: complete cases for THIS dataset's covariates only.
  md_b <- metadata[stats::complete.cases(metadata[, covars, drop = FALSE]), ,
                   drop = FALSE]
  dropped <- nrow(metadata) - nrow(md_b)
  if (dropped) message("  dropped ", dropped, " samples with missing covariates")

  # Re-check variation after the complete-case drop.
  keep_cov <- covars[vapply(covars, function(cv)
    length(unique(md_b[[cv]])) > 1, logical(1))]
  if (!length(keep_cov) || min(table(md_b$Condition)) < 1) {
    return(cbind(base, Covariates_in_model = paste(covars, collapse = "+"),
                 residual_df_adjusted = NA_integer_,
                 n_samples_adjusted = nrow(md_b),
                 Adjusted_model_status = "not estimable after complete-case filter"))
  }

  dge_b <- build_dge(counts_mat[, md_b$Sample, drop = FALSE], md_b)
  rhs <- paste(c("Condition", keep_cov), collapse = " + ")
  fb <- fit_model(dge_b, md_b, rhs)
  if (!isTRUE(fb$estimable)) {
    message("  MODEL NOT ESTIMABLE: ", rhs)
    return(cbind(base, Covariates_in_model = rhs,
                 residual_df_adjusted = NA_integer_,
                 n_samples_adjusted = nrow(md_b),
                 Adjusted_model_status = "not estimable (rank deficient / no residual df)"))
  }
  write.csv(fb$table[, c("Gene", "logFC", "AveExpr", "t", "P.Value",
                         "adj.P.Val", "SE")],
            file.path(de_dir, paste0(tag, "_adjusted.csv")), row.names = FALSE)
  message("  adjusted: ~ ", rhs, "  (residual df = ", fb$residual_df, ")")

  cbind(base, Covariates_in_model = rhs,
        residual_df_adjusted = fb$residual_df,
        n_samples_adjusted = nrow(md_b),
        Adjusted_model_status = "fitted")
}

summaries <- list()

# ------------------------------------------------------------- GSE67333 -----
# 8 subjects. age/sex/APOE deposited; APOE is 7/8 complete so it is excluded by
# the completeness rule. PMI, RIN, batch not deposited.
{
  ds <- "GSE67333"
  counts <- read_counts_csv(here::here("data", "processed", ds,
                                       "Edited_raw_counts.csv"))
  md <- as.data.frame(read_excel(here::here("data", "metadata", ds,
                                            "metadata.xlsx")))
  md$Sample <- md$identifier
  # ">90" is a GEO privacy cap; treated as 90 for modelling purposes.
  age  <- suppressWarnings(as.numeric(sub("^>", "", md$age)))
  sex  <- factor(tolower(trimws(md$gender)))
  apoe <- sub("^apoe genotype:\\s*", "", md$APOE)
  apoe[grepl("not available", apoe, ignore.case = TRUE)] <- NA
  pmi <- gse67333_pmi()
  meta <- data.frame(Sample = md$Sample,
                     Condition = ifelse(grepl("^LOAD", md$Sample), "AD", "Control"),
                     age = age, sex = sex, APOE = factor(apoe),
                     PMI = if (is.null(pmi)) NA_real_
                           else pmi$PMI[match(md$Sample, pmi$Sample)],
                     stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$Sample), , drop = FALSE]

  cov <- character()
  if (usable(ds, "age", meta$age))  cov <- c(cov, "age")
  if (usable(ds, "sex", meta$sex))  cov <- c(cov, "sex")
  # PMI is from the source publication, not the GEO deposit.
  if (usable(ds, "PMI", meta$PMI))  cov <- c(cov, "PMI")
  usable(ds, "RIN",   NULL, not_deposited = TRUE)
  usable(ds, "batch", NULL, not_deposited = TRUE)
  if (usable(ds, "APOE", meta$APOE)) cov <- c(cov, "APOE")

  summaries[[ds]] <- run_unit(ds, NA_character_, counts, meta, cov)
}

# ------------------------------------------------------------- GSE95587 -----
# SEX and AGE come from the curated workbook. GEO additionally deposits
# post_mortem_interval and rin_score for all 117 samples; those two fields are
# absent from metadata.xlsx, so they are joined in from the series matrix via
# the internal sample_id (column 3 of the workbook). APOE and batch: not deposited.
{
  ds <- "GSE95587"
  counts <- read_counts_csv(here::here("data", "processed", ds,
                                       "Edited_raw_counts.csv"))
  md <- as.data.frame(read_excel(here::here("data", "metadata", ds,
                                            "metadata.xlsx")))
  geo <- geo_characteristics(ds)
  j <- match(md[[3]], geo$`sample_id:ch1`)
  meta <- data.frame(Sample = md$IDENTIFIER,
                     Condition = ifelse(md$CONDITION == "AD", "AD", "Control"),
                     age = suppressWarnings(as.numeric(md$AGE)),
                     sex = factor(toupper(trimws(md$SEX))),
                     PMI = suppressWarnings(as.numeric(geo$`post_mortem_interval:ch1`[j])),
                     RIN = suppressWarnings(as.numeric(geo$`rin_score:ch1`[j])),
                     stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$Sample), , drop = FALSE]

  cov <- character()
  for (nm in c("age", "sex", "PMI", "RIN"))
    if (usable(ds, nm, meta[[nm]])) cov <- c(cov, nm)
  for (nm in c("batch", "APOE")) usable(ds, nm, NULL, TRUE)

  summaries[[ds]] <- run_unit(ds, NA_character_, counts, meta, cov)
}

# ------------------------------------------------------------ GSE159699 -----
# Gender, age at death, PMI (hr) and RIN all deposited and complete.
# Note: controls comprise young (HCY) and old (HCO) donors, so age is strongly
# associated with diagnosis in this dataset - see the findings summary.
{
  ds <- "GSE159699"
  counts <- read_counts_csv(here::here("data", "processed", ds,
                                       "EDITED_COUNTS.csv"))
  md <- as.data.frame(read_excel(here::here("data", "metadata", ds,
                                            "metadata.xlsx")))
  meta <- data.frame(Sample = md$identifier,
                     Condition = ifelse(grepl("^AD", md$identifier), "AD", "Control"),
                     age = suppressWarnings(as.numeric(md$`Age at death`)),
                     sex = factor(trimws(md$Gender)),
                     PMI = suppressWarnings(as.numeric(md$`PMI (hr)`)),
                     RIN = suppressWarnings(as.numeric(md$`RIN (RNA)`)),
                     stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$Sample), , drop = FALSE]

  cov <- character()
  for (nm in c("age", "sex", "PMI", "RIN"))
    if (usable(ds, nm, meta[[nm]])) cov <- c(cov, nm)
  for (nm in c("batch", "APOE")) usable(ds, nm, NULL, TRUE)

  summaries[[ds]] <- run_unit(ds, NA_character_, counts, meta, cov)
}

# ------------------------------------------------------------ GSE203206 -----
# SEX, AGE, RIN, APOE deposited. PMHOURS (PMI) is missing for 9/47 samples and
# so fails the completeness rule. No sequencing batch deposited.
{
  ds <- "GSE203206"
  counts_raw <- read.delim(
    here::here("data", "processed", ds,
               "GSE203206_Subramaniam.ADRC_brain.counts.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  colnames(counts_raw)[1] <- "Gene"
  counts_raw[is.na(counts_raw)] <- 0
  counts_raw$Gene <- make.unique(as.character(counts_raw$Gene))
  counts <- as.matrix(sapply(counts_raw[, -1, drop = FALSE], as.numeric))
  rownames(counts) <- counts_raw$Gene
  # Primary pipeline keeps AD (E*/L*) and control (C*) columns only.
  counts <- counts[, grepl("^(E|L|C)", colnames(counts)), drop = FALSE]

  md <- read.delim(here::here("data", "metadata", ds,
                              "GSE203206_Subramaniam.ADRC_brain.metadata.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  meta <- data.frame(
    Sample = md$SampleID,
    Condition = ifelse(grepl("^(E|L)", md$SampleID), "AD",
                       ifelse(grepl("^C", md$SampleID), "Control", NA)),
    age = suppressWarnings(as.numeric(md$AGE)),
    sex = factor(toupper(trimws(md$SEX))),
    PMI = suppressWarnings(as.numeric(md$PMHOURS)),
    RIN = suppressWarnings(as.numeric(md$RIN)),
    APOE = factor(trimws(as.character(md$APOE))),
    stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$Sample), , drop = FALSE]

  cov <- character()
  for (nm in c("age", "sex", "PMI", "RIN", "APOE"))
    if (usable(ds, nm, meta[[nm]])) cov <- c(cov, nm)
  usable(ds, "batch", NULL, TRUE)

  summaries[[ds]] <- run_unit(ds, NA_character_, counts, meta, cov)
}

# ------------------------------------------------------------ GSE278723 -----
# Five subregions from 28 retained donors. GEO deposits only tissue, subfield
# and group, but the source paper (Li C et al. 2024, Aging Dis, Supplementary
# Table 1) gives sex, age and PMI for all 50 donors; that table is parsed to
# data/metadata/GSE278723/GSE278723_supplementary_table1_donors.csv and joined
# here by donor ID. ABC and CAA scores are pathology staging, not covariates,
# so they are not adjusted for. RIN, batch and APOE are recorded nowhere.
#
# Condition labels come from the same Supplementary Table 1 join the corrected
# primary pipeline now uses (NC = 01-13, n = 13; AD = 14-28, n = 15), so Fit A
# continues to reproduce the deposited DE tables exactly.
{
  ds <- "GSE278723"
  donors <- read.csv(here::here("data", "metadata", ds,
                                "GSE278723_supplementary_table1_donors.csv"),
                     colClasses = c(Donor = "character"),
                     stringsAsFactors = FALSE)

  for (region in c("CA1", "CA2", "CA3", "CA4", "NX")) {
    counts <- read_counts_csv(here::here("data", "processed", ds,
                                         paste0("output_", region, ".csv")))
    donor <- sub("\\..*$", "", colnames(counts))
    j <- match(donor, donors$Donor)
    meta <- data.frame(
      Sample = colnames(counts),
      Donor = donor,
      Condition = ifelse(donors$Group[j] == "NC", "Control", "AD"),
      age = donors$Age[j],
      sex = factor(donors$Sex[j]),
      PMI = donors$PMI_h[j],
      stringsAsFactors = FALSE)

    cov <- character()
    for (nm in c("age", "sex", "PMI")) {
      ok <- usable(ds, nm, meta[[nm]])          # logged once per subregion
      if (ok) cov <- c(cov, nm)
    }
    for (nm in c("RIN", "batch", "APOE")) usable(ds, nm, NULL, TRUE)

    summaries[[paste0(ds, "_", region)]] <-
      run_unit(ds, region, counts, meta, cov)
  }
  # Collapse the repeated per-subregion availability rows to one per covariate.
  a <- do.call(rbind, avail)
  dup <- a$Dataset == ds & duplicated(paste(a$Dataset, a$Covariate))
  avail <- split(a[!dup, ], seq_len(sum(!dup)))
}

# ------------------------------------------------------------- outputs ------

model_summary <- do.call(rbind, summaries)
write.csv(model_summary, file.path(outdir, "model_summary.csv"), row.names = FALSE)

avail_long <- do.call(rbind, avail)
write.csv(avail_long, file.path(outdir, "covariate_availability_long.csv"),
          row.names = FALSE)
# Wide matrix: one row per dataset, one column per covariate (Table S1 feed).
avail_long$Cell <- ifelse(avail_long$Status == "Included", "Included",
                          ifelse(avail_long$Status == "Not deposited",
                                 "Not deposited",
                                 paste0("Recorded but not used (",
                                        avail_long$Reason, ")")))
avail_wide <- reshape(avail_long[, c("Dataset", "Covariate", "Cell")],
                      idvar = "Dataset", timevar = "Covariate",
                      direction = "wide")
names(avail_wide) <- sub("^Cell\\.", "", names(avail_wide))
avail_wide <- avail_wide[, c("Dataset",
                             intersect(c("age", "sex", "PMI", "RIN", "batch", "APOE"),
                                       names(avail_wide)))]
write.csv(avail_wide, file.path(outdir, "covariate_availability_matrix.csv"),
          row.names = FALSE)

# --- Fit A reproducibility check against the deposited primary DE tables -----
check_rows <- list()
deposited <- list(
  GSE67333  = here::here("results", "differential_expression", "GSE67333",  "full_DE_results.csv"),
  GSE95587  = here::here("results", "differential_expression", "GSE95587",  "full_DE_results.csv"),
  GSE159699 = here::here("results", "differential_expression", "GSE159699", "full_DE_results.csv"),
  GSE203206 = here::here("results", "differential_expression", "GSE203206", "full_DE_results.csv")
)
for (r in c("CA1", "CA2", "CA3", "CA4", "NX"))
  deposited[[paste0("GSE278723_", r)]] <-
    here::here("results", "differential_expression", "GSE278723", r, "full_DE_results.csv")

for (tag in names(deposited)) {
  new_path <- file.path(de_dir, paste0(tag, "_unadjusted.csv"))
  if (!file.exists(deposited[[tag]]) || !file.exists(new_path)) next
  old <- read.csv(deposited[[tag]], row.names = 1, check.names = FALSE)
  new <- read.csv(new_path, check.names = FALSE)
  common <- intersect(rownames(old), new$Gene)
  new <- new[match(common, new$Gene), ]
  check_rows[[tag]] <- data.frame(
    Unit = tag,
    n_genes_deposited = nrow(old), n_genes_refit = nrow(read.csv(new_path)),
    n_common = length(common),
    logFC_correlation = round(cor(old[common, "logFC"], new$logFC), 6),
    max_abs_logFC_diff = signif(max(abs(old[common, "logFC"] - new$logFC)), 4)
  )
}
fit_a_check <- do.call(rbind, check_rows)
write.csv(fit_a_check, file.path(outdir, "fitA_reproducibility_check.csv"),
          row.names = FALSE)
print(fit_a_check)

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_01.txt"))
message("\nWrote DE tables to ", de_dir)
