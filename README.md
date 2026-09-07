# Alzheimer disease transcriptomic meta-analysis

Reproducible preprocessing, differential-expression, meta-analysis, external
validation, and manuscript-figure workflows for five discovery datasets plus
GSE125583 and MSBB validation.

## Run order

Run from the repository root.

Install R dependencies once if needed:

```r
source("scripts/00_install_R_dependencies.R")
```

```r
# Both prespecified meta-analysis branches
source("scripts/03_meta_analysis/00_run_all_meta_analyses.R")

# Screening-sensitivity analyses (fold-change thresholds + summary tables)
source("scripts/03_meta_analysis/03_meta_analysis_logfc_threshold.R")

# Figures 1–4
source("scripts/05_figures/manuscript/figure1_consensus_meta_analysis/run_figure1.R")
source("scripts/05_figures/manuscript/figure2_enrichment_ppi/run_figure2.R")
source("scripts/05_figures/manuscript/shared_braak/run_figures3_and_4.R")
```

Figure 5 requires authorized local MSBB data:

```bash
export MSBB_EXPR=/authorized/path/MSBB_expression_data_ad_ctrl.csv
export MSBB_META=/authorized/path/metadata_with_diagnosis.csv
python3 scripts/05_figures/manuscript/figure5_external_validation/run_figure5.py
```

Final figures are written to `results/figures/main/`.

## Analysis branches

- `filtered_unadjusted`: dataset DEGs selected with unadjusted P < 0.05; this
  historical branch produced the 37-gene signature.
- `unfiltered_all_genes`: reviewer-requested analysis using all tested genes
  without a P-value or logFC prefilter.
- `fdr_screened`: FDR-screened branch (BH FDR <= 0.05 with the fold-change
  criterion held at 0.58).
- `logfc_gt_*`: fold-change screening sensitivity, holding unadjusted
  P <= 0.05 fixed and varying only the magnitude threshold.

Screening-sensitivity results for all of the above are collected into one
Panel A table in `results/meta_analysis/logfc_threshold_sensitivity/`; see the
README there.

## Data access and software

Public GEO provenance is documented under `data/raw/` and `data/metadata/`.
Large public count matrices and source metadata are not committed; local file
placement and availability are described in `data/DATA_AVAILABILITY.md`.
MSBB is controlled-access and intentionally absent. Never copy MSBB expression,
metadata, specimen IDs, or donor IDs into this repository. Python dependencies
are pinned in `requirements.txt`; R workflows write `sessionInfo()` alongside
outputs. Historical scripts with absolute paths are provenance copies, not
current runners.

## Documented limitations

- Figure 1 starts from an archived manually assembled five-dataset log2-CPM
  workbook; subsequent ComBat correction, z-scoring, and plotting are scripted.
- Figure 2 uses frozen original KEGG tables by default for stable offline
  reproduction; GO enrichment is recomputed locally.
- The Figure 4 hippocampal comparison has four low- and four high-Braak
  participants, so adjusted estimates have limited precision.
