# Meta-analysis workflow

There are **two prespecified analysis branches**, not two accidental copies:

| Branch | Input rule | Purpose | Output directory |
|---|---|---|---|
| 1. Filtered/unadjusted | Dataset DEGs selected using unadjusted `P.Value < 0.05` | Primary historical analysis that produced the 37-gene signature | `results/meta_analysis/filtered_unadjusted/` |
| 2. Unfiltered/all genes | All tested genes; no P-value or logFC prefilter | Reviewer-requested non-filtered sensitivity analysis | `results/meta_analysis/unfiltered_all_genes/` |
| 3. Fold-change threshold | Unadjusted `P <= 0.05` held fixed; `\|log2FC\|` varied | Reviewer-requested fold-change sensitivity analysis | `results/meta_analysis/logfc_gt_<cut>/` |

## Recommended command

Run both branches, in manuscript order, from the project root:

```r
source("scripts/03_meta_analysis/00_run_all_meta_analyses.R")
```

To run only one branch:

```r
source("scripts/03_meta_analysis/01_meta_analysis_filtered_unadjusted.R")
source("scripts/03_meta_analysis/02_meta_analysis_unfiltered.R")
```

## Branch 3 — fold-change threshold sensitivity

```r
source("scripts/03_meta_analysis/03_meta_analysis_logfc_threshold.R")
```

Runs the engine at `|log2FC| > 0.58, 1.00, 0.25` and `0.00`, all with unadjusted
`P <= 0.05` held fixed, driven from the unfiltered `full_DE_results.csv` tables.
`0.58` is a harness control (it must reproduce the primary analysis) and `0.00`
is a P-only reference used to attribute dropouts to the magnitude gate; only
`1.00` and `0.25` are reported specifications. Summary tables land in
`results/meta_analysis/logfc_threshold_sensitivity/`, which also collects the
two prespecified branches and the FDR-screened branch into one Panel A.

The branch sets `screen_lfc_cut` / `screen_p_cut` before sourcing the engine.
For branches 1 and 2 the config sets both to `NULL` and `stages/01` applies no
extra screen, so those branches are unchanged.

## Shared implementation

Both public branch runners call the internal `run_meta_analysis.R` engine. The
engine runs the same stages for either branch:

1. `00_config_and_inputs.R` — validates packages and selects branch inputs.
2. `stages/01_harmonize_inputs.R` — standardizes identifiers/statistics.
3. `stages/02_collapse_regions_and_qc.R` — evaluates and collapses GSE278723 regions.
4. `stages/03_stouffer_meta_analysis.R` — direction-aware Stouffer analysis.
5. `stages/04_effect_size_meta_analysis.R` — random-effects analysis, leave-one-out diagnostics, consensus calls, and heatmap preparation.
6. `stages/05_metavolcano_and_diagnostics.R` — MetaVolcanoR-compatible analysis and direction checks.

Only input selection changes between branches; the statistical implementation
is shared. The original monolithic scripts under `legacy/` are audit records,
not current runners.

Every run writes `sessionInfo.txt`. For an isolated test run:

```r
Sys.setenv(META_ANALYSIS_OUTPUT_DIR = tempfile("meta_test_"))
source("scripts/03_meta_analysis/01_meta_analysis_filtered_unadjusted.R")
Sys.unsetenv("META_ANALYSIS_OUTPUT_DIR")
```
