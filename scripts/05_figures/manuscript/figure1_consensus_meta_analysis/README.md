# Figure 1 workflow

The reproducible Figure 1 scope is the 37-gene heatmap plus the three
meta-analysis volcano panels. The historical Venn diagram is retained for
provenance but is not part of the reproducible runner.

The heatmap uses the manually assembled log2-CPM values for the 37 genes across
all five datasets, followed by ComBat correction and row z-score normalization.
That historical matrix is retained as a documented figure input.

The volcano panels are generated directly from the current filtered-unadjusted
meta-analysis CSV outputs; no manually assembled volcano workbook is required.

Run the complete heatmap-plus-volcano figure from the project root:

```r
source("scripts/05_figures/manuscript/figure1_consensus_meta_analysis/run_figure1.R")
```

The preparation step verifies the fixed 37-gene signature against effect-size,
Stouffer, and MetaVolcanoR outputs and writes method-specific CSVs under
`results/figures/intermediate/figure1`. Final volcano files are written under
`results/figures/main`.

The required heatmap workbook and phenotype/batch table are stored under
`data/figure_inputs/figure1/`. Final PNG and PDF files are written under
`results/figures/main/`. Historical scripts remain only for provenance.
