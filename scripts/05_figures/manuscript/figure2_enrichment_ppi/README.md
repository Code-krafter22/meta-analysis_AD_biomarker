# Figure 2 reproducible workflow

Figure 2 contains functional enrichment of the 37-gene signature and its STRING
protein-interaction network. Run from the project root:

```r
source("scripts/05_figures/manuscript/figure2_enrichment_ppi/run_figure2.R")
```

The signature is read from the 37-gene meta-effect table generated for Figure 1
and validated as 16 upregulated plus 21 downregulated genes. GO BP/MF/CC are
recomputed with the installed annotation packages. KEGG uses the frozen original
result tables by default so the figure is reproducible without a changing web
service. To deliberately refresh KEGG from its live service, run:

```r
Sys.setenv(USE_LIVE_KEGG = "true")
source("scripts/05_figures/manuscript/figure2_enrichment_ppi/run_figure2.R")
Sys.unsetenv("USE_LIVE_KEGG")
```

The exact downloaded STRING edge table is stored under
`data/figure_inputs/figure2/`. Tables and component panels are written under
`results/figures/intermediate/figure2/`; final PNG and PDF files are written
under `results/figures/main/`.
