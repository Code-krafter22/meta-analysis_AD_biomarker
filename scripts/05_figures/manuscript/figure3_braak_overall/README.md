# Figure 3 — overall Braak analysis

Figure 3 uses all 201 participants with recorded Braak stage across four
datasets. The source workbook uses numeric stages 0-6; the preparation code also
accepts Roman labels I-VI. Using AD cases alone is not possible for the intended
high-versus-low comparison because three regions contain no AD cases at stages
0-3.

- Panel A: signed fold change for Braak High (>3) versus Low (<=3). Asterisks
  represent BH-adjusted P-values from ANCOVA adjusted for region, age and sex.
- Panel B: odds ratios from logistic regression adjusted for region, age and
  sex. It intentionally contains no significance asterisks.

Run Figure 3 alone:

```r
source("scripts/05_figures/manuscript/figure3_braak_overall/run_figure3.R")
```

The final PNG and PDF are written to `results/figures/main/`.
