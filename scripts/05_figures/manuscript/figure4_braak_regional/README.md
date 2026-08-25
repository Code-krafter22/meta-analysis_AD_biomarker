# Figure 4 — regional Braak waterfall

Figure 4 uses the same 201 staged participants and shared ComBat preparation as
Figure 3. It shows signed fold change for Braak High (>3) versus Low (<=3) in
four regions. Asterisks are used only on the waterfall bars and represent
within-region BH-adjusted ANCOVA P-values; models adjust for age and sex.

Run Figure 4 alone:

```r
source("scripts/05_figures/manuscript/figure4_braak_regional/run_figure4.R")
```

The hippocampus contains only eight participants (four low and four high), so
its adjusted estimates have limited precision and should be interpreted
cautiously. This limitation is preserved rather than hidden.

The final PNG and PDF are written to `results/figures/main/`.
