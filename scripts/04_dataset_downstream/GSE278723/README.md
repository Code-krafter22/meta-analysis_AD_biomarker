# GSE278723 regional downstream analysis

This workflow starts from the five regional DEG analyses in
`results/differential_expression/GSE278723`. The input up/down lists were defined
using the study DEG criterion `abs(logFC) > 0.58` and unadjusted `P.Value < 0.05`.

It is a within-dataset regional analysis and is intentionally separate from the
cross-dataset meta-analysis in `scripts/03_meta_analysis`.

Run from the project root:

```r
source("scripts/04_dataset_downstream/GSE278723/run_all.R")
```

Or run the numbered scripts individually in order. Outputs are written under
`results/dataset_downstream/GSE278723`.

The enrichment input lists retain the unadjusted DEG definition. GO and KEGG
tables also report BH-adjusted enrichment p-values; this multiple-testing
adjustment applies to pathway terms, not to the original DEG selection.
