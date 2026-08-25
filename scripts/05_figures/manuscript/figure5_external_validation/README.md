# Figure 5 — external validation

Figure 5 compares donor-grouped five-fold validation in controlled-access MSBB
and public GSE125583. Elastic-net settings are fixed (`l1_ratio=0.5`, `C=1.0`)
and scaling is fit inside each training fold. MSBB samples are grouped by
`individualID`, preventing donor leakage. Confidence intervals use donor-level
bootstrap resampling.

```bash
export MSBB_EXPR=/authorized/path/MSBB_expression_data_ad_ctrl.csv
export MSBB_META=/authorized/path/metadata_with_diagnosis.csv
python3 scripts/05_figures/manuscript/figure5_external_validation/run_figure5.py
```

Restricted MSBB sources are never copied into the repository. Derived
predictions contain no specimen or donor identifiers. Actual CV assignments are
stored in the `fold` column; no fold AUC is hardcoded. Final PNG, PDF and TIFF
files are written to `results/figures/main/`.

Validated pooled out-of-fold AUCs are 0.784 for MSBB and 0.861 for GSE125583.
