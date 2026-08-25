# External validation

This stage validates the 37-gene signature in public GSE125583 and the
controlled-access Mount Sinai Brain Bank (MSBB) cohort.

## Public GSE125583

The project includes the public GEO count matrix, selected metadata, processed
37-gene matrix, and the historical preprocessing notebooks. GEO accession:
GSE125583.

Run the public grouped five-fold validation from the project root:

```bash
VALIDATION_MODE=grouped_cv VALIDATION_COHORT=GSE125583 \
python3 scripts/04_external_validation/05_grouped_validation.py
```

GSE125583 contains 36 of the 37 signature genes; COL25A1 is unavailable. The
script uses the available genes and reports the number used.

## Controlled-access MSBB

MSBB expression and participant metadata are not distributed with this project.
Authorized users must obtain them through the relevant AD Knowledge Portal data
access process and prepare these two local files:

- `MSBB_cqn_selected_37_genes.csv`: normalized expression, samples by genes (or
  genes by samples; the loader detects orientation).
- `metadata_with_diagnosis.csv`: must contain `specimenID`, `individualID`, and
  `Diagnosis`; optional `exclude` is honored.

Point to authorized local copies without placing them in the repository:

```bash
MSBB_EXPR=/authorized/path/MSBB_cqn_selected_37_genes.csv \
MSBB_META=/authorized/path/metadata_with_diagnosis.csv \
VALIDATION_MODE=grouped_cv VALIDATION_COHORT=MSBB \
python3 scripts/04_external_validation/05_grouped_validation.py
```

For cross-cohort validation set `VALIDATION_MODE=external` and optionally
`EXTERNAL_DIRECTION=MSBB->GSE125583` or `GSE125583->MSBB`.

The validated MSBB input contains 1,043 samples from 263 donors. Controlled
expression, metadata, specimen IDs, and donor IDs remain outside this project.
Derived prediction files contain only binary labels, probabilities, predicted
classes, and fold numbers; they contain no controlled identifiers.

To regenerate both grouped-CV cohorts and Figure 5 in one command:

```bash
export MSBB_EXPR=/authorized/path/MSBB_expression_data_ad_ctrl.csv
export MSBB_META=/authorized/path/metadata_with_diagnosis.csv
python3 scripts/05_figures/manuscript/figure5_external_validation/run_figure5.py
```
