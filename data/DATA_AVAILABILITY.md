# Data availability

Large public count matrices and source metadata are intentionally excluded from
Git. Dataset-specific provenance and download instructions remain under
`data/raw/<accession>/README.md` and `data/metadata/<accession>/README.md`.

After downloading the public GEO files, place them in the paths expected by the
preprocessing and differential-expression scripts. The large local inputs are
ignored by `.gitignore`, so they will not be staged accidentally.

Small derived figure inputs are retained under `data/figure_inputs/` because
they are necessary to reproduce the manuscript figures and are not raw count
matrices. The GSE125583 37-gene processed matrix is also retained because it is
a compact, signature-level external-validation input.

MSBB expression and participant metadata are controlled-access and must never
be placed in this repository. Authorized users supply local paths through
`MSBB_EXPR` and `MSBB_META`.
