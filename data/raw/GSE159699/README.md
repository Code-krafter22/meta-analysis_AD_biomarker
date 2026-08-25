# GSE159699

## Dataset source

The GSE159699 dataset was obtained from the NCBI Gene Expression Omnibus
(GEO):

- GEO accession: [GSE159699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699)
- Organism: *Homo sapiens*
- Experiment type: RNA-seq expression profiling
- Tissue: postmortem lateral temporal lobe
- Study groups: Alzheimer's disease, older controls, and younger controls

The GEO supplementary count table used for this project was downloaded and
decompressed:

- Compressed file: `GSE159699_summary_count.star.txt.gz`
- Decompressed file: `GSE159699_summary_count.star.txt`

GEO also provides an NCBI-generated GRCh38.p13 RNA-seq count matrix and the
Human Annotation Release 109.20190905 table. However, comparison of the files
confirmed that the processed input used in this project was created from the
GEO supplementary STAR summary count table, not from the NCBI-generated count
matrix. The supplementary table already contains gene symbols in its
`refGene` column, so no separate gene-annotation merge was applied.

Raw sequencing reads are available through the SRA link on the GEO record.

## Sample metadata

Selected sample metadata were obtained using GEOquery and saved at
`data/metadata/GSE159699/metadata.xlsx`. The workbook contains 30 samples and
selected demographic, neuropathological, RNA-quality, and tissue-composition
variables.

The cohort contains:

- 12 Alzheimer's disease samples: `AD-20`–`AD-31`
- 10 older control samples: `HCO-10`–`HCO-19`
- 8 younger control samples: `HCY-2`–`HCY-9`

The current differential-expression script combines the older and younger
control samples into one control group, producing an AD-versus-control contrast
of 12 AD samples versus 18 controls.

## Creation of the edited count table

`GSE159699_summary_count.star.txt` contains 27,130 gene rows and 30 sample count
columns. `EDITED_COUNTS.csv` was created by simplifying the source sample names:

- Samples ending in `-AD` were renamed using the `AD-<Sample ID>` convention.
- Samples ending in `-Old` were renamed using the `HCO-<Sample ID>` convention.
- Samples ending in `-Young` were renamed using the `HCY-<Sample ID>` convention.

The `refGene` identifiers, gene order, all 30 samples, and all raw count values
were retained unchanged. No gene annotation, filtering, normalization,
transformation, or aggregation was performed while creating the edited table.

This transformation is implemented in:

```text
scripts/01_preprocessing/03_GSE159699_preprocessing.R
```

To recreate the processed table, place the decompressed supplementary count
file in `data/raw/GSE159699/` and run the script from the repository root. The
downloaded count file does not need to be committed to the repository.

The resulting processed input is stored at:

```text
data/processed/GSE159699/EDITED_COUNTS.csv
```

## Differential expression

The processed count table is used by:

```text
scripts/02_differential_expression/03_GSE159699_DEG.Rmd
```
