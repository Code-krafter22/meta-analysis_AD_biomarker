# GSE95587

## Dataset source

The GSE95587 dataset was obtained from the NCBI Gene Expression Omnibus (GEO):

- GEO accession: [GSE95587](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95587)
- Organism: *Homo sapiens*
- Experiment type: RNA-seq expression profiling
- Tissue: postmortem fusiform gyrus
- Study groups: autopsy-confirmed Alzheimer's disease (AD) and neurologically
  normal age-matched controls

The NCBI-generated RNA-seq count matrix was downloaded using the GEO
**Download RNA-seq counts** option:

- Compressed download: `GSE95587_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- Local decompressed filename: `GSE95587_raw_counts.tsv`

The corresponding NCBI Homo sapiens GRCh38.p13 annotation table was also
downloaded:

- Compressed download: `Human.GRCh38.p13.annot.tsv.gz`
- Decompressed file: `Human.GRCh38.p13.annot.tsv`
- Annotation release: NCBI *Homo sapiens* Annotation Release 109.20190905

Raw sequencing reads are available through the SRA link on the GEO record.

## Sample metadata

Selected sample metadata were obtained using GEOquery and saved in
`data/metadata/GSE95587/metadata.xlsx`. The workbook contains 117 samples and
the selected fields needed to describe the cohort, including GEO accession,
internal sample ID, sex, age, condition, and Braak stage.

The analysis identifiers were assigned from the metadata as follows:

- 84 Alzheimer's disease samples: `AD1`–`AD84`
- 33 control samples: `HC1`–`HC33`

All 117 metadata samples are represented in the count matrix.

## Creation of the edited count table

The downloaded count matrix contains 39,376 rows, with NCBI GeneID values in
the first column and raw integer counts for 117 GEO samples.

`Edited_raw_counts.csv` was created as follows:

1. Each `GeneID` was matched to the `Symbol` column in
   `Human.GRCh38.p13.annot.tsv`.
2. The first column was changed from `GeneID` to `Symbol`.
3. The 117 GSM count columns were renamed to the corresponding `AD` or `HC`
   analysis identifiers in the GEOquery-derived metadata.
4. Gene order and all raw count values were retained unchanged.
5. No filtering, normalization, transformation, or aggregation was performed
   while creating the edited table.

This transformation is implemented in:

```text
scripts/01_preprocessing/02_GSE95587_preprocessing.R
```

To recreate the processed table, place the decompressed count and annotation
files in `data/raw/GSE95587/` and run the script from the repository root. The
script accepts either the official NCBI filename or the shorter local count
filename documented above. Raw downloads do not need to be committed.

The resulting processed input is stored at:

```text
data/processed/GSE95587/Edited_raw_counts.csv
```

## Differential expression

The processed count table is used by:

```text
scripts/02_differential_expression/02_GSE95587_DEG.Rmd
```
