# GSE67333

## Dataset source

The GSE67333 dataset was obtained from the NCBI Gene Expression Omnibus (GEO):

- GEO accession: [GSE67333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67333)
- Organism: *Homo sapiens*
- Experiment type: RNA-seq expression profiling
- Tissue: hippocampus
- Study groups: four late-onset Alzheimer's disease (LOAD) samples and four
  age-matched controls

The NCBI-generated RNA-seq count matrix was downloaded using the GEO
**Download RNA-seq counts** option:

- Compressed download: `GSE67333_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- Decompressed file: `GSE67333_raw_counts_GRCh38.p13_NCBI.tsv`

The corresponding NCBI Homo sapiens GRCh38.p13 annotation table was also
downloaded:

- Compressed download: `Human.GRCh38.p13.annot.tsv.gz`
- Decompressed file: `Human.GRCh38.p13.annot.tsv`

Raw sequencing reads are available through the SRA link on the GEO record.

## Sample metadata

Selected sample metadata were obtained using GEOquery and saved at
`data/metadata/GSE67333/metadata.xlsx`. The workbook contains the GSM accession,
analysis identifier, age, gender, APOE genotype, Braak stage, and tissue for
each sample.

The sample identifiers used in the analysis are:

| GEO sample | Analysis identifier | Group |
|---|---|---|
| GSM1644996 | LOAD1 | LOAD |
| GSM1644997 | LOAD2 | LOAD |
| GSM1644998 | LOAD3 | LOAD |
| GSM1644999 | LOAD4 | LOAD |
| GSM1645000 | CTRL1 | Control |
| GSM1645001 | CTRL2 | Control |
| GSM1645002 | CTRL3 | Control |
| GSM1645003 | CTRL4 | Control |

All samples are hippocampal tissue samples.

## Creation of the edited count table

The downloaded count matrix contains 39,376 rows, with NCBI GeneID values in
the first column and raw integer counts for the eight GEO samples.

`Edited_raw_counts.csv` was created as follows:

1. Each `GeneID` was matched to the `Symbol` column in
   `Human.GRCh38.p13.annot.tsv`.
2. The first column was changed from `GeneID` to `Symbol`.
3. The eight GSM count columns were renamed to `LOAD1`–`LOAD4` and
   `CTRL1`–`CTRL4` according to the GEOquery-derived metadata.
4. Gene order and all raw count values were retained unchanged.
5. No filtering, normalization, transformation, or aggregation was performed
   while creating the edited table.

This transformation is implemented in:

```text
scripts/01_preprocessing/01_GSE67333_preprocessing.R
```

To recreate the processed table, place the decompressed count and annotation
files in `data/raw/GSE67333/` and run the script from the repository root. The
raw downloads do not need to be committed to the repository.

The resulting processed input is stored at:

```text
data/processed/GSE67333/Edited_raw_counts.csv
```

## Differential expression

The processed count table is used by:

```text
scripts/02_differential_expression/01_GSE67333_DEG.Rmd
```
