# GSE95587 metadata

## Source

`metadata.xlsx` contains selected sample characteristics retrieved from the
GSE95587 GEO record using the Bioconductor `GEOquery` package.

- GEO accession: [GSE95587](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95587)
- Samples: 117
- Groups: 84 Alzheimer's disease and 33 controls
- Tissue: fusiform gyrus

## Selected fields

The workbook contains:

- `IDENTIFIER`: project analysis identifier (`AD1`–`AD84` or `HC1`–`HC33`)
- `GSM`: GEO sample accession
- `Unnamed: 2`: internal sample identifier reported by GEO
- `SEX`
- `AGE`
- `CONDITION`: `AD` or `CON`
- `BRAAK STAGING`

The analysis identifiers were used to rename the 117 GSM columns when creating
`data/processed/GSE95587/Edited_raw_counts.csv`.

The workbook is stored unchanged at:

```text
data/metadata/GSE95587/metadata.xlsx
```
