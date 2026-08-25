# GSE67333 metadata

## Source

`metadata.xlsx` contains selected sample characteristics retrieved from the
GSE67333 GEO record using the Bioconductor `GEOquery` package.

- GEO accession: [GSE67333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67333)
- Samples: 8
- Groups: 4 late-onset Alzheimer's disease (LOAD) and 4 controls
- Tissue: hippocampus

## Selected fields

The workbook contains:

- `GSM`: GEO sample accession
- `identifier`: analysis identifier (`LOAD1`–`LOAD4` or `CTRL1`–`CTRL4`)
- `age`
- `gender`
- `APOE`
- `BRAAK `: Braak-stage information as reported in GEO
- `tissue`

The analysis identifiers were used to rename the eight GSM columns when
creating `data/processed/GSE67333/Edited_raw_counts.csv`.

The workbook is stored unchanged at:

```text
data/metadata/GSE67333/metadata.xlsx
```
