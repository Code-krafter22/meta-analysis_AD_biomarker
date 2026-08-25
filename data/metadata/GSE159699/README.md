# GSE159699 metadata

## Source

`metadata.xlsx` contains selected sample characteristics retrieved from the
GSE159699 GEO record using the Bioconductor `GEOquery` package.

- GEO accession: [GSE159699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699)
- Samples: 30
- Groups: 12 Alzheimer's disease, 10 older controls, and 8 younger controls
- Tissue: lateral temporal lobe

## Selected fields

The workbook contains:

- `Sample ID`
- `identifier`: analysis identifier (`AD`, `HCO`, or `HCY` prefix)
- `Study Group`
- `Gender`
- `Age at death`
- `PMI (hr)`: postmortem interval
- `Braak`
- `Cerad`
- `RIN (RNA)`: RNA integrity number
- `neuron (%)`
- `Cause of Death`

The analysis identifiers were used to rename the 30 sample columns when
creating `data/processed/GSE159699/EDITED_COUNTS.csv`.

The workbook is stored unchanged at:

```text
data/metadata/GSE159699/metadata.xlsx
```
