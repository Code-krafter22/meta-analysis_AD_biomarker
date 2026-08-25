# GSE278723

## Dataset source

The GSE278723 dataset was downloaded from the NCBI Gene Expression Omnibus
(GEO):

- GEO accession: [GSE278723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278723)
- Organism: *Homo sapiens*
- Experiment type: RNA-seq expression profiling
- Tissue: postmortem brain tissue from the hippocampal-entorhinal system

The following supplementary archives were downloaded from GEO and decompressed:

- `GSE278723_Counts.txt.gz` → `GSE278723_Counts.txt`
- `GSE278723_sample_information.txt.gz` → `GSE278723_sample_information.txt`

The count matrix is not distributed through this repository and can be
downloaded from the GEO accession page above. The small sample-information file
is retained at:

```text
data/metadata/GSE278723/GSE278723_sample_information.txt
```

## Cohort selection

The analysis was designed to compare Alzheimer's disease (AD) without
documented cerebral amyloid angiopathy against neurologically normal controls
(NC). Based on `GSE278723_sample_information.txt`, the following participants
were retained:

- NC: participant IDs `01`–`14`
- AD: participant IDs `15`–`28`

The following participants were excluded because the planned contrast was AD
versus NC without CAA pathology:

- CAA: participant IDs `29`–`34`
- AD+CAA: participant IDs `35`–`50`

This exclusion prevents the CAA and combined AD+CAA groups from being included
in the AD-versus-control comparison.

## Anatomical regions

Counts were separated into five region-specific tables:

- CA1
- CA2
- CA3
- CA4
- Entorhinal cortex (EC)

The supplied count matrix uses the suffix `NX` for entorhinal cortex samples.
Therefore, `NX` and `EC` refer to the same anatomical region in this project.

## Sample counts after selection

| Region | NC | AD |
|---|---:|---:|
| CA1 | 12 | 14 |
| CA2 | 14 | 14 |
| CA3 | 14 | 14 |
| CA4 | 14 | 14 |
| NX/EC | 14 | 14 |

CA1 data are unavailable for NC participants `03` and `14`, explaining the
smaller NC sample count for that region.

## Preprocessing

The preprocessing script is located at:

```text
scripts/01_preprocessing/05_GSE278723_preprocessing.R
```

It reads the raw count matrix from:

```text
data/raw/GSE278723/GSE278723_Counts.txt
```

The script retains samples from the NC and AD groups, cleans column names, and
separates the count matrix by anatomical region.

## Processed outputs

The region-specific count tables are written to `data/processed/GSE278723/`:

```text
output_CA1.csv
output_CA2.csv
output_CA3.csv
output_CA4.csv
output_NX.csv
```

These files are the inputs to:

```text
scripts/02_differential_expression/05_GSE278723_DEG.Rmd
```
