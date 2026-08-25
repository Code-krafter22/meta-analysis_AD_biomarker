# GSE203206

## Dataset source

The GSE203206 dataset was downloaded from the NCBI Gene Expression Omnibus
(GEO):

- GEO accession: [GSE203206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203206)
- Organism: *Homo sapiens*
- Experiment type: RNA-seq expression profiling
- Tissue: postmortem occipital lobe, primary visual cortex (Brodmann area 17)

The following supplementary archives were downloaded from GEO and decompressed:

- `GSE203206_Subramaniam.ADRC_brain.counts.tsv.gz` →
  `GSE203206_Subramaniam.ADRC_brain.counts.tsv`
- `GSE203206_Subramaniam.ADRC_brain.metadata.tsv.gz` →
  `GSE203206_Subramaniam.ADRC_brain.metadata.tsv`

The GEO series page provides these processed supplementary files. Raw
sequencing data are available through the SRA link on the GEO record.

## Files used in this project

The downloaded count table is stored at:

```text
data/processed/GSE203206/GSE203206_Subramaniam.ADRC_brain.counts.tsv
```

The accompanying sample metadata are stored at:

```text
data/metadata/GSE203206/GSE203206_Subramaniam.ADRC_brain.metadata.tsv
```

The uncompressed source files do not need to be duplicated in `data/raw/`
because the count table is already a GEO-provided processed data product.

## Cohort and analysis groups

The metadata classify samples as `Alzheimers` or `Healthy`. The differential
expression script compares all available Alzheimer’s disease samples against
healthy controls. Sample names beginning with `E_` or `L_` are assigned to the
AD group, and sample names beginning with `C_` are assigned to the control
group.

In the supplied files:

- Metadata table: 47 samples (39 Alzheimer’s disease and 8 healthy controls)
- Count matrix: 46 samples (38 Alzheimer’s disease and 8 healthy controls)

`L_5017_OL_RNA_S87` appears in the metadata table but is absent from the count
matrix. Therefore, the differential expression analysis includes the 46
samples with available counts.

## Differential expression

The count table is used directly by:

```text
scripts/02_differential_expression/04_GSE203206_DEG.Rmd
```

No additional local preprocessing script was supplied for this dataset.
