# GSE278723 metadata

## Source

`GSE278723_sample_information.txt` is the decompressed sample-information
supplement distributed on the GSE278723 GEO series record.

- GEO accession: [GSE278723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278723)
- Participants: 50
- Diagnostic groups: NC, AD, CAA, and AD+CAA

## Cohort used in this project

- Retained NC participants: IDs `01`–`13` (n = 13)
- Retained AD participants: IDs `14`–`28` (n = 15)
- Excluded CAA participants: IDs `29`–`34`
- Excluded AD+CAA participants: IDs `35`–`50`

Group assignments follow Supplementary Table 1 of the source publication
(Li C et al. 2024, *Aging and Disease*, doi:10.14336/AD.2024.10530), whose
Supplementary Figure 1 caption independently states NC n = 13, AD n = 15,
CAA n = 6. The donor table is parsed to
`GSE278723_supplementary_table1_donors.csv` and joined by donor ID rather than
by column position, so subregions with donors missing are handled correctly.

Only NC and AD participants were retained for the planned AD-versus-control
analysis without CAA pathology. The count matrix uses `NX` for the entorhinal
cortex; `NX` and `EC` refer to the same region in this project.

The metadata file is stored unchanged at:

```text
data/metadata/GSE278723/GSE278723_sample_information.txt
```
