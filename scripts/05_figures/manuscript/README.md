# Manuscript figure provenance

This directory records the scripts traced from the manuscript figure documents
to the historical `alzheimer_meta/ad_meta/FIGURES` directory. The copies are
preserved as source evidence. Several still contain historical absolute Windows
paths and therefore must not be treated as fully portable runners until their
input data are added to this project and the paths are converted.

## Figure mapping

| Manuscript item | Final historical figure file | Source scripts copied here |
|---|---|---|
| Figure 1: consensus meta-analysis | `Figure-1_Meta-analysis of differentially expressed genes in Alzheimer's disease across five transcriptomic datasets.jpeg` (also `Figure_37_genes_new1.*`) | `figure1_consensus_meta_analysis/01_consensus_venn.R`, `02_consensus_heatmap.R`, `03_meta_volcanoes.R`, `04_assemble_figure1.R` |
| Figure 2: enrichment and PPI | `Figure-2_Functional enrichment of the dysregulated consensus genes across different AD brain regions..jpeg` (also `Figure_PPI_GO_KEGG.*`) | `figure2_enrichment_ppi/01_consensus_enrichment.R`, `02_ppi_network.R`; no final assembly script was found in the historical folder |
| Figure 3: overall Braak association | `Figure-3_Association of consensus gene expression with Braak stage across all brain regions..jpeg` (also `Figure_Combined_AB.*`) | shared `shared_braak/01_braak_analysis.R`, plus `figure3_braak_overall/02_braak_forest_plot.R` and `03_assemble_figure3.R` |
| Figure 4: region-stratified Braak differences | `Figre-4_Region-stratified gene expression differences between High and Low Braak groups.jpeg` (also `Final_Waterfall_Braak_by_region_FC_with_sig.*`) | shared `shared_braak/01_braak_analysis.R`; the regional waterfall section creates this panel |
| Figure 5: external validation | `Figure-5.pdf` (also `figure_external_validation.*`) | `figure5_external_validation/01_grouped_cv_figure.py`; the older alternative is retained as `legacy_external_validation_figure.py` |
| Supplementary Figure S2: Spearman–Braak correlation | `Final_Spearman_Braak_zscore.jpeg` | shared `shared_braak/01_braak_analysis.R`, section 10b. The older correlation-only implementation is retained under `supplementary_figure_s2/legacy_correlation_plot.R`. |

## Input gaps before these can be rerun

- Figure 1 requires the historical consensus expression workbook,
  phenotype/batch table, and volcano workbook, or equivalent files regenerated
  from the current filtered meta-analysis.
- Figure 2 requires the finalized up/down consensus gene lists and STRING edge
  table. The STRING TSV inputs exist historically but have not yet been copied
  into the reproducible data hierarchy.
- Figures 3, 4, and S2 require `BRAAK_DATA_4studies.xlsx` or a reproducibly
  generated equivalent.
- Figure 5 requires the saved prediction, coefficient, metric, and fold-level
  validation tables referenced by the Python script.

The next reproducibility pass should connect these scripts to current outputs,
remove embedded 37-gene vectors where possible, and write final artifacts to
`results/figures/`.

## Figure numbering note

The supplied Supplementary Figure legend assigns S1 panels as A=GSE203206,
B=GSE159699, C=GSE67333, D=GSE95587, and E=GSE278723. Any regenerated S1 must
either retain that ordering or update the legend at the same time.
