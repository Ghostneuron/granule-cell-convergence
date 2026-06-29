# GSE322785 Human H5 Epigenomic Targeting

Date built: 2026-06-26

## Purpose

This analysis maps the manuscript epigenomic target genes to gene-expression and nearby ATAC peak rows in the downloaded human adult cerebellar multiome H5 files.

## Scope

- Human H5 files analyzed: 3.
- Target genes represented in any human H5: 153/154.
- Per-sample target gene rows present: 459/462.
- Peak-target rows within gene body plus 100 kb: 9074.
- Selective manifest sample-feature rows: 9439 (2.08% of summed H5 matrix rows across analyzed files).

## Highest Peak-Count Genes

- `NRXN1`: 379 nearby peaks across human H5 files.
- `RBFOX3`: 292 nearby peaks across human H5 files.
- `ZBTB20`: 260 nearby peaks across human H5 files.
- `NFIA`: 259 nearby peaks across human H5 files.
- `TCF4`: 206 nearby peaks across human H5 files.
- `RTN1`: 163 nearby peaks across human H5 files.
- `ROBO2`: 151 nearby peaks across human H5 files.
- `NFIB`: 140 nearby peaks across human H5 files.

## Interpretation

The downloaded human cerebellar H5 files are usable for the epigenomic compatibility extension. They still lack manuscript-ready granule/Purkinje labels until cell-type annotation or label transfer is performed.

## Outputs

- H5 inventory: `Project/results/gse322785_human_h5_feature_inventory.tsv`
- Gene summary: `Project/results/gse322785_human_h5_epigenomic_gene_summary.tsv`
- Peak targets: `Project/results/gse322785_human_h5_epigenomic_peak_targets.tsv`
- Selective manifest: `Project/results/gse322785_human_h5_epigenomic_selective_manifest.tsv`
