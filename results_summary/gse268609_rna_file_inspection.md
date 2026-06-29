# GSE268609 RNA File Inspection

Date inspected: 2026-06-21

## Companion Files

- Features: 177931 total (36601 gene-expression rows, 141330 peak rows).
- Barcodes: 366175 total across 39 sample suffixes.
- Selected human-core genes present as gene-expression rows: 2083 / 2169.

## Matrix Header

- Header: `%%MatrixMarket matrix coordinate integer general`
- Dimensions: 177931 rows by 366175 columns.
- Non-zero entries: 1451618555.
- Matrix rows match feature count: True.
- Matrix columns match barcode count: True.
- Header read error: ``

## Extraction Decision

- Treat this as a primary candidate only through its RNA gene-expression rows for transcriptomic comparison.
- ATAC peak rows are valuable for later regulatory analysis, but they should not enter the first cross-dataset morphology-gene module projection.
- Full selected-gene extraction should be attempted only after confirming the non-zero count and available memory/disk, because the combined multiome matrix is much larger than the previous adult dentate anchor.

## Outputs

- Feature type summary: `Project/results/gse268609_feature_type_summary.tsv`
- Barcode suffix summary: `Project/results/gse268609_barcode_sample_suffix_summary.tsv`
- Selected gene presence: `Project/results/gse268609_selected_gene_presence.tsv`
- Matrix header: `Project/results/gse268609_matrix_header.tsv`
