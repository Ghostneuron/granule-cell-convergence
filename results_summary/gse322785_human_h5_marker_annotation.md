# GSE322785 Human H5 Selected Matrix and Marker Annotation

Date built: 2026-06-26

## Scope

- Human H5 files processed: 3.
- Barcodes processed: 2095603.
- Basic-QC barcodes: 24206.
- High/medium-confidence provisional marker calls: 6047.
- Selected feature columns per sample: up to 3249.
- Selected matrix nonzero entries across samples: 6384881.

## Top Provisional Calls Among Basic-QC Barcodes

- `ambiguous_neuronal`: 18056 QC barcodes.
- `purkinje_candidate`: 4054 QC barcodes.
- `cerebellar_granule_candidate`: 776 QC barcodes.
- `oligodendrocyte_candidate`: 537 QC barcodes.
- `astrocyte_bergmann_candidate`: 463 QC barcodes.
- `opc_candidate`: 105 QC barcodes.
- `ambiguous_non_neuronal_or_niche`: 75 QC barcodes.
- `vascular_candidate`: 56 QC barcodes.

## Interpretation

These are marker-score calls, not source-author cell labels. They are sufficient to prioritize label-transfer and selected count extraction, but any manuscript-level chromatin claim still needs verified clustering or label transfer.

## Outputs

- Summary: `Project/results/gse322785_human_h5_selected_matrix_summary.tsv`
- Marker panel coverage: `Project/results/gse322785_human_h5_marker_panel_coverage.tsv`
- Marker cell-type summary: `Project/results/gse322785_human_h5_marker_celltype_summary.tsv`
- High-confidence barcode table: `Project/results/gse322785_human_h5_marker_high_confidence_barcodes.tsv.gz`
- Per-donor selected matrices and metadata: `Project/processed/gse322785_human_h5_selected`
