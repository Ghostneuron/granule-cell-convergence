# GSE322785 Human H5 Cluster Validation

Date built: 2026-06-26

## Scope

- Basic-QC barcodes clustered: 24206.
- Donor-specific clusters: 36.
- Mean adjusted Rand index between provisional marker calls and clusters: 0.031.
- Mean normalized mutual information between provisional marker calls and clusters: 0.151.

## Focal Marker-Call Concentration

- `astrocyte_bergmann_candidate` in `H187`: 93 barcodes; dominant cluster `H187_C09` captures 88.2%.
- `cerebellar_granule_candidate` in `H187`: 104 barcodes; dominant cluster `H187_C11` captures 39.4%.
- `oligodendrocyte_candidate` in `H187`: 58 barcodes; dominant cluster `H187_C04` captures 89.7%.
- `purkinje_candidate` in `H187`: 1089 barcodes; dominant cluster `H187_C01` captures 25.3%.
- `astrocyte_bergmann_candidate` in `H628`: 126 barcodes; dominant cluster `H628_C11` captures 92.9%.
- `cerebellar_granule_candidate` in `H628`: 160 barcodes; dominant cluster `H628_C05` captures 30.6%.
- `oligodendrocyte_candidate` in `H628`: 74 barcodes; dominant cluster `H628_C11` captures 89.2%.
- `purkinje_candidate` in `H628`: 1570 barcodes; dominant cluster `H628_C05` captures 37.6%.
- `astrocyte_bergmann_candidate` in `H390`: 244 barcodes; dominant cluster `H390_C10` captures 73.4%.
- `cerebellar_granule_candidate` in `H390`: 512 barcodes; dominant cluster `H390_C09` captures 52.0%.
- `oligodendrocyte_candidate` in `H390`: 405 barcodes; dominant cluster `H390_C08` captures 64.7%.
- `purkinje_candidate` in `H390`: 1395 barcodes; dominant cluster `H390_C09` captures 48.0%.

## Interpretation

Selected-gene SVD/k-means clustering provides an internal validation layer for the provisional marker calls. This analysis can support prioritization of marker-group epigenomic contrasts, but it remains weaker than source-author taxonomy or full multimodal clustering.

## Outputs

- Barcode assignments: `Project/results/gse322785_human_h5_cluster_validation_barcode_assignments.tsv.gz`
- Cluster summary: `Project/results/gse322785_human_h5_cluster_validation_summary.tsv`
- Marker-call enrichment: `Project/results/gse322785_human_h5_cluster_validation_marker_call_enrichment.tsv`
- Marker support: `Project/results/gse322785_human_h5_cluster_validation_marker_support.tsv`
- Metrics: `Project/results/gse322785_human_h5_cluster_validation_metrics.tsv`
