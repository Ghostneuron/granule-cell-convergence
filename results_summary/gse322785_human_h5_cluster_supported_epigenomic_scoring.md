# GSE322785 Cluster-Supported Epigenomic Scoring

Date built: 2026-06-26

## Support Rule

Supported marker clusters required at least 20 barcodes of the marker call, cluster fraction at least 0.02, and at least 2.0-fold enrichment over the donor-level marker-call frequency. Ambiguous and low-information calls were excluded.

## Scope

- Supported cluster-marker rules: 19.
- Supported barcodes: 3550.
- Supported marker groups: 5.
- Module-score rows: 2630.
- Granule-versus-supported-comparator contrast rows: 2848.

## Supported Barcode Counts

- `purkinje_candidate`: 2235 supported barcodes.
- `oligodendrocyte_candidate`: 506 supported barcodes.
- `astrocyte_bergmann_candidate`: 378 supported barcodes.
- `cerebellar_granule_candidate`: 356 supported barcodes.
- `opc_candidate`: 75 supported barcodes.

## Top Pooled Granule-Positive Supported Contrasts

- `astrocyte_bergmann_candidate` / `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support`: delta 2.118.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `Tier 2 high-confidence wiring/synaptic executor`: delta 2.089.
- `oligodendrocyte_candidate` / `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support`: delta 2.008.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `downstream_synaptic_excitability`: delta 1.87.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `shared_postmitotic_granule_maturation`: delta 1.858.
- `astrocyte_bergmann_candidate` / `Gene Expression` / `ConstructionBalance` / `Tier 2 high-confidence wiring/synaptic executor`: delta 1.817.

## Interpretation

This stricter layer reduces false confidence by using only marker calls located in donor-specific clusters enriched for that same call. It is better suited for sensitivity analysis than the broader provisional marker-score layer, but it still does not replace source-author taxonomy or full multimodal clustering.

## Outputs

- Cluster support rules: `Project/results/gse322785_human_h5_cluster_supported_marker_rules.tsv`
- Supported barcodes: `Project/results/gse322785_human_h5_cluster_supported_marker_barcodes.tsv.gz`
- Supported feature scores: `Project/results/gse322785_human_h5_cluster_supported_epigenomic_feature_scores.tsv.gz`
- Supported module scores: `Project/results/gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv`
- Supported granule/comparator contrasts: `Project/results/gse322785_human_h5_cluster_supported_epigenomic_contrasts.tsv`
