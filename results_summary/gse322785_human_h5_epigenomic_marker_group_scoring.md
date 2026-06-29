# GSE322785 Human H5 Epigenomic Marker-Group Scoring

Date built: 2026-06-26

## Scope

- Marker groups scored: 26.
- Feature-score rows: 83071.
- Module-score rows: 4568.
- Granule-versus-comparator contrast rows: 6018.

## Top Pooled Granule-Positive Contrasts

- `oligodendrocyte_candidate` vs `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support`: delta 1.207.
- `oligodendrocyte_candidate` vs `Gene Expression` / `ConstructionBalance` / `Tier 2 high-confidence wiring/synaptic executor`: delta 1.2.
- `oligodendrocyte_candidate` vs `Gene Expression` / `ConstructionBalance` / `shared_postmitotic_granule_maturation`: delta 1.155.
- `astrocyte_bergmann_candidate` vs `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support`: delta 1.125.
- `oligodendrocyte_candidate` vs `Gene Expression` / `ConstructionBalance` / `downstream_synaptic_excitability`: delta 1.116.
- `astrocyte_bergmann_candidate` vs `Peaks` / `FatePolarity` / `shared_postmitotic_granule_maturation`: delta 1.038.
- `astrocyte_bergmann_candidate` vs `Peaks` / `EpigenomicCompatibility` / `dentate_fate_wnt_prox1`: delta 1.038.
- `astrocyte_bergmann_candidate` vs `Peaks` / `EpigenomicCompatibility` / `medial_pallium_dentate_lineage`: delta 1.038.

## Interpretation

This is a provisional marker-group scoring layer. It tests whether selected target genes and nearby ATAC peak features can be summarized by candidate cell groups, but it should not be interpreted as source-author cell-type annotation or causal chromatin evidence.

## Outputs

- Feature scores: `Project/results/gse322785_human_h5_epigenomic_marker_group_feature_scores.tsv.gz`
- Module scores: `Project/results/gse322785_human_h5_epigenomic_marker_group_module_scores.tsv`
- Granule/comparator contrasts: `Project/results/gse322785_human_h5_epigenomic_marker_group_contrasts.tsv`
