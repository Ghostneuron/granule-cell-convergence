# GSE322785 Epigenomic Robust Contrast Summary

Date built: 2026-06-26

## Rule

Robust-positive contrasts require both broad and cluster-supported pooled granule-minus-comparator mean delta scores >= 0.1. Strong robust-positive contrasts require both layers >= 0.25.

## Concordance Classes

- `cluster_supported_only_positive`: 122 contrasts.
- `discordant_strict_positive_broad_negative`: 14 contrasts.
- `missing_layer`: 372 contrasts.
- `robust_negative`: 7 contrasts.
- `robust_positive`: 166 contrasts.
- `robust_positive_strong`: 97 contrasts.
- `weak_or_neutral`: 152 contrasts.

## Robust-Positive Counts by Comparator

- `astrocyte_bergmann_candidate`: 142 robust-positive contrasts.
- `oligodendrocyte_candidate`: 121 robust-positive contrasts.

## Robust-Positive Counts by Term

- `Peaks` / `EpigenomicCompatibility`: 81 robust-positive contrasts.
- `Peaks` / `ConstructionBalance`: 61 robust-positive contrasts.
- `Peaks` / `FatePolarity`: 49 robust-positive contrasts.
- `Peaks` / `NicheSignal`: 24 robust-positive contrasts.
- `Gene Expression` / `EpigenomicCompatibility`: 21 robust-positive contrasts.
- `Gene Expression` / `ConstructionBalance`: 16 robust-positive contrasts.
- `Gene Expression` / `FatePolarity`: 7 robust-positive contrasts.
- `Gene Expression` / `NicheSignal`: 4 robust-positive contrasts.

## Top Robust-Positive Contrasts

- `oligodendrocyte_candidate` / `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support` / `gene_expression_or_marker`: broad 1.21, strict 2.01.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `Tier 2 high-confidence wiring/synaptic executor` / `gene_expression_or_marker`: broad 1.2, strict 2.09.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `shared_postmitotic_granule_maturation` / `gene_expression_or_marker`: broad 1.15, strict 1.86.
- `astrocyte_bergmann_candidate` / `Gene Expression` / `EpigenomicCompatibility` / `Tier 3 broad both-screen mechanism support` / `gene_expression_or_marker`: broad 1.12, strict 2.12.
- `oligodendrocyte_candidate` / `Gene Expression` / `ConstructionBalance` / `downstream_synaptic_excitability` / `gene_expression_or_marker`: broad 1.12, strict 1.87.
- `astrocyte_bergmann_candidate` / `Peaks` / `EpigenomicCompatibility` / `dentate_fate_wnt_prox1` / `gene_body_overlap`: broad 1.04, strict 1.7.
- `astrocyte_bergmann_candidate` / `Peaks` / `EpigenomicCompatibility` / `medial_pallium_dentate_lineage` / `gene_body_overlap`: broad 1.04, strict 1.7.
- `astrocyte_bergmann_candidate` / `Peaks` / `FatePolarity` / `shared_postmitotic_granule_maturation` / `gene_body_overlap`: broad 1.04, strict 1.7.

## Interpretation

This table identifies GSE322785 epigenomic signals that are stable to the stricter cluster-supported sensitivity filter. These robust contrasts are stronger candidates for discussion than broad-only or strict-only effects, while still remaining provisional until source taxonomy or full multimodal clustering is available.

## Outputs

- Full comparison: `Project/results/gse322785_epigenomic_broad_vs_cluster_supported_contrasts.tsv`
- Robust-positive contrasts: `Project/results/gse322785_epigenomic_robust_positive_contrasts.tsv`
- Summary: `Project/results/gse322785_epigenomic_robust_contrast_summary.tsv`
