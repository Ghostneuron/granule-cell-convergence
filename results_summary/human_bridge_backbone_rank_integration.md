# Human Bridge And Backbone Rank Integration

Date built: 2026-06-21

## Scope

This analysis integrates the existing refined mouse/human cerebellar-dentate backbone with the constructed human dentate/hippocampal bridge objects using within-sample rank metrics.

Raw module medians are retained, but the main integrated comparison uses rank metrics because the older backbone and newer human selected-gene bridge objects were not produced from one shared whole-transcriptome matrix.

Minimum unit size for statistics and plotting: 20 cells/spots.

## Eligible Unit Counts

- `backbone_refined` / `broad_neuronal_structural_warning`: 8 units, 17122 cells/spots.
- `backbone_refined` / `cerebellar_candidate`: 8 units, 36982 cells/spots.
- `backbone_refined` / `dentate_candidate`: 30 units, 11053 cells/spots.
- `backbone_refined` / `dentate_low_support`: 31 units, 1275 cells/spots.
- `backbone_refined` / `non_dentate_background`: 20 units, 14352 cells/spots.
- `backbone_refined` / `other_or_ambiguous`: 23 units, 61580 cells/spots.
- `gse268609_hippocampus_rna` / `broad_neuronal_structural_warning`: 39 units, 32081 cells/spots.
- `gse268609_hippocampus_rna` / `dentate_candidate`: 113 units, 32832 cells/spots.
- `gse268609_hippocampus_rna` / `other_or_ambiguous`: 39 units, 29474 cells/spots.
- `gse325391_adult_dg` / `dentate_candidate`: 48 units, 56808 cells/spots.
- `gse325391_adult_dg` / `non_dentate_background`: 8 units, 419 cells/spots.
- `human_core_tuned` / `broad_neuronal_structural_warning`: 25 units, 10520 cells/spots.
- `human_core_tuned` / `dentate_candidate`: 57 units, 67412 cells/spots.
- `human_core_tuned` / `dentate_low_support`: 19 units, 3074 cells/spots.
- `human_core_tuned` / `non_dentate_background`: 26 units, 37763 cells/spots.
- `human_core_tuned` / `other_or_ambiguous`: 52 units, 36805 cells/spots.

## Main Tests

- `gse325391_adult_dg_dentate_candidate_vs_gse268609_hippocampus_rna_dentate_candidate` / `identity_rank_contrast`: delta -0.3153 (n=48 vs 113, BH-adjusted p=5.93e-22).
- `dentate_candidate_vs_non_dentate_background` / `identity_rank_contrast`: delta 0.3423 (n=248 vs 54, BH-adjusted p=9.93e-22).
- `dentate_candidate_sign_test_greater_0.0` / `identity_rank_contrast`: delta 0.2553 (n=248 vs , BH-adjusted p=2.17e-17).
- `human_core_tuned_dentate_candidate_vs_gse268609_hippocampus_rna_dentate_candidate` / `structural_rank`: delta 0.2331 (n=57 vs 113, BH-adjusted p=2.68e-14).
- `dentate_candidate_sign_test_greater_0.5` / `structural_rank`: delta 0.0866 (n=248 vs , BH-adjusted p=1.14e-11).
- `dentate_candidate_vs_cerebellar_candidate` / `identity_rank_contrast`: delta 0.6472 (n=248 vs 8, BH-adjusted p=1.57e-10).
- `cerebellar_candidate_sign_test_less_0.0` / `identity_rank_contrast`: delta -0.3919 (n=8 vs , BH-adjusted p=0.00502).
- `cerebellar_candidate_sign_test_greater_0.5` / `structural_rank`: delta 0.1831 (n=8 vs , BH-adjusted p=0.0396).
- `dentate_candidate_vs_cerebellar_candidate` / `structural_rank`: delta -0.0965 (n=248 vs 8, BH-adjusted p=0.435).

## Source-Layer Median Signals

- `backbone_refined` / `cerebellar_candidate`: identity contrast -0.3919, structural rank 0.6831 (8 units; 36982 cells/spots).
- `backbone_refined` / `dentate_candidate`: identity contrast 0.3773, structural rank 0.6861 (30 units; 11053 cells/spots).
- `backbone_refined` / `non_dentate_background`: identity contrast -0.2433, structural rank 0.2527 (20 units; 14352 cells/spots).
- `gse268609_hippocampus_rna` / `dentate_candidate`: identity contrast 0.3074, structural rank 0.5487 (113 units; 32832 cells/spots).
- `gse325391_adult_dg` / `dentate_candidate`: identity contrast -0.0079, structural rank 0.3991 (48 units; 56808 cells/spots).
- `gse325391_adult_dg` / `non_dentate_background`: identity contrast -0.1892, structural rank 0.0744 (8 units; 419 cells/spots).
- `human_core_tuned` / `dentate_candidate`: identity contrast -0.0039, structural rank 0.7818 (57 units; 67412 cells/spots).
- `human_core_tuned` / `non_dentate_background`: identity contrast -0.0781, structural rank 0.2866 (26 units; 37763 cells/spots).

## Interpretation

- The integrated rank layer is a sanity-check bridge, not yet a replacement for a single harmonized object.
- A robust project signal would be: dentate candidate units show positive dentate-minus-cerebellar identity rank, cerebellar candidate units show negative identity rank contrast, and both groups show above-median structural-program rank.
- `GSE325391` and `GSE186538` provide the strongest source-aware human dentate anchors; their enriched DG composition makes them better anchors for structural/neurogenic state than for within-sample dentate-versus-cerebellar identity contrast.
- `GSE268609` broadens the human aging/AD context and shows useful projected dentate signal, but it remains projection-labeled until source taxonomy from the full Seurat object is added.

## Outputs

- Integrated units: `Project/results/human_bridge_backbone_rank_units.tsv`
- Integrated statistics: `Project/results/human_bridge_backbone_rank_statistics.tsv`
- Source-layer summary: `Project/results/human_bridge_backbone_rank_source_summary.tsv`
- Rank-unit plot: `Project/results/human_bridge_backbone_rank_units.png`
