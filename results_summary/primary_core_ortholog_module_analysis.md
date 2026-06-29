# Primary Core Ortholog-Aware Module Analysis

Date built: 2026-06-22

## Scope

This analysis freezes the strict 10-dataset primary core and tests the current module-level model: regional identity separation with structural-program convergence.

The ortholog-aware layer is currently marker-panel level, using curated same-root mouse-human symbols from the refined marker panel. It is not yet genome-wide ortholog-aware differential expression.

## Primary Core Composition

- `cerebellum`: 3 datasets.
- `human_dentate_hippocampus`: 3 datasets.
- `mouse_dentate`: 4 datasets.

## Eligible Module Units

- `cerebellum` / `broad_neuronal_structural_warning`: 6 units, 15415 cells/spots.
- `cerebellum` / `cerebellar_candidate`: 6 units, 36679 cells/spots.
- `cerebellum` / `other_or_ambiguous`: 6 units, 54983 cells/spots.
- `human_dentate_hippocampus` / `broad_neuronal_structural_warning`: 39 units, 32081 cells/spots.
- `human_dentate_hippocampus` / `dentate_candidate`: 167 units, 121644 cells/spots.
- `human_dentate_hippocampus` / `non_dentate_background`: 8 units, 419 cells/spots.
- `human_dentate_hippocampus` / `other_or_ambiguous`: 39 units, 29474 cells/spots.
- `mouse_dentate` / `dentate_candidate`: 29 units, 11015 cells/spots.
- `mouse_dentate` / `dentate_low_support`: 30 units, 1219 cells/spots.
- `mouse_dentate` / `non_dentate_background`: 20 units, 14352 cells/spots.
- `mouse_dentate` / `other_or_ambiguous`: 11 units, 395 cells/spots.

## Main Tests

- `human_dentate_branch_structural_rank_above_median` / `structural_rank`: median delta 0.0408; n=167 vs ; BH-adjusted p=0.00931.
- `primary_core_cerebellar_structural_rank_above_median` / `structural_rank`: median delta 0.2051; n=6 vs ; BH-adjusted p=0.0156.
- `primary_core_dentate_structural_rank_above_median` / `structural_rank`: median delta 0.0545; n=196 vs ; BH-adjusted p=0.000117.
- `primary_core_dentate_vs_cerebellar_identity_separation` / `identity_rank_contrast`: median delta 0.5356; n=196 vs 6; BH-adjusted p=2.95e-08.
- `primary_core_dentate_vs_cerebellar_structural_convergence` / `structural_rank`: median delta -0.1506; n=196 vs 6; BH-adjusted p=0.0105.

## Interpretation

- Leave-one-dataset-out identity separation is stable at BH-adjusted p<0.05 with positive dentate-minus-cerebellar delta.
- The strongest model remains: dentate and cerebellar granule candidates are identity-distinct, while structural/morphogenesis modules provide a shared elevated executor axis.
- The structural axis should not be framed as equal magnitude: in this strict-core rank layer, cerebellar candidates are higher than dentate candidates, while both candidate classes remain above the within-sample structural median.
- Human DG-enriched sources are most useful for anchoring human dentate state and structural programs; they can compress within-sample dentate-versus-cerebellar contrast.
- The next stricter step is genome-wide ortholog-aware pseudobulk or mixed-effect differential expression within the 10-dataset core.

## Outputs

- Ortholog marker map: `Project/results/primary_core_marker_panel_ortholog_map.tsv`
- Primary core module units: `Project/results/primary_core_integrated_module_units.tsv`
- Module statistics: `Project/results/primary_core_ortholog_module_statistics.tsv`
- Leave-one-dataset-out checks: `Project/results/primary_core_ortholog_module_leave_one_dataset_out.tsv`
- Primary core module plot: `Project/results/primary_core_identity_structural_modules.png`
