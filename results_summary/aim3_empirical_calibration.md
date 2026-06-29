# Aim 3 Empirical Sparse-Coding Calibration

Date built: 2026-06-24

## Purpose

This file adds an empirical fitting/calibration layer to the Aim 3 sparse expansion-coding model. It constrains the existing simulation grid with NeuroMorpho morphology and DANDI 000003 dentate granule activity/spatial-coding summaries.

## Empirical Targets

- `input_sampling_stems_cerebellar`: 4 (NeuroMorpho; direct lower-bound proxy for sparse input sampling/claws).
- `input_sampling_stems_dentate`: 1.5 (NeuroMorpho; direct lower-bound proxy for dentate primary dendritic input sampling).
- `branch_complexity_shared`: 20.5 (NeuroMorpho; morphology-complexity prior, not identical to model input_degree).
- `active_spatial_bin_fraction_dentate`: 0.5489 (DANDI 000003; upper-bound empirical activity support; not the same as instantaneous output active fraction).
- `spatial_information_dentate`: 0.78 (DANDI 000003; supports nontrivial information retention and spatial selectivity).
- `pv_far_minus_near_dentate`: 0.1229 (DANDI 000003; weak direct pattern-separation support).

## Calibration Objective

The empirical calibration score combines useful pattern separation, resource-adjusted useful score, information retention, distance-structure preservation, low collapse, sparse morphology/input-sampling plausibility, and a light DANDI activity proxy. A second resource-constrained calibration score upweights resource-adjusted performance and morphology/input-sampling plausibility. DANDI active spatial-bin fraction is treated as an upper-bound activity proxy, not as a direct instantaneous model output-active fraction.

## Best Grid Point

- Architecture: `intermediate`.
- Expansion ratio: 0.5.
- Input degree: 2.
- Output active fraction parameter: 0.05; observed 0.0625.
- Useful score: 0.094.
- Resource-adjusted useful score: 1.507.
- Empirical calibration score: 4.103.

## Best Resource-Constrained Grid Point

- Architecture: `intermediate`.
- Expansion ratio: 0.5.
- Input degree: 2.
- Output active fraction parameter: 0.05; observed 0.0625.
- Useful score: 0.094.
- Resource-adjusted useful score: 1.507.
- Resource-constrained calibration score: 5.268.

## Named Architecture Ranking

- Rank 1 `dense_expansion`: median calibration 0.712, resource-constrained rank 3 (median 0.249), useful 0.265, resource-adjusted 0.010, input degree 24.0, active fraction 0.150.
- Rank 2 `granule_like_sparse_expansion`: median calibration 0.584, resource-constrained rank 2 (median 0.456), useful 0.192, resource-adjusted 0.092, input degree 4.0, active fraction 0.051.
- Rank 3 `intermediate`: median calibration 0.538, resource-constrained rank 1 (median 0.568), useful 0.110, resource-adjusted 0.109, input degree 8.0, active fraction 0.051.
- Rank 4 `integrator_like_low_expansion`: median calibration 0.115, resource-constrained rank 4 (median -0.040), useful 0.121, resource-adjusted 0.053, input degree 24.0, active fraction 0.141.
- Rank 5 `excessively_sparse`: median calibration -0.703, resource-constrained rank 5 (median -0.579), useful 0.012, resource-adjusted 0.039, input degree 8.0, active fraction 0.011.

## Important Negative Control

The grid point that best matches DANDI active spatial-bin fraction alone is `dense_expansion` (input degree 32, observed output active fraction 0.297). This differs from the resource-constrained calibration optimum, reinforcing that spatial-bin activity should not be equated directly with the toy model's instantaneous output-active fraction.

## Interpretation

The best named architecture by raw empirical calibration is `dense_expansion`, whereas the best named architecture after resource and morphology constraints are emphasized is `intermediate`. This is the fitted version of the working-model claim: raw separation/information terms can favor dense expansion, but resource-constrained nontrivial expansion favors sparse granule-like designs; excessive sparsity loses useful information. NeuroMorpho directly supports limited-branch, compact input-sampling logic but also shows that dentate and cerebellar granule cells are not geometrically identical. DANDI supports nontrivial dentate granule spatial selectivity and weak-to-moderate population-vector separation, but it does not yet prove full behavioral pattern separation.

## Outputs

- Targets: `Project/results/aim3_empirical_calibration_targets.tsv`
- Calibrated grid: `Project/results/aim3_empirical_calibration_grid.tsv`
- Architecture summary: `Project/results/aim3_empirical_calibration_architecture_summary.tsv`
- Plot: `Project/results/aim3_empirical_calibration.png`
