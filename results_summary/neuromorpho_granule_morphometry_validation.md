# NeuroMorpho Granule Morphometry Validation

Date built: 2026-06-23

## Purpose

This is the first direct morphology validation layer for Aim 3. It uses all available cerebellar granule-cell hits from the current strict NeuroMorpho query and a reproducible species-stratified dentate granule-cell sample.

## Sampling

| Region | Metadata rows | Included rows |
|---|---:|---:|
| `cerebellum` | 62 | 62 |
| `dentate_gyrus` | 558 | 558 |

Included rows require dendrites in the reconstruction domain, moderate or complete dendritic integrity, and non-missing `n_stems`, `n_branch`, and `length` morphometry.

## Main Morphometry Result

| Metric | Dentate granule median (IQR) | Cerebellar granule median (IQR) | Interpretation |
|---|---:|---:|---|
| Primary stems | 1.50 (1.00-3.00) | 4.00 (3.00-5.00) | Different implementation: DG tends toward one primary dendritic tree, CB toward multiple short stems/claws. |
| Branches | 21.00 (13.00-32.00) | 20.00 (15.25-25.00) | Branch-count scale is closer than dendritic length. |
| Bifurcations | 10.00 (6.00-15.00) | 8.00 (5.25-10.00) | Useful proxy for limited branching complexity. |
| Dendritic length | 1352.94 (648.77-2836.80) | 274.72 (184.61-415.36) | DG has a much larger dendritic field in this sample. |
| Branch order | 5.00 (4.00-6.00) | 3.00 (2.00-4.00) | Branch order separates regional geometry. |

## Interpretation For The Project

The morphometry supports a refined version of the morphology hypothesis. Dentate and cerebellar granule cells are not geometrically identical. Instead, both are compact excitatory input-expansion neurons with constrained dendritic branching, but they implement input sampling differently: dentate granule cells through a larger dendritic tree, cerebellar granule cells through several short claw-like dendritic stems.

This means the Aim 3 computational model should not collapse morphology into one `input_degree` parameter. The better model has at least two anatomical knobs:

1. `primary_input_sampling_stems_or_claws`: lower-bound input-sampling geometry.
2. `dendritic_field_complexity`: branch count, bifurcation count, and dendritic length/scale.

## Caveats

- The dentate set is a species-stratified sample, not the full 9,672-cell NeuroMorpho dentate query.
- The cerebellar set is small and species-heterogeneous under the current strict query.
- Morphometric stems and branches are not direct synaptic input counts.
- Genotype, disease/condition, reconstruction completeness, and archive-specific methods should be filtered more strictly before final manuscript statistics.

## Outputs

- Metadata table: `Project/results/neuromorpho_granule_morphometry_metadata.tsv`
- Morphometry table: `Project/results/neuromorpho_granule_morphometry_sample.tsv`
- Summary table: `Project/results/neuromorpho_granule_morphometry_summary.tsv`
- DG-vs-CB comparison: `Project/results/neuromorpho_granule_morphometry_comparison.tsv`
- Aim 3 empirical priors: `Project/results/neuromorpho_aim3_input_degree_priors.tsv`
- Plot: `Project/results/neuromorpho_granule_morphometry_validation.png`
