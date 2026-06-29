# DANDI 000003 Spatial Pattern Pilot

Date built: 2026-06-23

## Purpose

This is the first position-linked pilot for Aim 3. It asks whether the already-downloaded DANDI 000003 session can support spatial coding and population-vector analyses for source-labeled dentate granule units.

## Position/Occupancy

- Position samples: 666029
- Position sampling rate: 39.0625 Hz
- Awake duration with finite position: 3317.0 s
- Awake-moving duration after speed filter: 2698.0 s
- Occupied spatial bins with >= 0.25 s occupancy: 205

Coordinates are treated as relative video-tracking coordinates. The NWB field labels the units as meters, but the observed values behave like tracking coordinates; no absolute anatomical scale is assumed.

## Granule-Cell Spatial Metrics

- Labeled granule units: 3
- Median granule spatial information: 0.3098 bits/spike
- Median granule spatial sparsity: 0.6273
- Median granule active spatial-bin fraction: 0.4634
- Median granule awake-moving rate: 0.6414 Hz
- Median granule max/mean spatial selectivity: 8.3017

## Population-Vector Check

- `granule_cell_labeled`: n_units=3, occupied_bins=125, rho(spatial, neural euclidean)=0.0942, far-minus-near euclidean=0.1842.
- `all_units`: n_units=23, occupied_bins=125, rho(spatial, neural euclidean)=0.1782, far-minus-near euclidean=3.9109.

## Interpretation

This pilot confirms that DANDI 000003 can support position-linked granule-cell activity analysis. It is still a feasibility layer: one session, three granule-labeled units, relative tracking coordinates, and no cross-session/task-condition model. The correct next step is to extend this workflow to more sessions and compute task-specific spatial information, active place-field fraction, and population-vector separation.

## Outputs

- Position summary: `Project/results/dandi_000003_pilot_position_summary.tsv`
- Unit spatial metrics: `Project/results/dandi_000003_pilot_spatial_unit_metrics.tsv`
- Cell-type spatial summary: `Project/results/dandi_000003_pilot_spatial_celltype_summary.tsv`
- Population-vector separation: `Project/results/dandi_000003_pilot_population_vector_separation.tsv`
- Plot: `Project/results/dandi_000003_pilot_granule_spatial_maps.png`
