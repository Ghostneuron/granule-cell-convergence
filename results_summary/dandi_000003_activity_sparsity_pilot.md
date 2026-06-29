# DANDI 000003 Activity/Sparsity Pilot

Date built: 2026-06-23

## Purpose

This layer prepares direct dentate activity validation for Aim 3 using DANDI 000003. The full archive is too large for blind download, so this step builds a complete asset manifest, chooses the smallest pilot NWB file, and analyzes local NWB files only when present.

## Asset Inventory

- Assets listed: 101
- Total archive size from manifest: 2559.25 GB
- Smallest asset: `sub-YutaMouse41/sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb` (4.66 GB)
- Pilot local file present: True

## Pilot Asset

| Asset | Subject | Session | Size GB | Local file |
|---|---|---|---:|---|
| `sub-YutaMouse41/sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb` | `YutaMouse41` | `YutaMouse41-150829` | 4.66 | `/Users/jili/Desktop/codex_play/Hippocanpus&Cerebellum/External_Data/DANDI/000003/sub-YutaMouse41/sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb` |

## Pilot Activity Result

- Units: 23
- Total spikes: 953531
- Recording duration: 17050.16 s
- Median unit firing rate: 1.1614 Hz
- Median unit active-bin fraction at 1 s bins: 0.4454
- Mean population active fraction per 1 s bin: 0.4705

Granule-cell-labeled pilot subset:
- Labeled granule units: 3
- Median granule-cell firing rate: 0.9131 Hz
- Median granule-cell active-bin fraction at 1 s bins: 0.3943

Behavior-state intervals in this pilot:
- `awake`: 5 intervals, 3317.0 s
- `nrem`: 13 intervals, 8329.0 s
- `rem`: 11 intervals, 913.0 s
- `transit`: 11 intervals, 125.0 s

Granule-cell state-dependent pooled rates:
- `awake`: 0.5674 Hz pooled unit rate; 1.000 active unit-interval fraction
- `nrem`: 1.3297 Hz pooled unit rate; 1.000 active unit-interval fraction
- `rem`: 0.4980 Hz pooled unit rate; 1.000 active unit-interval fraction
- `transit`: 1.0987 Hz pooled unit rate; 0.939 active unit-interval fraction

These are pilot values only. The source labels identify three units as granule cells, but broader conclusions require more sessions and source-paper validation of unit identity conventions.

## Caveats

- DANDI 000003 is excellent for direct activity validation, but the archive is 2.56 TB and the smallest NWB is 4.66 GB.
- Pilot extraction uses generic NWB `/units` HDF5 fields and does not yet classify units as granule cells versus mossy cells or other hippocampal units.
- Firing sparsity here is an electrophysiological activity measure, not a transcriptomic property.
- Pattern-separation metrics require position/task epochs and population-vector analysis after unit identity and behavior timestamps are verified.

## Outputs

- Asset manifest: `Project/results/dandi_000003_asset_manifest.tsv`
- Pilot download plan: `Project/results/dandi_000003_pilot_asset_plan.tsv`
- NWB structure table: `Project/results/dandi_000003_pilot_nwb_structure.tsv`
- Unit sparsity table: `Project/results/dandi_000003_pilot_unit_sparsity.tsv`
- Session sparsity summary: `Project/results/dandi_000003_pilot_session_sparsity_summary.tsv`
- Cell-type sparsity summary: `Project/results/dandi_000003_pilot_celltype_sparsity_summary.tsv`
- Behavior states: `Project/results/dandi_000003_pilot_behavior_states.tsv`
- State-unit firing table: `Project/results/dandi_000003_pilot_state_unit_firing.tsv`
- State-cell-type summary: `Project/results/dandi_000003_pilot_state_celltype_summary.tsv`
- Plot: `Project/results/dandi_000003_pilot_activity_sparsity.png`
