# DANDI 000003 Multi-Session Spatial Extension

Date built: 2026-06-24

## Purpose

This extends the one-session DANDI 000003 spatial pilot into a reusable multi-session workflow. It analyzes all locally available NWB files from the legacy internal data folder and the external-drive raw-data folder, then writes a ranked plan for downloading additional feasible sessions.

## Current Local Coverage

- Local NWB files detected: 6
- Local DANDI 000003 size analyzed or available: 39.82 GB
- Sessions successfully analyzed: 6
- Sessions with labeled granule units: 6
- Units analyzed: 124
- Labeled granule units analyzed: 26

## Pooled Granule Metrics

- Granule units: 26
- Median spatial information: 0.7800 bits/spike
- Median spatial sparsity: 0.4537
- Median active spatial-bin fraction: 0.5489
- Median awake-moving firing rate: 0.0627 Hz

## Granule Population-Vector Checks

- `YutaMouse41-150829`: n_units=3, occupied_bins=125, rho=0.0942, far-minus-near euclidean=0.1842.
- `YutaMouse37-150609`: n_units=5, occupied_bins=251, rho=0.1837, far-minus-near euclidean=1.2293.
- `YutaMouse42-151102`: n_units=2, occupied_bins=249, rho=0.0642, far-minus-near euclidean=0.0116.
- `YutaMouse55-160908`: n_units=7, occupied_bins=246, rho=0.1945, far-minus-near euclidean=0.0615.
- `YutaMouse55-160909`: n_units=2, occupied_bins=235, rho=0.0763, far-minus-near euclidean=0.0280.
- `YutaMouse55-160907`: n_units=7, occupied_bins=242, rho=0.1759, far-minus-near euclidean=0.2260.

## Next Size-Ranked Download Candidate

- Session: `YutaMouse44-151128`
- Subject: `YutaMouse44`
- Size: 7.04 GB
- Asset: `9a2b884e-620c-43cc-8071-31b9e201cadf`

## Targeted Download Priority

The size-ranked candidate is not necessarily the best biological next step. The targeted priority table favors expected granule-label yield per GB and subject-breadth tradeoffs.

- Targeted session: `YutaMouse55-160911`
- Targeted subject: `YutaMouse55`
- Targeted track: `yield_first`
- Targeted size: 9.01 GB
- Reason: same subject has 16 local granule units; 2-day gap from local session; 9.01 GB
- Priority table: `Project/results/dandi_000003_targeted_download_priority.tsv`

## Interpretation

This file is an expansion scaffold. With the current local sessions it provides a first multi-session spatial result; as more NWB files are downloaded, the same outputs will scale directly. A manuscript-level pattern-separation claim still needs source-paper unit validation and task/trajectory-specific comparisons.

## Outputs

- Download plan: `Project/results/dandi_000003_multisession_download_plan.tsv`
- Session summary: `Project/results/dandi_000003_multisession_session_summary.tsv`
- Unit spatial metrics: `Project/results/dandi_000003_multisession_spatial_unit_metrics.tsv`
- Cell-type by session: `Project/results/dandi_000003_multisession_spatial_celltype_by_session.tsv`
- Pooled cell-type summary: `Project/results/dandi_000003_multisession_spatial_celltype_pooled.tsv`
- Population-vector separation: `Project/results/dandi_000003_multisession_population_vector_separation.tsv`
- Plot: `Project/results/dandi_000003_multisession_spatial_summary.png`
