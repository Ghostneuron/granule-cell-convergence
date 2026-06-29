# External Data Drive Plan

Date updated: 2026-06-24

## External Download Root

Large raw-data downloads should now use:

`/Volumes/VV 2021 backup drive 01/Hippocampus_Cerebellum_downloads`

A workspace symlink points to the same location:

`Project/external_data`

## Directory Layout

- `raw/geo`: GEO/SRA/raw expression objects.
- `raw/dandi`: DANDI NWB files.
- `raw/allen`: Allen Institute raw or large intermediate files.
- `raw/neuromorpho`: NeuroMorpho raw exports.
- `cache`: resumable download caches and API caches.
- `derived`: external-drive derived files that are too large for the main workspace.
- `logs`: download manifests and processing logs.

## Use Policy

- Keep manuscript-facing summaries, figures, small TSV tables, and scripts in the main workspace.
- Put large raw downloads and large intermediate objects on the external drive.
- Record each new external download in a manifest before cleanup or deletion.
- Prefer targeted downloads that answer a figure or claim-specific question; avoid broad dataset hunting unless the manuscript review requires it.

## Current Rationale

The internal workspace has enough room for scripts and derived outputs, but not for broad raw GEO/DANDI/Allen exploration. The external drive has enough available capacity for optional raw-object rebuilds or one-off physiology downloads without destabilizing the working project folder.

## Current External Downloads

- `DANDI_000003_YutaMouse55_160907`: downloaded to `Project/external_data/raw/dandi/000003/sub-YutaMouse55/sub-YutaMouse55_ses-YutaMouse55-160907_behavior+ecephys.nwb`; HDF5/NWB validation passed; included in the refreshed DANDI multi-session spatial analysis and Aim 3 empirical calibration.
