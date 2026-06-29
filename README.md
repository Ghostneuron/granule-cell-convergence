# Granule-cell convergence analysis code release (v0.9-submission-prep)

This packet accompanies the manuscript **"Distinct Dentate and Cerebellar
Granule-Cell Lineages Converge through Niche and Circuit Constraints"**.

Packet kind: `github-ready repository directory`

Prepared: 2026-06-28

## Contents

- `scripts/`: analysis, curation, figure and manuscript-generation scripts.
- `config/`: marker panels and small configuration tables.
- `docs/`: manuscript-facing documentation, data-source inventory, figure plan and clean methods text.
- `results_summary/`: compact result summaries, plots and machine-readable outputs used to assemble figures/tables.
- `manuscript_outputs/final_figures/`: final main figure image files.
- `manuscript_outputs/supplementary_figures/`: final supplementary figure image files.
- `manuscript_outputs/supplementary_tables/`: ordered supplementary table packet and table archives.
- `included_file_manifest.tsv`: SHA-256 checksums for included files.
- `excluded_large_or_local_files.tsv`: local files deliberately excluded from the release.

## What is not included

Raw public datasets, large sparse matrices, DANDI NWB files, GEO H5/RDS
downloads, local render-QA files and project caches are not redistributed.
Public accessions and download sources are listed in `docs/downloaded_external_data_manifest.tsv`,
`docs/external_dataset_inventory.md`, and the supplementary table packet.

## Reproducibility note

The scripts were written as project-level workflows and assume public raw
data are downloaded or reconstructed according to the dataset manifests.
For manuscript review, the recommended reproducibility path is to inspect
the included scripts together with `manuscript_outputs/supplementary_tables/`
and `results_summary/`, which preserve the analysis products used for the
manuscript figures and tables.

## License status

A final license has not yet been selected. Before public GitHub/Zenodo
release, choose a code license (for example MIT/BSD-3-Clause/GPL) and a
data/documentation license if desired.

## Citation

Use the metadata in `CITATION.cff`. The current public records are:
GitHub `https://github.com/Ghostneuron/granule-cell-convergence` and Zenodo `https://doi.org/10.5281/zenodo.21018501`.
