# Storage Cleanup Log

Date: 2026-06-23

## Deleted Raw Or Cache Files

The following files were deleted to recover disk space after their derived analysis outputs were already written:

| File | Approx size | Reason |
|---|---:|---|
| `External_Data/GEO/GSE186538/GSE186538_Human_counts.mtx.gz` | 1.4 GB | DG granule candidate sparse subset already exported to `Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates`. |
| `External_Data/GEO/GSE268609/GSE268609_matrix.mtx.gz` | 4.3 GB | Selected-gene RNA object already exported to `Project/processed/gse268609_rna_selected`. |
| `External_Data/GEO/GSE325391/GSE325391_adultgc_filtered.RDS` | 3.3 GB | Selected sparse bridge already exported to `Project/processed/gse325391_adult_dg_selected`. |
| `External_Data/DANDI/000003/sub-YutaMouse37/sub-YutaMouse37_ses-YutaMouse37-150617_behavior+ecephys.nwb` | 6.2 GB | DANDI comparator session had 0 source-labeled granule units; six-session derived tables and plot were already generated. |
| `Project/scripts/__pycache__` | 1.1 MB | Regenerable Python bytecode cache. |
| `Project/results/mplconfig` | 124 KB | Regenerable Matplotlib cache. |

## Additional DNMT3L Cleanup

The following out-of-project raw/intermediate folders were also deleted after confirming that derived outputs, tables, figures, manuscripts, and project results were preserved:

| Folder | Approx size | Reason |
|---|---:|---|
| `/Users/jili/Desktop/codex_play/DNMT3L_Evolution/testis_transcriptome_h358/raw_balanced_pilot/raw_fastq` | 33 GB | Raw FASTQ files; quant/tables/figures/scripts were retained. |
| `/Users/jili/Desktop/codex_play/DNMT3L_Evolution/testis_transcriptome_h358/raw_balanced_pilot/sra_cache` | 20 GB | SRA cache files; regenerable from public accessions if needed. |
| `/Users/jili/Desktop/codex_play/DNMT3L_Evolution/testis_transcriptome_h358/raw_balanced_pilot/indices` | 1.9 GB | Regenerable alignment/quantification indices. |
| `/Users/jili/Desktop/codex_play/DNMT3L_Neuron/tmp` | 233 MB | Temporary working directory. |
| `/Users/jili/Desktop/codex_play/DNMT3L_Neuron/.matplotlib` | 124 KB | Regenerable Matplotlib cache. |

## Additional Antigravity DNMT3L Cleanup

The following regenerable infrastructure was deleted from `/Users/jili/Desktop/antigravity_play/DNMT3L_Evolution`:

| Folder | Approx size | Reason |
|---|---:|---|
| `PyRosetta4.Release.python39.mac.wheel` | 11 GB | PyRosetta wheel/download cache; regenerable if PyRosetta needs to be reinstalled. |
| `venv` | 767 MB | Project virtual environment; dependencies can be recreated. |
| `Revision/structural_h358_modeling/docking_haddock3/.venv_haddock3` | 351 MB | HADDOCK3 virtual environment; derived docking runs were retained. |
| `Revision/structural_h358_modeling/md_smd_h358/.venv_openmm` | 266 MB | OpenMM virtual environment; MD outputs were retained. |
| `Revision/structural_h358_modeling/foldx_h358_mutagenesis/.venv_pyrosetta39_x86` | 2.9 GB | PyRosetta virtual environment; FoldX/PyRosetta result folders were retained. |
| Selected `__pycache__` folders | <1 MB | Regenerable Python bytecode caches outside the deleted virtual environments. |

## Additional Desktop DNMT3L Neuron Manuscript Check

`/Users/jili/Desktop/DNMT3L Neuron Manuscript` was inspected and found to contain manuscript documents, final figures, main figures, extended data figures, supplementary figures, and supplementary tables rather than large raw/intermediate datasets. Four `.DS_Store` metadata files totaling <50 KB were removed; the manuscript content was otherwise left intact.

## Notes

- Current DANDI result tables remain a derived six-session snapshot.
- Five raw DANDI sessions with source-labeled granule units were kept for planned task/trajectory-specific spatial analyses.
- Rebuilding the deleted raw-source analyses from scratch would require re-downloading the listed GEO/DANDI source files.
- `Project/results/dandi_000003_asset_manifest.tsv`, `Project/results/gse268609_download_status.tsv`, `Project/results/gse325391_download_status.tsv`, and `Project/results/human_dentate_core_construction_status.tsv` were updated to mark the relevant raw files as no longer local.
- DNMT3L project Git histories, manuscripts, literature, compendia, derived results, quantification tables, figures, and scripts were left intact.
- No `/Users/jili/Desktop/antigravity_play/DNMT3L_Neuron` directory was found. The matching neuron folders found were `/Users/jili/Desktop/codex_play/DNMT3L_Neuron` and `/Users/jili/Desktop/DNMT3L Neuron Manuscript`; the `codex_play` neuron folder had already been cleaned conservatively.
- Antigravity DNMT3L structural docking/modeling run outputs were retained because they are analysis/result-bearing, not simple cache.
- `/Users/jili/Desktop/DNMT3L Neuron Manuscript` remains a small manuscript assembly folder (~175 MB) and was not a meaningful storage-recovery target.
