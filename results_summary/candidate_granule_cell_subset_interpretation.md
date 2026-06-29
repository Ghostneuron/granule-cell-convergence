# Candidate granule-cell subset interpretation

## What was added

I created a conservative candidate-cell classifier from the per-cell marker-module scores. This is not a final cell annotation; it is a practical triage layer for the next in silico comparison. The classifier uses within-sample ranks for dentate identity, cerebellar identity, shared neuronal score, and structural-program score, and it respects existing trusted metadata where available.

New files:

- `Project/scripts/classify_candidate_granule_cells.py`
- `Project/results/candidate_granule_cell_calls.tsv.gz`
- `Project/results/candidate_granule_cell_call_summary.tsv`

## Overall calls

The classifier produced 142,782 cell/spot calls:

- 30,797 `candidate_cerebellar_granule`
- 11,259 `candidate_dentate_granule`
- 14,352 `known_non_dentate_reference`
- 1,364 `dentate_like_low_support`
- 973 `organoid_granule_like_candidate`
- 22,000 `cerebellum_dentate_panel_high_warning`
- 62,037 `non_granule_or_ambiguous`

## Dentate reference set

`GSE104323` now uses its metadata directly. Only the curated dentate granule lineage groups are called as reference dentate granule cells:

- `GC-adult`: 2,613 cells
- `GC-juv`: 3,420 cells
- `Immature-GC`: 2,419 cells
- `Neuroblast`: 1,381 cells

This gives a clean local dentate reference set of 9,833 cells for comparing against candidate cerebellar granule-like cells.

Other dentate-enriched resources behave as expected. `GSE292261` shows strong postnatal DG signal, with high/medium candidate dentate fractions of 70.2% at P5, 48.7% at P7, 64.8% at P10, 52.7% at P15, and only 1.9% at P28 under the current stringent thresholds. The P28 drop likely reflects maturation and marker-panel/stringency effects rather than loss of dentate identity, so P28 should be retained but interpreted cautiously.

`GSE214309` nuclei remain dentate-biased across mature/immature and active/inactive groups. Candidate fractions are variable, but the low-support cells are still mostly dentate-like; this dataset is especially useful for activity-state analysis.

## Cerebellar candidate set

The mixed cerebellar datasets now yield usable candidate subsets instead of only whole-sample averages:

- `GSE122357` P0/P8 samples: about 29-31% candidate cerebellar granule-like cells.
- `GSE165657` human cerebellum aggregate: 19,882 candidate cerebellar granule-like cells, about 27.0% of the matrix.
- `GSE312658` mouse cerebellum: 1,839 Ctrl candidate cells, about 37.9%; 2,257 cKO candidate cells, about 32.5%.

These subsets are now much better inputs for testing whether cerebellar granule cells share structural/morphogenesis programs with dentate granule cells while preserving different identity signatures.

## Warning signal

A substantial `cerebellum_dentate_panel_high_warning` category appears in cerebellar and spatial sources. This is scientifically useful, not just a nuisance. It says the current dentate panel contains genes that are not dentate-specific enough for final cross-region claims, especially broad excitatory/neurodevelopmental genes. Before formal statistics, the panels should be split into:

- strict dentate identity genes
- strict cerebellar granule identity genes
- shared granule/excitatory neuronal genes
- downstream structural/morphogenesis executor genes

The spatial/proteomics-linked `GSE242688` remains especially cautionary: many spots show dentate-panel-high behavior despite the cerebellar source. Treat this dataset as validation/context only until spatial spot identities are better resolved.

## Next analytical move

The next best step is to run the core comparison on candidate sets only:

1. Use `candidate_dentate_granule` and `candidate_cerebellar_granule` calls as the working comparison groups.
2. Recompute panel summaries on these subsets.
3. Quantify whether region-identity contrast separates the two groups while structural-program scores overlap or converge.
4. Repeat after refining identity panels to remove broad developmental/excitatory genes.
5. Then move from this lightweight marker-score approach into full Seurat/SCE normalization and ortholog-aware module scoring.

## Current conclusion

This stage strengthens the project. We now have not only evidence of a shared structural program, but also a practical candidate-cell subset strategy for separating cerebellar granule-like cells from mixed cerebellar matrices. The key manuscript claim is becoming more concrete: similar morphology likely reflects convergence on downstream structural/wiring programs, not shared regional identity.
