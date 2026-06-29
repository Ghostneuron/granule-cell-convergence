# Dataset-level candidate-program interpretation

## What was added

I added a dataset-aware analysis layer that collapses the refined per-cell/per-spot candidate calls into dataset, sample, group, region, platform, species, and analysis-class units. This is an anti-pseudoreplication check: it asks whether the dentate/cerebellar pattern remains visible after summarizing groups rather than treating every cell as an independent replicate.

New files:

- `Project/scripts/dataset_level_candidate_program_analysis.py`
- `Project/results/refined_dataset_level_granule_program_units.tsv`
- `Project/results/refined_dataset_level_granule_program_statistics.tsv`
- `Project/results/refined_dataset_level_granule_program_leave_one_dataset_out.tsv`
- `Project/results/refined_dataset_level_granule_program_identity_structural_units.png`

## Unit counts

The refined unit table contains:

- 41 dentate candidate units
- 8 cerebellar candidate units
- 20 known non-dentate reference units
- 42 other or ambiguous units
- 8 cerebellum warning units
- 38 dentate low-support units
- 2 organoid granule-like units

The candidate units come from multiple datasets rather than a single source. Dentate candidate units are represented in `GSE104323`, `GSE214309`, `GSE214905`, `GSE292261`, and `GSE95752`; cerebellar candidate units are represented in `GSE122357`, `GSE165657`, `GSE242688`, and `GSE312658`.

## Main result

The regional identity contrast is robust at the unit level.

- Dentate candidate versus cerebellar candidate identity contrast: median 0.649 versus -0.237; difference 0.886; Mann-Whitney p = 9.74e-06.
- Dentate candidates versus known non-dentate references: median 0.649 versus 0; Mann-Whitney p = 2.76e-10.
- Cerebellar candidates versus ambiguous cells/spots for cerebellar identity: median 0.367 versus 0.063; Mann-Whitney p = 5.68e-05.
- Sign test: 41/41 dentate candidate units have positive dentate-minus-cerebellar identity contrast, p = 4.55e-13.
- Sign test: 8/8 cerebellar candidate units have negative dentate-minus-cerebellar identity contrast, p = 0.00391.

This strengthens the central interpretation: dentate and cerebellar granule-cell candidates are not being merged by a generic granule marker set. They retain distinct regional identity signatures.

## Structural-program result

The structural/morphogenesis signal behaves as expected for a convergence hypothesis.

- Dentate candidate units: median structural rank 0.623.
- Cerebellar candidate units: median structural rank 0.683.
- Dentate versus cerebellar structural rank is not significantly different, Mann-Whitney p = 0.760.
- 33/41 dentate candidate units are above the within-sample median structural rank, p = 5.61e-05.
- 7/8 cerebellar candidate units are above the within-sample median structural rank, p = 0.0352.

This supports partial convergence on structural/wiring programs. It does not mean the structural module uniquely identifies granule cells. Some non-dentate developmental neuronal reference units also show high structural rank, which is biologically plausible because axon guidance, cytoskeletal remodeling, and synapse genes are reused by many developing neurons.

## Leave-one-dataset-out check

The identity result is not driven by one dataset. After excluding each dataset one at a time, the dentate-minus-cerebellar identity difference remains positive in every case, ranging from 0.825 to 1.746. The structural-rank medians for both dentate and cerebellar candidates also remain above 0.5 in every leave-one-dataset-out run.

This is an important project milestone. It means the working model is not just a single-study artifact.

## Interpretation for the project

The project remains a strong "yes." The refined, dataset-aware result supports a publishable conceptual frame:

Dentate and cerebellar granule cells preserve distinct regional identity programs, but converge partially on downstream neuronal structural and wiring modules. This convergence is a plausible molecular route to similar compact granule-cell morphology despite different developmental origins and circuit contexts.

The novelty is not simply saying that both are called granule cells. The stronger and more modern question is whether similar morphology reflects reuse of a shared morphogenesis toolkit downstream of different lineage-specifying programs. The current data support that question well enough to justify a focused in silico project.

## Caveats

This remains a first-pass, marker-panel analysis. Important limitations remain:

- Units from the same dataset are not fully independent.
- Platform, species, age, and tissue-preparation effects are still present.
- Cerebellar candidate units are fewer than dentate candidate units.
- The module scores are targeted marker summaries, not yet full integrated normalization.
- Structural modules should be interpreted as shared developmental/morphogenetic programs, not granule-specific identity markers.

These caveats are manageable. They define the next analysis layer rather than weakening the project.

## Recommended strategy

Use this refined dataset-level result as the foundation for the next phase.

1. Build object-level Seurat/SCE or AnnData files for the most informative datasets.
2. Harmonize gene symbols and add ortholog mapping before mouse, rat, and human comparisons.
3. Curate cell-type annotations from source metadata, marker calls, and clustering, especially for cerebellar granule cells.
4. Compute normalized module scores or rank-based scores within each dataset, with control genes.
5. Use pseudobulk or mixed-effect models with dataset/sample as the replicate structure.
6. Split analyses by developmental stage, because structural convergence is likely strongest during maturation, migration, axon extension, dendrite remodeling, and synapse formation.
7. Add regulatory-network analysis to distinguish upstream regional identity regulators from downstream morphology executors.
8. Use the result to write either a review with new computational synthesis, or an original analysis paper if the object-level model confirms the unit-level result.

The current evidence is strong enough to proceed toward a manuscript-shaped analysis.
