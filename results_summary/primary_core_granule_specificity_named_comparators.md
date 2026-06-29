# Granule Specificity Named-Comparator Analysis

Date built: 2026-06-22

## Purpose

This first-pass analysis asks whether the mechanism-axis similarity between dentate and cerebellar granule cells is higher than in named non-granule neuronal comparators. It uses only local datasets with explicit source labels for the relevant named groups.

## Local Named Comparators Used

- `GSE104323`: dentate granule-lineage groups (`GC-adult`, `GC-juv`, `Immature-GC`, `Neuroblast`) versus pyramidal comparators (`CA3-Pyr`, `Immature-Pyr`).
- `GSE122357`: cerebellar `Granule cells` and `Granule precursor` versus `Purkinje cells` across P0, P8a, and P8b.

The analysis scores Tier 1-4 mechanism-axis genes, then ranks source groups within each sample and axis. A positive specificity result requires both dentate granule > pyramidal comparator and cerebellar granule > Purkinje comparator.

## Scope

- GSE104323 source groups scored: 24.
- GSE122357 source groups scored: 8.
- Mechanism axes tested: 4.
- Specificity calls: not_granule_specific_vs_named_comparators: 3, cerebellar_granule_enriched_but_not_dentate_vs_pyramidal: 1.

## Axis Results

- Developmental regulatory: cerebellar_granule_enriched_but_not_dentate_vs_pyramidal; dentate-vs-pyramidal delta -0.031, cerebellar-vs-Purkinje delta 0.625, granule specificity index -0.031.
- Neurite/cytoskeleton: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.333, cerebellar-vs-Purkinje delta -0.062, granule specificity index -0.333.
- Axon guidance/adhesion: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.531, cerebellar-vs-Purkinje delta -0.250, granule specificity index -0.604.
- Synaptic/excitability: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.062, cerebellar-vs-Purkinje delta -0.188, granule specificity index -0.188.

## Interpretation

- This analysis directly addresses the user's specificity question for the named local comparators that are currently available.
- A pathway is treated as granule-enriched only if both granule branches exceed their named local comparator.
- If an axis fails this test, it may still be biologically important, but it should be phrased as a broader neuronal maturation or morphology pathway rather than a unique granule-cell pathway.
- This is still a module-level analysis. A stronger next step would add more explicit pyramidal/Purkinje datasets or Allen expression matrices and test gene-level specificity with raw-count models.

## Outputs

- Unit table: `Project/results/primary_core_granule_specificity_named_comparator_units.tsv`
- Axis summary: `Project/results/primary_core_granule_specificity_named_comparator_axis_summary.tsv`
- Gene coverage: `Project/results/primary_core_granule_specificity_named_comparator_gene_coverage.tsv`
- Plot: `Project/results/primary_core_granule_specificity_named_comparators.png`
