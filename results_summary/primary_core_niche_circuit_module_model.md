# Niche/Fate Versus Circuit/Morphology Module Model

Date built: 2026-06-22

## Purpose

This analysis asks whether dentate and cerebellar granule-cell similarity is better explained by shared upstream stem-cell niche/fate programs or by convergent downstream circuit, morphology, and maturation programs.

## Module Families

- Cerebellar fate/rhombic-lip/SHH: 13 genes; cerebellar-specific upstream fate and progenitor-expansion logic.
- Dentate fate/WNT/PROX1: 17 genes; dentate-specific hippocampal granule fate and WNT-associated identity logic.
- Shared neurogenic niche/progenitor state: 21 genes; general neurogenic stem/progenitor and local niche signaling state.
- Downstream neurite/morphology: 21 genes; convergent neurite, axon, adhesion, and morphology implementation.
- Downstream synaptic/excitability: 20 genes; convergent sparse-coding, synaptic maturation, and excitability implementation.

## Formal-Core Result

- Genes scored in the formal model: 92/92.
- Median upstream/niche convergence delta: -0.500.
- Median downstream circuit/morphology convergence delta: 0.500.
- Mann-Whitney test, downstream greater than upstream/niche: p=1.56e-09.

Top downstream modules by formal convergence:
- Downstream neurite/morphology: median convergence 0.500, 14/21 genes shared-positive in at least one screen.
- Downstream synaptic/excitability: median convergence 0.500, 17/20 genes shared-positive in at least one screen.

Best upstream/niche modules by formal convergence:
- Dentate fate/WNT/PROX1: median convergence 0.000, 5/17 genes shared-positive in at least one screen.
- Cerebellar fate/rhombic-lip/SHH: median convergence -0.250, 2/13 genes shared-positive in at least one screen.

## Named-Comparator Result

- Comparator units scored: 240.
- Direct comparator rule: dentate granule-lineage groups must exceed pyramidal labels, and cerebellar granule-lineage groups must exceed Purkinje labels.
- Cerebellar fate/rhombic-lip/SHH: cerebellar_granule_enriched_only; dentate-vs-pyramidal delta 0.000, cerebellar-vs-Purkinje delta 0.562.
- Dentate fate/WNT/PROX1: dentate_granule_enriched_only; dentate-vs-pyramidal delta 0.063, cerebellar-vs-Purkinje delta -0.125.
- Shared neurogenic niche/progenitor state: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.260, cerebellar-vs-Purkinje delta -0.375.
- Downstream neurite/morphology: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.229, cerebellar-vs-Purkinje delta -0.188.
- Downstream synaptic/excitability: not_granule_specific_vs_named_comparators; dentate-vs-pyramidal delta -0.104, cerebellar-vs-Purkinje delta -0.188.

## Interpretation

- The current evidence favors a convergence model: upstream fate modules are branch-specific or mixed, while downstream morphology/excitability modules carry the strongest formal shared-convergence signal.
- The named-comparator layer remains a constraint: several downstream modules are not uniquely granule-specific versus pyramidal and Purkinje comparators, so the claim should be circuit-convergence rather than pathway uniqueness.
- This analysis can be strengthened later by adding developmental time-resolved trajectory inference and perturbation/lineage-tracing evidence from the literature.

## Outputs

- Gene-set table: `Project/results/primary_core_niche_circuit_module_gene_sets.tsv`
- Formal gene scores: `Project/results/primary_core_niche_circuit_module_formal_gene_scores.tsv`
- Formal module summary: `Project/results/primary_core_niche_circuit_module_formal_summary.tsv`
- Named-comparator units: `Project/results/primary_core_niche_circuit_module_named_comparator_units.tsv`
- Named-comparator summary: `Project/results/primary_core_niche_circuit_module_named_comparator_summary.tsv`
- Plot: `Project/results/primary_core_niche_circuit_module_model.png`
