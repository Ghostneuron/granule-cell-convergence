# Transcriptomic Configuration Model

Date built: 2026-06-22

## Purpose

This analysis tests whether granule cells are distinguished by a configuration of shared neuronal construction programs rather than by unique single pathways.

## Configuration Scores

- Downstream construction balance = mean(neurite/morphology, synaptic/excitability) - neurogenic niche/progenitor rank.
- Regional fate balance = branch-matched fate rank - branch-opposed fate rank.
- Combined configuration score = downstream construction balance + regional fate balance.

## Main Result

- Named granule-versus-comparator contrasts tested: 4.
- Contrasts with positive combined configuration score: 4/4.
- Wilcoxon p, combined configuration greater than comparator: 0.0625.
- Wilcoxon p, construction-over-niche balance greater than comparator: 0.438.
- Wilcoxon p, regional fate balance greater than comparator: 0.125.

Median role-level combined configuration scores:
- Dentate granule: 1.052.
- Pyramidal comparator: 0.938.
- Cerebellar granule: 0.406.
- Purkinje comparator: -0.312.

Sample-level contrasts:
- Dentate granule vs pyramidal comparator (10X_all_cells): construction delta 0.115, fate-polarity delta -0.083, combined configuration delta 0.115.
- Cerebellar granule vs Purkinje comparator (P0): construction delta -0.281, fate-polarity delta 0.375, combined configuration delta 0.094.
- Cerebellar granule vs Purkinje comparator (P8a): construction delta 0.156, fate-polarity delta 0.562, combined configuration delta 0.719.
- Cerebellar granule vs Purkinje comparator (P8b): construction delta 0.188, fate-polarity delta 0.812, combined configuration delta 1.000.

## Interpretation

- This supports the user's idea that morphology may be encoded in the transcriptome as a configuration of shared neuronal construction programs, not as a unique granule-only gene list.
- The configuration signal is stronger than single downstream-module specificity because it accounts for module balance: construction over progenitor/niche state plus correct regional fate polarity.
- This remains transcriptomic inference. Protein localization, local translation, post-translational control, activity, and physical circuit constraints still require external validation.

## Outputs

- Configuration units: `Project/results/primary_core_transcriptomic_configuration_units.tsv`
- Role summary: `Project/results/primary_core_transcriptomic_configuration_role_summary.tsv`
- Contrast summary: `Project/results/primary_core_transcriptomic_configuration_contrasts.tsv`
- Plot: `Project/results/primary_core_transcriptomic_configuration_model.png`
