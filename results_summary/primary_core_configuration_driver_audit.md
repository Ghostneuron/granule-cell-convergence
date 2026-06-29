# Configuration Driver Audit

Date built: 2026-06-22

## Purpose

This audit decomposes the transcriptomic configuration score into downstream construction balance and regional fate polarity. It asks whether the assembly-plan signal is morphology-weighted, identity-weighted, or both.

## Driver Classes

- All contrasts: 56/67 configuration-positive.
- Both components positive: 28.
- Fate-driven positive: 27.
- Construction-driven positive: 1.
- Configuration not positive: 11.

Primary-core candidate-background layer:
- Configuration-positive: 52/63.
- Both components positive: 26.
- Fate-driven positive: 26.
- Construction-driven positive: 0.
- Median construction delta: 0.000.
- Median fate-polarity delta: 0.417.
- Median combined configuration delta: 0.417.

Local named-comparator layer:
- Configuration-positive: 4/4.
- Both components positive: 2.
- Fate-driven positive: 1.
- Construction-driven positive: 1.

## Base Module Delta Summary

- Cerebellar fate: median candidate-background delta 0.208; 44/63 positive.
- Dentate fate: median candidate-background delta 0.500; 49/63 positive.
- Niche/progenitor: median candidate-background delta 0.250; 49/63 positive.
- Neurite/morphology: median candidate-background delta 0.500; 38/63 positive.
- Synaptic/excitability: median candidate-background delta 0.500; 45/63 positive.
- Downstream mean: median candidate-background delta 0.500; 44/63 positive.
- Matched fate: median candidate-background delta 0.500; 61/63 positive.
- Opposed fate: median candidate-background delta 0.042; 32/63 positive.

## Top Downstream Configuration Genes

- `GPM6A` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `ROBO2` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `DCC` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `CADM3` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `STMN2` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `STMN3` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `DPYSL2` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `MAP1B` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `BASP1` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `CFL1` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `RTN1` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.
- `RTN3` (Downstream neurite/morphology): convergence 0.500, branch bias 0.000, pattern shared_positive_both_screens.

## Interpretation

- The broad primary-core configuration signal is real, but it is strongly identity-coupled: regional fate polarity is positive in most primary contrasts, while downstream construction balance is more selective.
- This supports the phrase 'transcriptomic assembly configuration' more than a pure morphology-only transcriptomic code.
- The most defensible claim is that distinct regional fate programs place granule cells into a permissive context, while shared downstream neurite/synaptic machinery contributes the morphology implementation layer.
- A stronger morphology-specific test should add explicit pyramidal/Purkinje comparator labels in more datasets, morphology-linked datasets, or spatial/proteomic localization.

## Outputs

- Driver contrast table: `Project/results/primary_core_configuration_driver_audit_contrasts.tsv`
- Module delta table: `Project/results/primary_core_configuration_driver_audit_module_deltas.tsv`
- Summary table: `Project/results/primary_core_configuration_driver_audit_summary.tsv`
- Gene priority table: `Project/results/primary_core_configuration_driver_audit_gene_priorities.tsv`
- Plot: `Project/results/primary_core_configuration_driver_audit.png`
