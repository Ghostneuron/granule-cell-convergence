# Primary-Core Transcriptomic Configuration Validation

Date built: 2026-06-24

## Purpose

This analysis broadens the local transcriptomic configuration test from named pyramidal/Purkinje comparators to primary-core candidate-versus-background pseudobulk contrasts.

## Inputs

- Full MGI one-to-one ortholog matrix expression layer.
- 2,169-gene selected-feature primary-core pseudobulk layer.
- The five niche/circuit module gene sets used in the local configuration model.

## Scale

- Configuration class units: 210 across 10 datasets.
- Candidate-versus-background contrasts: 63 across 7 datasets.
- Expression layers represented: 2.

## Main Result

- Combined configuration positive contrasts: 52/63.
- Median candidate-background combined configuration delta: 0.417.
- Sign-test p for positive combined configuration deltas: 8.37e-08.
- Wilcoxon p for combined configuration greater than zero: 4.89e-08.

Layer-level summary:
- full_mgi_ortholog_matrix: 8/8 positive, median delta 1.500, Wilcoxon p=0.00391.
- selected_feature_matrix: 44/55 positive, median delta 0.333, Wilcoxon p=2.96e-06.

Branch-level summary:
- cerebellar: 12/12 positive, median delta 1.000, Wilcoxon p=0.000244.
- dentate: 40/51 positive, median delta 0.250, Wilcoxon p=2.18e-05.

Minimum module gene coverage across expression layers:
- Shared neurogenic niche/progenitor state: 4 genes present.
- Cerebellar fate/rhombic-lip/SHH: 6 genes present.
- Dentate fate/WNT/PROX1: 11 genes present.
- Downstream synaptic/excitability: 20 genes present.
- Downstream neurite/morphology: 21 genes present.

## Interpretation

- This strengthens the configuration hypothesis beyond the four named local contrasts: candidate granule classes usually show higher combined configuration scores than local backgrounds across the primary-core pseudobulk layers.
- The test is broader but less specific than the named-comparator audit because most datasets do not preserve explicit pyramidal or Purkinje comparator labels.
- The result supports manuscript language that morphology is partly encoded as transcriptomic module balance, while final geometry still requires spatial, proteomic, activity, and lineage/perturbation validation.

## Outputs

- Configuration units: `Project/results/primary_core_transcriptomic_configuration_primary_units.tsv.gz`
- Contrast table: `Project/results/primary_core_transcriptomic_configuration_primary_contrasts.tsv`
- Summary table: `Project/results/primary_core_transcriptomic_configuration_primary_summary.tsv`
- Coverage table: `Project/results/primary_core_transcriptomic_configuration_primary_coverage.tsv`
- Plot label key: `Project/results/primary_core_transcriptomic_configuration_primary_validation_label_key.tsv`
- Plot: `Project/results/primary_core_transcriptomic_configuration_primary_validation.png`
