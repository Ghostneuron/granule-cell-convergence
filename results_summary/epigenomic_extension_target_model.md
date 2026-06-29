# Epigenomic Extension Target Model

Date built: 2026-06-26

## Purpose

This selective-download scaffold defines the regulatory targets and model terms needed to add a matched or comparable scATAC/multiome/methylome layer to the granule-cell convergence model.

## Target Summary

- Target rows: 252.
- Unique genes: 154.
- Unique genes present as `GSE268609` gene-expression features: 153/154.
- Full `GSE268609` matrix status: not present locally.
- Unique target genes present in downloaded human `GSE322785` H5 files: 153/154.
- Per-sample target-gene rows present in downloaded human `GSE322785` H5 files: 459/462.
- Human `GSE322785` peak-target rows: 9074.
- Human `GSE322785` selective manifest rows: 9533.
- Human `GSE322785` barcodes processed for selected count extraction: 2095603.
- Human `GSE322785` basic-QC barcodes: 24206.
- Human `GSE322785` high/medium-confidence provisional marker calls: 6047.
- Human `GSE322785` selected matrix nonzero entries: 6384881.
- Human `GSE322785` provisional marker-group module-score rows: 4568.
- Human `GSE322785` selected-gene cluster validation mean ARI: 0.031.
- Human `GSE322785` selected-gene cluster validation mean NMI: 0.151.
- Human `GSE322785` cluster-supported barcodes: 3550.
- Human `GSE322785` cluster-supported marker groups: 5.
- Human `GSE322785` cluster-supported granule-candidate barcodes: 356.
- Human `GSE322785` cluster-supported module-score rows: 2630.
- Human `GSE322785` robust-positive broad-vs-supported contrasts: 263.
- Human `GSE322785` strong robust-positive broad-vs-supported contrasts: 97.

## Integrative Model

`GranuleDesign_i = beta0 + betaF FatePolarity_i + betaC ConstructionBalance_i + betaT Stage_i + betaN NicheSignal_i + betaE EpigenomicCompatibility_i + betaM Morphology_i + betaA Activity_i + betaR CircuitConstraint_i + interactions + random effects + error_i`.

The new epigenomic term is `EpigenomicCompatibility`, estimated from promoter/gene-body accessibility, linked enhancer accessibility, methylation/accessibility support, and motif deviation near fate, construction, and candidate genes.

## Top Resource Route

- `GSE268609` (high): Use the existing peak-target manifest for selective extraction of promoter/gene-body and nearby enhancer accessibility around fate, niche, and construction modules
- `GSE322785` (high): Use existing selected matrices and provisional marker-group scores to prioritize regulatory modules, then refine granule/Purkinje/glial labels with clustering or source-taxonomy transfer before final chromatin claims
- `GSE245367` (medium_high): Use as mouse hippocampal regulatory validation; prioritize gene-linked methylation/accessibility around conservative seed and construction genes
- `GSE118987` (medium): Use as older hippocampal ATAC sanity check for candidate regulatory-element openness and motif enrichment

## Outputs

- Regulatory targets: `Project/results/epigenomic_extension_regulatory_targets.tsv`
- Model-term specification: `Project/results/integrative_granule_model_term_specification.tsv`
- Target summary: `Project/results/epigenomic_extension_target_summary.tsv`
