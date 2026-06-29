# Epigenomic Extension Strategy

Date built: 2026-06-26

## Bottom line

Matched or comparable scATAC/multiome data can be added, but it should be treated as an epigenomic extension layer rather than a replacement for the current strict 10-dataset transcriptomic core.

The strongest immediate route is:

1. Use `GSE268609` as the human hippocampus/dentate multiome extension, because it is already represented in the core through RNA/selected-gene projection and its feature table supports selective promoter/gene-body and nearby-peak targeting.
2. Use the downloaded human `GSE322785` H5 files as the cerebellar multiome counterpart, because they provide adult human cerebellar cortex gene-expression and ATAC peak rows for the same regulatory target set; selected matrices and provisional marker-group scores are now available.
3. Use mouse hippocampal chromatin resources such as `GSE245367` and `GSE118987` as supporting validation, not as strict dentate/cerebellar matched-core replacements.

## Model term to add

The integrative model can add an epigenomic compatibility term:

`E_i = accessibility_or_methylation_support_i`

where `E_i` can be estimated from promoter/gene-body/gene-linked enhancer accessibility, chromVAR motif deviation, or methylation/accessibility support around curated gene modules.

The model then becomes:

`GranuleDesign_i = beta0 + betaF FatePolarity_i + betaC ConstructionBalance_i + betaT Stage_i + betaN NicheSignal_i + betaE EpigenomicCompatibility_i + betaM Morphology_i + betaA Activity_i + betaR CircuitConstraint_i + interactions + random effects + error_i`.

## Practical scoring plan

- For matched multiome datasets, transfer or use source labels to define granule, pyramidal/Purkinje, and local niche comparators.
- Link ATAC peaks to genes using promoter windows and, where feasible, gene-activity or peak-to-gene correlation.
- Score the existing modules: cerebellar fate/rhombic-lip, dentate fate/WNT/PROX1, neurogenic niche/progenitor, neurite/morphology, synaptic/excitability, and the six-gene conservative seed set.
- Compute `EpigenomicCompatibility` as module-level accessibility support in candidate granule populations minus named comparators.
- Use motif scores for candidate regulators, especially `NFIA`, `NFIB`, `RFX3`, `ATOH1`, `PROX1`, `NEUROD1`, `ZIC1/2`, and `PAX6`.

## Recommended manuscript status

This extension should be described as a feature-level epigenomic scaffold with selective peak/gene target manifests plus provisional `GSE322785` marker-group count scoring, selected-gene cluster validation, stricter cluster-supported sensitivity scoring, and broad-versus-strict robust contrast classification, not as a completed chromatin result, unless source-label quality is verified and matched dentate count-level scoring is completed. It should not change the current core dataset count yet.

Candidate resources are listed in `Project/results/epigenomic_extension_candidate_resources.tsv` and `Project/manuscript/Supplementary tables/Table_S54_epigenomic_extension_candidate_resources.tsv`; GSE268609 and GSE322785 selective feature manifests are listed in `Table_S58` through `Table_S60` and `Table_S65` through `Table_S68`; provisional `GSE322785` selected-matrix, marker-group scores, cluster-validation outputs, cluster-supported sensitivity scores, and robust contrast summaries are listed in `Table_S69` through `Table_S88`.
