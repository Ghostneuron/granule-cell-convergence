# Hierarchical Integrative Granule-Cell Model

## Purpose

This model formalizes the current project as a hierarchical evidence synthesis. The available evidence is not measured on one matched per-cell table: transcriptomic terms are sample/contrast-level, stage and niche terms are branch-level, epigenomic terms are provisional marker-group contrasts, and morphology/activity/resource terms are external calibration layers. The model therefore keeps each observation level explicit.

## Ideal Matched-Data Model

If matched transcriptome, epigenome, morphology, activity, and circuit measurements were available for the same units, the target model would be:

```text
GranuleDesign_i = beta0
  + betaF * FatePolarity_i
  + betaC * ConstructionBalance_i
  + betaI * FatePolarity_i:ConstructionBalance_i
  + betaT * StageWindow_branch/stage(i)
  + betaN * NicheSignal_branch(i)
  + betaE * EpigenomicCompatibility_resource(i)
  + betaM * MorphologySparseSampling_branch(i)
  + betaA * ActivitySparsity_branch(i)
  + betaR * CircuitResourceConstraint
  + random effects for dataset, species, source layer, and assay
  + error_i
```

## Implemented Current Model

The current implementation uses weighted hierarchical evidence units:

```text
S_level,term = sum_j(weight_j * normalized_score_j) / sum_j(weight_j)
S_component = weighted average of term-level S values within each component
```

Positive scores support the term's predicted contribution to granule-like configuration or calibration. Negative scores indicate branch-specific opposition, absent support, or a result that argues against a simple shared direction. The component scores are term-balanced so that a large table with many rows does not automatically dominate a smaller but relevant evidence layer.

## Outputs

- `hierarchical_integrative_model_terms.tsv`: 10 model terms and levels.
- `hierarchical_integrative_model_evidence_units.tsv`: 560 evidence units.
- `hierarchical_integrative_model_layer_summary.tsv`: level-by-term evidence summaries.
- `hierarchical_integrative_model_branch_summary.tsv`: branch-level term summaries.
- `hierarchical_integrative_model_component_scores.tsv`: component-level scores.

## Current Summary

- Direct transcriptomic configuration: term-balanced score 0.348 (moderate_support).
- Epigenomic extension: term-balanced score 0.328 (moderate_support); this remains provisional.
- Overall hierarchical model: term-balanced score 0.456 (moderate_support).

The model supports the manuscript interpretation that dentate and cerebellar granule cells are not the same recent lineage or a single gene-barcode-defined cell type. Instead, distinct regional lineages appear to converge on a granule-like design through branch-specific fate polarity, shared downstream construction modules, stage/niche gating, regulatory compatibility, and resource-constrained sparse-expansion architecture.
