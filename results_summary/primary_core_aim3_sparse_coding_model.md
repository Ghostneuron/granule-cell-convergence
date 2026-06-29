# Sparse Expansion-Coding Model

Date built: 2026-06-22

## Purpose

This analysis tests whether a granule-like architecture - high expansion, sparse input sampling, and sparse output activity - can improve pattern separation while penalizing information loss.

## Model Design

- Input patterns: 144 correlated binary patterns over 64 input channels.
- Pattern generation: 18 prototypes, input active fraction 0.12, noise rate 0.08.
- Random sparse projections define output populations.
- Parameters varied: expansion ratio, input degree, and output active fraction.
- Main score: near-pair Jaccard/overlap separation gain multiplied by information retention, output entropy, and an activity-balance penalty.

This score rewards useful separation of similar input patterns but penalizes collapse, excessive sparsity, and loss of input-distance structure. Jaccard distance is used for the main separation metric because normalized Hamming distance underestimates separation quality in very sparse binary codes.

## Main Results

- Best grid point: expansion ratio 16.0, input degree 16, output active fraction 0.2; useful score 0.318.
- Best resource-adjusted nontrivial expansion grid point: expansion ratio 4.0, input degree 2, output active fraction 0.05; resource-adjusted score 0.394.
- Best named architecture: dense_expansion_high_activity with median useful score 0.284.
- Best resource-adjusted named architecture: balanced_granule_like with median resource-adjusted score 0.122.
- Cerebellar granule-like architecture: useful score 0.229, near-pair separation gain 1.172, retention 0.192.
- Dentate granule-like architecture: useful score 0.081, near-pair separation gain 1.192, retention 0.088.
- Pyramidal/integrator-like architecture: useful score 0.164.
- Purkinje/integrator-like architecture: useful score 0.078.
- Excessive sparsity: useful score 0.004, collapse rate 0.000.

Parameter-zone summary:
- dense_expansion: median useful 0.265, near-pair separation 1.076, retention 0.418, collapse 0.000, resource-adjusted 0.010.
- granule_like_sparse_expansion: median useful 0.192, near-pair separation 1.167, retention 0.165, collapse 0.000, resource-adjusted 0.092.
- integrator_like_low_expansion: median useful 0.121, near-pair separation 1.067, retention 0.197, collapse 0.000, resource-adjusted 0.053.
- intermediate: median useful 0.110, near-pair separation 1.186, retention 0.112, collapse 0.000, resource-adjusted 0.109.
- excessively_sparse: median useful 0.012, near-pair separation 1.250, retention 0.020, collapse 0.005, resource-adjusted 0.039.

## Link To Transcriptomic Results

The model supports the computational plausibility of a sparse-expansion convergence model: sparse expansion architectures can separate similar input patterns better than low-expansion integrator-like architectures, but only when sparsity is balanced enough to preserve information.

The transcriptomic data connect to this model at the parameter-implementation level rather than as direct performance measurements:
- expansion_ratio: downstream neurite/morphology module plus granule-cell identity (downstream neurite/morphology median formal convergence 0.5). Caveat: cell number and packing are anatomical parameters, not directly inferred from expression.
- input_degree: axon guidance, adhesion, neurite, and cytoskeleton genes (14/21 neurite/morphology genes shared-positive in at least one screen). Caveat: actual dendrite number and synaptic input count require morphology or connectomics.
- output_active_fraction: synaptic/excitability module (downstream synaptic/excitability median formal convergence 0.5). Caveat: activity sparsity requires electrophysiology or calcium-imaging validation.
- architecture_configuration: identity-coupled transcriptomic configuration score (52/63 primary-core configuration-positive contrasts; median configuration delta 0.4166666666666667). Caveat: configuration evidence is identity-coupled and does not by itself prove morphology or performance.

## Interpretation

- Granule-like expansion plus sparse sampling is a plausible convergent solution for useful pattern separation in this conceptual computational validation.
- Dense expansion with high activity can score higher in raw computational terms, but it is much more expensive in wiring/activity load; resource-adjusted scoring favors sparse expansion designs.
- The model also warns against a simplistic 'more sparse is always better' story, because excessive sparsity causes information collapse.
- The transcriptomic evidence supports construction and excitability modules that could implement the parameters, but direct morphology, connectomics, or activity data would be needed to prove the parameter values in vivo.
- This fits the broader project model: morphology similarity is likely constrained by circuit computation, while the implementation remains identity-coupled and region-specific.

## Outputs

- Parameter grid: `Project/results/primary_core_aim3_sparse_coding_parameter_grid.tsv`
- Named architecture summary: `Project/results/primary_core_aim3_sparse_coding_architecture_summary.tsv`
- Transcriptomic parameter mapping: `Project/results/primary_core_aim3_transcriptomic_parameter_mapping.tsv`
- Plot: `Project/results/primary_core_aim3_sparse_coding_model.png`
