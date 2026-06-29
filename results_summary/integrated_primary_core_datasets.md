# Integrated Primary Core Datasets

Date built: 2026-06-22

## Definition

This is the strict 10-dataset primary core for the manuscript-scale in silico project. It keeps the analysis balanced across mouse dentate, cerebellum, and human dentate/hippocampus.

`GSE185277` and `GSE185553` remain important human construction scaffolds for label tuning and marker validation, but they are not counted in this strict 10-dataset primary core.

## The 10 Primary Core Datasets

| Dataset | Branch | Primary role | Why included |
|---|---|---|---|
| `GSE104323` | mouse dentate | Primary mouse dentate reference | Best annotated mouse dentate granule lineage plus non-dentate controls. |
| `GSE95752` | mouse dentate | Dentate maturation validation | Independent dentate single-cell maturation resource. |
| `GSE292261` | mouse dentate | Postnatal dentate developmental validation | Usable stage metadata and gene symbols after curation. |
| `GSE214309` | mouse dentate | Adult/activity-state dentate validation | Mature/immature and activity-state dentate labels. |
| `GSE122357` | cerebellum | Primary mouse cerebellar developmental comparison | Best current cerebellar developmental granule-cell anchor. |
| `GSE165657` | cerebellum | Primary human cerebellar validation | Large human cerebellum aggregate for cross-species cerebellar validation. |
| `GSE312658` | cerebellum | Mouse cerebellar perturbation validation | Control/cKO dataset for perturbation sensitivity of granule programs. |
| `GSE186538` | human dentate/hippocampus | Human DG taxonomy anchor | Explicit human DG GC PROX1-labeled subset. |
| `GSE325391` | human dentate/hippocampus | Primary adult human dentate anchor | Direct adult human dentate granule-cell nuclei. |
| `GSE268609` | human dentate/hippocampus | Broader human hippocampal aging/AD RNA expansion | Large human hippocampal RNA branch with projected DG/imGC labels and diagnosis-aware module tests. |

## Current Interpretation

This is the cleanest primary core for the paper:

- Mouse dentate branch: 4 datasets.
- Cerebellar branch: 3 datasets.
- Human dentate/hippocampal branch: 3 datasets.

The rank-level integrated analysis already connects this primary core plus the constructed human bridge:

- `Project/results/human_bridge_backbone_rank_units.tsv`
- `Project/results/human_bridge_backbone_rank_statistics.tsv`
- `Project/results/human_bridge_backbone_rank_source_summary.tsv`
- `Project/results/human_bridge_backbone_rank_units.png`

The key conclusion is that dentate and cerebellar candidates remain regionally identity-separated, while structural/morphogenesis-related programs behave as the convergence axis.

## Formal Module Analysis

The strict 10-dataset core now has a first marker-panel ortholog-aware module analysis:

- `Project/results/primary_core_marker_panel_ortholog_map.tsv`
- `Project/results/primary_core_integrated_module_units.tsv`
- `Project/results/primary_core_ortholog_module_statistics.tsv`
- `Project/results/primary_core_ortholog_module_leave_one_dataset_out.tsv`
- `Project/results/primary_core_identity_structural_modules.png`
- `Project/results/primary_core_ortholog_module_analysis.md`

Current formal result: identity separation is robust, and both dentate and cerebellar candidates are above the structural-program median. The structural signal should be framed as a shared elevated executor axis, not equal magnitude across regions.

## Candidate-Gene Pseudobulk Analysis

The strict 10-dataset core now has a focused candidate-gene pseudobulk screen:

- `Project/results/primary_core_candidate_gene_pseudobulk_expression.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_hits.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_effects.png`
- `Project/results/primary_core_candidate_gene_pseudobulk_analysis.md`

Current candidate-gene result: all 10 primary datasets contribute expression evidence for the 67-gene shortlist; 14 structural-executor genes have positive candidate dentate and candidate cerebellar rank deltas, and 12 pass the exploratory BH<0.2 criterion in both branches. The strongest shared executor signals include `CFL1`, `GAP43`, `ROBO2`, `STMN2`, `STMN3`, `CDK5R1`, `DPYSL2`, `EPHB2`, `L1CAM`, `MAPT`, `DPYSL3`, and `ELAVL4`.

## Expanded Selected-Gene Screen

The strict 10-dataset core now also has a 2,169-gene selected-feature pseudobulk screen and mechanism triage:

- `Project/results/primary_core_expanded_gene_pseudobulk_expression.tsv.gz`
- `Project/results/primary_core_expanded_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_branch_specific.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.png`
- `Project/results/primary_core_expanded_gene_pseudobulk_analysis.md`
- `Project/results/primary_core_expanded_gene_mechanism_triage.tsv`
- `Project/results/primary_core_expanded_gene_mechanism_triage.md`

Current expanded result: 537,070 pseudobulk expression rows, 2,169 genes tested, 1,279 shared-positive rank genes, 579 shared-positive genes passing BH<0.10 in both branches, and 71 mechanism figure candidates after triage. The strongest recovered structural genes are `STMN2`, `GAP43`, `ROBO2`, `CFL1`, and `STMN3`; strong new candidates include `GPM6A`, `MAP1B`, `TCF4`, `NFIA`, `NFIB`, `NRXN1`, `BASP1`, `CADM1`, and `DCC`.

## Full-Matrix Same-Symbol Screen

The strict primary core now has a full-matrix same-symbol pseudobulk screen plus a cross-screen mechanism consensus:

- `Project/results/primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_coverage.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_statistics.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_branch_specific.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.png`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_analysis.md`
- `Project/results/primary_core_genomewide_symbol_mechanism_triage.tsv`
- `Project/results/primary_core_genomewide_symbol_mechanism_triage.md`
- `Project/results/primary_core_cross_screen_mechanism_consensus.tsv`
- `Project/results/primary_core_cross_screen_mechanism_consensus.md`

Current full-matrix result: 33,939 target symbols, 628,339 pseudobulk expression rows, 21,253 genes tested, 6,440 shared-positive genes, and 283 shared-positive genes passing BH<0.10 in both branches. The cross-screen consensus identifies 24 figure-level mechanism candidates that survive both selected-feature and full-matrix screens, led by `GPM6A`, `NFIA`, `NFIB`, `PPP3CA`, `CAMTA1`, `MAPK1`, `STXBP1`, `CALM2`, `CACNA2D1`, and `SYNPR`.

## Dataset-Aware Consensus Validation

The 24 cross-screen consensus candidates now have a dataset/sample-aware robustness check:

- `Project/results/primary_core_consensus_candidate_dataset_deltas.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation_heatmap.png`
- `Project/results/primary_core_consensus_candidate_dataset_validation.md`

Current validation result: 1,608 dataset/sample/gene branch-delta rows support a six-gene strongest shortlist robust across all available selected-feature/full-matrix and dentate/cerebellar tests: `GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, and `RFX3`. Eight additional genes are robust in 3/4 tests: `CACNA2D1`, `GABRB3`, `GRIN2B`, `KCND2`, `KCNJ3`, `KCNJ6`, `PPP3CA`, and `STXBP5L`.

## MGI Ortholog Meta-Model

The strict core now has an MGI-filtered ortholog-aware meta-model:

- `Project/results/primary_core_mgi_ortholog_meta_model_map.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_meta_model_branch_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_top_hits.png`
- `Project/results/primary_core_mgi_ortholog_meta_model.md`

Current ortholog result: 16,245 strict same-symbol one-to-one human-mouse pairs, 241,690 dataset/sample unit deltas, 1,304 shared ortholog-filtered hits, and 36 mechanism-prioritized hits. The six dataset-robust consensus genes all remain supported in this strict MGI frame: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.

## Expanded Ortholog Screen

The MGI ortholog screen now includes non-identical mouse/human one-to-one symbols:

- `Project/results/primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_full_matrix_statistics.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model.md`

Current expanded result: 17,611 MGI one-to-one target genes, 508,917 full-matrix pseudobulk rows, 6,413 shared-positive ortholog genes, 246 shared-positive non-identical-symbol genes, and 1,370 dataset-aware expanded meta-model shared hits. The manuscript mechanism shortlist remains led by `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.

## Formal Rank-Meta Validation

The expanded MGI ortholog model now has a formal dataset-level rank-meta validation:

- `Project/results/primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_formal_rank_branch_tests.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_model_long.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.png`
- `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`

Current formal result: 116,013 dataset-level deltas, 36,303 branch tests, 1,370 formal shared hits, and 158 formal hits supported in both selected and full-matrix screens. The six-gene seed set (`GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`) remains the safest mechanism tier, with `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`, and `ROBO2` as strong second-tier candidates.

The manuscript-facing tier packet is:

- `Project/results/primary_core_manuscript_candidate_tiers.tsv`
- `Project/results/primary_core_manuscript_candidate_tiers.md`

## Mechanism Axis Model

The manuscript candidate tiers now have a biological interpretation layer:

- `Project/results/primary_core_mechanism_axis_gene_table.tsv`
- `Project/results/primary_core_mechanism_axis_summary.tsv`
- `Project/results/primary_core_mechanism_axis_branch_summary.tsv`
- `Project/results/primary_core_mechanism_axis_model.png`
- `Project/results/primary_core_mechanism_axis_model.md`

Current model: the central hypothesis is a layered shared toolkit, not shared regional identity. Tier 1 anchors developmental regulatory control (`NFIB`, `NFIA`, `RFX3`), neurite/cytoskeleton morphogenesis (`GPM6A`), and synaptic/excitability maturation (`KCNK1`, `GABRA2`). Tier 2 adds `ROBO2`-linked guidance plus synaptic, calcium, potassium, glutamate, and vesicle-release support.

## Named-Comparator Specificity Audit

The mechanism axes have now been tested against explicit named non-granule comparators in the local primary datasets:

- `GSE104323`: dentate granule-lineage groups versus `CA3-Pyr` and `Immature-Pyr`.
- `GSE122357`: cerebellar granule-lineage groups versus `Purkinje cells`.

Main outputs:

- `Project/results/primary_core_granule_specificity_named_comparator_units.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_axis_summary.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_gene_coverage.tsv`
- `Project/results/primary_core_granule_specificity_named_comparators.png`
- `Project/results/primary_core_granule_specificity_named_comparators.md`

Current specificity result: the developmental regulatory axis is cerebellar granule-enriched versus Purkinje cells but not dentate granule-enriched versus pyramidal comparators. The neurite/cytoskeleton, axon guidance/adhesion, and synaptic/excitability axes do not pass strict named-comparator specificity. This supports a convergence model rather than a claim that the pathways themselves are uniquely granule-specific.

## Niche/Fate Versus Circuit/Morphology Model

The primary core now has a hypothesis-driven model that separates upstream fate/niche logic from downstream circuit and morphology implementation.

- Upstream modules: cerebellar fate/rhombic-lip/SHH, dentate fate/WNT/PROX1, and shared neurogenic niche/progenitor state.
- Downstream modules: neurite/morphology and synaptic/excitability.
- Formal result: downstream modules have median convergence delta 0.500, while upstream/niche modules have median convergence delta -0.500.
- Comparator result: upstream fate modules are branch-specific, and downstream modules are formally convergent but not uniquely granule-specific versus named pyramidal/Purkinje comparators.

Main outputs:

- `Project/results/primary_core_niche_circuit_module_gene_sets.tsv`
- `Project/results/primary_core_niche_circuit_module_formal_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_named_comparator_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_model.png`
- `Project/results/primary_core_niche_circuit_module_model.md`

## Transcriptomic Configuration Model

The primary core now has a local named-comparator test of the transcriptomic assembly-plan hypothesis.

- Construction-over-niche balance: downstream neurite/synaptic rank minus neurogenic niche/progenitor rank.
- Regional fate polarity: branch-matched fate rank minus branch-opposed fate rank.
- Combined configuration score: construction-over-niche balance plus regional fate polarity.
- Result: combined configuration score is positive in 4/4 named granule-versus-comparator contrasts, but the small-n Wilcoxon p is 0.0625.

Main outputs:

- `Project/results/primary_core_transcriptomic_configuration_units.tsv`
- `Project/results/primary_core_transcriptomic_configuration_role_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_model.png`
- `Project/results/primary_core_transcriptomic_configuration_model.md`

## Primary-Core Configuration Validation

The configuration model has now been broadened across the primary-core pseudobulk layers.

- Configuration class units: 210 across all 10 primary datasets.
- Candidate-versus-background contrasts: 63 across 7 datasets.
- Positive combined configuration contrasts: 52/63.
- Median candidate-background configuration delta: 0.417.
- Wilcoxon p for configuration delta greater than zero: 4.89e-08.
- This is broader than the named-comparator audit, but less specific because not every dataset has explicit pyramidal/Purkinje labels.

Main outputs:

- `Project/results/primary_core_transcriptomic_configuration_primary_units.tsv.gz`
- `Project/results/primary_core_transcriptomic_configuration_primary_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_coverage.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.png`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`

## Configuration Driver Audit

The configuration evidence has been decomposed into downstream construction balance versus regional fate polarity.

- Total contrasts audited: 67.
- Configuration-positive: 56/67.
- Both components positive: 28.
- Fate-driven positive: 27.
- Construction-driven positive: 1.
- Primary-core candidate-background layer: 52/63 positive, split between 26 both-component and 26 fate-driven positives.

Interpretation: the transcriptomic assembly-plan model is supported, but it is identity-coupled. The safest claim is not a pure morphology-only transcriptomic code, but a regional-fate context plus selective downstream construction module balance.

Main outputs:

- `Project/results/primary_core_configuration_driver_audit_contrasts.tsv`
- `Project/results/primary_core_configuration_driver_audit_module_deltas.tsv`
- `Project/results/primary_core_configuration_driver_audit_summary.tsv`
- `Project/results/primary_core_configuration_driver_audit_gene_priorities.tsv`
- `Project/results/primary_core_configuration_driver_audit.png`
- `Project/results/primary_core_configuration_driver_audit.md`

## Output

- Table: `Project/results/integrated_primary_core_datasets.tsv`
