# Revised core and supporting dataset tiers

Date revised: 2026-06-22

## Current Local Core Datasets

These are the datasets that are already local and suitable for object-level analysis:

- Dentate gyrus: `GSE104323`, `GSE95752`, `GSE292261`, `GSE214309`
- Cerebellum: `GSE122357`, `GSE165657`, `GSE312658`

`GSE292261` and `GSE214309` are no longer merely secondary. They are promoted to core current validation datasets because their annotation/gene-identifier issues are manageable and their biology is directly aligned with the project.

## Constructed Human Dentate/Hippocampal Core Extension

The earlier local core was asymmetric: it had a human cerebellum anchor (`GSE165657`) but no constructed human dentate or hippocampal granule-cell counterpart. That gap is now corrected at the bridge-analysis level. The current human dentate/hippocampal extension has sparse objects, QC, marker triage, GEO specimen/donor metadata, a normalized reduced object, tuned labels, first-pass dataset-aware module tests, and rank-level integration with the refined dentate/cerebellar backbone.

Current staged build order and status:

1. `GSE185277`: small, fast human hippocampal/dentate imGC scaffold; raw archive downloaded, sparse libraries built, QC-harmonized, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
2. `GSE185553`: companion broader human hippocampal reference; raw archive downloaded, sparse libraries built, QC-harmonized, marker-scored, pooled specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
3. `GSE186538`: human DG taxonomy and cross-species reference; human files downloaded and 32,067 DG granule-cell candidates extracted as a sparse, marker-validated, donor-curated subset and included as the tuned normalized DG anchor.
4. `GSE325391`: primary modern adult human dentate acquisition; adult RDS downloaded, inspected, converted to a selected sparse bridge, and projected into the tuned human-core label convention.
5. `GSE268609`: broader human hippocampal RNA/multiome branch; RNA matrix/barcodes/features are now downloaded, selected-gene bridged, and projected into human-core labels, while ATAC/full RDS remain deferred.

## Integrated Core Evidence

The refined local backbone and the human dentate/hippocampal extension now have a shared rank-level summary layer:

- `Project/results/integrated_primary_core_datasets.tsv`: strict 10-dataset primary core table.
- `Project/results/integrated_primary_core_datasets.md`: readable summary of the primary core definition.

- `Project/results/human_bridge_backbone_rank_units.tsv`: integrated unit table.
- `Project/results/human_bridge_backbone_rank_statistics.tsv`: identity-separation and structural-rank tests.
- `Project/results/human_bridge_backbone_rank_source_summary.tsv`: source-layer medians for project narration and figure planning.
- `Project/results/human_bridge_backbone_rank_units.png`: current integrated rank plot.

This analysis supports the core project premise: dentate and cerebellar candidate cells remain distinct by regional identity rank, while both show elevated structural/morphogenesis-related rank signal. The result is strong enough to justify manuscript-style figure building, with the caveat that a full cross-species harmonized object would still be needed for final gene-level modeling.

The first human-bridge candidate-gene packet is also built:

- `Project/results/human_bridge_marker_gene_summary.tsv`: gene/source/class summaries for 67 refined panel genes.
- `Project/results/human_bridge_candidate_gene_packet.tsv`: first regional-identity versus shared-structural-executor shortlist.
- `Project/results/human_bridge_structural_executor_candidates.png`: plot of top structural/executor candidates.
- `Project/results/human_bridge_candidate_gene_packet.md`: interpretation and caveats.

The strict 10-dataset core now also has a first formal marker-panel ortholog-aware module analysis:

- `Project/results/primary_core_marker_panel_ortholog_map.tsv`: mouse-human marker-panel symbol map.
- `Project/results/primary_core_integrated_module_units.tsv`: strict primary-core module units.
- `Project/results/primary_core_ortholog_module_statistics.tsv`: primary module tests.
- `Project/results/primary_core_ortholog_module_leave_one_dataset_out.tsv`: robustness checks.
- `Project/results/primary_core_identity_structural_modules.png`: primary-core identity/structural module plot.
- `Project/results/primary_core_ortholog_module_analysis.md`: interpretation and next step.

The strict 10-dataset core also now has a focused candidate-gene pseudobulk screen:

- `Project/results/primary_core_candidate_gene_pseudobulk_expression.tsv`: 16,719 pseudobulk expression rows.
- `Project/results/primary_core_candidate_gene_pseudobulk_coverage.tsv`: source coverage for the 67 candidate genes.
- `Project/results/primary_core_candidate_gene_pseudobulk_statistics.tsv`: dentate/cerebellar candidate-versus-background tests.
- `Project/results/primary_core_candidate_gene_pseudobulk_hits.tsv`: ranked structural-executor support table.
- `Project/results/primary_core_candidate_gene_pseudobulk_effects.png`: shared structural-executor effect plot.
- `Project/results/primary_core_candidate_gene_pseudobulk_analysis.md`: interpretation and next step.

The strict 10-dataset core now also has an expanded selected-gene pseudobulk discovery and triage layer:

- `Project/results/primary_core_expanded_gene_pseudobulk_expression.tsv.gz`: 537,070 pseudobulk expression rows.
- `Project/results/primary_core_expanded_gene_pseudobulk_coverage.tsv`: source coverage for 2,169 selected genes.
- `Project/results/primary_core_expanded_gene_pseudobulk_statistics.tsv`: broad gene statistics.
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.tsv`: 1,279 shared-positive rank hits.
- `Project/results/primary_core_expanded_gene_pseudobulk_branch_specific.tsv`: branch-biased genes.
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.png`: expanded shared-hit plot.
- `Project/results/primary_core_expanded_gene_mechanism_triage.tsv`: mechanism-class triage.
- `Project/results/primary_core_expanded_gene_mechanism_triage.md`: interpretation.

The strict primary core now also has a full-matrix same-symbol pseudobulk and consensus layer:

- `Project/results/primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz`: 628,339 pseudobulk expression rows.
- `Project/results/primary_core_genomewide_symbol_pseudobulk_coverage.tsv`: full-matrix symbol coverage.
- `Project/results/primary_core_genomewide_symbol_pseudobulk_statistics.tsv`: 21,253 same-symbol gene statistics.
- `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.tsv`: 6,440 shared-positive hits.
- `Project/results/primary_core_genomewide_symbol_pseudobulk_branch_specific.tsv`: branch-biased same-symbol genes.
- `Project/results/primary_core_genomewide_symbol_mechanism_triage.tsv`: full-matrix mechanism triage.
- `Project/results/primary_core_cross_screen_mechanism_consensus.tsv`: consensus between selected-feature and full-matrix screens.
- `Project/results/primary_core_cross_screen_mechanism_consensus.md`: interpretation.

The strict primary core now also has a dataset-aware validation layer for the cross-screen consensus candidates:

- `Project/results/primary_core_consensus_candidate_dataset_deltas.tsv`: 1,608 dataset/sample/gene branch-delta rows.
- `Project/results/primary_core_consensus_candidate_dataset_validation.tsv`: robustness summary for 24 consensus candidates across selected-feature/full-matrix screens and dentate/cerebellar branches.
- `Project/results/primary_core_consensus_candidate_dataset_validation_heatmap.png`: visual summary of dataset-aware branch support.
- `Project/results/primary_core_consensus_candidate_dataset_validation.md`: interpretation.

The strongest current figure-level mechanism candidates are the six genes robust across all available screen/branch tests: `GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, and `RFX3`. The next tier robust in 3/4 tests is `CACNA2D1`, `GABRB3`, `GRIN2B`, `KCND2`, `KCNJ3`, `KCNJ6`, `PPP3CA`, and `STXBP5L`.

The strict primary core now also has an MGI-filtered ortholog-aware meta-model:

- `External_Data/Orthology/HOM_MouseHumanSequence.rpt`: official MGI human-mouse homology report.
- `Project/results/primary_core_mgi_ortholog_meta_model_map.tsv`: 17,611 one-to-one human-mouse pairs and 16,245 strict same-symbol one-to-one pairs.
- `Project/results/primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz`: 241,690 dataset/sample/gene branch-delta rows.
- `Project/results/primary_core_mgi_ortholog_meta_model_branch_summary.tsv`: 33,829 branch summary rows.
- `Project/results/primary_core_mgi_ortholog_meta_model_gene_summary.tsv`: 15,345 gene summary rows.
- `Project/results/primary_core_mgi_ortholog_meta_model_shared_hits.tsv`: 1,304 shared ortholog-filtered hits.
- `Project/results/primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv`: 36 mechanism-prioritized shared hits.
- `Project/results/primary_core_mgi_ortholog_meta_model.md`: interpretation and caveats.

This is now the strongest ortholog-aware evidence layer, with the caveat that non-identical one-to-one orthologs still need a mouse-symbol-aware extraction pass.

The non-identical ortholog extraction pass is now complete:

- `Project/results/primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`: 508,917 full-matrix ortholog pseudobulk rows.
- `Project/results/primary_core_mgi_ortholog_full_matrix_statistics.tsv`: 16,704 gene statistics.
- `Project/results/primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv`: 246 shared-positive non-identical-symbol hits.
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv`: 16,708 expanded meta-model gene summaries.
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_shared_hits.tsv`: 1,370 shared hits.
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv`: 36 mechanism-prioritized hits.
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model.md`: current interpretation.

The formal rank-meta validation pass is now also complete:

- `Project/results/primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz`: 116,013 dataset-level deltas.
- `Project/results/primary_core_mgi_ortholog_formal_rank_branch_tests.tsv`: 36,303 formal branch tests.
- `Project/results/primary_core_mgi_ortholog_formal_rank_gene_summary.tsv`: 16,708 gene summaries.
- `Project/results/primary_core_mgi_ortholog_formal_rank_shared_hits.tsv`: 1,370 formal shared hits.
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv`: 36 mechanism-prioritized genes with formal rank-meta and model summaries.
- `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`: current interpretation.

The current caveat is now narrower: the rank-based formal MGI layer is complete, but raw-count/object-level DE and external validation would still strengthen final claim-level support.

The formal candidate tiers have also been translated into a mechanism-axis model:

- `Project/results/primary_core_manuscript_candidate_tiers.tsv`: 66 tiered manuscript candidate genes.
- `Project/results/primary_core_mechanism_axis_gene_table.tsv`: candidate genes assigned to mechanism axes.
- `Project/results/primary_core_mechanism_axis_summary.tsv`: axis-level counts and tier composition.
- `Project/results/primary_core_mechanism_axis_branch_summary.tsv`: branch support by axis.
- `Project/results/primary_core_mechanism_axis_model.md`: biological interpretation.

This is the current paper-argument layer: the central claim should be a layered shared toolkit across developmental regulatory control, neurite/cytoskeleton morphogenesis, axon guidance/adhesion, and synaptic/excitability maturation.

The named-comparator specificity audit now adds an important boundary condition:

- `Project/results/primary_core_granule_specificity_named_comparator_axis_summary.tsv`: axis-level comparison against `GSE104323` pyramidal labels and `GSE122357` Purkinje labels.
- `Project/results/primary_core_granule_specificity_named_comparators.md`: interpretation and outputs.

The audit supports convergence, not pathway uniqueness: developmental regulatory control is cerebellar granule-enriched versus Purkinje cells but not dentate granule-enriched versus pyramidal comparators, and the other three Tier 1-4 mechanism axes do not pass strict named-comparator specificity.

The niche/circuit module model now refines the causal interpretation:

- `Project/results/primary_core_niche_circuit_module_formal_summary.tsv`: formal comparison of cerebellar fate, dentate fate, shared neurogenic niche, downstream neurite/morphology, and downstream synaptic/excitability modules.
- `Project/results/primary_core_niche_circuit_module_model.md`: interpretation and outputs.

The model supports distinct upstream fate logic with stronger downstream convergence: median formal convergence is 0.500 for downstream circuit/morphology modules versus -0.500 for upstream/niche modules.

The transcriptomic configuration model adds a module-balance test:

- `Project/results/primary_core_transcriptomic_configuration_contrasts.tsv`: four local granule-versus-comparator contrasts.
- `Project/results/primary_core_transcriptomic_configuration_model.md`: interpretation and outputs.

The combined configuration score is positive in all four named contrasts, supporting the idea that transcriptomic configuration may encode morphology better than unique pathway membership, although the current local test remains small-n.

The broader primary-core configuration validation strengthens this layer:

- `Project/results/primary_core_transcriptomic_configuration_primary_summary.tsv`: 52/63 positive candidate-background contrasts across 7 contrastable datasets.
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`: interpretation and outputs.

This broader layer supports the assembly-plan model across the core, while retaining the caveat that most datasets use local background classes rather than explicit pyramidal/Purkinje labels.

The configuration driver audit adds the claim-safety layer:

- `Project/results/primary_core_configuration_driver_audit_summary.tsv`: driver-class counts and module-delta summaries.
- `Project/results/primary_core_configuration_driver_audit_gene_priorities.tsv`: top downstream genes supporting the assembly-plan layer.
- `Project/results/primary_core_configuration_driver_audit.md`: interpretation and outputs.

The audit shows that the configuration signal is strong but identity-coupled: 56/67 contrasts are configuration-positive, with 28 both-component positives, 27 fate-driven positives, and only 1 construction-driven positive.

`GSE233363` remains highly useful, but it is a mouse DG aging/neurogenic-lineage validation dataset. It should follow the human dentate build rather than substitute for it.

## Supporting And Context Datasets

Use these for validation, disease/spatial context, or marker refinement, not as the main discovery backbone:

- `GSE214905`: small dentate patch-seq physiology/projection validation.
- `GSE242688`: spatial/proteomics-linked validation; analyze separately from cell-level statistics.
- `spatial_DG_lifespan`: human DG lifespan spatial context from the Cell Reports 2025 paper already in `Literature`.
- `GSE198323`: Nature 2022 human hippocampal immature-neuron AD/disease context.
- `GSE216877`: epilepsy disease-context validation.
- `GSE317381`: human DG spatial annotation-transfer context.
- `GSE150153`: organoid context.
- Allen WMB taxonomy: reference atlas for labels and marker refinement.

## Practical Rule

The main manuscript claim should not rely on a human cerebellum dataset without a human dentate/hippocampal comparator. The staged human dentate branch now exists as a normalized reduced object with tuned human DG/imGC, broader hippocampal, background labels, a direct adult dentate anchor from `GSE325391`, a broader aging/AD hippocampal RNA branch from `GSE268609`, an integrated rank comparison to the dentate/cerebellar backbone, a dataset-aware consensus-candidate validation layer, an expanded MGI-filtered ortholog meta-model, and a formal rank-meta validation layer. Spatial, patch-seq, disease, and organoid resources should support or qualify that claim unless they are converted into harmonized objects with appropriate dataset-specific statistics.
