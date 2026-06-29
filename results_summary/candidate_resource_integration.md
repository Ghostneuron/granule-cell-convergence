# Candidate resource integration

Date checked: 2026-06-23

## Decision

The newly inspected resources are useful, and the human dentate/hippocampal core expansion now includes a direct adult dentate anchor plus a broader human hippocampal aging/AD RNA branch from `GSE268609`. These human bridge objects have also been integrated with the refined dentate/cerebellar backbone at a rank-summary level.

## Resource groups

### human_dentate_core_construction
- `GSE185277`: human imGC marker scaffold and lifespan neurogenesis framing (core_build_first_human_dentate_reference; raw_downloaded;sparse_objects_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests).
- `GSE185553`: broader human hippocampal context for GSE185277 (core_build_first_human_hippocampus_reference; raw_downloaded;sparse_objects_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests).
- `GSE186538`: human DG taxonomy and cross-species reference (core_build_first_human_dg_taxonomy; human_files_downloaded;sparse_subset_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests).
- `GSE325391`: primary modern adult human dentate comparator (core_build_first_human_dentate_primary; adult_rds_downloaded;geo_metadata_curated;rds_inspected;selected_sparse_bridge_built;human_core_label_projected).
- `GSE268609`: broader human hippocampal RNA/multiome expansion (core_build_first_human_hippocampus_multiome; source_listing_saved;geo_metadata_curated;rna_matrix_barcodes_features_downloaded;rna_selected_npz_built;human_core_label_projected;diagnosis_module_tests_complete;atac_and_full_rds_deferred).
- `primary_core_consensus_candidate_dataset_validation`: dataset-aware validation of the 24 cross-screen consensus mechanism candidates, with six genes robust across all available screen/branch tests (`GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, `RFX3`).
- `primary_core_mgi_ortholog_meta_model`: MGI-filtered ortholog-aware meta-model over selected-feature and full-matrix pseudobulk screens, with 16,245 strict same-symbol one-to-one pairs, 1,304 shared hits, and 36 mechanism-prioritized shared hits.
- `primary_core_mgi_ortholog_expanded_meta_model`: expanded MGI one-to-one meta-model including non-identical mouse/human symbols, with 17,611 one-to-one targets, 1,370 shared hits, 64 non-identical-symbol shared hits, and 36 mechanism-prioritized shared hits.
- `primary_core_mgi_ortholog_formal_rank_model`: formal dataset-level rank-meta validation of the expanded MGI one-to-one ortholog model, with 116,013 dataset-level deltas, 36,303 branch tests, 1,370 formal shared hits, and 36 mechanism-prioritized genes modeled.
- `primary_core_manuscript_candidate_tiers`: concise manuscript-facing candidate packet with 66 tiered genes, including 6 Tier 1 core convergent-program genes and 9 Tier 2 high-confidence wiring/synaptic/executor genes.
- `primary_core_mechanism_axis_model`: biological interpretation layer that organizes the formal candidate tiers into developmental regulatory control, neurite/cytoskeleton morphogenesis, axon guidance/adhesion, synaptic/excitability maturation, and exploratory ortholog-completeness axes.
- `primary_core_granule_specificity_named_comparators`: named-comparator specificity audit that scores Tier 1-4 mechanism axes in `GSE104323` dentate granule-lineage versus pyramidal labels and `GSE122357` cerebellar granule-lineage versus Purkinje labels.
- `primary_core_niche_circuit_module_model`: hypothesis-driven module model separating cerebellar fate, dentate fate, shared neurogenic niche/progenitor state, downstream neurite/morphology, and downstream synaptic/excitability programs.
- `primary_core_transcriptomic_configuration_model`: module-balance test of the transcriptomic assembly-plan hypothesis, asking whether construction-over-niche balance plus regional fate polarity separates granule cells from named pyramidal/Purkinje comparators.
- `primary_core_transcriptomic_configuration_primary_validation`: broader primary-core validation of the transcriptomic assembly-plan score across full MGI and selected-feature pseudobulk layers.
- `primary_core_configuration_driver_audit`: claim-safety audit that decomposes the configuration score into downstream construction balance and regional fate polarity.
- `manuscript_planning_packet`: current manuscript scaffold with title options, central claim, abstract draft, result spine, figure plan, and claim/evidence/caveat table.
- `specific_aims_completion_audit`: aim-by-aim status check against the starting three-aim design, marking all three aims as first-pass computationally complete with optional validation upgrades.
- `primary_core_aim2_niche_pathway_model`: targeted Aim 2 pathway-readiness and ligand-receptor audit for TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP, Reelin, Semaphorin, SHH, WNT, FGF, and Notch pathways.
- `primary_core_aim2b_stage_resolved_tgf_bdnf`: stage-aware Aim 2 refinement testing TGF-beta/BDNF and stop/permissive module readiness across adult dentate lineage states, postnatal dentate stages, postnatal cerebellar candidates, and adult dentate activity/maturation states.
- `cerebellar_conditioned_medium_secretome_candidates`: sequencing-derived secreted/ligand candidate screen for possible cerebellar conditioned-medium anti-proliferative or pro-differentiation factors beyond the 2005 TGF-beta2/BDNF anchors.
- `primary_core_2005_endpoint_pseudotime_audit`: stage/module pseudotime audit of the 2005 paper endpoints, including proliferation, p21/p27, differentiation, TGF/SMAD, BDNF/ERK, apoptosis/survival, and secreted stop-factor modules.
- `primary_core_2005_endpoint_graph_pseudotime`: first cell-level graph pseudotime overlay for the 2005 paper endpoint modules across GSE104323, GSE292261, and GSE122357.
- `primary_core_full_transcriptome_diffusion`: full primary-core HVG diffusion/pseudotime layer across all 10 strict datasets, with 56,892 cells/nuclei in trajectory tables and a Fig1-5 impact audit.
- `primary_core_aim3_sparse_coding_model`: targeted Aim 3 sparse expansion-coding model that links expansion, sparse input sampling, and output sparsity to useful pattern separation and maps transcriptomic modules onto these computational parameters.
- `secondary_phys_morph_validation_layer`: supporting physiology/morphology validation layer for Aim 3, integrating Allen Cell Types as comparator ephys calibration, DANDI 000003 as the main dentate activity resource, DANDI 000165 as hippocampal network support, and NeuroMorpho as direct granule-cell morphometry evidence.
- `neuromorpho_granule_morphometry_validation`: first direct morphology validation, comparing 558 dentate granule reconstructions with all 62 cerebellar granule reconstructions available under the strict NeuroMorpho query; branch counts are close, while primary stems and dendritic field scale differ strongly.
- `dandi_000003_activity_sparsity_pilot`: first direct activity pilot from DANDI 000003, downloading the smallest 4.66 GB NWB session and extracting unit/cell-type/state firing summaries; three source-labeled granule units show median firing rate 0.9131 Hz and median active 1 s-bin fraction 0.3943.
- `dandi_000003_spatial_pattern_pilot`: first position-linked DANDI 000003 pilot, using awake-moving position samples from the downloaded NWB session; three source-labeled granule units show median spatial information 0.3098 bits/spike, spatial sparsity 0.6273, active spatial-bin fraction 0.4634, and a feasible but underpowered population-vector check.
- `dandi_000003_multisession_spatial_extension`: six-session expansion of the position-linked DANDI workflow, now combining five legacy local NWB files with the external-drive `YutaMouse55-160907` download. It analyzes 124 units across `YutaMouse41-150829`, `YutaMouse37-150609`, `YutaMouse42-151102`, `YutaMouse55-160908`, `YutaMouse55-160909`, and `YutaMouse55-160907`. Direct granule-labeled evidence comes from all six sessions, totaling 26 granule units with pooled median spatial information 0.7800 bits/spike, active spatial-bin fraction 0.5489, and granule population-vector far-minus-near Euclidean separation of 0.1842, 1.2293, 0.0116, 0.0615, 0.0280, and 0.2260.
- `dandi_000003_targeted_download_priority`: yield/breadth download-priority layer for DANDI 000003. After the `YutaMouse55-160907` external-drive follow-up, the next yield-first candidate is `YutaMouse55-160911` (9.01 GB), while the smallest size-ranked new-subject option is `YutaMouse44-151128` (7.04 GB). The previously downloaded/deleted `YutaMouse37-150617` file is excluded because it had 0 source-labeled granule units.
- `manuscript_figure_assembly`: revised manuscript figure packet, with six composite figures exported as PNG/PDF plus a manifest and assembly report; Fig5 now includes the completed stage-windowed pseudotime layer and Fig6 integrates the Fig1/2/4/5 flow charts into one working model.

### nature_neuroscience_2025_mouse_dg_aging
- `GSE233363`: mouse DG aging, neurogenic-lineage maturation, niche inflammation, and spatial validation (core_validation_after_human_dentate; source_listing_saved;not_downloaded).

### cell_reports_2025_human_dg_spatial_lifespan
- `spatial_DG_lifespan`: human DG spatial/lifespan validation for granule-cell maturation, ECM, and neuroinflammation modules (supporting_spatial_context; article_pdf_present;zenodo_metadata_saved;data_not_downloaded).

### nature_2022_human_hippocampal_immature_neurons
- `GSE198323`: AD-related human immature-neuron disease context; keep separate from healthy/reference analyses (supporting_disease_reference; source_listing_saved;not_downloaded).

## Practical use

1. Use the tuned human-core labels and dataset-aware module tests as the convention for the human dentate/hippocampal branch.
2. Treat the `GSE186538` DG GC subset as the marker-validated human DG taxonomy anchor.
3. Treat `GSE325391` as the primary modern adult human dentate anchor and `GSE268609` as the broader human hippocampal aging/AD RNA expansion with projected, not source-taxonomy, labels.
4. Use the integrated rank layer (`human_bridge_backbone_rank_*`) as the current comparative bridge for manuscript-scale figure planning.
5. Use `human_bridge_candidate_gene_packet_*` as the first mechanism shortlist for regional identity versus shared structural-executor genes.
6. Use `primary_core_ortholog_module_*` as the first formal marker-panel module result for the strict 10-dataset core.
7. Use `primary_core_candidate_gene_pseudobulk_*` as the first DE-adjacent candidate-gene evidence layer for shared structural executors versus regional identity/wiring genes.
8. Use `primary_core_expanded_gene_pseudobulk_*` and `primary_core_expanded_gene_mechanism_triage.*` as the broad selected-gene discovery layer before whole-transcriptome DE.
9. Use `primary_core_genomewide_symbol_pseudobulk_*`, `primary_core_genomewide_symbol_mechanism_triage.*`, and `primary_core_cross_screen_mechanism_consensus.*` as the current full-matrix same-symbol discovery and consensus layer.
10. Use `primary_core_consensus_candidate_dataset_validation.*` as the current strongest dataset-aware mechanism shortlist for figure planning.
11. Use `primary_core_mgi_ortholog_meta_model.*` as the current strongest ortholog-aware evidence layer, prioritizing the six dataset-robust consensus genes plus `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`, and `ROBO2`.
12. Use `primary_core_mgi_ortholog_full_matrix_*` and `primary_core_mgi_ortholog_expanded_meta_model.*` as the complete rank-based MGI ortholog evidence layer.
13. Use `primary_core_mgi_ortholog_formal_rank_model.*` as the current strongest statistical evidence layer for manuscript mechanism claims; it retains the six-gene seed set (`GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`) and formalizes the broader shared-hit table.
14. Use `primary_core_manuscript_candidate_tiers.*` for Results drafting and figure planning.
15. Use `primary_core_mechanism_axis_model.*` as the current biological model layer for the paper argument.
16. Use `primary_core_granule_specificity_named_comparators.*` to temper uniqueness claims: the current axes support shared/convergent granule-cell programs, but they do not behave as uniquely granule-specific pathways versus the named pyramidal and Purkinje comparators tested here.
17. Use `primary_core_niche_circuit_module_model.*` as the current answer to the niche/circuit-constraint question: upstream fate/niche modules are branch-specific or mixed, whereas downstream neurite/morphology and synaptic/excitability modules show the strongest formal convergence.
18. Use `primary_core_transcriptomic_configuration_model.*` as the current transcriptomic-configuration layer: the combined module-balance score is positive in all four available named granule-versus-comparator contrasts, but it remains a small-n local test.
19. Use `primary_core_transcriptomic_configuration_primary_validation.*` as the broader primary-core validation layer: 52/63 candidate-background contrasts are positive, median delta 0.417, Wilcoxon p=4.89e-08.
20. Use `primary_core_configuration_driver_audit.*` to phrase the configuration result carefully: the broad signal is strong but identity-coupled, with selective downstream construction-balance support.
21. Use `manuscript_planning_packet.*`, `manuscript_claim_evidence_caveat_table.tsv`, and `manuscript_figure_plan.tsv` as the current drafting scaffold and claim-strength guardrail.
22. Use `primary_core_aim2_niche_pathway_model.*` as the Aim 2 completion layer: it supports context-dependent niche/maturation signaling, with dentate TGF-beta/BDNF enrichment and cerebellar SHH enrichment, but not a simple cerebellar-biased TGF-beta/BDNF stop-signaling program.
23. Use `primary_core_aim2b_stage_tgf_bdnf.*` as the Aim 2 stage refinement: it supports TGF-beta/BDNF as a timing/state-sensitive maturation-readiness axis rather than a universal granule-cell stop switch.
24. Use `cerebellar_conditioned_medium_secretome_candidates.*` as a validation-prioritization list for the 2005 conditioned-medium question; top inferred classes are WNT antagonists, BMP/GDF/activin factors, Reelin/guidance cues, and matricellular/context factors.
25. Use `primary_core_2005_endpoint_pseudotime_audit.*` and `primary_core_2005_endpoint_graph_pseudotime.*` as the first trajectory evidence for the 2005 paper endpoints; the key conclusion is stage/trajectory dependence, not simple monotonic age dependence.
26. Use `primary_core_full_transcriptome_diffusion.*` as the stronger primary-core trajectory refinement: ordered anchors validate the pseudotime direction, Fig1-4 are refined rather than overturned, and Fig5 should gain an explicit pseudotime/stage-window layer.
27. Use `primary_core_aim3_sparse_coding_model.*` as the Aim 3 completion layer: raw computational score can favor dense expansion, but resource-constrained nontrivial expansion favors sparse granule-like designs.
28. Use `secondary_phys_morph_validation_layer.*` as the concrete Aim 3 validation-resource map: NeuroMorpho can quantify dendrite/stem/branch morphology, DANDI 000003 can support dentate firing sparsity and behavior-linked separation metrics, and Allen Cell Types should be used as intrinsic-ephys comparator calibration rather than direct granule-cell primary evidence.
29. Use `neuromorpho_granule_morphometry_validation.*` as the first direct morphology validation: the data support convergent compact/limited-branch architecture, not literal identical dendritic geometry. The Aim 3 model should split morphology into primary stems/claws and dendritic-field complexity.
30. Use `dandi_000003_activity_sparsity_pilot.*` as the first direct dentate activity validation: it establishes feasible NWB unit extraction and labeled granule-cell sparsity summaries, but pattern-separation claims still need position/task-specific analysis across more sessions.
31. Use `dandi_000003_spatial_pattern_pilot.*` as the first position-linked feasibility layer: it supports spatial coding analysis in DANDI, with granule-cell spatial information and active spatial-bin metrics, but pattern-separation claims still require multi-session task/trajectory comparisons.
32. Use `dandi_000003_multisession_spatial_extension.*` as the current stronger DANDI spatial layer: it extends the pilot to six analyzed sessions, with 26 granule-labeled units across six granule-containing sessions, provides the reusable download plan, and should be expanded further only for a figure-specific or review-specific need.
33. Use `specific_aims_completion_audit.md` to keep the original project design honest: all three aims now have first-pass computational completion, with optional spatial sender-receiver and direct morphology/connectomics/activity data extraction upgrades.
34. Use `manuscript_figures/*` as the current complete draft figure packet for manuscript writing; Fig5 has been revised around stage-windowed pseudotime, Fig6 now provides the integrated working model, and journal sizing/typography still need polish before submission.
35. Move next into manuscript text drafting and optional external validation with `GSE233363`, `spatial_DG_lifespan`, and `GSE198323`, plus raw-count/object-level validation if stronger claim support is needed.
