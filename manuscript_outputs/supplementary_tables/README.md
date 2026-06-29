# Supplementary Tables

Ordered manuscript support tables copied from `Project/results`.

These files are copies, not replacements; the original analysis outputs remain in `Project/results`.

## Submission Strategy

The 103 numbered tables are retained for reproducibility, but they are triaged so reviewers do not have to treat every table as equally central.

- Tables S1-S40 and S102-S103 are the reader-facing primary supplementary tables tied to the main figures.
- Tables S41-S53 and S89-S93 are reader-facing supplementary controls tied to supplementary figures.
- Tables S54-S88 and S94-S101 are best treated as a machine-readable Supplementary Data archive and provenance packet.

## Packet Summary

| Tables | Packet | Tier | Topic | Recommended use |
|---|---|---|---|---|
| Tables S1-S2 | Reader-facing primary supplementary tables | Primary result | Dataset frame and curation | Directly supports Figure 1 and the primary-core dataset design. |
| Tables S3-S8 | Reader-facing primary supplementary tables | Primary result | Candidate discovery and tiering | Directly supports Figure 2 and the ortholog-aware candidate-tiering pipeline. |
| Tables S9-S18 | Reader-facing primary supplementary tables | Primary result | Comparator and configuration tests | Directly supports Figures 3 and 4, including named comparators and configuration scores. |
| Tables S19-S28 | Reader-facing primary supplementary tables | Primary result | Stage-window and signaling reinterpretation | Directly supports Figure 5 pathway, stage-window, conditioned-medium, and 2005-paper comparison analyses. |
| Tables S29-S37 | Reader-facing primary supplementary tables | Primary result | Morphology, physiology and sparse-coding calibration | Directly supports Figure 5 morphology/activity validation and model calibration. |
| Tables S38-S40 | Reader-facing primary supplementary tables | Primary result | Sender-receiver predictions | Directly supports Figure 6 ligand-receptor sender-receiver predictions. |
| Tables S41-S53 | Reader-facing supplementary controls | Supporting result | Stage split, granule specificity and developmental-origin controls | Supports Supplementary Figures S1-S4 and the argument against a recent common granule progenitor. |
| Tables S54-S88 | Supplementary Data archive | Archive/provenance | Regulatory-compatibility extension | Retain as machine-readable regulatory-resource, feature-manifest, and sensitivity-analysis support rather than as primary evidence. |
| Tables S89-S93 | Reader-facing supplementary controls | Supporting result | Hierarchical integrative evidence model | Supports Supplementary Figure S5 and summarizes the term-balanced evidence synthesis. |
| Tables S94-S100 | Supplementary Data archive | Archive/provenance | Public perturbation audit | Retain as hypothesis-generating perturbation-resource triage and module-shift support for Supplementary Figure S6. |
| Table S101 | Supplementary Data archive | Provenance map | Result-level evidence anchors | Maps each Results section to source analysis files and the corresponding supplementary tables. |
| Tables S102-S103 | Reader-facing primary supplementary tables | Primary result | Hypothesis-support matrix | Directly supports Figure 7 and the evidence-weighted H1/H2/H3 comparison. |

## Detailed Index


### Reader-facing primary supplementary tables


#### Dataset frame and curation

- Table S1: `Table_S1_primary_core_dataset_frame.tsv` - Primary-core dataset frame (Figure 1; dataset design; Primary result).
- Table S2: `Table_S2_supporting_dataset_tiers.tsv` - Supporting dataset tiers (Figure 1; dataset design; Primary result).

#### Candidate discovery and tiering

- Table S3: `Table_S3_candidate_gene_tiers.tsv` - Candidate gene tiers (Figure 2; Primary result).
- Table S4: `Table_S4_formal_rank_gene_summary.tsv` - Formal ortholog rank-meta gene summary (Figure 2; candidate discovery; Primary result).
- Table S5: `Table_S5_formal_rank_shared_hits.tsv` - Formal shared hits (Figure 2; candidate discovery; Primary result).
- Table S6: `Table_S6_formal_rank_mechanism_hits.tsv` - Formal mechanism-axis hits (Figure 2; Primary result).
- Table S7: `Table_S7_mechanism_axis_gene_table.tsv` - Mechanism-axis gene table (Figure 2; Primary result).
- Table S8: `Table_S8_mechanism_axis_summary.tsv` - Mechanism-axis summary (Figure 2; Primary result).

#### Comparator and configuration tests

- Table S9: `Table_S9_named_comparator_axis_summary.tsv` - Named-comparator mechanism-axis summary (Figure 3a; Primary result).
- Table S10: `Table_S10_named_comparator_axis_units.tsv` - Named-comparator mechanism-axis units (Figure 3a; Primary result).
- Table S11: `Table_S11_niche_circuit_gene_sets.tsv` - Niche/circuit module gene sets (Figure 3b; Primary result).
- Table S12: `Table_S12_niche_circuit_formal_summary.tsv` - Niche/circuit formal summary (Figure 3b; Primary result).
- Table S13: `Table_S13_niche_circuit_named_comparator_summary.tsv` - Niche/circuit named-comparator summary (Figure 3b; Primary result).
- Table S14: `Table_S14_configuration_local_contrasts.tsv` - Transcriptomic configuration local contrasts (Figure 4b; Primary result).
- Table S15: `Table_S15_configuration_primary_summary.tsv` - Primary-core configuration summary (Figure 4c; Primary result).
- Table S16: `Table_S16_configuration_primary_contrasts.tsv` - Primary-core configuration contrasts (Figure 4c; Primary result).
- Table S17: `Table_S17_configuration_driver_summary.tsv` - Configuration driver summary (Figure 4d; Primary result).
- Table S18: `Table_S18_configuration_driver_contrasts.tsv` - Configuration driver contrasts (Figure 4d; Primary result).

#### Stage-window and signaling reinterpretation

- Table S19: `Table_S19_aim2_pathway_gene_sets.tsv` - Pathway-readiness gene sets (Figure 5a; Primary result).
- Table S20: `Table_S20_aim2_pathway_summary.tsv` - Pathway-readiness summary (Figure 5a; Primary result).
- Table S21: `Table_S21_aim2_pathway_contrasts.tsv` - Pathway-readiness contrasts (Figure 5a; Primary result).
- Table S22: `Table_S22_aim2_ligand_receptor_summary.tsv` - Ligand-receptor readiness summary (Figure 5a; Primary result).
- Table S23: `Table_S23_tgf_bdnf_stage_summary.tsv` - Stage-resolved TGF-beta/BDNF summary (Figure 5b; Primary result).
- Table S24: `Table_S24_tgf_bdnf_stage_transitions.tsv` - Stage-resolved TGF-beta/BDNF transitions (Figure 5b; Primary result).
- Table S25: `Table_S25_stage_window_branch_summary.tsv` - Stage-window model branch summary (Figure 5b; Primary result).
- Table S26: `Table_S26_stage_window_coefficients.tsv` - Stage-window model coefficients (Figure 5b; Primary result).
- Table S27: `Table_S27_conditioned_medium_secretome_candidates.tsv` - Conditioned-medium secretome candidates (Figure 5; Discussion; Primary result).
- Table S28: `Table_S28_paper_2005_support_revision.tsv` - 2005 paper support/revision table (Figure 5; Discussion; Primary result).

#### Morphology, physiology and sparse-coding calibration

- Table S29: `Table_S29_neuromorpho_summary.tsv` - NeuroMorpho summary (Figure 5c; Primary result).
- Table S30: `Table_S30_neuromorpho_comparison.tsv` - NeuroMorpho comparison (Figure 5c; Primary result).
- Table S31: `Table_S31_dandi_session_summary.tsv` - DANDI multisession summary (Figure 5c; Primary result).
- Table S32: `Table_S32_dandi_unit_metrics.tsv` - DANDI spatial unit metrics (Figure 5c; Primary result).
- Table S33: `Table_S33_dandi_population_vector.tsv` - DANDI population-vector separation (Figure 5c; Primary result).
- Table S34: `Table_S34_sparse_coding_architecture_summary.tsv` - Sparse-coding architecture summary (Figure 5c; Primary result).
- Table S35: `Table_S35_sparse_coding_parameter_grid.tsv` - Sparse-coding parameter grid (Figure 5c; Primary result).
- Table S36: `Table_S36_empirical_calibration_architecture_summary.tsv` - Empirical calibration architecture summary (Figure 5c; Primary result).
- Table S37: `Table_S37_empirical_calibration_grid.tsv` - Empirical calibration grid (Figure 5c; Primary result).

#### Sender-receiver predictions

- Table S38: `Table_S38_sender_receiver_summary.tsv` - Sender-receiver summary (Figure 6; Primary result).
- Table S39: `Table_S39_sender_receiver_top_predictions.tsv` - Sender-receiver top predictions (Figure 6; Primary result).
- Table S40: `Table_S40_sender_receiver_predictions.tsv.gz` - Sender-receiver prediction units (Figure 6; Primary result).

### Reader-facing supplementary controls


#### Stage split, granule specificity and developmental-origin controls

- Table S41: `Table_S41_stage_split_module_branch_summary.tsv` - Stage-split module branch summary (Supplementary Figure S1; Supporting result).
- Table S42: `Table_S42_stage_split_similarity.tsv` - Stage-split similarity (Supplementary Figure S1; Supporting result).
- Table S43: `Table_S43_stage_split_transitions.tsv` - Stage-split transitions (Supplementary Figure S1; Supporting result).
- Table S44: `Table_S44_stage_split_group_calls.tsv` - Stage-split group calls (Supplementary Figure S1; Supporting result).
- Table S45: `Table_S45_granule_special_top_candidates.tsv` - Granule-special top candidates (Supplementary Figure S2; Supporting result).
- Table S46: `Table_S46_granule_special_summary.tsv` - Granule-special genome-wide summary (Supplementary Figure S2; Supporting result).
- Table S47: `Table_S47_regional_origin_timing_state_summary.tsv` - Regional-origin timing state summary (Supplementary Figure S3; Supporting result).
- Table S48: `Table_S48_regional_origin_timing_metrics.tsv` - Regional-origin timing metrics (Supplementary Figure S3; Supporting result).
- Table S49: `Table_S49_regional_origin_timing_gene_units.tsv` - Regional-origin timing gene units (Supplementary Figure S3; Supporting result).
- Table S50: `Table_S50_developmental_origin_gene_sets.tsv` - Developmental-origin gene sets (Supplementary Figure S4; Supporting result).
- Table S51: `Table_S51_developmental_origin_state_summary.tsv` - Developmental-origin state summary (Supplementary Figure S4; Supporting result).
- Table S52: `Table_S52_developmental_origin_branch_metrics.tsv` - Developmental-origin branch metrics (Supplementary Figure S4; Supporting result).
- Table S53: `Table_S53_developmental_origin_units.tsv` - Developmental-origin module units (Supplementary Figure S4; Supporting result).

### Supplementary Data archive


#### Regulatory-compatibility extension

- Table S54: `Table_S54_epigenomic_extension_candidate_resources.tsv` - Epigenomic extension candidate resources (Future epigenomic extension; Archive/provenance).
- Table S55: `Table_S55_epigenomic_regulatory_targets.tsv` - Epigenomic regulatory targets (Future epigenomic extension; Archive/provenance).
- Table S56: `Table_S56_integrative_model_term_specification.tsv` - Integrative model term specification (Future epigenomic extension; Archive/provenance).
- Table S57: `Table_S57_epigenomic_target_summary.tsv` - Epigenomic target summary (Future epigenomic extension; Archive/provenance).
- Table S58: `Table_S58_gse268609_epigenomic_peak_targets.tsv` - GSE268609 epigenomic peak targets (Future epigenomic extension; Archive/provenance).
- Table S59: `Table_S59_gse268609_epigenomic_peak_gene_summary.tsv` - GSE268609 epigenomic peak gene summary (Future epigenomic extension; Archive/provenance).
- Table S60: `Table_S60_gse268609_epigenomic_selective_extraction_manifest.tsv` - GSE268609 epigenomic selective extraction manifest (Future epigenomic extension; Archive/provenance).
- Table S61: `Table_S61_gse322785_cerebellar_multiome_sample_metadata.tsv` - GSE322785 cerebellar multiome sample metadata (Future epigenomic extension; Archive/provenance).
- Table S62: `Table_S62_gse322785_cerebellar_multiome_file_inventory.tsv` - GSE322785 cerebellar multiome file inventory (Future epigenomic extension; Archive/provenance).
- Table S63: `Table_S63_gse322785_cerebellar_multiome_donor_summary.tsv` - GSE322785 cerebellar multiome donor summary (Future epigenomic extension; Archive/provenance).
- Table S64: `Table_S64_gse322785_cerebellar_multiome_download_plan.tsv` - GSE322785 cerebellar multiome download plan (Future epigenomic extension; Archive/provenance).
- Table S65: `Table_S65_gse322785_human_h5_feature_inventory.tsv` - GSE322785 human H5 feature inventory (Future epigenomic extension; Archive/provenance).
- Table S66: `Table_S66_gse322785_human_h5_epigenomic_gene_summary.tsv` - GSE322785 human H5 epigenomic gene summary (Future epigenomic extension; Archive/provenance).
- Table S67: `Table_S67_gse322785_human_h5_epigenomic_peak_targets.tsv` - GSE322785 human H5 epigenomic peak targets (Future epigenomic extension; Archive/provenance).
- Table S68: `Table_S68_gse322785_human_h5_epigenomic_selective_manifest.tsv` - GSE322785 human H5 epigenomic selective manifest (Future epigenomic extension; Archive/provenance).
- Table S69: `Table_S69_gse322785_human_h5_selected_matrix_summary.tsv` - GSE322785 human H5 selected matrix summary (Future epigenomic extension; Archive/provenance).
- Table S70: `Table_S70_gse322785_human_h5_marker_panel_coverage.tsv` - GSE322785 human H5 marker panel coverage (Future epigenomic extension; Archive/provenance).
- Table S71: `Table_S71_gse322785_human_h5_marker_celltype_summary.tsv` - GSE322785 human H5 marker cell-type summary (Future epigenomic extension; Archive/provenance).
- Table S72: `Table_S72_gse322785_human_h5_marker_high_confidence_barcodes.tsv.gz` - GSE322785 human H5 high-confidence marker barcodes (Future epigenomic extension; Archive/provenance).
- Table S73: `Table_S73_gse322785_human_h5_epigenomic_marker_group_feature_scores.tsv.gz` - GSE322785 marker-group feature scores (Future epigenomic extension; Archive/provenance).
- Table S74: `Table_S74_gse322785_human_h5_epigenomic_marker_group_module_scores.tsv` - GSE322785 marker-group module scores (Future epigenomic extension; Archive/provenance).
- Table S75: `Table_S75_gse322785_human_h5_epigenomic_marker_group_contrasts.tsv` - GSE322785 marker-group epigenomic contrasts (Future epigenomic extension; Archive/provenance).
- Table S76: `Table_S76_gse322785_human_h5_cluster_validation_barcode_assignments.tsv.gz` - GSE322785 cluster-validation barcode assignments (Future epigenomic extension; Archive/provenance).
- Table S77: `Table_S77_gse322785_human_h5_cluster_validation_summary.tsv` - GSE322785 cluster-validation summary (Future epigenomic extension; Archive/provenance).
- Table S78: `Table_S78_gse322785_human_h5_cluster_validation_marker_call_enrichment.tsv` - GSE322785 cluster marker-call enrichment (Future epigenomic extension; Archive/provenance).
- Table S79: `Table_S79_gse322785_human_h5_cluster_validation_marker_support.tsv` - GSE322785 cluster marker support (Future epigenomic extension; Archive/provenance).
- Table S80: `Table_S80_gse322785_human_h5_cluster_validation_metrics.tsv` - GSE322785 cluster-validation metrics (Future epigenomic extension; Archive/provenance).
- Table S81: `Table_S81_gse322785_human_h5_cluster_supported_marker_rules.tsv` - GSE322785 cluster-supported marker rules (Future epigenomic extension; Archive/provenance).
- Table S82: `Table_S82_gse322785_human_h5_cluster_supported_marker_barcodes.tsv.gz` - GSE322785 cluster-supported marker barcodes (Future epigenomic extension; Archive/provenance).
- Table S83: `Table_S83_gse322785_human_h5_cluster_supported_epigenomic_feature_scores.tsv.gz` - GSE322785 cluster-supported feature scores (Future epigenomic extension; Archive/provenance).
- Table S84: `Table_S84_gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv` - GSE322785 cluster-supported module scores (Future epigenomic extension; Archive/provenance).
- Table S85: `Table_S85_gse322785_human_h5_cluster_supported_epigenomic_contrasts.tsv` - GSE322785 cluster-supported epigenomic contrasts (Future epigenomic extension; Archive/provenance).
- Table S86: `Table_S86_gse322785_epigenomic_broad_vs_cluster_supported_contrasts.tsv` - GSE322785 broad-versus-supported contrast comparison (Future epigenomic extension; Archive/provenance).
- Table S87: `Table_S87_gse322785_epigenomic_robust_positive_contrasts.tsv` - GSE322785 robust-positive epigenomic contrasts (Future epigenomic extension; Archive/provenance).
- Table S88: `Table_S88_gse322785_epigenomic_robust_contrast_summary.tsv` - GSE322785 robust contrast summary (Future epigenomic extension; Archive/provenance).

### Reader-facing supplementary controls


#### Hierarchical integrative evidence model

- Table S89: `Table_S89_hierarchical_integrative_model_terms.tsv` - Hierarchical integrative model terms (Integrative model; Supporting result).
- Table S90: `Table_S90_hierarchical_integrative_model_evidence_units.tsv` - Hierarchical integrative model evidence units (Integrative model; Supporting result).
- Table S91: `Table_S91_hierarchical_integrative_model_layer_summary.tsv` - Hierarchical integrative model layer summary (Integrative model; Supporting result).
- Table S92: `Table_S92_hierarchical_integrative_model_branch_summary.tsv` - Hierarchical integrative model branch summary (Integrative model; Supporting result).
- Table S93: `Table_S93_hierarchical_integrative_model_component_scores.tsv` - Hierarchical integrative model component scores (Integrative model; Supporting result).

### Supplementary Data archive


#### Public perturbation audit

- Table S94: `Table_S94_causal_perturbation_dataset_triage.tsv` - Public perturbation dataset triage (Causal follow-up resources; Archive/provenance).
- Table S95: `Table_S95_causal_perturbation_node_summary.tsv` - Public perturbation node summary (Causal follow-up resources; Archive/provenance).
- Table S96: `Table_S96_causal_perturbation_processing_status.tsv` - Public perturbation processing status (Causal follow-up resources; Archive/provenance).
- Table S97: `Table_S97_causal_perturbation_module_catalog.tsv` - Public perturbation module catalog (Causal follow-up resources; Archive/provenance).
- Table S98: `Table_S98_causal_perturbation_module_shift_gene_effects.tsv` - Public perturbation gene-level module effects (Causal follow-up resources; Archive/provenance).
- Table S99: `Table_S99_causal_perturbation_module_shift_summary.tsv` - Public perturbation contrast-module shifts (Causal follow-up resources; Archive/provenance).
- Table S100: `Table_S100_causal_perturbation_module_shift_node_summary.tsv` - Public perturbation node-module shifts (Causal follow-up resources; Archive/provenance).

#### Result-level evidence anchors

- Table S101: `Table_S101_result_evidence_anchor_provenance.tsv` - Result-level evidence anchor provenance (Analysis provenance; Provenance map).

### Reader-facing primary supplementary tables


#### Hypothesis-support matrix

- Table S102: `Table_S102_hypothesis_support_score_matrix.tsv` - Hypothesis support coefficient matrix (Figure 7; Primary result).
- Table S103: `Table_S103_hypothesis_support_scores.tsv` - Hypothesis support scores (Figure 7; Primary result).

## Prebuilt Upload Packets

- `Supplementary_Tables_Primary_Reader_Facing_S1_S40_S102_S103.zip`: main-figure, reader-facing result tables.
- `Supplementary_Tables_Supplementary_Controls_S41_S53_S89_S93.zip`: supplementary-figure control tables.
- `Supplementary_Data_Archive_S54_S88_S94_S101.zip`: machine-readable regulatory, perturbation, audit and provenance archive.

The ZIP files are convenience copies only. The individual numbered tables remain the authoritative files.

Machine-readable table index: `Supplementary_Table_Index.tsv`
Submission grouping guide: `Supplementary_Table_Submission_Guide.tsv`
