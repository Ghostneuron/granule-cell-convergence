# Suggested Run Order

This project was assembled as script-level workflows. Exact reruns require
public raw datasets or local selected matrices reconstructed from GEO,
Allen-related resources, NeuroMorpho and DANDI as documented in the
manifests.

1. Dataset discovery and metadata curation:
   `curate_secondary_and_human_candidates.py`,
   `integrate_candidate_resources.py`,
   `prioritize_datasets_for_next_phase.py`.
2. Human bridge/core construction and label validation:
   `inspect_human_seed_archives.py`,
   `curate_human_core_geo_sample_metadata.py`,
   `build_human_seed_sparse_objects.py`,
   `qc_harmonize_human_core_sparse_objects.py`,
   `build_human_core_normalized_reduced_object.py`,
   `tune_human_core_labels_and_test_modules.py`,
   `validate_human_core_marker_programs.py`.
3. Primary-core annotation, pseudobulk and ortholog rank-meta analyses:
   `classify_candidate_granule_cells.py`,
   `build_primary_core_*pseudobulk.py`,
   `build_primary_core_mgi_ortholog_*meta_model.py`,
   `build_primary_core_mgi_ortholog_formal_rank_model.py`.
4. Candidate tiering, comparator and configuration analyses:
   `build_primary_core_manuscript_candidate_packet.py`,
   `build_primary_core_granule_specificity_named_comparators.py`,
   `build_primary_core_transcriptomic_configuration_model.py`,
   `build_primary_core_transcriptomic_configuration_primary_validation.py`,
   `build_primary_core_configuration_driver_audit.py`.
5. Stage, pathway, ligand-receptor and conditioned-medium analyses:
   `build_primary_core_aim2_niche_pathway_model.py`,
   `build_primary_core_aim2b_stage_resolved_tgf_bdnf.py`,
   `fit_aim2_stage_pseudotime_model.py`,
   `build_aim2_sender_receiver_ligand_receptor.py`,
   `build_cerebellar_conditioned_medium_secretome_candidates.py`.
6. External validation and sparse-coding model:
   `build_neuromorpho_granule_morphometry_validation.py`,
   `prioritize_dandi_000003_targeted_downloads.py`,
   `build_dandi_000003_*`,
   `build_primary_core_aim3_sparse_coding_model.py`,
   `calibrate_aim3_empirical_sparse_model.py`.
7. Regulatory, hierarchical and hypothesis-support summaries:
   `build_epigenomic_extension_targets.py`,
   `build_gse322785_*`,
   `build_hierarchical_integrative_granule_model.py`,
   `build_hypothesis_support_score_matrix.py`,
   `build_fig7_hypothesis_support_figure.py`.
8. Manuscript table/figure assembly:
   `organize_supplementary_tables.py`,
   `assemble_manuscript_figures.py`,
   `build_development_submission_manuscript.py`,
   `build_development_formatted_submission_package.py`.
