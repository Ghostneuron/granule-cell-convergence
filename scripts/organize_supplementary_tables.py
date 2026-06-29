#!/usr/bin/env python3
"""Assemble ordered supplementary tables for the manuscript packet."""

from __future__ import annotations

import csv
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
OUT_DIR = ROOT / "Project/manuscript/Supplementary tables"
INDEX = OUT_DIR / "Supplementary_Table_Index.tsv"
GUIDE = OUT_DIR / "Supplementary_Table_Submission_Guide.tsv"
README = OUT_DIR / "README.md"
GENERATED_PROVENANCE_SOURCE = "generated:result_evidence_anchor_provenance.tsv"

ZIP_PACKETS = [
    (
        "Supplementary_Tables_Primary_Reader_Facing_S1_S40_S102_S103.zip",
        lambda row: row["reviewer_priority"] == "1_read_first",
    ),
    (
        "Supplementary_Tables_Supplementary_Controls_S41_S53_S89_S93.zip",
        lambda row: row["submission_packet"] == "Reader-facing supplementary controls",
    ),
    (
        "Supplementary_Data_Archive_S54_S88_S94_S101.zip",
        lambda row: row["submission_packet"] == "Supplementary Data archive",
    ),
]


SUBMISSION_GROUPS = [
    {
        "start": 1,
        "end": 2,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Dataset frame and curation",
        "recommended_use": "Directly supports Figure 1 and the primary-core dataset design.",
    },
    {
        "start": 3,
        "end": 8,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Candidate discovery and tiering",
        "recommended_use": "Directly supports Figure 2 and the ortholog-aware candidate-tiering pipeline.",
    },
    {
        "start": 9,
        "end": 18,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Comparator and configuration tests",
        "recommended_use": "Directly supports Figures 3 and 4, including named comparators and configuration scores.",
    },
    {
        "start": 19,
        "end": 28,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Stage-window and signaling reinterpretation",
        "recommended_use": "Directly supports Figure 5 pathway, stage-window, conditioned-medium, and 2005-paper comparison analyses.",
    },
    {
        "start": 29,
        "end": 37,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Morphology, physiology and sparse-coding calibration",
        "recommended_use": "Directly supports Figure 5 morphology/activity validation and model calibration.",
    },
    {
        "start": 38,
        "end": 40,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Sender-receiver predictions",
        "recommended_use": "Directly supports Figure 6 ligand-receptor sender-receiver predictions.",
    },
    {
        "start": 41,
        "end": 53,
        "submission_packet": "Reader-facing supplementary controls",
        "evidence_tier": "Supporting result",
        "reviewer_priority": "2_read_for_supplementary_figures",
        "topic_group": "Stage split, granule specificity and developmental-origin controls",
        "recommended_use": "Supports Supplementary Figures S1-S4 and the argument against a recent common granule progenitor.",
    },
    {
        "start": 54,
        "end": 88,
        "submission_packet": "Supplementary Data archive",
        "evidence_tier": "Archive/provenance",
        "reviewer_priority": "3_machine_readable_archive",
        "topic_group": "Regulatory-compatibility extension",
        "recommended_use": "Retain as machine-readable regulatory-resource, feature-manifest, and sensitivity-analysis support rather than as primary evidence.",
    },
    {
        "start": 89,
        "end": 93,
        "submission_packet": "Reader-facing supplementary controls",
        "evidence_tier": "Supporting result",
        "reviewer_priority": "2_read_for_supplementary_figures",
        "topic_group": "Hierarchical integrative evidence model",
        "recommended_use": "Supports Supplementary Figure S5 and summarizes the term-balanced evidence synthesis.",
    },
    {
        "start": 94,
        "end": 100,
        "submission_packet": "Supplementary Data archive",
        "evidence_tier": "Archive/provenance",
        "reviewer_priority": "3_machine_readable_archive",
        "topic_group": "Public perturbation audit",
        "recommended_use": "Retain as hypothesis-generating perturbation-resource triage and module-shift support for Supplementary Figure S6.",
    },
    {
        "start": 101,
        "end": 101,
        "submission_packet": "Supplementary Data archive",
        "evidence_tier": "Provenance map",
        "reviewer_priority": "3_machine_readable_archive",
        "topic_group": "Result-level evidence anchors",
        "recommended_use": "Maps each Results section to source analysis files and the corresponding supplementary tables.",
    },
    {
        "start": 102,
        "end": 103,
        "submission_packet": "Reader-facing primary supplementary tables",
        "evidence_tier": "Primary result",
        "reviewer_priority": "1_read_first",
        "topic_group": "Hypothesis-support matrix",
        "recommended_use": "Directly supports Figure 7 and the evidence-weighted H1/H2/H3 comparison.",
    },
]


PROVENANCE_ROWS = [
    {
        "result_id": "Result 1",
        "manuscript_section": "A strict 10-dataset primary core frames the granule-cell convergence problem",
        "claim_focus": "Primary-core dataset frame, strict versus supporting dataset separation, and branch-level interpretation",
        "provenance_files": "Project/results/integrated_primary_core_datasets.md; Project/results/revised_core_supporting_dataset_tiers.md",
        "corresponding_supplementary_tables": "Table_S1; Table_S2",
        "notes": "Documents the core dataset inclusion frame and supporting/scaffold resource separation.",
    },
    {
        "result_id": "Result 2",
        "manuscript_section": "Ortholog-aware rank-meta modeling identifies a conservative shared candidate set",
        "claim_focus": "Formal MGI rank-meta discovery, manuscript candidate tiers, and mechanism-axis organization",
        "provenance_files": "Project/results/primary_core_mgi_ortholog_formal_rank_model.md; Project/results/primary_core_manuscript_candidate_tiers.md; Project/results/primary_core_mechanism_axis_model.md",
        "corresponding_supplementary_tables": "Table_S3; Table_S4; Table_S5; Table_S6; Table_S7; Table_S8",
        "notes": "Documents the discovery-to-tiering route from formal shared hits to Tier 1/Tier 2 mechanism candidates.",
    },
    {
        "result_id": "Result 3",
        "manuscript_section": "Shared signals are constrained by named-comparator specificity",
        "claim_focus": "Named-comparator specificity, upstream versus downstream module separation, regional-origin timing, and developmental-origin control",
        "provenance_files": "Project/results/primary_core_granule_specificity_named_comparators.md; Project/results/primary_core_niche_circuit_module_model.md; Project/results/primary_core_stage_split_granule_comparison.md; Project/results/regional_origin_shared_toolkit_timing.md; Project/results/developmental_origin_divergence_audit.md",
        "corresponding_supplementary_tables": "Table_S9; Table_S10; Table_S11; Table_S12; Table_S13; Table_S41; Table_S42; Table_S43; Table_S44; Table_S45; Table_S46; Table_S47; Table_S48; Table_S49; Table_S50; Table_S51; Table_S52; Table_S53",
        "notes": "Documents the evidence separating shared downstream construction modules from broad neuronal or recent-lineage explanations.",
    },
    {
        "result_id": "Result 4",
        "manuscript_section": "Candidate granule populations show an identity-coupled transcriptomic assembly configuration",
        "claim_focus": "Transcriptomic configuration score, primary-core validation, and driver decomposition",
        "provenance_files": "Project/results/primary_core_transcriptomic_configuration_model.md; Project/results/primary_core_transcriptomic_configuration_primary_validation.md; Project/results/primary_core_configuration_driver_audit.md",
        "corresponding_supplementary_tables": "Table_S14; Table_S15; Table_S16; Table_S17; Table_S18",
        "notes": "Documents local and primary-core configuration contrasts and whether positive configuration is driven by construction, fate polarity, or both.",
    },
    {
        "result_id": "Result 5",
        "manuscript_section": "Stage-windowed niche signaling and sparse-coding constraints refine the convergence model",
        "claim_focus": "Pathway readiness, stage-window modeling, conditioned-medium reinterpretation, morphology/activity validation, and sparse-coding calibration",
        "provenance_files": "Project/results/primary_core_aim2_niche_pathway_model.md; Project/results/primary_core_aim2b_stage_tgf_bdnf.md; Project/results/aim2_stage_window_model.md; Project/results/cerebellar_conditioned_medium_secretome_candidates.md; Project/results/neuromorpho_granule_morphometry_validation.md; Project/results/dandi_000003_multisession_spatial_extension.md; Project/results/primary_core_aim3_sparse_coding_model.md; Project/results/aim3_empirical_calibration.md",
        "corresponding_supplementary_tables": "Table_S19; Table_S20; Table_S21; Table_S22; Table_S23; Table_S24; Table_S25; Table_S26; Table_S27; Table_S28; Table_S29; Table_S30; Table_S31; Table_S32; Table_S33; Table_S34; Table_S35; Table_S36; Table_S37",
        "notes": "Documents the stage-dependent pathway model and the external morphology, activity, and computational calibration layers.",
    },
    {
        "result_id": "Result 6",
        "manuscript_section": "Focused sender-receiver ligand-receptor prediction nominates testable niche cues",
        "claim_focus": "Directional sender-receiver predictions, hierarchical evidence summary, and public perturbation module-shift support",
        "provenance_files": "Project/results/aim2_sender_receiver_lr.md; Project/results/aim2_sender_receiver_lr_top_predictions.tsv; Project/results/aim2_sender_receiver_lr_predictions.tsv.gz; Project/results/hierarchical_integrative_model.md; Project/results/causal_perturbation_dataset_triage.md; Project/results/causal_perturbation_module_shift.md",
        "corresponding_supplementary_tables": "Table_S38; Table_S39; Table_S40; Table_S89; Table_S90; Table_S91; Table_S92; Table_S93; Table_S94; Table_S95; Table_S96; Table_S97; Table_S98; Table_S99; Table_S100",
        "notes": "Documents expression-supported sender-receiver hypotheses and the supplementary integrative/perturbation support layers.",
    },
    {
        "result_id": "Result 7",
        "manuscript_section": "Evidence-weighted hypothesis comparison favors an integrated convergence model",
        "claim_focus": "Term-level evidence scores, hypothesis coefficients, and support-index summaries for H1, H2, H3, and H2+H3",
        "provenance_files": "Project/results/hierarchical_integrative_model.md; Project/results/hypothesis_support_score_model.md; Project/results/hypothesis_support_score_matrix.tsv; Project/results/hypothesis_support_scores.tsv",
        "corresponding_supplementary_tables": "Table_S89; Table_S90; Table_S91; Table_S92; Table_S93; Table_S102; Table_S103",
        "notes": "Documents the manuscript-defined evidence-alignment calculation used for Fig. 7.",
    },
]


TABLES = [
    (
        "primary_core_dataset_frame",
        "Primary-core dataset frame",
        "Figure 1; dataset design",
        "Strict 10-dataset manuscript core with branch, species, and role annotations.",
        "integrated_primary_core_datasets.tsv",
    ),
    (
        "supporting_dataset_tiers",
        "Supporting dataset tiers",
        "Figure 1; dataset design",
        "Secondary/scaffold resources retained outside the strict manuscript core.",
        "revised_core_supporting_dataset_tiers.tsv",
    ),
    (
        "candidate_gene_tiers",
        "Candidate gene tiers",
        "Figure 2",
        "Tiered candidate genes used for manuscript interpretation.",
        "primary_core_manuscript_candidate_tiers.tsv",
    ),
    (
        "formal_rank_gene_summary",
        "Formal ortholog rank-meta gene summary",
        "Figure 2; candidate discovery",
        "Gene-level formal MGI ortholog rank-meta summary across primary-core layers.",
        "primary_core_mgi_ortholog_formal_rank_gene_summary.tsv",
    ),
    (
        "formal_rank_shared_hits",
        "Formal shared hits",
        "Figure 2; candidate discovery",
        "Genes passing formal shared-hit criteria in the ortholog rank-meta model.",
        "primary_core_mgi_ortholog_formal_rank_shared_hits.tsv",
    ),
    (
        "formal_rank_mechanism_hits",
        "Formal mechanism-axis hits",
        "Figure 2",
        "Formal model hits annotated by mechanism axis.",
        "primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv",
    ),
    (
        "mechanism_axis_gene_table",
        "Mechanism-axis gene table",
        "Figure 2",
        "Gene-to-axis assignments for developmental, synaptic, neurite, and guidance axes.",
        "primary_core_mechanism_axis_gene_table.tsv",
    ),
    (
        "mechanism_axis_summary",
        "Mechanism-axis summary",
        "Figure 2",
        "Summary statistics for prioritized mechanism axes.",
        "primary_core_mechanism_axis_summary.tsv",
    ),
    (
        "named_comparator_axis_summary",
        "Named-comparator mechanism-axis summary",
        "Figure 3a",
        "Mechanism-axis ranks in dentate granule, pyramidal, cerebellar granule, and Purkinje comparators.",
        "primary_core_granule_specificity_named_comparator_axis_summary.tsv",
    ),
    (
        "named_comparator_axis_units",
        "Named-comparator mechanism-axis units",
        "Figure 3a",
        "Source group-level units underlying named-comparator mechanism-axis ranks.",
        "primary_core_granule_specificity_named_comparator_units.tsv",
    ),
    (
        "niche_circuit_gene_sets",
        "Niche/circuit module gene sets",
        "Figure 3b",
        "Curated upstream fate, neurogenic niche, and downstream construction module genes.",
        "primary_core_niche_circuit_module_gene_sets.tsv",
    ),
    (
        "niche_circuit_formal_summary",
        "Niche/circuit formal summary",
        "Figure 3b",
        "Formal-core module convergence summary for upstream versus downstream modules.",
        "primary_core_niche_circuit_module_formal_summary.tsv",
    ),
    (
        "niche_circuit_named_comparator_summary",
        "Niche/circuit named-comparator summary",
        "Figure 3b",
        "Local named-comparator specificity summary for niche/circuit modules.",
        "primary_core_niche_circuit_module_named_comparator_summary.tsv",
    ),
    (
        "configuration_local_contrasts",
        "Transcriptomic configuration local contrasts",
        "Figure 4b",
        "Local named-comparator configuration contrasts.",
        "primary_core_transcriptomic_configuration_contrasts.tsv",
    ),
    (
        "configuration_primary_summary",
        "Primary-core configuration summary",
        "Figure 4c",
        "Primary-core configuration score summary.",
        "primary_core_transcriptomic_configuration_primary_summary.tsv",
    ),
    (
        "configuration_primary_contrasts",
        "Primary-core configuration contrasts",
        "Figure 4c",
        "Candidate-background configuration deltas across full-MGI and selected-feature layers.",
        "primary_core_transcriptomic_configuration_primary_contrasts.tsv",
    ),
    (
        "configuration_driver_summary",
        "Configuration driver summary",
        "Figure 4d",
        "Driver decomposition summary for configuration-positive contrasts.",
        "primary_core_configuration_driver_audit_summary.tsv",
    ),
    (
        "configuration_driver_contrasts",
        "Configuration driver contrasts",
        "Figure 4d",
        "Contrast-level driver classes for construction and fate-polarity components.",
        "primary_core_configuration_driver_audit_contrasts.tsv",
    ),
    (
        "aim2_pathway_gene_sets",
        "Pathway-readiness gene sets",
        "Figure 5a",
        "Curated pathway/readiness gene sets used for TGF-beta/BDNF and niche analyses.",
        "primary_core_aim2_pathway_gene_sets.tsv",
    ),
    (
        "aim2_pathway_summary",
        "Pathway-readiness summary",
        "Figure 5a",
        "Branch-level pathway/readiness summary.",
        "primary_core_aim2_pathway_summary.tsv",
    ),
    (
        "aim2_pathway_contrasts",
        "Pathway-readiness contrasts",
        "Figure 5a",
        "Candidate-background contrasts for pathway/readiness modules.",
        "primary_core_aim2_pathway_contrasts.tsv",
    ),
    (
        "aim2_ligand_receptor_summary",
        "Ligand-receptor readiness summary",
        "Figure 5a",
        "Nondirectional ligand-receptor readiness summaries from the pathway-readiness analysis.",
        "primary_core_aim2_ligand_receptor_summary.tsv",
    ),
    (
        "tgf_bdnf_stage_summary",
        "Stage-resolved TGF-beta/BDNF summary",
        "Figure 5b",
        "Stage-resolved TGF-beta/BDNF and related pathway summaries.",
        "primary_core_aim2b_stage_tgf_bdnf_summary.tsv",
    ),
    (
        "tgf_bdnf_stage_transitions",
        "Stage-resolved TGF-beta/BDNF transitions",
        "Figure 5b",
        "Transition summaries for TGF-beta/BDNF-related stage effects.",
        "primary_core_aim2b_stage_tgf_bdnf_transitions.tsv",
    ),
    (
        "stage_window_branch_summary",
        "Stage-window model branch summary",
        "Figure 5b",
        "Branch summaries from the fitted stage-window model.",
        "aim2_stage_window_model_branch_summary.tsv",
    ),
    (
        "stage_window_coefficients",
        "Stage-window model coefficients",
        "Figure 5b",
        "Regression coefficients for fitted stage-window models.",
        "aim2_stage_window_model_coefficients.tsv",
    ),
    (
        "conditioned_medium_secretome_candidates",
        "Conditioned-medium secretome candidates",
        "Figure 5; Discussion",
        "Prioritized non-TGF-beta/BDNF secretome candidates from cerebellar conditioned-medium analysis.",
        "cerebellar_conditioned_medium_secretome_ranked_candidates.tsv",
    ),
    (
        "paper_2005_support_revision",
        "2005 paper support/revision table",
        "Figure 5; Discussion",
        "How current sequencing analyses support, revise, or do not directly test the 2005 paper results.",
        "paper_2005_support_revision_table.tsv",
    ),
    (
        "neuromorpho_summary",
        "NeuroMorpho summary",
        "Figure 5c",
        "Dendritic morphology summary for public dentate and cerebellar granule-cell reconstructions.",
        "neuromorpho_granule_morphometry_summary.tsv",
    ),
    (
        "neuromorpho_comparison",
        "NeuroMorpho comparison",
        "Figure 5c",
        "Comparison statistics for dentate versus cerebellar granule-cell morphometry.",
        "neuromorpho_granule_morphometry_comparison.tsv",
    ),
    (
        "dandi_session_summary",
        "DANDI multisession summary",
        "Figure 5c",
        "Session-level DANDI dentate activity-validation summary.",
        "dandi_000003_multisession_session_summary.tsv",
    ),
    (
        "dandi_unit_metrics",
        "DANDI spatial unit metrics",
        "Figure 5c",
        "Unit-level spatial information, sparsity, and firing metrics.",
        "dandi_000003_multisession_spatial_unit_metrics.tsv",
    ),
    (
        "dandi_population_vector",
        "DANDI population-vector separation",
        "Figure 5c",
        "Population-vector near/far separation checks from DANDI dentate sessions.",
        "dandi_000003_multisession_population_vector_separation.tsv",
    ),
    (
        "sparse_coding_architecture_summary",
        "Sparse-coding architecture summary",
        "Figure 5c",
        "Named architecture summaries from the sparse expansion-coding simulation.",
        "primary_core_aim3_sparse_coding_architecture_summary.tsv",
    ),
    (
        "sparse_coding_parameter_grid",
        "Sparse-coding parameter grid",
        "Figure 5c",
        "Full sparse expansion-coding simulation parameter grid.",
        "primary_core_aim3_sparse_coding_parameter_grid.tsv",
    ),
    (
        "empirical_calibration_architecture_summary",
        "Empirical calibration architecture summary",
        "Figure 5c",
        "Architecture-level empirical calibration against morphology/activity constraints.",
        "aim3_empirical_calibration_architecture_summary.tsv",
    ),
    (
        "empirical_calibration_grid",
        "Empirical calibration grid",
        "Figure 5c",
        "Full empirical calibration grid for sparse-coding model fitting.",
        "aim3_empirical_calibration_grid.tsv",
    ),
    (
        "sender_receiver_summary",
        "Sender-receiver summary",
        "Figure 6",
        "Focused directional ligand-receptor sender-receiver summary.",
        "aim2_sender_receiver_lr_summary.tsv",
    ),
    (
        "sender_receiver_top_predictions",
        "Sender-receiver top predictions",
        "Figure 6",
        "Top expression-supported sender-receiver ligand-receptor predictions.",
        "aim2_sender_receiver_lr_top_predictions.tsv",
    ),
    (
        "sender_receiver_predictions",
        "Sender-receiver prediction units",
        "Figure 6",
        "Full expression-supported sender-receiver ligand-receptor prediction table.",
        "aim2_sender_receiver_lr_predictions.tsv.gz",
    ),
    (
        "stage_split_module_branch_summary",
        "Stage-split module branch summary",
        "Supplementary Figure S1",
        "Branch/stage module summaries for immature-versus-mature granule-cell comparison.",
        "primary_core_stage_split_granule_module_branch_summary.tsv",
    ),
    (
        "stage_split_similarity",
        "Stage-split similarity",
        "Supplementary Figure S1",
        "Cross-branch similarity scores for immature and mature/maturing bins.",
        "primary_core_stage_split_granule_stage_similarity.tsv",
    ),
    (
        "stage_split_transitions",
        "Stage-split transitions",
        "Supplementary Figure S1",
        "Mature-minus-immature similarity transition summary.",
        "primary_core_stage_split_granule_stage_transitions.tsv",
    ),
    (
        "stage_split_group_calls",
        "Stage-split group calls",
        "Supplementary Figure S1",
        "Dataset/group stage calls used in the stage-split analysis.",
        "primary_core_stage_split_granule_group_calls.tsv",
    ),
    (
        "granule_special_top_candidates",
        "Granule-special top candidates",
        "Supplementary Figure S2",
        "Top granule-enriched candidates relative to named pyramidal and Purkinje comparators.",
        "primary_core_granule_special_gene_named_comparator_top_candidates.tsv",
    ),
    (
        "granule_special_summary",
        "Granule-special genome-wide summary",
        "Supplementary Figure S2",
        "Genome-wide named-comparator screen summary for granule-special genes.",
        "primary_core_granule_special_gene_named_comparator_summary.tsv",
    ),
    (
        "regional_origin_timing_state_summary",
        "Regional-origin timing state summary",
        "Supplementary Figure S3",
        "Ordered state-level summary for regional-origin/shared-toolkit timing analysis.",
        "regional_origin_shared_toolkit_timing_state_summary.tsv",
    ),
    (
        "regional_origin_timing_metrics",
        "Regional-origin timing metrics",
        "Supplementary Figure S3",
        "Early-to-late metric shifts for regional-origin/shared-toolkit timing analysis.",
        "regional_origin_shared_toolkit_timing_metrics.tsv",
    ),
    (
        "regional_origin_timing_gene_units",
        "Regional-origin timing gene units",
        "Supplementary Figure S3",
        "Gene-level timing units for `NFIA`, `NEUROD1`, `RBFOX3`, and `HMGN2`.",
        "regional_origin_shared_toolkit_timing_gene_units.tsv",
    ),
    (
        "developmental_origin_gene_sets",
        "Developmental-origin gene sets",
        "Supplementary Figure S4",
        "Curated marker modules for deep origin, regional divergence, and later toolkit reuse.",
        "developmental_origin_divergence_audit_gene_sets.tsv",
    ),
    (
        "developmental_origin_state_summary",
        "Developmental-origin state summary",
        "Supplementary Figure S4",
        "Ordered state-level module ranks for the developmental-origin audit.",
        "developmental_origin_divergence_audit_state_summary.tsv",
    ),
    (
        "developmental_origin_branch_metrics",
        "Developmental-origin branch metrics",
        "Supplementary Figure S4",
        "Branch-level early-to-late metrics for the developmental-origin audit.",
        "developmental_origin_divergence_audit_branch_metrics.tsv",
    ),
    (
        "developmental_origin_units",
        "Developmental-origin module units",
        "Supplementary Figure S4",
        "Group-level module scoring units for the developmental-origin audit.",
        "developmental_origin_divergence_audit_units.tsv",
    ),
    (
        "epigenomic_extension_candidate_resources",
        "Epigenomic extension candidate resources",
        "Future epigenomic extension",
        "Candidate matched or comparable scATAC/multiome/spatial epigenomic resources for adding a regulatory-compatibility layer.",
        "epigenomic_extension_candidate_resources.tsv",
    ),
    (
        "epigenomic_regulatory_targets",
        "Epigenomic regulatory targets",
        "Future epigenomic extension",
        "Regulatory target genes/modules for promoter, enhancer, methylation, and motif support in matched or comparable epigenomic datasets.",
        "epigenomic_extension_regulatory_targets.tsv",
    ),
    (
        "integrative_model_term_specification",
        "Integrative model term specification",
        "Future epigenomic extension",
        "Single-equation model terms, data sources, current fitting status, and recommended next actions.",
        "integrative_granule_model_term_specification.tsv",
    ),
    (
        "epigenomic_target_summary",
        "Epigenomic target summary",
        "Future epigenomic extension",
        "Compact summary of epigenomic target rows by source, model term, priority, and GSE268609 feature coverage.",
        "epigenomic_extension_target_summary.tsv",
    ),
    (
        "gse268609_epigenomic_peak_targets",
        "GSE268609 epigenomic peak targets",
        "Future epigenomic extension",
        "Nearby ATAC peak feature rows mapped to regulatory target genes within gene body plus 100 kb.",
        "gse268609_epigenomic_peak_targets.tsv",
    ),
    (
        "gse268609_epigenomic_peak_gene_summary",
        "GSE268609 epigenomic peak gene summary",
        "Future epigenomic extension",
        "Gene-level nearby ATAC peak counts for future selective GSE268609 epigenomic extraction.",
        "gse268609_epigenomic_peak_gene_summary.tsv",
    ),
    (
        "gse268609_epigenomic_selective_extraction_manifest",
        "GSE268609 epigenomic selective extraction manifest",
        "Future epigenomic extension",
        "Target gene-expression and nearby ATAC peak feature rows for selective matrix extraction.",
        "gse268609_epigenomic_selective_extraction_manifest.tsv",
    ),
    (
        "gse322785_cerebellar_multiome_sample_metadata",
        "GSE322785 cerebellar multiome sample metadata",
        "Future epigenomic extension",
        "Parsed GEO sample metadata for the adult primate cerebellar cortex multiome resource.",
        "gse322785_cerebellar_multiome_sample_metadata.tsv",
    ),
    (
        "gse322785_cerebellar_multiome_file_inventory",
        "GSE322785 cerebellar multiome file inventory",
        "Future epigenomic extension",
        "GEO supplementary-file inventory for GSE322785 H5 and ATAC fragment resources.",
        "gse322785_cerebellar_multiome_file_inventory.tsv",
    ),
    (
        "gse322785_cerebellar_multiome_donor_summary",
        "GSE322785 cerebellar multiome donor summary",
        "Future epigenomic extension",
        "Donor-level summary of available GSE322785 cerebellar multiome files by species.",
        "gse322785_cerebellar_multiome_donor_summary.tsv",
    ),
    (
        "gse322785_cerebellar_multiome_download_plan",
        "GSE322785 cerebellar multiome download plan",
        "Future epigenomic extension",
        "Tiered download plan separating compact H5 matrices from larger deferred ATAC fragment files.",
        "gse322785_cerebellar_multiome_download_plan.tsv",
    ),
    (
        "gse322785_human_h5_feature_inventory",
        "GSE322785 human H5 feature inventory",
        "Future epigenomic extension",
        "Validated feature and barcode counts for downloaded human cerebellar multiome H5 files.",
        "gse322785_human_h5_feature_inventory.tsv",
    ),
    (
        "gse322785_human_h5_epigenomic_gene_summary",
        "GSE322785 human H5 epigenomic gene summary",
        "Future epigenomic extension",
        "Target-gene feature coverage and nearby ATAC peak counts in human cerebellar H5 files.",
        "gse322785_human_h5_epigenomic_gene_summary.tsv",
    ),
    (
        "gse322785_human_h5_epigenomic_peak_targets",
        "GSE322785 human H5 epigenomic peak targets",
        "Future epigenomic extension",
        "Nearby ATAC peak feature rows mapped to regulatory target genes in downloaded human cerebellar H5 files.",
        "gse322785_human_h5_epigenomic_peak_targets.tsv",
    ),
    (
        "gse322785_human_h5_epigenomic_selective_manifest",
        "GSE322785 human H5 epigenomic selective manifest",
        "Future epigenomic extension",
        "Target gene-expression and nearby ATAC peak feature rows for selective extraction from GSE322785 human H5 matrices.",
        "gse322785_human_h5_epigenomic_selective_manifest.tsv",
    ),
    (
        "gse322785_human_h5_selected_matrix_summary",
        "GSE322785 human H5 selected matrix summary",
        "Future epigenomic extension",
        "Per-donor summary of selected RNA/ATAC matrix extraction and provisional marker-call coverage.",
        "gse322785_human_h5_selected_matrix_summary.tsv",
    ),
    (
        "gse322785_human_h5_marker_panel_coverage",
        "GSE322785 human H5 marker panel coverage",
        "Future epigenomic extension",
        "Marker-gene panel coverage in downloaded human cerebellar H5 feature tables.",
        "gse322785_human_h5_marker_panel_coverage.tsv",
    ),
    (
        "gse322785_human_h5_marker_celltype_summary",
        "GSE322785 human H5 marker cell-type summary",
        "Future epigenomic extension",
        "Provisional marker-call counts by donor, confidence class, and basic-QC status.",
        "gse322785_human_h5_marker_celltype_summary.tsv",
    ),
    (
        "gse322785_human_h5_marker_high_confidence_barcodes",
        "GSE322785 human H5 high-confidence marker barcodes",
        "Future epigenomic extension",
        "Barcode-level high/medium-confidence provisional marker calls for selected downstream scoring.",
        "gse322785_human_h5_marker_high_confidence_barcodes.tsv.gz",
    ),
    (
        "gse322785_human_h5_epigenomic_marker_group_feature_scores",
        "GSE322785 marker-group feature scores",
        "Future epigenomic extension",
        "Selected target-gene and ATAC-peak feature scores by provisional marker group.",
        "gse322785_human_h5_epigenomic_marker_group_feature_scores.tsv.gz",
    ),
    (
        "gse322785_human_h5_epigenomic_marker_group_module_scores",
        "GSE322785 marker-group module scores",
        "Future epigenomic extension",
        "Module-level selected RNA/ATAC scores by provisional marker group.",
        "gse322785_human_h5_epigenomic_marker_group_module_scores.tsv",
    ),
    (
        "gse322785_human_h5_epigenomic_marker_group_contrasts",
        "GSE322785 marker-group epigenomic contrasts",
        "Future epigenomic extension",
        "Granule-candidate versus Purkinje/glial/ambiguous-neuronal provisional marker-group contrasts.",
        "gse322785_human_h5_epigenomic_marker_group_contrasts.tsv",
    ),
    (
        "gse322785_human_h5_cluster_validation_barcode_assignments",
        "GSE322785 cluster-validation barcode assignments",
        "Future epigenomic extension",
        "Basic-QC barcode-level selected-gene SVD/k-means cluster assignments and provisional marker calls.",
        "gse322785_human_h5_cluster_validation_barcode_assignments.tsv.gz",
    ),
    (
        "gse322785_human_h5_cluster_validation_summary",
        "GSE322785 cluster-validation summary",
        "Future epigenomic extension",
        "Donor-specific selected-gene cluster summaries with marker composition and median marker scores.",
        "gse322785_human_h5_cluster_validation_summary.tsv",
    ),
    (
        "gse322785_human_h5_cluster_validation_marker_call_enrichment",
        "GSE322785 cluster marker-call enrichment",
        "Future epigenomic extension",
        "Cluster-by-marker-call enrichment table for provisional GSE322785 labels.",
        "gse322785_human_h5_cluster_validation_marker_call_enrichment.tsv",
    ),
    (
        "gse322785_human_h5_cluster_validation_marker_support",
        "GSE322785 cluster marker support",
        "Future epigenomic extension",
        "Dominant-cluster concentration summary for each provisional marker call.",
        "gse322785_human_h5_cluster_validation_marker_support.tsv",
    ),
    (
        "gse322785_human_h5_cluster_validation_metrics",
        "GSE322785 cluster-validation metrics",
        "Future epigenomic extension",
        "Per-donor adjusted Rand index and normalized mutual information between selected-gene clusters and provisional marker calls.",
        "gse322785_human_h5_cluster_validation_metrics.tsv",
    ),
    (
        "gse322785_human_h5_cluster_supported_marker_rules",
        "GSE322785 cluster-supported marker rules",
        "Future epigenomic extension",
        "Cluster-marker enrichment rules used to define stricter supported provisional marker groups.",
        "gse322785_human_h5_cluster_supported_marker_rules.tsv",
    ),
    (
        "gse322785_human_h5_cluster_supported_marker_barcodes",
        "GSE322785 cluster-supported marker barcodes",
        "Future epigenomic extension",
        "Barcode-level cluster-supported provisional marker calls used for stricter sensitivity scoring.",
        "gse322785_human_h5_cluster_supported_marker_barcodes.tsv.gz",
    ),
    (
        "gse322785_human_h5_cluster_supported_epigenomic_feature_scores",
        "GSE322785 cluster-supported feature scores",
        "Future epigenomic extension",
        "Selected target-gene and ATAC-peak feature scores after restricting to cluster-supported provisional marker groups.",
        "gse322785_human_h5_cluster_supported_epigenomic_feature_scores.tsv.gz",
    ),
    (
        "gse322785_human_h5_cluster_supported_epigenomic_module_scores",
        "GSE322785 cluster-supported module scores",
        "Future epigenomic extension",
        "Module-level selected RNA/ATAC scores after restricting to cluster-supported provisional marker groups.",
        "gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv",
    ),
    (
        "gse322785_human_h5_cluster_supported_epigenomic_contrasts",
        "GSE322785 cluster-supported epigenomic contrasts",
        "Future epigenomic extension",
        "Granule-candidate versus supported comparator contrasts after cluster-enrichment filtering.",
        "gse322785_human_h5_cluster_supported_epigenomic_contrasts.tsv",
    ),
    (
        "gse322785_epigenomic_broad_vs_cluster_supported_contrasts",
        "GSE322785 broad-versus-supported contrast comparison",
        "Future epigenomic extension",
        "Concordance comparison between broad provisional and cluster-supported GSE322785 epigenomic contrasts.",
        "gse322785_epigenomic_broad_vs_cluster_supported_contrasts.tsv",
    ),
    (
        "gse322785_epigenomic_robust_positive_contrasts",
        "GSE322785 robust-positive epigenomic contrasts",
        "Future epigenomic extension",
        "Contrasts with positive granule-candidate deltas in both broad provisional and cluster-supported sensitivity layers.",
        "gse322785_epigenomic_robust_positive_contrasts.tsv",
    ),
    (
        "gse322785_epigenomic_robust_contrast_summary",
        "GSE322785 robust contrast summary",
        "Future epigenomic extension",
        "Summary of robust, layer-specific, discordant, and weak GSE322785 epigenomic contrast classes.",
        "gse322785_epigenomic_robust_contrast_summary.tsv",
    ),
    (
        "hierarchical_integrative_model_terms",
        "Hierarchical integrative model terms",
        "Integrative model",
        "Term definitions, hierarchy levels, evidence quality classes, and caveats for the hierarchical integrative granule-cell model.",
        "hierarchical_integrative_model_terms.tsv",
    ),
    (
        "hierarchical_integrative_model_evidence_units",
        "Hierarchical integrative model evidence units",
        "Integrative model",
        "Evidence-unit table combining direct transcriptomic, stage/niche, epigenomic, morphology/activity, and resource-constraint layers.",
        "hierarchical_integrative_model_evidence_units.tsv",
    ),
    (
        "hierarchical_integrative_model_layer_summary",
        "Hierarchical integrative model layer summary",
        "Integrative model",
        "Level-by-term weighted summaries for the hierarchical integrative evidence model.",
        "hierarchical_integrative_model_layer_summary.tsv",
    ),
    (
        "hierarchical_integrative_model_branch_summary",
        "Hierarchical integrative model branch summary",
        "Integrative model",
        "Branch-level term summaries for dentate, cerebellar, cross-branch, and simulation evidence.",
        "hierarchical_integrative_model_branch_summary.tsv",
    ),
    (
        "hierarchical_integrative_model_component_scores",
        "Hierarchical integrative model component scores",
        "Integrative model",
        "Component-level term-balanced scores and caveats for the hierarchical integrative model.",
        "hierarchical_integrative_model_component_scores.tsv",
    ),
    (
        "causal_perturbation_dataset_triage",
        "Public perturbation dataset triage",
        "Causal follow-up resources",
        "Curated public NFIA, BDNF/TrkB, TGF-beta/SMAD, SHH, RBFOX3, and HMGN2 perturbation or pathway-response datasets for causal follow-up planning.",
        "causal_perturbation_dataset_triage.tsv",
    ),
    (
        "causal_perturbation_node_summary",
        "Public perturbation node summary",
        "Causal follow-up resources",
        "Node-level summary of public perturbation-resource strength, recommended model role, and key caveats.",
        "causal_perturbation_node_summary.tsv",
    ),
    (
        "causal_perturbation_processing_status",
        "Public perturbation processing status",
        "Causal follow-up resources",
        "Processing status for downloaded perturbation datasets used in module-shift scoring.",
        "causal_perturbation_processing_status.tsv",
    ),
    (
        "causal_perturbation_module_catalog",
        "Public perturbation module catalog",
        "Causal follow-up resources",
        "Curated fate, niche, pathway, neurite, synaptic, and excitability modules used to score perturbation module shifts.",
        "causal_perturbation_module_catalog.tsv",
    ),
    (
        "causal_perturbation_module_shift_gene_effects",
        "Public perturbation gene-level module effects",
        "Causal follow-up resources",
        "Gene-level perturbation effects for module genes recovered from public RBFOX3/RBFOX1, BDNF/NTRK2, and SHH/PTCH/Norrin datasets.",
        "causal_perturbation_module_shift_gene_effects.tsv",
    ),
    (
        "causal_perturbation_module_shift_summary",
        "Public perturbation contrast-module shifts",
        "Causal follow-up resources",
        "Contrast-level module-shift summaries with module coverage, signed shift scores, context-match calls, and top shifted genes.",
        "causal_perturbation_module_shift_summary.tsv",
    ),
    (
        "causal_perturbation_module_shift_node_summary",
        "Public perturbation node-module shifts",
        "Causal follow-up resources",
        "Node-level summaries of matched-context and off-context perturbation-sensitive modules.",
        "causal_perturbation_module_shift_node_summary.tsv",
    ),
    (
        "result_evidence_anchor_provenance",
        "Result-level evidence anchor provenance",
        "Analysis provenance",
        "Result-level mapping from each Results section to source analysis files and corresponding supplementary tables.",
        GENERATED_PROVENANCE_SOURCE,
    ),
    (
        "hypothesis_support_score_matrix",
        "Hypothesis support coefficient matrix",
        "Figure 7",
        "Evidence-term scores and H1/H2/H3 prediction coefficients used for the Fig. 7 support-index calculation.",
        "hypothesis_support_score_matrix.tsv",
    ),
    (
        "hypothesis_support_scores",
        "Hypothesis support scores",
        "Figure 7",
        "Final support-index scores and domain summaries for H1, H2, H3, and the integrated H2+H3 synthesis.",
        "hypothesis_support_scores.tsv",
    ),
]


def destination_name(index: int, slug: str, source: str) -> str:
    suffix = ".tsv.gz" if source.endswith(".tsv.gz") else ".tsv"
    return f"Table_S{index}_{slug}{suffix}"


def cleanup_generated_outputs() -> None:
    for path in OUT_DIR.iterdir():
        if path.name.startswith("Table_S") and (path.name.endswith(".tsv") or path.name.endswith(".tsv.gz")):
            path.unlink()
        elif path.name.startswith("Supplementary_") and path.name.endswith(".zip"):
            path.unlink()


def submission_metadata(index: int) -> dict[str, str]:
    for group in SUBMISSION_GROUPS:
        if group["start"] <= index <= group["end"]:
            return {k: str(v) for k, v in group.items() if k not in {"start", "end"}}
    raise ValueError(f"No submission metadata for table {index}")


def table_range(group: dict[str, object]) -> str:
    start = int(group["start"])
    end = int(group["end"])
    if start == end:
        return f"Table S{start}"
    return f"Tables S{start}-S{end}"


def write_result_evidence_anchor_provenance() -> Path:
    path = OUT_DIR / "Table_S101_result_evidence_anchor_provenance.tsv"
    fieldnames = list(PROVENANCE_ROWS[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(PROVENANCE_ROWS)
    return path


def source_path(source_name: str) -> Path:
    if source_name == GENERATED_PROVENANCE_SOURCE:
        return write_result_evidence_anchor_provenance()
    return RESULTS / source_name


def write_submission_guide() -> None:
    fieldnames = [
        "table_range",
        "submission_packet",
        "evidence_tier",
        "reviewer_priority",
        "topic_group",
        "recommended_use",
    ]
    with GUIDE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for group in SUBMISSION_GROUPS:
            writer.writerow(
                {
                    "table_range": table_range(group),
                    "submission_packet": group["submission_packet"],
                    "evidence_tier": group["evidence_tier"],
                    "reviewer_priority": group["reviewer_priority"],
                    "topic_group": group["topic_group"],
                    "recommended_use": group["recommended_use"],
                }
            )


def write_zip_packets(rows: list[dict[str, object]]) -> list[Path]:
    common_files = [README, INDEX, GUIDE]
    written = []
    for zip_name, selector in ZIP_PACKETS:
        zip_path = OUT_DIR / zip_name
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for common_file in common_files:
                archive.write(common_file, arcname=common_file.name)
            for row in rows:
                if selector(row):
                    table_path = ROOT / str(row["copied_file"])
                    archive.write(table_path, arcname=table_path.name)
        written.append(zip_path)
    return written


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cleanup_generated_outputs()
    rows = []
    for i, (slug, title, manuscript_location, description, source_name) in enumerate(TABLES, start=1):
        source = source_path(source_name)
        if not source.exists():
            raise FileNotFoundError(source)
        dest_name = destination_name(i, slug, source_name)
        dest = OUT_DIR / dest_name
        if source.resolve() != dest.resolve():
            shutil.copy2(source, dest)
        metadata = submission_metadata(i)
        rows.append(
            {
                "table_id": f"Table S{i}",
                "file_name": dest_name,
                "title": title,
                **metadata,
                "manuscript_location": manuscript_location,
                "description": description,
                "source_file": str(source.relative_to(ROOT)),
                "copied_file": str(dest.relative_to(ROOT)),
                "file_size_bytes": dest.stat().st_size,
            }
        )

    with INDEX.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_submission_guide()

    readme_lines = [
        "# Supplementary Tables",
        "",
        "Ordered manuscript support tables copied from `Project/results`.",
        "",
        "These files are copies, not replacements; the original analysis outputs remain in `Project/results`.",
        "",
        "## Submission Strategy",
        "",
        "The 103 numbered tables are retained for reproducibility, but they are triaged so reviewers do not have to treat every table as equally central.",
        "",
        "- Tables S1-S40 and S102-S103 are the reader-facing primary supplementary tables tied to the main figures.",
        "- Tables S41-S53 and S89-S93 are reader-facing supplementary controls tied to supplementary figures.",
        "- Tables S54-S88 and S94-S101 are best treated as a machine-readable Supplementary Data archive and provenance packet.",
        "",
        "## Packet Summary",
        "",
        "| Tables | Packet | Tier | Topic | Recommended use |",
        "|---|---|---|---|---|",
    ]
    for group in SUBMISSION_GROUPS:
        readme_lines.append(
            "| {table_range} | {submission_packet} | {evidence_tier} | {topic_group} | {recommended_use} |".format(
                table_range=table_range(group),
                submission_packet=group["submission_packet"],
                evidence_tier=group["evidence_tier"],
                topic_group=group["topic_group"],
                recommended_use=group["recommended_use"],
            )
        )
    readme_lines.extend(
        [
            "",
            "## Detailed Index",
            "",
        ]
    )
    current_packet = None
    current_topic = None
    for row in rows:
        if row["submission_packet"] != current_packet:
            current_packet = row["submission_packet"]
            readme_lines.extend(["", f"### {current_packet}", ""])
            current_topic = None
        if row["topic_group"] != current_topic:
            current_topic = row["topic_group"]
            readme_lines.extend(["", f"#### {current_topic}", ""])
        readme_lines.append(
            f"- {row['table_id']}: `{row['file_name']}` - {row['title']} ({row['manuscript_location']}; {row['evidence_tier']})."
        )
    readme_lines.extend(
        [
            "",
            "## Prebuilt Upload Packets",
            "",
            "- `Supplementary_Tables_Primary_Reader_Facing_S1_S40_S102_S103.zip`: main-figure, reader-facing result tables.",
            "- `Supplementary_Tables_Supplementary_Controls_S41_S53_S89_S93.zip`: supplementary-figure control tables.",
            "- `Supplementary_Data_Archive_S54_S88_S94_S101.zip`: machine-readable regulatory, perturbation, audit and provenance archive.",
            "",
            "The ZIP files are convenience copies only. The individual numbered tables remain the authoritative files.",
            "",
            f"Machine-readable table index: `{INDEX.name}`",
            f"Submission grouping guide: `{GUIDE.name}`",
            "",
        ]
    )
    README.write_text("\n".join(readme_lines))
    written_zips = write_zip_packets(rows)

    print(f"Wrote {len(rows)} supplementary tables to {OUT_DIR.relative_to(ROOT)}")
    print(f"Wrote {INDEX.relative_to(ROOT)}")
    print(f"Wrote {GUIDE.relative_to(ROOT)}")
    print(f"Wrote {README.relative_to(ROOT)}")
    for zip_path in written_zips:
        print(f"Wrote {zip_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
