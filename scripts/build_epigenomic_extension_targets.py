#!/usr/bin/env python3
"""Build epigenomic-extension target and model-term tables.

This is a no-heavy-download scaffold. It uses already curated transcriptomic
modules and candidate tiers to define which genes/modules should be tested once
ATAC, multiome, methylome, or spatial epigenomic matrices are processed.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

MODULES = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
ORIGIN_MODULES = RESULTS / "developmental_origin_divergence_audit_gene_sets.tsv"
CANDIDATE_TIERS = RESULTS / "primary_core_manuscript_candidate_tiers.tsv"
RESOURCE_TABLE = RESULTS / "epigenomic_extension_candidate_resources.tsv"
GSE268609_FEATURES = EXTERNAL / "GEO/GSE268609/GSE268609_features.tsv.gz"
GSE268609_MATRIX = EXTERNAL / "GEO/GSE268609/GSE268609_matrix.mtx.gz"
GSE322785_GENE_SUMMARY = RESULTS / "gse322785_human_h5_epigenomic_gene_summary.tsv"
GSE322785_PEAK_TARGETS = RESULTS / "gse322785_human_h5_epigenomic_peak_targets.tsv"
GSE322785_MANIFEST = RESULTS / "gse322785_human_h5_epigenomic_selective_manifest.tsv"
GSE322785_SELECTED_SUMMARY = RESULTS / "gse322785_human_h5_selected_matrix_summary.tsv"
GSE322785_MARKER_SCORING = RESULTS / "gse322785_human_h5_epigenomic_marker_group_module_scores.tsv"
GSE322785_CLUSTER_METRICS = RESULTS / "gse322785_human_h5_cluster_validation_metrics.tsv"
GSE322785_SUPPORTED_BARCODES = RESULTS / "gse322785_human_h5_cluster_supported_marker_barcodes.tsv.gz"
GSE322785_SUPPORTED_MODULE = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv"
GSE322785_ROBUST = RESULTS / "gse322785_epigenomic_robust_positive_contrasts.tsv"

OUT_TARGETS = RESULTS / "epigenomic_extension_regulatory_targets.tsv"
OUT_MODEL_TERMS = RESULTS / "integrative_granule_model_term_specification.tsv"
OUT_SUMMARY = RESULTS / "epigenomic_extension_target_summary.tsv"
OUT_MD = RESULTS / "epigenomic_extension_target_model.md"


MODEL_TERM_BY_MODULE_FAMILY = {
    "upstream_region_fate": "FatePolarity",
    "downstream_circuit_morphology": "ConstructionBalance",
    "shared_niche_state": "NicheSignal",
    "deep_origin": "FatePolarity",
    "regional_origin": "FatePolarity",
    "postmitotic_construction": "ConstructionBalance",
}

TF_OR_REGULATORY_GENES = {
    "ATOH1",
    "ASCL1",
    "BCL11B",
    "EOMES",
    "EMX1",
    "EMX2",
    "EN1",
    "EN2",
    "FOXG1",
    "GLI1",
    "GLI2",
    "HES1",
    "HES5",
    "LEF1",
    "LHX2",
    "MEIS1",
    "MYCN",
    "NEUROD1",
    "NEUROD2",
    "NFIA",
    "NFIB",
    "OTX2",
    "PAX6",
    "PROX1",
    "RFX3",
    "RFX7",
    "SOX1",
    "SOX2",
    "SOX3",
    "TCF7L2",
    "TBR1",
    "ZBTB20",
    "ZIC1",
    "ZIC2",
    "ZIC3",
}


def read_gse268609_gene_features() -> dict[str, dict[str, object]]:
    if not GSE268609_FEATURES.exists():
        return {}

    records: dict[str, dict[str, object]] = {}
    with gzip.open(GSE268609_FEATURES, "rt") as handle:
        for idx, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            feature_id, name, feature_type, chrom, start, end = parts[:6]
            if feature_type != "Gene Expression":
                continue
            gene = name.upper()
            if gene in records:
                records[gene]["n_matching_gene_rows"] += 1
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                start_i = None
                end_i = None
            if start_i is not None and end_i is not None:
                promoter_gene_window = f"{chrom}:{max(1, start_i - 2000)}-{end_i + 2000}"
            else:
                promoter_gene_window = ""
            records[gene] = {
                "gse268609_feature_id": feature_id,
                "gse268609_feature_name": name,
                "gse268609_feature_row_1based": idx,
                "gse268609_chrom": chrom,
                "gse268609_start": start_i,
                "gse268609_end": end_i,
                "gse268609_promoter_gene_window_plusminus_2kb": promoter_gene_window,
                "n_matching_gene_rows": 1,
            }
    return records


def add_feature_info(row: dict[str, object], feature_map: dict[str, dict[str, object]]) -> dict[str, object]:
    info = feature_map.get(str(row["gene"]).upper())
    if info is None:
        row.update(
            {
                "present_as_gse268609_gene_expression_feature": False,
                "gse268609_feature_id": "",
                "gse268609_feature_name": "",
                "gse268609_feature_row_1based": "",
                "gse268609_chrom": "",
                "gse268609_start": "",
                "gse268609_end": "",
                "gse268609_promoter_gene_window_plusminus_2kb": "",
                "n_matching_gene_rows": 0,
            }
        )
    else:
        row["present_as_gse268609_gene_expression_feature"] = True
        row.update(info)
    return row


def build_targets() -> pd.DataFrame:
    feature_map = read_gse268609_gene_features()
    rows: list[dict[str, object]] = []

    modules = pd.read_csv(MODULES, sep="\t")
    for rec in modules.to_dict("records"):
        family = rec["module_family"]
        model_term = MODEL_TERM_BY_MODULE_FAMILY.get(family, "EpigenomicCompatibility")
        gene = str(rec["gene"]).upper()
        row = {
            "target_source": "niche_circuit_module",
            "target_set": rec["module_id"],
            "target_label": rec["module_label"],
            "model_term_supported": model_term,
            "gene": gene,
            "mouse_symbol": rec.get("default_mouse_symbol", ""),
            "priority": "high" if model_term in {"FatePolarity", "ConstructionBalance", "NicheSignal"} else "medium",
            "is_tf_or_regulatory_candidate": gene in TF_OR_REGULATORY_GENES,
            "recommended_epigenomic_signal": "module-level promoter/gene-body accessibility plus linked enhancer accessibility",
            "rationale": rec.get("hypothesis_role", ""),
            "current_evidence_source": str(MODULES.relative_to(ROOT)),
        }
        rows.append(add_feature_info(row, feature_map))

    origin = pd.read_csv(ORIGIN_MODULES, sep="\t")
    for rec in origin.to_dict("records"):
        family = rec["module_family"]
        model_term = MODEL_TERM_BY_MODULE_FAMILY.get(family, "EpigenomicCompatibility")
        gene = str(rec["gene"]).upper()
        row = {
            "target_source": "developmental_origin_module",
            "target_set": rec["module_id"],
            "target_label": rec["module_label"],
            "model_term_supported": model_term,
            "gene": gene,
            "mouse_symbol": rec.get("default_mouse_symbol", ""),
            "priority": "medium_high" if model_term in {"FatePolarity", "ConstructionBalance"} else "medium",
            "is_tf_or_regulatory_candidate": gene in TF_OR_REGULATORY_GENES,
            "recommended_epigenomic_signal": "stage-linked accessibility or methylation support near lineage and construction genes",
            "rationale": rec.get("interpretation_role", ""),
            "current_evidence_source": str(ORIGIN_MODULES.relative_to(ROOT)),
        }
        rows.append(add_feature_info(row, feature_map))

    tiers = pd.read_csv(CANDIDATE_TIERS, sep="\t")
    for rec in tiers.to_dict("records"):
        gene = str(rec["gene"]).upper()
        tier_rank = int(rec["tier_rank"])
        if tier_rank <= 1:
            priority = "highest"
        elif tier_rank <= 2:
            priority = "high"
        elif tier_rank <= 3:
            priority = "medium_high"
        else:
            priority = "medium"
        row = {
            "target_source": "manuscript_candidate_tier",
            "target_set": str(rec["manuscript_tier"]),
            "target_label": str(rec["manuscript_tier"]),
            "model_term_supported": "ConstructionBalance",
            "gene": gene,
            "mouse_symbol": rec.get("mouse_symbol", ""),
            "priority": priority,
            "is_tf_or_regulatory_candidate": gene in TF_OR_REGULATORY_GENES
            or "regulatory" in str(rec.get("mechanism_class", "")),
            "recommended_epigenomic_signal": "candidate gene promoter/gene-body accessibility; linked enhancer accessibility; motif support for regulatory candidates",
            "rationale": rec.get("suggested_role", ""),
            "current_evidence_source": str(CANDIDATE_TIERS.relative_to(ROOT)),
        }
        rows.append(add_feature_info(row, feature_map))

    targets = pd.DataFrame(rows)
    targets = targets.sort_values(
        ["priority", "target_source", "target_set", "gene"],
        key=lambda col: col.map(
            {
                "highest": 0,
                "high": 1,
                "medium_high": 2,
                "medium": 3,
            }
        )
        if col.name == "priority"
        else col,
    )
    return targets


def build_model_terms() -> pd.DataFrame:
    rows = [
        {
            "symbol": "Y",
            "term": "GranuleDesign",
            "meaning": "response variable: granule-like transcriptomic configuration or design score",
            "current_data_source": "primary_core_transcriptomic_configuration_primary_contrasts.tsv; primary_core_configuration_driver_audit_contrasts.tsv",
            "current_status": "fitted in current manuscript as transcriptomic configuration delta",
            "recommended_next_action": "Use as the dependent variable for an integrative mixed model.",
        },
        {
            "symbol": "F",
            "term": "FatePolarity",
            "meaning": "branch-matched regional fate rank minus opposed fate rank",
            "current_data_source": "primary_core_transcriptomic_configuration_primary_contrasts.tsv",
            "current_status": "direct transcriptomic term already available",
            "recommended_next_action": "Retain as core predictor and test interaction with construction balance.",
        },
        {
            "symbol": "C",
            "term": "ConstructionBalance",
            "meaning": "downstream neurite/synapse/excitability rank relative to progenitor/niche state",
            "current_data_source": "primary_core_transcriptomic_configuration_primary_contrasts.tsv; primary_core_niche_circuit_module_model.md",
            "current_status": "direct transcriptomic term already available",
            "recommended_next_action": "Retain as core predictor and test interaction with circuit/resource terms.",
        },
        {
            "symbol": "T",
            "term": "Stage/Pseudotime",
            "meaning": "normalized developmental stage or diffusion/pseudotime axis",
            "current_data_source": "aim2_stage_window_model_group_fits.tsv; primary_core_full_transcriptome_diffusion_pseudotime_scatter.png",
            "current_status": "available as fitted stage-window and diffusion support, but not paired with every contrast row",
            "recommended_next_action": "Use dataset/stage-level mapping where available; otherwise include as a sensitivity layer.",
        },
        {
            "symbol": "N",
            "term": "NicheSignal",
            "meaning": "pathway-readiness or sender-receiver ligand-receptor support",
            "current_data_source": "primary_core_aim2_niche_pathway_model.md; aim2_sender_receiver_lr_summary.tsv",
            "current_status": "available as pathway and directional niche evidence",
            "recommended_next_action": "Use as branch/dataset-level predictor; avoid causal wording without perturbation.",
        },
        {
            "symbol": "E",
            "term": "EpigenomicCompatibility",
            "meaning": "accessibility, methylation, or motif support near fate/construction/candidate genes",
            "current_data_source": "epigenomic_extension_regulatory_targets.tsv; epigenomic_extension_candidate_resources.tsv; gse268609_epigenomic_peak_targets.tsv; gse322785_human_h5_epigenomic_peak_targets.tsv; gse322785_human_h5_epigenomic_marker_group_module_scores.tsv; gse322785_human_h5_cluster_validation_metrics.tsv; gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv; gse322785_epigenomic_robust_positive_contrasts.tsv",
            "current_status": "feature/peak target manifests available for GSE268609; downloaded human GSE322785 H5 files have selected count matrices, provisional marker-group RNA/ATAC scoring, selected-gene cluster validation, stricter cluster-supported sensitivity scoring, and broad-vs-strict robust contrast classification; not yet source-label verified",
            "recommended_next_action": "Verify GSE322785 labels with clustering or source taxonomy and restore GSE268609 full matrix for matched dentate peak-count scoring.",
        },
        {
            "symbol": "M",
            "term": "MorphologySparseSampling",
            "meaning": "external morphology support for compact input-sampling architecture",
            "current_data_source": "neuromorpho_granule_morphometry_summary.tsv; aim3_empirical_calibration_grid.tsv",
            "current_status": "external calibration/proxy, not matched to transcriptome rows",
            "recommended_next_action": "Keep as calibration term or branch-level prior.",
        },
        {
            "symbol": "A",
            "term": "ActivitySparsity",
            "meaning": "activity sparsity, spatial selectivity, and population-vector support",
            "current_data_source": "dandi_000003_multisession_spatial_celltype_pooled.tsv; dandi_000003_multisession_population_vector_separation.tsv",
            "current_status": "dentate external calibration/proxy, not cerebellar physiology",
            "recommended_next_action": "Use as asymmetric dentate support and avoid claiming direct cerebellar validation.",
        },
        {
            "symbol": "R",
            "term": "CircuitResourceConstraint",
            "meaning": "expansion ratio, input degree, output sparsity, and resource-adjusted useful score",
            "current_data_source": "primary_core_aim3_sparse_coding_model.md; aim3_empirical_calibration_grid.tsv",
            "current_status": "simulation and empirical calibration term",
            "recommended_next_action": "Keep as computational constraint term, not as measured developmental cause.",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    targets = build_targets()
    targets.to_csv(OUT_TARGETS, sep="\t", index=False)

    terms = build_model_terms()
    terms.to_csv(OUT_MODEL_TERMS, sep="\t", index=False)

    summary = (
        targets.groupby(["target_source", "model_term_supported", "priority"], dropna=False)
        .agg(
            n_rows=("gene", "size"),
            n_unique_genes=("gene", "nunique"),
            n_gse268609_gene_expression_present=("present_as_gse268609_gene_expression_feature", "sum"),
            n_regulatory_candidates=("is_tf_or_regulatory_candidate", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    resources = pd.read_csv(RESOURCE_TABLE, sep="\t") if RESOURCE_TABLE.exists() else pd.DataFrame()
    top_resources = []
    if not resources.empty:
        for rec in resources.head(4).to_dict("records"):
            top_resources.append(
                f"- `{rec['dataset_or_resource']}` ({rec['priority']}): {rec['recommended_use']}"
            )

    n_unique = targets["gene"].nunique()
    n_present = int(
        targets.drop_duplicates("gene")["present_as_gse268609_gene_expression_feature"].sum()
    )
    matrix_status = "present" if GSE268609_MATRIX.exists() else "not present locally"
    if GSE322785_GENE_SUMMARY.exists():
        gse322785_gene_summary = pd.read_csv(GSE322785_GENE_SUMMARY, sep="\t")
        gse322785_unique = gse322785_gene_summary.loc[
            gse322785_gene_summary["present_in_h5"].astype(bool), "gene"
        ].nunique()
        gse322785_gene_rows = int(gse322785_gene_summary["present_in_h5"].sum())
        gse322785_gene_tests = len(gse322785_gene_summary)
    else:
        gse322785_unique = 0
        gse322785_gene_rows = 0
        gse322785_gene_tests = 0

    gse322785_peak_rows = (
        sum(1 for _ in GSE322785_PEAK_TARGETS.open()) - 1
        if GSE322785_PEAK_TARGETS.exists()
        else 0
    )
    gse322785_manifest_rows = (
        sum(1 for _ in GSE322785_MANIFEST.open()) - 1
        if GSE322785_MANIFEST.exists()
        else 0
    )
    if GSE322785_SELECTED_SUMMARY.exists():
        gse322785_selected = pd.read_csv(GSE322785_SELECTED_SUMMARY, sep="\t")
        gse322785_processed_barcodes = int(gse322785_selected["n_processed_barcodes"].sum())
        gse322785_basic_qc_barcodes = int(gse322785_selected["n_basic_qc_barcodes"].sum())
        gse322785_marker_calls = int(gse322785_selected["n_high_or_medium_confidence_calls"].sum())
        gse322785_selected_nnz = int(gse322785_selected["selected_matrix_nnz"].sum())
    else:
        gse322785_processed_barcodes = 0
        gse322785_basic_qc_barcodes = 0
        gse322785_marker_calls = 0
        gse322785_selected_nnz = 0
    gse322785_marker_module_rows = (
        sum(1 for _ in GSE322785_MARKER_SCORING.open()) - 1
        if GSE322785_MARKER_SCORING.exists()
        else 0
    )
    if GSE322785_CLUSTER_METRICS.exists():
        cluster_metrics = pd.read_csv(GSE322785_CLUSTER_METRICS, sep="\t")
        gse322785_mean_ari = float(cluster_metrics["adjusted_rand_marker_call_vs_cluster"].mean())
        gse322785_mean_nmi = float(cluster_metrics["normalized_mutual_info_marker_call_vs_cluster"].mean())
    else:
        gse322785_mean_ari = 0.0
        gse322785_mean_nmi = 0.0
    if GSE322785_SUPPORTED_BARCODES.exists():
        supported_barcodes = pd.read_csv(GSE322785_SUPPORTED_BARCODES, sep="\t")
        gse322785_supported_barcodes = len(supported_barcodes)
        gse322785_supported_groups = supported_barcodes["marker_call"].nunique()
        gse322785_supported_granule = int(
            supported_barcodes["marker_call"].eq("cerebellar_granule_candidate").sum()
        )
    else:
        gse322785_supported_barcodes = 0
        gse322785_supported_groups = 0
        gse322785_supported_granule = 0
    gse322785_supported_module_rows = (
        sum(1 for _ in GSE322785_SUPPORTED_MODULE.open()) - 1
        if GSE322785_SUPPORTED_MODULE.exists()
        else 0
    )
    if GSE322785_ROBUST.exists():
        robust = pd.read_csv(GSE322785_ROBUST, sep="\t")
        gse322785_robust_positive = len(robust)
        gse322785_robust_strong = int(robust["concordance_class"].eq("robust_positive_strong").sum())
    else:
        gse322785_robust_positive = 0
        gse322785_robust_strong = 0

    lines = [
        "# Epigenomic Extension Target Model",
        "",
        "Date built: 2026-06-26",
        "",
        "## Purpose",
        "",
        "This selective-download scaffold defines the regulatory targets and model terms needed to add a matched or comparable scATAC/multiome/methylome layer to the granule-cell convergence model.",
        "",
        "## Target Summary",
        "",
        f"- Target rows: {len(targets)}.",
        f"- Unique genes: {n_unique}.",
        f"- Unique genes present as `GSE268609` gene-expression features: {n_present}/{n_unique}.",
        f"- Full `GSE268609` matrix status: {matrix_status}.",
        f"- Unique target genes present in downloaded human `GSE322785` H5 files: {gse322785_unique}/{n_unique}.",
        f"- Per-sample target-gene rows present in downloaded human `GSE322785` H5 files: {gse322785_gene_rows}/{gse322785_gene_tests}.",
        f"- Human `GSE322785` peak-target rows: {gse322785_peak_rows}.",
        f"- Human `GSE322785` selective manifest rows: {gse322785_manifest_rows}.",
        f"- Human `GSE322785` barcodes processed for selected count extraction: {gse322785_processed_barcodes}.",
        f"- Human `GSE322785` basic-QC barcodes: {gse322785_basic_qc_barcodes}.",
        f"- Human `GSE322785` high/medium-confidence provisional marker calls: {gse322785_marker_calls}.",
        f"- Human `GSE322785` selected matrix nonzero entries: {gse322785_selected_nnz}.",
        f"- Human `GSE322785` provisional marker-group module-score rows: {gse322785_marker_module_rows}.",
        f"- Human `GSE322785` selected-gene cluster validation mean ARI: {gse322785_mean_ari:.3f}.",
        f"- Human `GSE322785` selected-gene cluster validation mean NMI: {gse322785_mean_nmi:.3f}.",
        f"- Human `GSE322785` cluster-supported barcodes: {gse322785_supported_barcodes}.",
        f"- Human `GSE322785` cluster-supported marker groups: {gse322785_supported_groups}.",
        f"- Human `GSE322785` cluster-supported granule-candidate barcodes: {gse322785_supported_granule}.",
        f"- Human `GSE322785` cluster-supported module-score rows: {gse322785_supported_module_rows}.",
        f"- Human `GSE322785` robust-positive broad-vs-supported contrasts: {gse322785_robust_positive}.",
        f"- Human `GSE322785` strong robust-positive broad-vs-supported contrasts: {gse322785_robust_strong}.",
        "",
        "## Integrative Model",
        "",
        "`GranuleDesign_i = beta0 + betaF FatePolarity_i + betaC ConstructionBalance_i + betaT Stage_i + betaN NicheSignal_i + betaE EpigenomicCompatibility_i + betaM Morphology_i + betaA Activity_i + betaR CircuitConstraint_i + interactions + random effects + error_i`.",
        "",
        "The new epigenomic term is `EpigenomicCompatibility`, estimated from promoter/gene-body accessibility, linked enhancer accessibility, methylation/accessibility support, and motif deviation near fate, construction, and candidate genes.",
        "",
        "## Top Resource Route",
        "",
        *(top_resources or ["- No resource table found."]),
        "",
        "## Outputs",
        "",
        f"- Regulatory targets: `{OUT_TARGETS.relative_to(ROOT)}`",
        f"- Model-term specification: `{OUT_MODEL_TERMS.relative_to(ROOT)}`",
        f"- Target summary: `{OUT_SUMMARY.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
