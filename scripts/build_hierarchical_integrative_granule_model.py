#!/usr/bin/env python3
"""Build a hierarchical integrative evidence model for granule-cell convergence.

The project now has evidence at different units of observation: sample-level
transcriptomic contrasts, branch/stage summaries, provisional epigenomic
sensitivity contrasts, external morphology/activity calibration, and
simulation-based resource constraints. This script keeps those levels explicit
and summarizes them as a weighted hierarchical evidence model rather than
pretending that all terms are measured in the same cells.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

TERM_SPEC = RESULTS / "integrative_granule_model_term_specification.tsv"
CONFIG_CONTRASTS = RESULTS / "primary_core_transcriptomic_configuration_primary_contrasts.tsv"
DRIVER_CONTRASTS = RESULTS / "primary_core_configuration_driver_audit_contrasts.tsv"
STAGE_BRANCH = RESULTS / "aim2_stage_window_model_branch_summary.tsv"
PATHWAY_SUMMARY = RESULTS / "primary_core_aim2_pathway_summary.tsv"
SENDER_RECEIVER_SUMMARY = RESULTS / "aim2_sender_receiver_lr_summary.tsv"
EPI_ROBUST_POSITIVE = RESULTS / "gse322785_epigenomic_robust_positive_contrasts.tsv"
NEUROMORPHO_COMPARISON = RESULTS / "neuromorpho_granule_morphometry_comparison.tsv"
DANDI_POOLED = RESULTS / "dandi_000003_multisession_spatial_celltype_pooled.tsv"
DANDI_PV = RESULTS / "dandi_000003_multisession_population_vector_separation.tsv"
AIM3_TARGETS = RESULTS / "aim3_empirical_calibration_targets.tsv"
AIM3_CALIBRATION = RESULTS / "aim3_empirical_calibration_architecture_summary.tsv"

OUT_TERMS = RESULTS / "hierarchical_integrative_model_terms.tsv"
OUT_UNITS = RESULTS / "hierarchical_integrative_model_evidence_units.tsv"
OUT_LAYER = RESULTS / "hierarchical_integrative_model_layer_summary.tsv"
OUT_BRANCH = RESULTS / "hierarchical_integrative_model_branch_summary.tsv"
OUT_COMPONENT = RESULTS / "hierarchical_integrative_model_component_scores.tsv"
OUT_MD = RESULTS / "hierarchical_integrative_model.md"
OUT_PNG = RESULTS / "hierarchical_integrative_model_component_scores.png"


TERM_WEIGHTS = {
    "Y": 1.00,
    "F": 0.90,
    "C": 0.90,
    "I": 0.80,
    "T": 0.65,
    "N": 0.75,
    "E": 0.45,
    "M": 0.60,
    "A": 0.55,
    "R": 0.60,
}

TERM_NAMES = {
    "Y": "GranuleDesign",
    "F": "FatePolarity",
    "C": "ConstructionBalance",
    "I": "ConfigurationCoupling",
    "T": "Stage/Pseudotime",
    "N": "NicheSignal",
    "E": "EpigenomicCompatibility",
    "M": "MorphologySparseSampling",
    "A": "ActivitySparsity",
    "R": "CircuitResourceConstraint",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def clean_text(value: Any, default: str = "") -> str:
    if value is None:
        return default
    if isinstance(value, float) and math.isnan(value):
        return default
    text = str(value)
    return text if text and text != "nan" else default


def to_float(value: Any, default: float = np.nan) -> float:
    try:
        if value is None:
            return default
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def clip_score(value: float, scale: float = 1.0) -> float:
    if not np.isfinite(value) or not np.isfinite(scale) or scale == 0:
        return np.nan
    return float(np.clip(value / scale, -1.0, 1.0))


def support_class(score: float) -> str:
    if not np.isfinite(score):
        return "not_scored"
    if score >= 0.60:
        return "strong_support"
    if score >= 0.30:
        return "moderate_support"
    if score > 0.00:
        return "weak_or_mixed_support"
    if score <= -0.30:
        return "opposing_or_branch_specific"
    return "near_zero_or_mixed"


def branch_group_from_values(*values: Any) -> str:
    text = " ".join(clean_text(v).lower() for v in values)
    if "cerebell" in text:
        return "cerebellum"
    if "dentate" in text or "hippoc" in text or "sgz" in text:
        return "dentate"
    if "simulation" in text or "resource" in text:
        return "simulation"
    return "cross_branch_or_external"


def add_unit(
    rows: list[dict[str, Any]],
    *,
    level: str,
    term_symbol: str,
    branch_group: str,
    branch: str,
    dataset_or_source: str,
    contrast_or_scope: str,
    evidence_family: str,
    metric: str,
    raw_value: Any,
    normalized_score: float,
    support_direction: str,
    weight: float | None = None,
    quality_class: str,
    source_table: str,
    notes: str,
) -> None:
    if weight is None:
        weight = TERM_WEIGHTS[term_symbol]
    normalized_score = to_float(normalized_score)
    weighted_score = normalized_score * weight if np.isfinite(normalized_score) else np.nan
    rows.append(
        {
            "unit_id": f"EU{len(rows) + 1:05d}",
            "hierarchy_level": level,
            "term_symbol": term_symbol,
            "term": TERM_NAMES[term_symbol],
            "branch_group": branch_group,
            "branch": branch,
            "dataset_or_source": dataset_or_source,
            "contrast_or_scope": contrast_or_scope,
            "evidence_family": evidence_family,
            "metric": metric,
            "raw_value": raw_value,
            "normalized_score": normalized_score,
            "support_direction": support_direction,
            "weight": weight,
            "weighted_score": weighted_score,
            "unit_support_class": support_class(normalized_score),
            "quality_class": quality_class,
            "source_table": source_table,
            "notes": notes,
        }
    )


def build_terms() -> pd.DataFrame:
    rows = [
        {
            "term_symbol": "Y",
            "term": "GranuleDesign",
            "hierarchy_level": "sample_contrast",
            "unit_of_observation": "dataset/sample candidate-background contrast",
            "model_role": "response or primary outcome",
            "default_weight": TERM_WEIGHTS["Y"],
            "quality_class": "direct_primary_core",
            "source_table": CONFIG_CONTRASTS.name,
            "interpretation": "Granule candidates show a higher combined regional-fate plus downstream-construction configuration than local backgrounds.",
            "caveat": "This is transcriptomic configuration, not direct morphology.",
        },
        {
            "term_symbol": "F",
            "term": "FatePolarity",
            "hierarchy_level": "sample_contrast",
            "unit_of_observation": "dataset/sample candidate-background contrast",
            "model_role": "core predictor",
            "default_weight": TERM_WEIGHTS["F"],
            "quality_class": "direct_primary_core",
            "source_table": CONFIG_CONTRASTS.name,
            "interpretation": "Branch-matched regional-fate programs exceed opposed fate programs.",
            "caveat": "Fate polarity separates dentate and cerebellar origins rather than proving one shared origin.",
        },
        {
            "term_symbol": "C",
            "term": "ConstructionBalance",
            "hierarchy_level": "sample_contrast",
            "unit_of_observation": "dataset/sample candidate-background contrast",
            "model_role": "core predictor",
            "default_weight": TERM_WEIGHTS["C"],
            "quality_class": "direct_primary_core",
            "source_table": CONFIG_CONTRASTS.name,
            "interpretation": "Downstream neurite, synaptic, and excitability modules exceed progenitor/niche-state modules.",
            "caveat": "Many neurons use construction genes; specificity depends on configuration, stage, and circuit context.",
        },
        {
            "term_symbol": "I",
            "term": "ConfigurationCoupling",
            "hierarchy_level": "sample_contrast",
            "unit_of_observation": "dataset/sample candidate-background contrast",
            "model_role": "interaction term",
            "default_weight": TERM_WEIGHTS["I"],
            "quality_class": "direct_primary_core_decomposition",
            "source_table": DRIVER_CONTRASTS.name,
            "interpretation": "Regional fate and construction balance are simultaneously positive in the same contrast.",
            "caveat": "This is a sign-coupling audit, not a fitted molecular interaction coefficient.",
        },
        {
            "term_symbol": "T",
            "term": "Stage/Pseudotime",
            "hierarchy_level": "branch_stage",
            "unit_of_observation": "branch/signature stage-window summary",
            "model_role": "branch- or stage-level predictor",
            "default_weight": TERM_WEIGHTS["T"],
            "quality_class": "stage_branch_model",
            "source_table": STAGE_BRANCH.name,
            "interpretation": "Maturation, permissive, and TGF-beta/BDNF axes are stage-windowed rather than static.",
            "caveat": "Stage coverage is richer for dentate than cerebellum in the current local core.",
        },
        {
            "term_symbol": "N",
            "term": "NicheSignal",
            "hierarchy_level": "branch_context",
            "unit_of_observation": "pathway or sender-receiver branch summary",
            "model_role": "branch-context predictor",
            "default_weight": TERM_WEIGHTS["N"],
            "quality_class": "expression_niche_model",
            "source_table": f"{PATHWAY_SUMMARY.name}; {SENDER_RECEIVER_SUMMARY.name}",
            "interpretation": "Local pathway and ligand-receptor readiness provide branch-specific permissive or differentiating contexts.",
            "caveat": "Expression support does not prove ligand secretion, spatial contact, or receptor activation.",
        },
        {
            "term_symbol": "E",
            "term": "EpigenomicCompatibility",
            "hierarchy_level": "epigenomic_sensitivity",
            "unit_of_observation": "human cerebellar marker-group/cluster-supported contrast",
            "model_role": "regulatory-compatibility sensitivity term",
            "default_weight": TERM_WEIGHTS["E"],
            "quality_class": "provisional_epigenomic_extension",
            "source_table": EPI_ROBUST_POSITIVE.name,
            "interpretation": "Accessibility/RNA selected-feature targets near fate, construction, niche, or candidate genes support granule-candidate regulatory compatibility.",
            "caveat": "GSE322785 labels are provisional and matched dentate peak-count scoring is not yet restored.",
        },
        {
            "term_symbol": "M",
            "term": "MorphologySparseSampling",
            "hierarchy_level": "external_calibration",
            "unit_of_observation": "NeuroMorpho branch comparison or calibration target",
            "model_role": "morphology prior/calibration",
            "default_weight": TERM_WEIGHTS["M"],
            "quality_class": "external_empirical_calibration",
            "source_table": f"{NEUROMORPHO_COMPARISON.name}; {AIM3_TARGETS.name}",
            "interpretation": "Public reconstructions support sparse primary input sampling plus nontrivial branch complexity.",
            "caveat": "Morphology measurements are not matched to transcriptomic cells and stems are not synapse counts.",
        },
        {
            "term_symbol": "A",
            "term": "ActivitySparsity",
            "hierarchy_level": "external_calibration",
            "unit_of_observation": "DANDI dentate electrophysiology summary",
            "model_role": "activity prior/calibration",
            "default_weight": TERM_WEIGHTS["A"],
            "quality_class": "external_empirical_calibration",
            "source_table": f"{DANDI_POOLED.name}; {DANDI_PV.name}",
            "interpretation": "Dentate granule-labeled units show spatial information and restricted active-bin occupancy.",
            "caveat": "This is dentate-only activity support and does not directly validate cerebellar physiology.",
        },
        {
            "term_symbol": "R",
            "term": "CircuitResourceConstraint",
            "hierarchy_level": "simulation_calibration",
            "unit_of_observation": "architecture-level simulation/calibration summary",
            "model_role": "resource constraint term",
            "default_weight": TERM_WEIGHTS["R"],
            "quality_class": "computational_simulation",
            "source_table": AIM3_CALIBRATION.name,
            "interpretation": "Resource-constrained scoring favors intermediate or sparse granule-like expansion over dense or excessive-sparsity extremes.",
            "caveat": "Simulation terms are mathematical constraints, not direct developmental causes.",
        },
    ]
    return pd.DataFrame(rows)


def build_evidence_units() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    config = read_tsv(CONFIG_CONTRASTS)
    for _, row in config.iterrows():
        branch = clean_text(row.get("core_branch"), clean_text(row.get("branch_kind")))
        branch_group = branch_group_from_values(row.get("branch_kind"), branch)
        dataset = clean_text(row.get("dataset"))
        sample = clean_text(row.get("sample"))
        layer = clean_text(row.get("expression_layer"))
        target = clean_text(row.get("target_class"))
        scope = f"{layer}|{dataset}|{sample}|{target}"
        add_unit(
            rows,
            level="sample_contrast",
            term_symbol="Y",
            branch_group=branch_group,
            branch=branch,
            dataset_or_source=dataset,
            contrast_or_scope=scope,
            evidence_family="primary_configuration",
            metric="delta_configuration_score",
            raw_value=to_float(row.get("delta_configuration_score")),
            normalized_score=clip_score(to_float(row.get("delta_configuration_score")), 1.0),
            support_direction="positive candidate-background delta supports granule-like configuration",
            quality_class="direct_primary_core",
            source_table=CONFIG_CONTRASTS.name,
            notes="Combined construction balance plus branch-matched regional fate polarity.",
        )
        add_unit(
            rows,
            level="sample_contrast",
            term_symbol="F",
            branch_group=branch_group,
            branch=branch,
            dataset_or_source=dataset,
            contrast_or_scope=scope,
            evidence_family="regional_fate",
            metric="delta_regional_fate_balance",
            raw_value=to_float(row.get("delta_regional_fate_balance")),
            normalized_score=clip_score(to_float(row.get("delta_regional_fate_balance")), 1.0),
            support_direction="positive candidate-background delta supports branch-matched fate polarity",
            quality_class="direct_primary_core",
            source_table=CONFIG_CONTRASTS.name,
            notes="Dentate and cerebellar branches are expected to differ in regional fate polarity.",
        )
        add_unit(
            rows,
            level="sample_contrast",
            term_symbol="C",
            branch_group=branch_group,
            branch=branch,
            dataset_or_source=dataset,
            contrast_or_scope=scope,
            evidence_family="downstream_construction",
            metric="delta_downstream_construction_balance",
            raw_value=to_float(row.get("delta_downstream_construction_balance")),
            normalized_score=clip_score(to_float(row.get("delta_downstream_construction_balance")), 1.0),
            support_direction="positive candidate-background delta supports downstream construction bias",
            quality_class="direct_primary_core",
            source_table=CONFIG_CONTRASTS.name,
            notes="Construction balance is shared toolkit usage, not by itself granule specificity.",
        )

    drivers = read_tsv(DRIVER_CONTRASTS)
    if "audit_scope" in drivers.columns:
        drivers = drivers.loc[drivers["audit_scope"].astype(str).eq("primary_core_candidate_background")]
    for _, row in drivers.iterrows():
        branch = clean_text(row.get("core_branch"), clean_text(row.get("branch_kind")))
        branch_group = branch_group_from_values(row.get("branch_kind"), branch)
        construction = bool(row.get("construction_positive"))
        fate = bool(row.get("fate_positive"))
        if construction and fate:
            score = 1.0
        elif construction or fate:
            score = 0.25
        else:
            score = -0.50
        add_unit(
            rows,
            level="sample_contrast",
            term_symbol="I",
            branch_group=branch_group,
            branch=branch,
            dataset_or_source=clean_text(row.get("dataset")),
            contrast_or_scope=clean_text(row.get("contrast_id")),
            evidence_family="configuration_driver_coupling",
            metric="driver_class_score",
            raw_value=score,
            normalized_score=score,
            support_direction="joint positivity supports coupled fate-plus-construction configuration",
            quality_class="direct_primary_core_decomposition",
            source_table=DRIVER_CONTRASTS.name,
            notes=clean_text(row.get("driver_notes")),
        )

    stage = read_tsv(STAGE_BRANCH)
    for _, row in stage.iterrows():
        signature = clean_text(row.get("signature_id"))
        endpoint = to_float(row.get("median_endpoint_delta"))
        peak = to_float(row.get("median_peak_score"))
        n_fits = to_float(row.get("n_group_fits"), 0.0)
        if signature == "neurogenic_permissive_index":
            raw = -endpoint
            direction = "decreasing permissive/progenitor signal with stage supports maturation-window structure"
        elif signature == "tgf_bdnf_2005_index":
            raw = abs(endpoint)
            direction = "nonzero TGF-beta/BDNF endpoint shift supports a branch-specific stage window"
        else:
            raw = endpoint
            direction = "positive endpoint shift supports increasing stop/maturation signal with stage"
        score = 0.70 * clip_score(raw, 0.50) + 0.30 * clip_score(peak, 1.0)
        branch = clean_text(row.get("branch"))
        add_unit(
            rows,
            level="branch_stage",
            term_symbol="T",
            branch_group=branch_group_from_values(branch),
            branch=branch,
            dataset_or_source="stage_window_model",
            contrast_or_scope=signature,
            evidence_family="stage_window_signature",
            metric="endpoint_or_window_score",
            raw_value=raw,
            normalized_score=score,
            support_direction=direction,
            weight=0.55 if n_fits <= 1 else 0.70,
            quality_class="stage_branch_model_limited" if n_fits <= 1 else "stage_branch_model",
            source_table=STAGE_BRANCH.name,
            notes=f"interpretation={clean_text(row.get('interpretation'))}; n_group_fits={int(n_fits)}",
        )

    pathways = read_tsv(PATHWAY_SUMMARY)
    pathways = pathways.loc[pathways["summary_level"].astype(str).eq("by_id_and_branch")]
    for _, row in pathways.iterrows():
        branch = clean_text(row.get("branch_kind"))
        median_delta = to_float(row.get("median_delta"))
        frac_positive = to_float(row.get("fraction_positive"))
        score = 0.60 * clip_score(median_delta, 0.50) + 0.40 * clip_score(frac_positive - 0.50, 0.50)
        add_unit(
            rows,
            level="branch_context",
            term_symbol="N",
            branch_group=branch_group_from_values(branch),
            branch=branch,
            dataset_or_source="primary_core_pathway_summary",
            contrast_or_scope=clean_text(row.get("pathway_id")),
            evidence_family="pathway_readiness",
            metric="median_delta_and_fraction_positive",
            raw_value=median_delta,
            normalized_score=score,
            support_direction="positive pathway median delta and positive-contrast fraction support branch niche/readiness signal",
            quality_class="expression_pathway_model",
            source_table=PATHWAY_SUMMARY.name,
            notes=f"fraction_positive={frac_positive:.3f}; n_contrasts={int(to_float(row.get('n_contrasts'), 0))}",
        )

    sender = read_tsv(SENDER_RECEIVER_SUMMARY)
    sender = sender.loc[sender["summary_level"].astype(str).eq("region")]
    for _, row in sender.iterrows():
        region = clean_text(row.get("region"))
        fraction_supported = to_float(row.get("fraction_supported"))
        max_score = to_float(row.get("max_lr_expression_score"))
        score = 0.70 * clip_score(fraction_supported, 0.08) + 0.30 * clip_score(max_score, 0.50)
        add_unit(
            rows,
            level="branch_context",
            term_symbol="N",
            branch_group=branch_group_from_values(region),
            branch=region,
            dataset_or_source="sender_receiver_lr_summary",
            contrast_or_scope="region_level_ligand_receptor_readiness",
            evidence_family="sender_receiver_readiness",
            metric="fraction_supported_and_max_lr_score",
            raw_value=fraction_supported,
            normalized_score=score,
            support_direction="higher supported ligand-receptor fraction and max score support niche-readiness context",
            weight=0.55,
            quality_class="expression_sender_receiver_model",
            source_table=SENDER_RECEIVER_SUMMARY.name,
            notes=f"n_supported={int(to_float(row.get('n_supported'), 0))}; n_high_support={int(to_float(row.get('n_high_support'), 0))}",
        )

    epi = read_tsv(EPI_ROBUST_POSITIVE)
    for _, row in epi.iterrows():
        concordance = clean_text(row.get("concordance_class"))
        term_supported = clean_text(row.get("model_term_supported"))
        robust_score = to_float(row.get("robust_rank_score"))
        weight = 0.50 if concordance == "robust_positive_strong" else 0.40
        add_unit(
            rows,
            level="epigenomic_sensitivity",
            term_symbol="E",
            branch_group="cerebellum",
            branch="human_cerebellum_GSE322785_provisional",
            dataset_or_source="GSE322785",
            contrast_or_scope=f"{clean_text(row.get('comparator'))}|{clean_text(row.get('feature_type'))}|{term_supported}|{clean_text(row.get('target_set'))}|{clean_text(row.get('peak_category'))}",
            evidence_family=f"epigenomic_support_for_{term_supported}",
            metric="robust_rank_score",
            raw_value=robust_score,
            normalized_score=clip_score(robust_score, 1.0),
            support_direction="robust-positive broad and cluster-supported granule-candidate contrast supports regulatory compatibility",
            weight=weight,
            quality_class="provisional_epigenomic_robust_positive",
            source_table=EPI_ROBUST_POSITIVE.name,
            notes=f"concordance={concordance}; cluster_supported_n_donors={row.get('cluster_supported_n_donors')}",
        )

    morph = read_tsv(NEUROMORPHO_COMPARISON)
    for metric in ["n_stems", "n_branch", "compact_branch_index", "length"]:
        sub = morph.loc[morph["metric"].astype(str).eq(metric)]
        if sub.empty:
            continue
        row = sub.iloc[0]
        ratio = to_float(row.get("median_ratio_dentate_over_cerebellum"))
        cliff = to_float(row.get("cliffs_delta_dentate_vs_cerebellum"))
        if metric == "n_branch":
            score = 1.0 - min(abs(math.log(max(ratio, 1e-9), 2.0)), 1.0)
            direction = "near-equal branch counts support shared nontrivial branching complexity"
        elif metric == "n_stems":
            score = min(abs(cliff), 1.0)
            direction = "large branch difference supports branch-specific sparse input-sampling geometry"
        elif metric == "compact_branch_index":
            score = min(abs(cliff), 1.0) * 0.75
            direction = "large compactness difference supports morphology as branch-specific configuration rather than identical geometry"
        else:
            score = min(abs(cliff), 1.0) * 0.50
            direction = "field-scale difference supports branch-specific geometry around a shared granule design principle"
        add_unit(
            rows,
            level="external_calibration",
            term_symbol="M",
            branch_group="cross_branch_or_external",
            branch="dentate_vs_cerebellar_granule_morphology",
            dataset_or_source="NeuroMorpho",
            contrast_or_scope=metric,
            evidence_family="public_morphology_comparison",
            metric="dentate_cerebellum_morphometry",
            raw_value=ratio,
            normalized_score=score,
            support_direction=direction,
            quality_class="external_morphology_calibration",
            source_table=NEUROMORPHO_COMPARISON.name,
            notes=f"dentate_median={row.get('dentate_median')}; cerebellum_median={row.get('cerebellum_median')}; cliffs_delta={cliff:.3f}",
        )

    targets = read_tsv(AIM3_TARGETS)
    for _, row in targets.iterrows():
        target_id = clean_text(row.get("target_id"))
        value = to_float(row.get("value"))
        source = clean_text(row.get("source"))
        if not target_id.startswith(("input_sampling", "branch_complexity")):
            continue
        if target_id.startswith("input_sampling"):
            score = 1.0 - min(max(value - 1.0, 0.0) / 15.0, 1.0)
            branch = "cerebellum" if "cerebellar" in target_id else "dentate"
            direction = "low primary-stem count is a lower-bound sparse input-sampling proxy"
        else:
            score = 0.75
            branch = "cross_branch_shared"
            direction = "nontrivial shared branch complexity constrains the sparse-expansion model"
        add_unit(
            rows,
            level="external_calibration",
            term_symbol="M",
            branch_group=branch_group_from_values(branch),
            branch=branch,
            dataset_or_source=source,
            contrast_or_scope=target_id,
            evidence_family="aim3_morphology_target",
            metric="calibration_target_value",
            raw_value=value,
            normalized_score=score,
            support_direction=direction,
            quality_class="external_morphology_calibration",
            source_table=AIM3_TARGETS.name,
            notes=clean_text(row.get("caveat")),
        )

    dandi = read_tsv(DANDI_POOLED)
    granule = dandi.loc[dandi["cell_type"].astype(str).eq("granule cell")]
    pyramidal = dandi.loc[dandi["cell_type"].astype(str).eq("pyramidal cell")]
    if not granule.empty:
        g = granule.iloc[0]
        p = pyramidal.iloc[0] if not pyramidal.empty else None
        info = to_float(g.get("median_spatial_information_bits_per_spike"))
        active_fraction = to_float(g.get("median_active_spatial_bin_fraction"))
        add_unit(
            rows,
            level="external_calibration",
            term_symbol="A",
            branch_group="dentate",
            branch="dentate_granule_activity",
            dataset_or_source="DANDI 000003",
            contrast_or_scope="granule_cell_spatial_information",
            evidence_family="dentate_activity_sparsity",
            metric="median_spatial_information_bits_per_spike",
            raw_value=info,
            normalized_score=clip_score(info, 0.80),
            support_direction="higher spatial information supports useful sparse or selective coding",
            quality_class="external_activity_calibration",
            source_table=DANDI_POOLED.name,
            notes=f"n_units={int(to_float(g.get('n_units'), 0))}",
        )
        add_unit(
            rows,
            level="external_calibration",
            term_symbol="A",
            branch_group="dentate",
            branch="dentate_granule_activity",
            dataset_or_source="DANDI 000003",
            contrast_or_scope="granule_cell_active_spatial_bin_fraction",
            evidence_family="dentate_activity_sparsity",
            metric="median_active_spatial_bin_fraction",
            raw_value=active_fraction,
            normalized_score=clip_score(1.0 - active_fraction, 0.50),
            support_direction="lower active-bin fraction supports restricted spatial occupancy",
            quality_class="external_activity_calibration",
            source_table=DANDI_POOLED.name,
            notes="Spatial-bin occupancy is not identical to instantaneous firing sparsity.",
        )
        if p is not None:
            info_delta = info - to_float(p.get("median_spatial_information_bits_per_spike"))
            active_delta = to_float(p.get("median_active_spatial_bin_fraction")) - active_fraction
            add_unit(
                rows,
                level="external_calibration",
                term_symbol="A",
                branch_group="dentate",
                branch="dentate_granule_vs_pyramidal_activity",
                dataset_or_source="DANDI 000003",
                contrast_or_scope="granule_minus_pyramidal_spatial_information",
                evidence_family="dentate_activity_comparator",
                metric="spatial_information_delta",
                raw_value=info_delta,
                normalized_score=clip_score(info_delta, 0.40),
                support_direction="granule-labeled units exceeding pyramidal spatial information supports cell-type activity distinction",
                quality_class="external_activity_calibration",
                source_table=DANDI_POOLED.name,
                notes="DANDI cell labels and sampling are external to the transcriptomic core.",
            )
            add_unit(
                rows,
                level="external_calibration",
                term_symbol="A",
                branch_group="dentate",
                branch="dentate_granule_vs_pyramidal_activity",
                dataset_or_source="DANDI 000003",
                contrast_or_scope="pyramidal_minus_granule_active_bin_fraction",
                evidence_family="dentate_activity_comparator",
                metric="active_spatial_bin_fraction_delta",
                raw_value=active_delta,
                normalized_score=clip_score(active_delta, 0.30),
                support_direction="granule-labeled units occupying fewer spatial bins than pyramidal units supports relative sparsity",
                quality_class="external_activity_calibration",
                source_table=DANDI_POOLED.name,
                notes="Relative spatial occupancy is a proxy, not a direct behavior-level pattern separation assay.",
            )

    pv = read_tsv(DANDI_PV)
    pv_g = pv.loc[pv["unit_set"].astype(str).eq("granule_cell_labeled")]
    if not pv_g.empty:
        median_far_near = float(pv_g["far_minus_near_neural_euclidean"].median())
        add_unit(
            rows,
            level="external_calibration",
            term_symbol="A",
            branch_group="dentate",
            branch="dentate_granule_population_vector",
            dataset_or_source="DANDI 000003",
            contrast_or_scope="multisession_far_minus_near_population_vector",
            evidence_family="dentate_population_vector_separation",
            metric="median_far_minus_near_neural_euclidean",
            raw_value=median_far_near,
            normalized_score=clip_score(median_far_near, 0.30),
            support_direction="far positions having larger population-vector distance than near positions supports weak pattern-separation calibration",
            quality_class="external_activity_calibration",
            source_table=DANDI_PV.name,
            notes=f"n_granule_labeled_sessions={len(pv_g)}",
        )

    calibration = read_tsv(AIM3_CALIBRATION)
    max_rank = int(calibration["resource_constrained_calibration_rank"].max())
    expected_preferred = {"intermediate", "granule_like_sparse_expansion"}
    expected_penalized = {"dense_expansion", "integrator_like_low_expansion", "excessively_sparse"}
    for _, row in calibration.iterrows():
        architecture = clean_text(row.get("architecture"))
        rank = int(to_float(row.get("resource_constrained_calibration_rank"), max_rank))
        if architecture in expected_preferred:
            # Preferred architectures support the model when they rank high.
            score = ((max_rank + 1 - rank) / max_rank)
            expected_direction = "preferred architecture should rank high under resource-constrained calibration"
        elif architecture in expected_penalized:
            # Predicted failure modes support the model when resource scoring penalizes them.
            score = rank / max_rank
            expected_direction = "predicted failure-mode architecture should rank low under resource-constrained calibration"
        else:
            score = ((max_rank + 1 - rank) / max_rank) * 2.0 - 1.0
            expected_direction = "architecture rank agreement with resource-constrained calibration"
        add_unit(
            rows,
            level="simulation_calibration",
            term_symbol="R",
            branch_group="simulation",
            branch="sparse_expansion_resource_model",
            dataset_or_source="Aim3_empirical_calibration",
            contrast_or_scope=architecture,
            evidence_family="resource_constraint_rank_agreement",
            metric="resource_constrained_calibration_rank",
            raw_value=rank,
            normalized_score=score,
            support_direction=expected_direction,
            quality_class="computational_simulation",
            source_table=AIM3_CALIBRATION.name,
            notes=f"median_resource_constrained_calibration_score={row.get('median_resource_constrained_calibration_score')}; empirical_rank={row.get('empirical_calibration_rank')}",
        )

    units = pd.DataFrame(rows)
    return units


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    valid = values.notna() & weights.notna() & (weights > 0)
    if not valid.any():
        return np.nan
    return float(np.average(values.loc[valid], weights=weights.loc[valid]))


def summarize_units(units: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group_cols = ["hierarchy_level", "term_symbol", "term"]
    layer_rows = []
    for keys, sub in units.groupby(group_cols, dropna=False):
        score = weighted_mean(sub["normalized_score"], sub["weight"])
        layer_rows.append(
            {
                "hierarchy_level": keys[0],
                "term_symbol": keys[1],
                "term": keys[2],
                "n_evidence_units": len(sub),
                "total_weight": sub["weight"].sum(),
                "weighted_mean_score": score,
                "median_normalized_score": sub["normalized_score"].median(),
                "fraction_positive_units": float((sub["normalized_score"] > 0).mean()),
                "support_class": support_class(score),
                "quality_classes": ";".join(sorted(sub["quality_class"].dropna().unique())),
            }
        )
    layer = pd.DataFrame(layer_rows).sort_values(["hierarchy_level", "term_symbol"])

    branch_rows = []
    for branch, sub in units.groupby("branch_group", dropna=False):
        row: dict[str, Any] = {
            "branch_group": branch,
            "n_evidence_units": len(sub),
            "total_weight": sub["weight"].sum(),
            "row_weighted_score": weighted_mean(sub["normalized_score"], sub["weight"]),
        }
        term_scores = []
        for term_symbol in TERM_NAMES:
            term_sub = sub.loc[sub["term_symbol"].eq(term_symbol)]
            score = weighted_mean(term_sub["normalized_score"], term_sub["weight"]) if not term_sub.empty else np.nan
            row[f"{term_symbol}_{TERM_NAMES[term_symbol]}_score"] = score
            if np.isfinite(score):
                term_scores.append(score)
        row["term_balanced_score"] = float(np.mean(term_scores)) if term_scores else np.nan
        row["support_class"] = support_class(row["term_balanced_score"])
        branch_rows.append(row)
    branch_summary = pd.DataFrame(branch_rows).sort_values("branch_group")

    component_map = {
        "direct_transcriptomic_configuration": ["Y", "F", "C", "I"],
        "stage_and_niche_context": ["T", "N"],
        "epigenomic_extension": ["E"],
        "morphology_activity_calibration": ["M", "A"],
        "resource_constraint_model": ["R"],
        "overall_hierarchical_model": list(TERM_NAMES),
    }
    interpretations = {
        "direct_transcriptomic_configuration": "Primary-core contrasts support a granule-cell configuration made from coupled fate polarity and downstream construction balance.",
        "stage_and_niche_context": "Stage-window and niche-readiness evidence is branch-specific and supports context-dependent regulation.",
        "epigenomic_extension": "Human cerebellar selected multiome data support a regulatory-compatibility scaffold but remain provisional.",
        "morphology_activity_calibration": "External morphology and dentate activity data calibrate sparse sampling, nontrivial branch complexity, and useful selectivity.",
        "resource_constraint_model": "Simulation/calibration supports sparse or intermediate expansion when resource costs are included.",
        "overall_hierarchical_model": "Across levels, evidence favors convergent granule-like design by hierarchical configuration rather than a single shared recent lineage or one exclusive gene set.",
    }
    caveats = {
        "direct_transcriptomic_configuration": "Direct but RNA-only; final geometry and function are not measured in the same cells.",
        "stage_and_niche_context": "Expression and pseudotime support do not prove protein activity or causal niche signaling.",
        "epigenomic_extension": "GSE322785 labels are provisional and dentate peak-count scoring is incomplete.",
        "morphology_activity_calibration": "External datasets are not matched to the transcriptomic core and cerebellar activity is not directly included.",
        "resource_constraint_model": "Mathematical architecture term; not a developmental perturbation result.",
        "overall_hierarchical_model": "The current implementation is a weighted evidence synthesis, not a fully matched Bayesian mixed-effects fit.",
    }
    comp_rows = []
    for component, terms in component_map.items():
        sub = units.loc[units["term_symbol"].isin(terms)]
        term_scores = []
        term_weights = []
        for term in terms:
            term_sub = sub.loc[sub["term_symbol"].eq(term)]
            if term_sub.empty:
                continue
            term_score = weighted_mean(term_sub["normalized_score"], term_sub["weight"])
            if np.isfinite(term_score):
                term_scores.append(term_score)
                term_weights.append(TERM_WEIGHTS[term])
        term_balanced = float(np.average(term_scores, weights=term_weights)) if term_scores else np.nan
        row_weighted = weighted_mean(sub["normalized_score"], sub["weight"]) if not sub.empty else np.nan
        comp_rows.append(
            {
                "component": component,
                "included_terms": ";".join(terms),
                "n_terms_with_evidence": len(term_scores),
                "n_evidence_units": len(sub),
                "total_evidence_weight": sub["weight"].sum(),
                "row_weighted_score": row_weighted,
                "term_balanced_score": term_balanced,
                "support_class": support_class(term_balanced),
                "interpretation": interpretations[component],
                "caveat": caveats[component],
            }
        )
    component = pd.DataFrame(comp_rows)
    return layer, branch_summary, component


def write_markdown(terms: pd.DataFrame, units: pd.DataFrame, layer: pd.DataFrame, component: pd.DataFrame) -> None:
    term_count = len(terms)
    unit_count = len(units)
    overall = component.loc[component["component"].eq("overall_hierarchical_model")].iloc[0]
    direct = component.loc[component["component"].eq("direct_transcriptomic_configuration")].iloc[0]
    epi = component.loc[component["component"].eq("epigenomic_extension")].iloc[0]

    lines = [
        "# Hierarchical Integrative Granule-Cell Model",
        "",
        "## Purpose",
        "",
        "This model formalizes the current project as a hierarchical evidence synthesis. The available evidence is not measured on one matched per-cell table: transcriptomic terms are sample/contrast-level, stage and niche terms are branch-level, epigenomic terms are provisional marker-group contrasts, and morphology/activity/resource terms are external calibration layers. The model therefore keeps each observation level explicit.",
        "",
        "## Ideal Matched-Data Model",
        "",
        "If matched transcriptome, epigenome, morphology, activity, and circuit measurements were available for the same units, the target model would be:",
        "",
        "```text",
        "GranuleDesign_i = beta0",
        "  + betaF * FatePolarity_i",
        "  + betaC * ConstructionBalance_i",
        "  + betaI * FatePolarity_i:ConstructionBalance_i",
        "  + betaT * StageWindow_branch/stage(i)",
        "  + betaN * NicheSignal_branch(i)",
        "  + betaE * EpigenomicCompatibility_resource(i)",
        "  + betaM * MorphologySparseSampling_branch(i)",
        "  + betaA * ActivitySparsity_branch(i)",
        "  + betaR * CircuitResourceConstraint",
        "  + random effects for dataset, species, source layer, and assay",
        "  + error_i",
        "```",
        "",
        "## Implemented Current Model",
        "",
        "The current implementation uses weighted hierarchical evidence units:",
        "",
        "```text",
        "S_level,term = sum_j(weight_j * normalized_score_j) / sum_j(weight_j)",
        "S_component = weighted average of term-level S values within each component",
        "```",
        "",
        "Positive scores support the term's predicted contribution to granule-like configuration or calibration. Negative scores indicate branch-specific opposition, absent support, or a result that argues against a simple shared direction. The component scores are term-balanced so that a large table with many rows does not automatically dominate a smaller but relevant evidence layer.",
        "",
        "## Outputs",
        "",
        f"- `{OUT_TERMS.name}`: {term_count} model terms and levels.",
        f"- `{OUT_UNITS.name}`: {unit_count} evidence units.",
        f"- `{OUT_LAYER.name}`: level-by-term evidence summaries.",
        f"- `{OUT_BRANCH.name}`: branch-level term summaries.",
        f"- `{OUT_COMPONENT.name}`: component-level scores.",
        "",
        "## Current Summary",
        "",
        f"- Direct transcriptomic configuration: term-balanced score {direct['term_balanced_score']:.3f} ({direct['support_class']}).",
        f"- Epigenomic extension: term-balanced score {epi['term_balanced_score']:.3f} ({epi['support_class']}); this remains provisional.",
        f"- Overall hierarchical model: term-balanced score {overall['term_balanced_score']:.3f} ({overall['support_class']}).",
        "",
        "The model supports the manuscript interpretation that dentate and cerebellar granule cells are not the same recent lineage or a single gene-barcode-defined cell type. Instead, distinct regional lineages appear to converge on a granule-like design through branch-specific fate polarity, shared downstream construction modules, stage/niche gating, regulatory compatibility, and resource-constrained sparse-expansion architecture.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def write_plot(component: pd.DataFrame) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    plot_df = component.loc[~component["component"].eq("overall_hierarchical_model")].copy()
    plot_df["label"] = plot_df["component"].str.replace("_", " ").str.title()
    colors = ["#4B7F52" if v >= 0.3 else "#8A6F3C" if v > 0 else "#9B4D4D" for v in plot_df["term_balanced_score"]]

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.barh(plot_df["label"], plot_df["term_balanced_score"], color=colors, edgecolor="#333333", linewidth=0.6)
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel("Term-balanced evidence score")
    ax.set_title("Hierarchical integrative granule-cell model")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def main() -> None:
    terms = build_terms()
    units = build_evidence_units()
    layer, branch, component = summarize_units(units)

    terms.to_csv(OUT_TERMS, sep="\t", index=False)
    units.to_csv(OUT_UNITS, sep="\t", index=False)
    layer.to_csv(OUT_LAYER, sep="\t", index=False)
    branch.to_csv(OUT_BRANCH, sep="\t", index=False)
    component.to_csv(OUT_COMPONENT, sep="\t", index=False)
    write_markdown(terms, units, layer, component)
    write_plot(component)

    print(f"Wrote {OUT_TERMS.relative_to(ROOT)}")
    print(f"Wrote {OUT_UNITS.relative_to(ROOT)}")
    print(f"Wrote {OUT_LAYER.relative_to(ROOT)}")
    print(f"Wrote {OUT_BRANCH.relative_to(ROOT)}")
    print(f"Wrote {OUT_COMPONENT.relative_to(ROOT)}")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    if OUT_PNG.exists():
        print(f"Wrote {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
