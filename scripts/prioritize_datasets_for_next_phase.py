#!/usr/bin/env python3
"""Prioritize datasets for the next phase of the granule-cell project."""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
ALLEN = ROOT / "External_Data" / "Allen_Institute"

UNITS_IN = Path(os.environ.get("DATASET_LEVEL_UNITS_INPUT", RESULTS / "refined_dataset_level_granule_program_units.tsv"))
AUDIT_IN = Path(os.environ.get("MATRIX_AUDIT_INPUT", RESULTS / "matrix_dimension_audit.tsv"))
ALLEN_CLUSTER_IN = Path(os.environ.get("ALLEN_CLUSTER_INPUT", ALLEN / "WMB_granule_relevant_cluster_annotation.tsv"))
ALLEN_SUBCLASS_IN = Path(os.environ.get("ALLEN_SUBCLASS_INPUT", ALLEN / "WMB_granule_relevant_subclass_annotation.tsv"))

OUTPUT_PREFIX = os.environ.get("PRIORITY_PREFIX", "next_phase")
SAMPLE_PRIORITY_OUT = RESULTS / f"{OUTPUT_PREFIX}_dataset_sample_priority.tsv"
DATASET_PRIORITY_OUT = RESULTS / f"{OUTPUT_PREFIX}_dataset_priority_summary.tsv"
ALLEN_REFERENCE_OUT = RESULTS / f"{OUTPUT_PREFIX}_allen_reference_summary.tsv"
FIG_OUT = RESULTS / f"{OUTPUT_PREFIX}_dataset_priority_map.png"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


CORE_DENTATE = {"GSE104323", "GSE95752", "GSE214309", "GSE292261"}
SUPPORT_DENTATE = {"GSE214905"}
CORE_CEREBELLAR = {"GSE122357", "GSE165657", "GSE312658"}
SUPPORT_CEREBELLAR = {"GSE242688"}
CONTEXT_DATASETS = {"GSE150153"}
PROMOTED_PRIMARY_VALIDATION = {"GSE214309", "GSE292261"}

CLASS_COLS = [
    "dentate_candidate",
    "cerebellar_candidate",
    "known_non_dentate_reference",
    "cerebellum_warning",
    "dentate_low_support",
    "organoid_granule_like",
    "other_or_ambiguous",
]


def sample_key(value: object) -> str:
    text = str(value)
    return re.sub(r"^GSM\d+_", "", text)


def first_non_null(values: pd.Series) -> str:
    vals = [str(v) for v in values.dropna().tolist() if str(v) and str(v) != "nan"]
    return vals[0] if vals else ""


def joined_unique(values: pd.Series) -> str:
    vals = sorted({str(v) for v in values.dropna().tolist() if str(v) and str(v) != "nan"})
    return ";".join(vals)


def median_for_class(units: pd.DataFrame, cls: str, metric: str) -> pd.DataFrame:
    col = f"{metric}_median"
    if col not in units.columns:
        return pd.DataFrame()
    sub = units.loc[units["analysis_class"] == cls].copy()
    if sub.empty:
        return pd.DataFrame()
    out = (
        sub.groupby(["dataset", "sample_key"], dropna=False)[col]
        .median()
        .reset_index()
        .rename(columns={col: f"{cls}_{metric}_unit_median"})
    )
    return out


def summarize_units(units: pd.DataFrame) -> pd.DataFrame:
    units = units.copy()
    units["sample_key"] = units["sample"].map(sample_key)
    units["n_cells_or_spots"] = pd.to_numeric(units["n_cells_or_spots"], errors="coerce").fillna(0)
    for col in ["n_high_confidence", "n_medium_confidence", "n_reference_confidence"]:
        units[col] = pd.to_numeric(units[col], errors="coerce").fillna(0)

    index_cols = ["dataset", "sample_key"]
    base = (
        units.groupby(index_cols, dropna=False)
        .agg(
            sample=("sample", first_non_null),
            region=("region", joined_unique),
            platform=("platform", joined_unique),
            species=("species", joined_unique),
            total_cells_or_spots=("n_cells_or_spots", "sum"),
            high_confidence_calls=("n_high_confidence", "sum"),
            medium_confidence_calls=("n_medium_confidence", "sum"),
            reference_calls=("n_reference_confidence", "sum"),
            n_analysis_units=("analysis_class", "size"),
        )
        .reset_index()
    )

    count_pivot = (
        units.pivot_table(
            index=index_cols,
            columns="analysis_class",
            values="n_cells_or_spots",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
    )
    count_pivot.columns.name = None
    for col in CLASS_COLS:
        if col not in count_pivot.columns:
            count_pivot[col] = 0
    count_pivot = count_pivot.rename(columns={col: f"{col}_cells_or_spots" for col in CLASS_COLS})

    out = base.merge(count_pivot, on=index_cols, how="left")
    for cls in ["dentate_candidate", "cerebellar_candidate", "cerebellum_warning", "organoid_granule_like"]:
        for metric in ["identity_contrast", "dentate_identity", "cerebellar_identity", "structural_rank"]:
            med = median_for_class(units, cls, metric)
            if not med.empty:
                out = out.merge(med, on=index_cols, how="left")

    out["candidate_cells_or_spots"] = out["dentate_candidate_cells_or_spots"] + out["cerebellar_candidate_cells_or_spots"]
    out["candidate_fraction"] = np.where(
        out["total_cells_or_spots"] > 0,
        out["candidate_cells_or_spots"] / out["total_cells_or_spots"],
        np.nan,
    )
    return out


def summarize_audit(audit: pd.DataFrame) -> pd.DataFrame:
    audit = audit.copy()
    audit["sample_key"] = audit["sample"].map(sample_key)
    for col in ["features_or_genes", "observations", "nonzero_entries"]:
        audit[col] = pd.to_numeric(audit[col], errors="coerce")
    return (
        audit.groupby(["dataset", "sample_key"], dropna=False)
        .agg(
            matrix_sample=("sample", first_non_null),
            matrix_format=("format", joined_unique),
            matrix_source_path=("source_path", joined_unique),
            features_or_genes=("features_or_genes", "max"),
            observations=("observations", "sum"),
            nonzero_entries=("nonzero_entries", "sum"),
            matrix_notes=("notes", joined_unique),
        )
        .reset_index()
    )


def infer_role(row: pd.Series) -> str:
    dataset = row["dataset"]
    region = row.get("region", "")
    has_dentate = row.get("dentate_candidate_cells_or_spots", 0) > 0
    has_cerebellar = row.get("cerebellar_candidate_cells_or_spots", 0) > 0

    if dataset in CORE_DENTATE and has_dentate:
        if dataset == "GSE104323":
            return "core dentate annotated reference and non-dentate controls"
        if dataset == "GSE95752":
            return "core dentate single-cell maturation validation"
        if dataset == "GSE214309":
            return "core dentate adult/activity-state validation"
        if dataset == "GSE292261":
            return "core dentate postnatal developmental trajectory"
    if dataset in SUPPORT_DENTATE and has_dentate:
        return "supporting dentate patch-seq physiology-linked validation"
    if dataset in CORE_CEREBELLAR and has_cerebellar:
        if dataset == "GSE122357":
            return "core cerebellar postnatal developmental comparison"
        if dataset == "GSE165657":
            return "core cerebellar large atlas-style validation"
        if dataset == "GSE312658":
            return "core cerebellar perturbation validation"
    if dataset in SUPPORT_CEREBELLAR and has_cerebellar:
        return "supporting cerebellar spatial/proteomics-linked validation"
    if dataset in CONTEXT_DATASETS or "organoid" in str(region).lower():
        return "context organoid granule-like comparison"
    return "supporting background or control"


def priority_score(row: pd.Series) -> float:
    score = 0.0
    dentate = row.get("dentate_candidate_cells_or_spots", 0)
    cerebellar = row.get("cerebellar_candidate_cells_or_spots", 0)
    candidate = dentate + cerebellar
    total = row.get("total_cells_or_spots", 0)
    observations = row.get("observations", np.nan)
    role = row.get("recommended_role", "")
    platform = str(row.get("platform", "")).lower() + ";" + str(row.get("matrix_format", "")).lower()
    region = str(row.get("region", "")).lower()

    if dentate > 0:
        score += 3.0
    if cerebellar > 0:
        score += 3.0
    if candidate >= 1000:
        score += 2.0
    elif candidate >= 100:
        score += 1.0
    if total >= 1000:
        score += 1.0
    if row.get("high_confidence_calls", 0) > 0 or row.get("reference_calls", 0) > 0:
        score += 1.0
    if row.get("other_or_ambiguous_cells_or_spots", 0) > 0 or row.get("known_non_dentate_reference_cells_or_spots", 0) > 0:
        score += 1.0
    if "10x" in platform or "single_cell" in platform:
        score += 1.0
    if "smart" in platform or "c1" in platform:
        score += 0.5
    if "core" in role:
        score += 1.0
    if "supporting" in role:
        score += 0.5
    if "spatial" in platform:
        score -= 1.0
    if "organoid" in region or "context organoid" in role:
        score -= 1.0
    if pd.notna(observations) and observations < 100:
        score -= 1.0
    return max(score, 0.0)


def priority_tier(row: pd.Series) -> str:
    role = row.get("recommended_role", "")
    score = row.get("priority_score", 0)
    if row.get("dataset", "") in PROMOTED_PRIMARY_VALIDATION and "core dentate" in role:
        return "core"
    if "context organoid" in role:
        return "context"
    if "spatial" in role or "patch-seq" in role:
        return "supporting"
    if score >= 8:
        return "core"
    if score >= 6:
        return "priority"
    if score >= 4:
        return "supporting"
    return "context"


def main_recommendation(row: pd.Series) -> str:
    tier = row["priority_tier"]
    role = row["recommended_role"]
    if tier == "core":
        return f"Use for primary object-level analysis: {role}."
    if tier == "priority":
        return f"Use in primary or secondary object-level analysis after annotation cleanup: {role}."
    if tier == "supporting":
        return f"Use as validation/supporting evidence, not as the main discovery set: {role}."
    return f"Use for biological context or review synthesis: {role}."


def caveats(row: pd.Series) -> str:
    notes = []
    dataset = row["dataset"]
    if row.get("cerebellar_candidate_cells_or_spots", 0) > 0 and row.get("cerebellum_warning_cells_or_spots", 0) > 0:
        notes.append("contains cerebellum warning calls; keep strict identity and shared structural modules separate")
    if dataset == "GSE242688":
        notes.append("spatial/proteomics-linked spots are not single cells")
    if dataset == "GSE150153":
        notes.append("organoid metadata were previously flagged as unmatched")
    if dataset == "GSE214905":
        notes.append("small patch-seq sample; best for targeted validation")
    if row.get("observations", np.nan) < 100:
        notes.append("small observation count")
    if not notes:
        notes.append("suitable for next-phase object-level curation")
    return "; ".join(notes)


def add_priority_fields(samples: pd.DataFrame) -> pd.DataFrame:
    samples = samples.copy()
    samples["recommended_role"] = samples.apply(infer_role, axis=1)
    samples["priority_score"] = samples.apply(priority_score, axis=1)
    samples["priority_tier"] = samples.apply(priority_tier, axis=1)
    samples["main_recommendation"] = samples.apply(main_recommendation, axis=1)
    samples["caveats"] = samples.apply(caveats, axis=1)
    tier_order = {"core": 0, "priority": 1, "supporting": 2, "context": 3}
    samples["_tier_rank"] = samples["priority_tier"].map(tier_order).fillna(9)
    samples = samples.sort_values(["_tier_rank", "priority_score", "dataset", "sample"], ascending=[True, False, True, True])
    return samples.drop(columns=["_tier_rank"])


def summarize_datasets(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, sub in samples.groupby("dataset", dropna=False):
        roles = joined_unique(sub["recommended_role"])
        tiers = joined_unique(sub["priority_tier"])
        rows.append(
            {
                "dataset": dataset,
                "n_samples": len(sub),
                "region": joined_unique(sub["region"]),
                "platform": joined_unique(sub["platform"]),
                "species": joined_unique(sub["species"]),
                "total_cells_or_spots": int(sub["total_cells_or_spots"].sum()),
                "observations_from_audit": int(sub["observations"].sum()) if sub["observations"].notna().any() else "",
                "dentate_candidate_cells_or_spots": int(sub["dentate_candidate_cells_or_spots"].sum()),
                "cerebellar_candidate_cells_or_spots": int(sub["cerebellar_candidate_cells_or_spots"].sum()),
                "cerebellum_warning_cells_or_spots": int(sub["cerebellum_warning_cells_or_spots"].sum()),
                "best_priority_score": sub["priority_score"].max(),
                "priority_tiers": tiers,
                "recommended_roles": roles,
                "main_use": main_recommendation(sub.sort_values("priority_score", ascending=False).iloc[0]),
                "caveats": joined_unique(sub["caveats"]),
            }
        )
    out = pd.DataFrame(rows)
    tier_order = {"core": 0, "priority": 1, "supporting": 2, "context": 3}
    out["_tier_rank"] = out["priority_tiers"].map(
        lambda value: min(tier_order.get(tier, 9) for tier in str(value).split(";") if tier)
    )
    return out.sort_values(["_tier_rank", "best_priority_score", "dataset"], ascending=[True, False, True]).drop(
        columns=["_tier_rank"]
    )


def summarize_allen_reference() -> pd.DataFrame:
    if not ALLEN_CLUSTER_IN.exists():
        return pd.DataFrame()
    clusters = pd.read_csv(ALLEN_CLUSTER_IN, sep="\t")
    subclasses = pd.read_csv(ALLEN_SUBCLASS_IN, sep="\t") if ALLEN_SUBCLASS_IN.exists() else pd.DataFrame()
    for col in ["v3.size", "v2.size", "multiome.size"]:
        if col in clusters.columns:
            clusters[col] = pd.to_numeric(clusters[col], errors="coerce").fillna(0).astype(int)
    rows = []
    for subclass, sub in clusters.groupby("subclass_label", dropna=False):
        sub_info = subclasses.loc[subclasses.get("subclass_label", pd.Series(dtype=str)) == subclass]
        markers = first_non_null(sub_info["subclass.markers.combo"]) if not sub_info.empty and "subclass.markers.combo" in sub_info else ""
        tf_markers = first_non_null(sub_info["subclass.tf.markers.combo"]) if not sub_info.empty and "subclass.tf.markers.combo" in sub_info else ""
        rows.append(
            {
                "source": "Allen_Institute_WMB_taxonomy",
                "subclass_label": subclass,
                "n_clusters": len(sub),
                "n_supertypes": sub["supertype_label"].nunique() if "supertype_label" in sub else "",
                "total_v3_size": int(sub["v3.size"].sum()) if "v3.size" in sub else "",
                "total_v2_size": int(sub["v2.size"].sum()) if "v2.size" in sub else "",
                "markers": markers,
                "tf_markers": tf_markers,
                "recommended_use": "reference atlas for cell-label validation, marker refinement, and figure annotation",
            }
        )
    return pd.DataFrame(rows).sort_values(["subclass_label"])


def plot_priority(samples: pd.DataFrame) -> None:
    order = samples.sort_values(["priority_score", "candidate_cells_or_spots"], ascending=[True, True]).copy()
    labels = order["dataset"] + ":" + order["sample"].astype(str)
    colors = {
        "core": "#2a9d8f",
        "priority": "#457b9d",
        "supporting": "#d08c60",
        "context": "#8d99ae",
    }
    fig_h = max(4.5, 0.34 * len(order) + 1.5)
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    ax.barh(labels, order["priority_score"], color=[colors.get(t, "#8d99ae") for t in order["priority_tier"]])
    ax.set_xlabel("Next-phase priority score")
    ax.set_ylabel("")
    ax.set_title("Granule-cell project dataset/sample prioritization")
    ax.grid(axis="x", linewidth=0.4, color="#d9d9d9", alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=color, label=tier) for tier, color in colors.items()]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=180)


def main() -> None:
    units = pd.read_csv(UNITS_IN, sep="\t")
    audit = pd.read_csv(AUDIT_IN, sep="\t")
    samples = summarize_units(units)
    audit_summary = summarize_audit(audit)
    samples = samples.merge(audit_summary, on=["dataset", "sample_key"], how="left")
    samples = add_priority_fields(samples)
    datasets = summarize_datasets(samples)
    allen = summarize_allen_reference()

    samples.to_csv(SAMPLE_PRIORITY_OUT, sep="\t", index=False, float_format="%.6g")
    datasets.to_csv(DATASET_PRIORITY_OUT, sep="\t", index=False)
    if not allen.empty:
        allen.to_csv(ALLEN_REFERENCE_OUT, sep="\t", index=False)
    plot_priority(samples)

    print(f"Wrote {SAMPLE_PRIORITY_OUT}")
    print(f"Wrote {DATASET_PRIORITY_OUT}")
    if not allen.empty:
        print(f"Wrote {ALLEN_REFERENCE_OUT}")
    print(f"Wrote {FIG_OUT}")


if __name__ == "__main__":
    main()
