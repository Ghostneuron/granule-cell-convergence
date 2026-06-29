#!/usr/bin/env python3
"""Integrate human dentate bridge units with the existing dentate/cerebellar backbone.

The older backbone and the newer human bridge objects were scored through
different count matrices, so this script uses within-sample rank metrics for the
main cross-branch comparison. Raw module medians are retained for inspection but
not used as the primary combined statistic.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

BACKBONE_UNITS = RESULTS / "refined_dataset_level_granule_program_units.tsv"
HUMAN_CORE_LABELS = RESULTS / "human_core_tuned_labels.tsv.gz"
GSE325391_LABELS = RESULTS / "gse325391_human_core_label_projection.tsv.gz"
GSE268609_LABELS = RESULTS / "gse268609_human_core_label_projection.tsv.gz"

OUT_UNITS = RESULTS / "human_bridge_backbone_rank_units.tsv"
OUT_STATS = RESULTS / "human_bridge_backbone_rank_statistics.tsv"
OUT_SOURCE_SUMMARY = RESULTS / "human_bridge_backbone_rank_source_summary.tsv"
OUT_PLOT = RESULTS / "human_bridge_backbone_rank_units.png"
OUT_MD = RESULTS / "human_bridge_backbone_rank_integration.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


RAW_MODULES = [
    "dentate_identity",
    "cerebellar_identity",
    "shared_granule_neuronal",
    "morphogenesis_cytoskeleton",
    "axon_guidance_synapse",
    "structural_program_mean",
    "background_max",
]

RANK_MODULES = [
    "dentate_rank",
    "cerebellar_rank",
    "shared_rank",
    "structural_rank",
    "identity_rank_contrast",
]

MIN_UNIT_CELLS = 20


def percentile_rank(series: pd.Series) -> pd.Series:
    if len(series) <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def broad_class(analysis_class: str) -> str:
    if analysis_class in {
        "backbone_dentate_candidate",
        "human_dentate_anchor",
        "human_dentate_candidate",
        "human_immature_neurogenic_candidate",
        "human_projected_dentate_candidate",
        "human_projected_immature_neurogenic_candidate",
    }:
        return "dentate_candidate"
    if analysis_class == "backbone_cerebellar_candidate":
        return "cerebellar_candidate"
    if analysis_class in {"backbone_non_dentate_reference", "human_non_neuronal_background"}:
        return "non_dentate_background"
    if analysis_class in {"backbone_cerebellum_warning", "human_broad_neuronal_structural_warning"}:
        return "broad_neuronal_structural_warning"
    if analysis_class in {"backbone_dentate_low_support", "human_dentate_low_support"}:
        return "dentate_low_support"
    return "other_or_ambiguous"


def classify_human_core(label: str) -> str:
    if label == "curated_human_dg_gc_anchor":
        return "human_dentate_anchor"
    if label == "human_dg_like_high_confidence":
        return "human_dentate_candidate"
    if label == "human_dg_like_candidate":
        return "human_dentate_low_support"
    if label == "immature_neurogenic_candidate":
        return "human_immature_neurogenic_candidate"
    if label == "immature_neurogenic_candidate_low_support":
        return "human_dentate_low_support"
    if label == "non_neuronal_background":
        return "human_non_neuronal_background"
    if label == "broad_neuronal_structural_warning":
        return "human_broad_neuronal_structural_warning"
    return "human_other_or_ambiguous"


def classify_gse325391(label: str) -> str:
    if label in {"adult_human_dg_mature_anchor", "adult_human_dg_differentiating_anchor"}:
        return "human_dentate_anchor"
    if label == "adult_dg_background_warning":
        return "human_non_neuronal_background"
    if label == "adult_dg_doublet_flag":
        return "human_other_or_ambiguous"
    return "human_other_or_ambiguous"


def classify_gse268609(label: str) -> str:
    if label == "human_dg_like_high_confidence":
        return "human_projected_dentate_candidate"
    if label == "human_dg_like_candidate":
        return "human_projected_dentate_candidate"
    if label == "immature_neurogenic_candidate":
        return "human_projected_immature_neurogenic_candidate"
    if label == "non_neuronal_background":
        return "human_non_neuronal_background"
    if label == "broad_neuronal_structural_warning":
        return "human_broad_neuronal_structural_warning"
    return "human_other_or_ambiguous"


def load_backbone_units() -> pd.DataFrame:
    df = pd.read_csv(BACKBONE_UNITS, sep="\t")
    class_map = {
        "dentate_candidate": "backbone_dentate_candidate",
        "cerebellar_candidate": "backbone_cerebellar_candidate",
        "known_non_dentate_reference": "backbone_non_dentate_reference",
        "cerebellum_warning": "backbone_cerebellum_warning",
        "dentate_low_support": "backbone_dentate_low_support",
        "organoid_granule_like": "backbone_organoid_granule_like",
        "other_or_ambiguous": "backbone_other_or_ambiguous",
    }
    out = pd.DataFrame(
        {
            "source_layer": "backbone_refined",
            "dataset": df["dataset"],
            "species": df["species"],
            "region": df["region"],
            "platform": df["platform"],
            "sample": df["sample"],
            "sample_context": df["region"],
            "label": df["group"],
            "analysis_class": df["analysis_class"].map(class_map).fillna("backbone_other_or_ambiguous"),
            "n_cells_or_spots": df["n_cells_or_spots"],
        }
    )
    for metric in RAW_MODULES:
        col = f"{metric}_median"
        out[f"{metric}_median"] = pd.to_numeric(df[col], errors="coerce") if col in df else np.nan
    out["dentate_rank_median"] = pd.to_numeric(df["dentate_rank_median"], errors="coerce")
    out["cerebellar_rank_median"] = pd.to_numeric(df["cerebellar_rank_median"], errors="coerce")
    out["shared_rank_median"] = pd.to_numeric(df["shared_rank_median"], errors="coerce")
    out["structural_rank_median"] = pd.to_numeric(df["structural_rank_median"], errors="coerce")
    out["identity_rank_contrast_median"] = out["dentate_rank_median"] - out["cerebellar_rank_median"]
    out["identity_rank_contrast_mean"] = out["identity_rank_contrast_median"]
    out["broad_class"] = out["analysis_class"].map(broad_class)
    out["unit_id"] = (
        out["source_layer"]
        + "|"
        + out["dataset"].astype(str)
        + "|"
        + out["sample"].astype(str)
        + "|"
        + out["label"].astype(str)
        + "|"
        + out["analysis_class"].astype(str)
    )
    return out


def standardize_human_core() -> pd.DataFrame:
    df = pd.read_csv(HUMAN_CORE_LABELS, sep="\t", low_memory=False)
    out = pd.DataFrame(
        {
            "source_layer": "human_core_tuned",
            "dataset": df["dataset"],
            "species": "human",
            "region": "hippocampus_dentate_gyrus",
            "platform": "snRNA_seq_selected_bridge",
            "sample": df["replicate_unit"].astype(str),
            "sample_context": df["geo_age"].fillna("").astype(str),
            "cell_id": df["cell_id"],
            "label": df["tuned_label"],
            "analysis_class": df["tuned_label"].map(classify_human_core),
        }
    )
    for metric in RAW_MODULES:
        source = f"norm_{metric}" if metric != "background_max" else "norm_background_max"
        if source in df:
            out[metric] = pd.to_numeric(df[source], errors="coerce")
    return out


def standardize_gse325391() -> pd.DataFrame:
    df = pd.read_csv(GSE325391_LABELS, sep="\t", low_memory=False)
    include = df["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    df = df.loc[include].copy()
    out = pd.DataFrame(
        {
            "source_layer": "gse325391_adult_dg",
            "dataset": "GSE325391",
            "species": "human",
            "region": "dentate_gyrus",
            "platform": "single_nucleus_selected_bridge",
            "sample": df["sample"].astype(str),
            "sample_context": df["group"].astype(str),
            "cell_id": df["cell_id"],
            "label": df["source_anchor_label"],
            "analysis_class": df["source_anchor_label"].map(classify_gse325391),
        }
    )
    for metric in RAW_MODULES:
        source = f"norm_{metric}" if metric != "background_max" else "norm_background_max"
        if source in df:
            out[metric] = pd.to_numeric(df[source], errors="coerce")
    return out


def standardize_gse268609() -> pd.DataFrame:
    df = pd.read_csv(GSE268609_LABELS, sep="\t", low_memory=False, dtype={"sample_id": str})
    include = df["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    df = df.loc[include].copy()
    out = pd.DataFrame(
        {
            "source_layer": "gse268609_hippocampus_rna",
            "dataset": "GSE268609",
            "species": "human",
            "region": "hippocampus_dentate_gyrus",
            "platform": "snRNA_multiome_RNA_selected_bridge",
            "sample": df["sample_id"].astype(str),
            "sample_context": df["diagnosis"].astype(str),
            "cell_id": df["cell_id"],
            "label": df["projected_label"],
            "analysis_class": df["projected_label"].map(classify_gse268609),
        }
    )
    for metric in RAW_MODULES:
        source = f"norm_{metric}" if metric != "background_max" else "norm_background_max"
        if source in df:
            out[metric] = pd.to_numeric(df[source], errors="coerce")
    return out


def summarize_human_units(cells: pd.DataFrame) -> pd.DataFrame:
    cells = cells.copy()
    for metric in ["dentate_identity", "cerebellar_identity", "shared_granule_neuronal", "structural_program_mean"]:
        cells[f"{metric}_rank"] = cells.groupby(["source_layer", "dataset", "sample"], dropna=False)[metric].transform(percentile_rank)
    cells["identity_rank_contrast"] = cells["dentate_identity_rank"] - cells["cerebellar_identity_rank"]
    cells["broad_class"] = cells["analysis_class"].map(broad_class)

    group_cols = [
        "source_layer",
        "dataset",
        "species",
        "region",
        "platform",
        "sample",
        "sample_context",
        "label",
        "analysis_class",
        "broad_class",
    ]
    rows = []
    numeric_cols = RAW_MODULES + [
        "dentate_identity_rank",
        "cerebellar_identity_rank",
        "shared_granule_neuronal_rank",
        "structural_program_mean_rank",
        "identity_rank_contrast",
    ]
    rename = {
        "dentate_identity_rank": "dentate_rank",
        "cerebellar_identity_rank": "cerebellar_rank",
        "shared_granule_neuronal_rank": "shared_rank",
        "structural_program_mean_rank": "structural_rank",
    }
    for key, sub in cells.groupby(group_cols, dropna=False):
        row = dict(zip(group_cols, key))
        row["n_cells_or_spots"] = len(sub)
        for col in numeric_cols:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()
            out_name = rename.get(col, col)
            row[f"{out_name}_mean"] = vals.mean() if len(vals) else np.nan
            row[f"{out_name}_median"] = vals.median() if len(vals) else np.nan
            row[f"{out_name}_q25"] = vals.quantile(0.25) if len(vals) else np.nan
            row[f"{out_name}_q75"] = vals.quantile(0.75) if len(vals) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    out["unit_id"] = (
        out["source_layer"]
        + "|"
        + out["dataset"].astype(str)
        + "|"
        + out["sample"].astype(str)
        + "|"
        + out["label"].astype(str)
        + "|"
        + out["analysis_class"].astype(str)
    )
    return out


def load_all_units() -> pd.DataFrame:
    backbone = load_backbone_units()
    human_cells = pd.concat(
        [standardize_human_core(), standardize_gse325391(), standardize_gse268609()],
        ignore_index=True,
        sort=False,
    )
    human_units = summarize_human_units(human_cells)
    all_cols = sorted(set(backbone.columns) | set(human_units.columns))
    return pd.concat([backbone.reindex(columns=all_cols), human_units.reindex(columns=all_cols)], ignore_index=True, sort=False)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    idx = p_values.index[valid].to_numpy()
    if len(idx) == 0:
        return out
    order = p_values.loc[idx].to_numpy().argsort()
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    m = len(p_sorted)
    adjusted = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def compare(units: pd.DataFrame, class_a: str, class_b: str, metric: str, note: str) -> dict[str, object]:
    col = f"{metric}_median"
    a = units.loc[units["broad_class"].eq(class_a), col].dropna().to_numpy(dtype=float)
    b = units.loc[units["broad_class"].eq(class_b), col].dropna().to_numpy(dtype=float)
    if len(a) == 0 or len(b) == 0:
        p_value = np.nan
    else:
        p_value = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return {
        "comparison": f"{class_a}_vs_{class_b}",
        "metric": metric,
        "n_units_a": len(a),
        "n_units_b": len(b),
        "median_a": float(np.median(a)) if len(a) else np.nan,
        "median_b": float(np.median(b)) if len(b) else np.nan,
        "delta_a_minus_b": float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan,
        "p_value": p_value,
        "note": note,
    }


def sign_test(units: pd.DataFrame, broad: str, metric: str, direction: str, threshold: float, note: str) -> dict[str, object]:
    vals = units.loc[units["broad_class"].eq(broad), f"{metric}_median"].dropna().to_numpy(dtype=float)
    if direction == "greater":
        successes = int(np.sum(vals > threshold))
    elif direction == "less":
        successes = int(np.sum(vals < threshold))
    else:
        raise ValueError(direction)
    p_value = stats.binomtest(successes, len(vals), 0.5, alternative="greater").pvalue if len(vals) else np.nan
    return {
        "comparison": f"{broad}_sign_test_{direction}_{threshold}",
        "metric": metric,
        "n_units_a": len(vals),
        "n_units_b": "",
        "median_a": float(np.median(vals)) if len(vals) else np.nan,
        "median_b": threshold,
        "delta_a_minus_b": float(np.median(vals) - threshold) if len(vals) else np.nan,
        "p_value": p_value,
        "note": f"{note}; successes={successes}/{len(vals)}",
    }


def origin_compare(units: pd.DataFrame, origin_a: str, origin_b: str, broad: str, metric: str, note: str) -> dict[str, object]:
    col = f"{metric}_median"
    a = units.loc[units["source_layer"].eq(origin_a) & units["broad_class"].eq(broad), col].dropna().to_numpy(dtype=float)
    b = units.loc[units["source_layer"].eq(origin_b) & units["broad_class"].eq(broad), col].dropna().to_numpy(dtype=float)
    p_value = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) and len(b) else np.nan
    return {
        "comparison": f"{origin_a}_{broad}_vs_{origin_b}_{broad}",
        "metric": metric,
        "n_units_a": len(a),
        "n_units_b": len(b),
        "median_a": float(np.median(a)) if len(a) else np.nan,
        "median_b": float(np.median(b)) if len(b) else np.nan,
        "delta_a_minus_b": float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan,
        "p_value": p_value,
        "note": note,
    }


def compute_stats(units: pd.DataFrame) -> pd.DataFrame:
    eligible = units.loc[units["n_cells_or_spots"] >= MIN_UNIT_CELLS].copy()
    rows = [
        compare(
            eligible,
            "dentate_candidate",
            "cerebellar_candidate",
            "identity_rank_contrast",
            "Rank-based identity contrast across the integrated backbone plus human bridge.",
        ),
        compare(
            eligible,
            "dentate_candidate",
            "non_dentate_background",
            "identity_rank_contrast",
            "Dentate units should separate from non-dentate/background units.",
        ),
        compare(
            eligible,
            "dentate_candidate",
            "cerebellar_candidate",
            "structural_rank",
            "Structural rank is expected to show convergence rather than identity separation.",
        ),
        sign_test(
            eligible,
            "dentate_candidate",
            "identity_rank_contrast",
            "greater",
            0.0,
            "Dentate candidate units should rank higher for dentate than cerebellar identity within their sample.",
        ),
        sign_test(
            eligible,
            "cerebellar_candidate",
            "identity_rank_contrast",
            "less",
            0.0,
            "Cerebellar candidate units should rank higher for cerebellar than dentate identity within their sample.",
        ),
        sign_test(
            eligible,
            "dentate_candidate",
            "structural_rank",
            "greater",
            0.5,
            "Dentate candidate units above within-sample structural median support structural-program activity.",
        ),
        sign_test(
            eligible,
            "cerebellar_candidate",
            "structural_rank",
            "greater",
            0.5,
            "Cerebellar candidate units above within-sample structural median support structural-program activity.",
        ),
        origin_compare(
            eligible,
            "gse325391_adult_dg",
            "gse268609_hippocampus_rna",
            "dentate_candidate",
            "identity_rank_contrast",
            "Source-anchored adult DG units versus projected broader hippocampal RNA units.",
        ),
        origin_compare(
            eligible,
            "human_core_tuned",
            "gse268609_hippocampus_rna",
            "dentate_candidate",
            "structural_rank",
            "Tuned human-core dentate units versus projected GSE268609 dentate units.",
        ),
    ]
    out = pd.DataFrame(rows)
    out["p_adj_bh"] = bh_adjust(out["p_value"])
    out["min_unit_cells"] = MIN_UNIT_CELLS
    return out


def compute_source_summary(units: pd.DataFrame) -> pd.DataFrame:
    eligible = units.loc[units["n_cells_or_spots"] >= MIN_UNIT_CELLS].copy()
    summary = (
        eligible.groupby(["source_layer", "broad_class"], dropna=False)
        .agg(
            n_units=("unit_id", "size"),
            n_cells_or_spots=("n_cells_or_spots", "sum"),
            median_identity_rank_contrast=("identity_rank_contrast_median", "median"),
            median_structural_rank=("structural_rank_median", "median"),
            median_dentate_rank=("dentate_rank_median", "median"),
            median_cerebellar_rank=("cerebellar_rank_median", "median"),
        )
        .reset_index()
        .sort_values(["source_layer", "broad_class"])
    )
    return summary


def plot_units(units: pd.DataFrame) -> None:
    plot_df = units.loc[units["n_cells_or_spots"] >= MIN_UNIT_CELLS].copy()
    plot_df = plot_df.loc[
        plot_df["broad_class"].isin(
            [
                "dentate_candidate",
                "cerebellar_candidate",
                "non_dentate_background",
                "broad_neuronal_structural_warning",
                "dentate_low_support",
                "other_or_ambiguous",
            ]
        )
    ]
    palette = {
        "dentate_candidate": "#168b7a",
        "cerebellar_candidate": "#6d3bbd",
        "non_dentate_background": "#7f8790",
        "broad_neuronal_structural_warning": "#c97938",
        "dentate_low_support": "#60a99f",
        "other_or_ambiguous": "#c2c2c2",
    }
    markers = {
        "backbone_refined": "o",
        "human_core_tuned": "s",
        "gse325391_adult_dg": "D",
        "gse268609_hippocampus_rna": "^",
    }
    fig, ax = plt.subplots(figsize=(10.6, 6.2))
    for (broad, layer), sub in plot_df.groupby(["broad_class", "source_layer"], dropna=False):
        ax.scatter(
            sub["identity_rank_contrast_median"],
            sub["structural_rank_median"],
            s=18 + np.clip(np.sqrt(sub["n_cells_or_spots"].astype(float)), 2, 85),
            c=palette.get(broad, "#999999"),
            marker=markers.get(layer, "o"),
            alpha=0.78,
            edgecolors="white",
            linewidths=0.5,
            label=f"{broad} | {layer}",
        )
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.axhline(0.5, color="#777777", linewidth=0.9, linestyle=":")
    ax.set_xlabel("Identity rank contrast (dentate - cerebellar)")
    ax.set_ylabel("Structural-program rank")
    ax.set_title("Integrated granule-program rank units")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.grid(True, linewidth=0.4, color="#d9d9d9", alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 0.74, 1])
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_md(units: pd.DataFrame, stats_df: pd.DataFrame, source_summary: pd.DataFrame) -> None:
    eligible = units.loc[units["n_cells_or_spots"] >= MIN_UNIT_CELLS].copy()
    counts = (
        eligible.groupby(["source_layer", "broad_class"], dropna=False)
        .agg(n_units=("unit_id", "size"), n_cells_or_spots=("n_cells_or_spots", "sum"))
        .reset_index()
        .sort_values(["source_layer", "broad_class"])
    )
    lines = [
        "# Human Bridge And Backbone Rank Integration",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This analysis integrates the existing refined mouse/human cerebellar-dentate backbone with the constructed human dentate/hippocampal bridge objects using within-sample rank metrics.",
        "",
        "Raw module medians are retained, but the main integrated comparison uses rank metrics because the older backbone and newer human selected-gene bridge objects were not produced from one shared whole-transcriptome matrix.",
        "",
        f"Minimum unit size for statistics and plotting: {MIN_UNIT_CELLS} cells/spots.",
        "",
        "## Eligible Unit Counts",
        "",
    ]
    for _, row in counts.iterrows():
        lines.append(
            f"- `{row['source_layer']}` / `{row['broad_class']}`: "
            f"{int(row['n_units'])} units, {int(row['n_cells_or_spots'])} cells/spots."
        )
    lines.extend(["", "## Main Tests", ""])
    for _, row in stats_df.sort_values(["p_adj_bh", "comparison"]).iterrows():
        lines.append(
            f"- `{row['comparison']}` / `{row['metric']}`: delta {row['delta_a_minus_b']:.4f} "
            f"(n={row['n_units_a']} vs {row['n_units_b']}, BH-adjusted p={row['p_adj_bh']:.3g})."
        )
    lines.extend(["", "## Source-Layer Median Signals", ""])
    signal_rows = source_summary.loc[
        source_summary["broad_class"].isin(["dentate_candidate", "cerebellar_candidate", "non_dentate_background"])
    ].copy()
    for _, row in signal_rows.iterrows():
        lines.append(
            f"- `{row['source_layer']}` / `{row['broad_class']}`: "
            f"identity contrast {row['median_identity_rank_contrast']:.4f}, "
            f"structural rank {row['median_structural_rank']:.4f} "
            f"({int(row['n_units'])} units; {int(row['n_cells_or_spots'])} cells/spots)."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The integrated rank layer is a sanity-check bridge, not yet a replacement for a single harmonized object.",
            "- A robust project signal would be: dentate candidate units show positive dentate-minus-cerebellar identity rank, cerebellar candidate units show negative identity rank contrast, and both groups show above-median structural-program rank.",
            "- `GSE325391` and `GSE186538` provide the strongest source-aware human dentate anchors; their enriched DG composition makes them better anchors for structural/neurogenic state than for within-sample dentate-versus-cerebellar identity contrast.",
            "- `GSE268609` broadens the human aging/AD context and shows useful projected dentate signal, but it remains projection-labeled until source taxonomy from the full Seurat object is added.",
            "",
            "## Outputs",
            "",
            f"- Integrated units: `{OUT_UNITS.relative_to(ROOT)}`",
            f"- Integrated statistics: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Source-layer summary: `{OUT_SOURCE_SUMMARY.relative_to(ROOT)}`",
            f"- Rank-unit plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    units = load_all_units()
    units.to_csv(OUT_UNITS, sep="\t", index=False, float_format="%.6g")
    stats_df = compute_stats(units)
    stats_df.to_csv(OUT_STATS, sep="\t", index=False, float_format="%.6g")
    source_summary = compute_source_summary(units)
    source_summary.to_csv(OUT_SOURCE_SUMMARY, sep="\t", index=False, float_format="%.6g")
    plot_units(units)
    write_md(units, stats_df, source_summary)

    print(f"Wrote {OUT_UNITS}")
    print(f"Wrote {OUT_STATS}")
    print(f"Wrote {OUT_SOURCE_SUMMARY}")
    print(f"Wrote {OUT_PLOT}")
    print(f"Wrote {OUT_MD}")
    print(f"units={len(units)}; eligible={(units['n_cells_or_spots'] >= MIN_UNIT_CELLS).sum()}")


if __name__ == "__main__":
    main()
