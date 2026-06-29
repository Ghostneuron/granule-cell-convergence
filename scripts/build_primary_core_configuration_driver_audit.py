#!/usr/bin/env python3
"""Audit what drives the transcriptomic configuration score.

The primary-core configuration validation is strong overall, but the combined
score has two components:

1. downstream construction balance
2. regional fate polarity

This audit decomposes every contrast to determine whether the signal is driven
by construction, fate, or both. It is meant to keep the manuscript claim honest:
configuration support is not automatically morphology-only support.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

PRIMARY_UNITS = RESULTS / "primary_core_transcriptomic_configuration_primary_units.tsv.gz"
PRIMARY_CONTRASTS = RESULTS / "primary_core_transcriptomic_configuration_primary_contrasts.tsv"
LOCAL_CONTRASTS = RESULTS / "primary_core_transcriptomic_configuration_contrasts.tsv"
FORMAL_GENE_SCORES = RESULTS / "primary_core_niche_circuit_module_formal_gene_scores.tsv"

OUT_DRIVER_CONTRASTS = RESULTS / "primary_core_configuration_driver_audit_contrasts.tsv"
OUT_MODULE_DELTAS = RESULTS / "primary_core_configuration_driver_audit_module_deltas.tsv"
OUT_SUMMARY = RESULTS / "primary_core_configuration_driver_audit_summary.tsv"
OUT_GENE_PRIORITIES = RESULTS / "primary_core_configuration_driver_audit_gene_priorities.tsv"
OUT_PLOT = RESULTS / "primary_core_configuration_driver_audit.png"
OUT_MD = RESULTS / "primary_core_configuration_driver_audit.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE_SCORE_COLS = [
    "cerebellar_fate_rank",
    "dentate_fate_rank",
    "niche_rank",
    "neurite_morphology_rank",
    "synaptic_excitability_rank",
    "downstream_mean_rank",
    "branch_matched_fate_rank",
    "branch_opposed_fate_rank",
]

MODULE_LABELS = {
    "cerebellar_fate_rank": "Cerebellar fate",
    "dentate_fate_rank": "Dentate fate",
    "niche_rank": "Niche/progenitor",
    "neurite_morphology_rank": "Neurite/morphology",
    "synaptic_excitability_rank": "Synaptic/excitability",
    "downstream_mean_rank": "Downstream mean",
    "branch_matched_fate_rank": "Matched fate",
    "branch_opposed_fate_rank": "Opposed fate",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def median_or_nan(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def driver_class(row: pd.Series) -> str:
    construction = float(row["delta_downstream_construction_balance"])
    fate = float(row["delta_regional_fate_balance"])
    config = float(row["delta_configuration_score"])
    if not np.isfinite(config) or config <= 0:
        return "configuration_not_positive"
    construction_pos = np.isfinite(construction) and construction > 0
    fate_pos = np.isfinite(fate) and fate > 0
    if construction_pos and fate_pos:
        return "both_components_positive"
    if construction_pos and not fate_pos:
        return "construction_driven_positive"
    if fate_pos and not construction_pos:
        return "fate_driven_positive"
    return "weak_mixed_positive"


def read_and_classify_contrasts() -> pd.DataFrame:
    primary = pd.read_csv(PRIMARY_CONTRASTS, sep="\t")
    primary["audit_scope"] = "primary_core_candidate_background"
    primary["contrast_id"] = (
        primary["audit_scope"]
        + "|"
        + primary["expression_layer"].astype(str)
        + "|"
        + primary["dataset"].astype(str)
        + "|"
        + primary["sample"].astype(str)
        + "|"
        + primary["target_class"].astype(str)
    )

    local = pd.read_csv(LOCAL_CONTRASTS, sep="\t")
    local = local.rename(
        columns={
            "granule_role": "target_class",
            "comparator_role": "background_classes",
            "granule_median_downstream_construction_balance": "candidate_median_downstream_construction_balance",
            "comparator_median_downstream_construction_balance": "background_median_downstream_construction_balance",
            "granule_median_regional_fate_balance": "candidate_median_regional_fate_balance",
            "comparator_median_regional_fate_balance": "background_median_regional_fate_balance",
            "granule_median_configuration_score": "candidate_median_configuration_score",
            "comparator_median_configuration_score": "background_median_configuration_score",
        }
    )
    local["expression_layer"] = "local_named_comparator_modules"
    local["branch_kind"] = local["branch"]
    local["core_branch"] = np.where(local["branch"].eq("dentate"), "mouse_dentate", "cerebellum")
    local["source_layer"] = "source_group_module_scores"
    local["expression_scope"] = "local_named_comparator"
    local["audit_scope"] = "local_named_comparator"
    local["configuration_positive"] = local["delta_configuration_score"].gt(0)
    local["contrast_id"] = (
        local["audit_scope"]
        + "|"
        + local["dataset"].astype(str)
        + "|"
        + local["sample"].astype(str)
        + "|"
        + local["target_class"].astype(str)
        + "_vs_"
        + local["background_classes"].astype(str)
    )
    keep_cols = [
        "audit_scope",
        "contrast_id",
        "expression_layer",
        "dataset",
        "core_branch",
        "branch_kind",
        "sample",
        "source_layer",
        "expression_scope",
        "target_class",
        "background_classes",
        "candidate_median_downstream_construction_balance",
        "background_median_downstream_construction_balance",
        "delta_downstream_construction_balance",
        "candidate_median_regional_fate_balance",
        "background_median_regional_fate_balance",
        "delta_regional_fate_balance",
        "candidate_median_configuration_score",
        "background_median_configuration_score",
        "delta_configuration_score",
        "configuration_positive",
    ]
    contrasts = pd.concat([primary[keep_cols], local[keep_cols]], ignore_index=True)
    contrasts["driver_class"] = contrasts.apply(driver_class, axis=1)
    contrasts["construction_positive"] = contrasts["delta_downstream_construction_balance"].gt(0)
    contrasts["fate_positive"] = contrasts["delta_regional_fate_balance"].gt(0)
    contrasts["driver_notes"] = np.select(
        [
            contrasts["driver_class"].eq("both_components_positive"),
            contrasts["driver_class"].eq("fate_driven_positive"),
            contrasts["driver_class"].eq("construction_driven_positive"),
            contrasts["driver_class"].eq("configuration_not_positive"),
        ],
        [
            "construction and regional fate both support the configuration",
            "configuration is positive mainly because regional fate polarity is positive",
            "configuration is positive mainly because construction balance is positive",
            "combined configuration is not positive",
        ],
        default="weak or mixed positive configuration",
    )
    contrasts.to_csv(OUT_DRIVER_CONTRASTS, sep="\t", index=False)
    return contrasts


def build_primary_module_deltas() -> pd.DataFrame:
    units = pd.read_csv(PRIMARY_UNITS, sep="\t")
    contrasts = pd.read_csv(PRIMARY_CONTRASTS, sep="\t")
    rows: list[dict[str, object]] = []
    key_cols = ["expression_layer", "dataset", "core_branch", "sample", "source_layer", "expression_scope"]
    for _, contrast in contrasts.iterrows():
        sub = units.copy()
        for col in key_cols:
            sub = sub.loc[sub[col].astype(str).eq(str(contrast[col]))]
        if sub.empty:
            continue
        candidate = sub.loc[sub["broad_class"].eq(contrast["target_class"])]
        background_classes = str(contrast["background_classes"]).split(",")
        background = sub.loc[sub["broad_class"].astype(str).isin(background_classes)]
        if candidate.empty or background.empty:
            continue
        for score_col in MODULE_SCORE_COLS:
            candidate_median = median_or_nan(candidate[score_col])
            background_median = median_or_nan(background[score_col])
            rows.append(
                {
                    "contrast_id": (
                        "primary_core_candidate_background"
                        + "|"
                        + str(contrast["expression_layer"])
                        + "|"
                        + str(contrast["dataset"])
                        + "|"
                        + str(contrast["sample"])
                        + "|"
                        + str(contrast["target_class"])
                    ),
                    "audit_scope": "primary_core_candidate_background",
                    "expression_layer": contrast["expression_layer"],
                    "dataset": contrast["dataset"],
                    "core_branch": contrast["core_branch"],
                    "branch_kind": contrast["branch_kind"],
                    "sample": contrast["sample"],
                    "target_class": contrast["target_class"],
                    "background_classes": contrast["background_classes"],
                    "module_score": score_col,
                    "module_label": MODULE_LABELS.get(score_col, score_col),
                    "candidate_median": candidate_median,
                    "background_median": background_median,
                    "delta_candidate_minus_background": candidate_median - background_median
                    if np.isfinite(candidate_median) and np.isfinite(background_median)
                    else np.nan,
                }
            )
    module_deltas = pd.DataFrame(rows)
    module_deltas.to_csv(OUT_MODULE_DELTAS, sep="\t", index=False)
    return module_deltas


def binom_and_wilcoxon(values: pd.Series) -> tuple[int, int, float, float, float]:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    n = int(arr.size)
    if n == 0:
        return 0, 0, np.nan, np.nan, np.nan
    positives = int(np.sum(arr > 0))
    median = float(np.median(arr))
    sign_p = float(stats.binomtest(positives, n=n, p=0.5, alternative="greater").pvalue)
    wilcox = float(stats.wilcoxon(arr, alternative="greater").pvalue) if n > 1 and np.any(arr != 0) else np.nan
    return positives, n, median, sign_p, wilcox


def summarize_drivers(contrasts: pd.DataFrame, module_deltas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_specs = [
        ("all", []),
        ("by_audit_scope", ["audit_scope"]),
        ("by_primary_expression_layer", ["audit_scope", "expression_layer"]),
        ("by_primary_branch", ["audit_scope", "branch_kind"]),
        ("by_primary_layer_branch", ["audit_scope", "expression_layer", "branch_kind"]),
    ]
    for summary_level, cols in group_specs:
        if not cols:
            iterator = [((), contrasts)]
        else:
            iterator = contrasts.groupby(cols, sort=False)
        for keys, sub in iterator:
            if not isinstance(keys, tuple):
                keys = (keys,)
            row: dict[str, object] = {
                "summary_level": summary_level,
                "n_contrasts": int(len(sub)),
                "n_configuration_positive": int(sub["delta_configuration_score"].gt(0).sum()),
                "n_both_components_positive": int(sub["driver_class"].eq("both_components_positive").sum()),
                "n_fate_driven_positive": int(sub["driver_class"].eq("fate_driven_positive").sum()),
                "n_construction_driven_positive": int(sub["driver_class"].eq("construction_driven_positive").sum()),
                "n_configuration_not_positive": int(sub["driver_class"].eq("configuration_not_positive").sum()),
            }
            for col, value in zip(cols, keys):
                row[col] = value
            for score in [
                "delta_downstream_construction_balance",
                "delta_regional_fate_balance",
                "delta_configuration_score",
            ]:
                positives, n, median, sign_p, wilcox = binom_and_wilcoxon(sub[score])
                row[f"n_positive_{score}"] = positives
                row[f"n_tested_{score}"] = n
                row[f"median_{score}"] = median
                row[f"sign_p_{score}"] = sign_p
                row[f"wilcoxon_p_{score}"] = wilcox
            rows.append(row)

    module_rows: list[dict[str, object]] = []
    for module_score, sub in module_deltas.groupby("module_score", sort=False):
        positives, n, median, sign_p, wilcox = binom_and_wilcoxon(sub["delta_candidate_minus_background"])
        module_rows.append(
            {
                "summary_level": "primary_module_delta",
                "module_score": module_score,
                "module_label": MODULE_LABELS.get(module_score, module_score),
                "n_contrasts": n,
                "n_positive_delta": positives,
                "median_delta_candidate_minus_background": median,
                "sign_p_delta_positive": sign_p,
                "wilcoxon_p_delta_positive": wilcox,
            }
        )

    summary = pd.concat([pd.DataFrame(rows), pd.DataFrame(module_rows)], ignore_index=True, sort=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    return summary


def build_gene_priorities() -> pd.DataFrame:
    genes = pd.read_csv(FORMAL_GENE_SCORES, sep="\t")
    present = genes.loc[genes["present_in_formal_model"].astype(bool)].copy()
    present["module_family_priority"] = np.select(
        [
            present["module_family"].eq("downstream_circuit_morphology"),
            present["module_family"].eq("upstream_region_fate"),
            present["module_family"].eq("shared_niche_state"),
        ],
        [0, 1, 2],
        default=3,
    )
    present["driver_priority_score"] = (
        present["overall_convergence_delta"].fillna(-1)
        + present["shared_positive_both_screens"].astype(float)
        + 0.5 * present["shared_positive_any_screen"].astype(float)
        - 0.25 * present["mean_branch_bias"].fillna(0)
    )
    out = present.sort_values(
        ["module_family_priority", "driver_priority_score", "overall_convergence_delta"],
        ascending=[True, False, False],
    )
    keep = [
        "module_id",
        "module_label",
        "module_family",
        "gene",
        "canonical_gene",
        "formal_mouse_symbol",
        "overall_convergence_delta",
        "mean_branch_bias",
        "shared_positive_any_screen",
        "shared_positive_both_screens",
        "branch_pattern",
        "formal_rank_tier",
        "formal_rank_priority_score",
        "driver_priority_score",
    ]
    out[keep].to_csv(OUT_GENE_PRIORITIES, sep="\t", index=False)
    return out[keep]


def plot_audit(contrasts: pd.DataFrame, module_deltas: pd.DataFrame) -> None:
    primary = contrasts.loc[contrasts["audit_scope"].eq("primary_core_candidate_background")].copy()
    driver_counts = (
        contrasts.groupby(["audit_scope", "driver_class"], sort=False).size().reset_index(name="n")
    )
    module_plot = module_deltas.groupby(["module_score", "module_label"], sort=False)[
        "delta_candidate_minus_background"
    ].median().reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.2), constrained_layout=True)

    ax = axes[0]
    class_order = [
        "both_components_positive",
        "fate_driven_positive",
        "construction_driven_positive",
        "weak_mixed_positive",
        "configuration_not_positive",
    ]
    colors = {
        "both_components_positive": "#2f7f8f",
        "fate_driven_positive": "#7b6d8d",
        "construction_driven_positive": "#8a9a5b",
        "weak_mixed_positive": "#b48a4a",
        "configuration_not_positive": "#8a4f4f",
    }
    bottom_by_scope: dict[str, float] = {}
    scopes = list(driver_counts["audit_scope"].drop_duplicates())
    x = np.arange(len(scopes))
    bottoms = np.zeros(len(scopes))
    for cls in class_order:
        vals = []
        for scope in scopes:
            val = driver_counts.loc[
                driver_counts["audit_scope"].eq(scope) & driver_counts["driver_class"].eq(cls), "n"
            ]
            vals.append(int(val.iloc[0]) if len(val) else 0)
        ax.bar(x, vals, bottom=bottoms, label=cls.replace("_", " "), color=colors[cls])
        bottoms += np.asarray(vals)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", "\n") for s in scopes], fontsize=8)
    ax.set_ylabel("Contrasts")
    ax.set_title("Driver class counts")
    ax.legend(
        frameon=True,
        fontsize=7,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        facecolor="white",
        edgecolor="none",
        framealpha=0.92,
        borderpad=0.35,
        handlelength=1.8,
    )

    ax = axes[1]
    primary = primary.sort_values(["expression_layer", "branch_kind", "dataset", "sample"])
    y = np.arange(len(primary))
    ax.barh(y - 0.18, primary["delta_downstream_construction_balance"], height=0.34, color="#2f7f8f", label="construction")
    ax.barh(y + 0.18, primary["delta_regional_fate_balance"], height=0.34, color="#7b6d8d", label="fate polarity")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y[:: max(1, len(y) // 18)])
    ax.set_yticklabels((primary["dataset"] + " " + primary["sample"].astype(str)).iloc[:: max(1, len(y) // 18)], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Candidate-background delta")
    ax.set_title("Primary component deltas")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    ax = axes[2]
    module_plot = module_plot.loc[module_plot["module_score"].isin(MODULE_SCORE_COLS[:5])]
    ax.barh(module_plot["module_label"], module_plot["delta_candidate_minus_background"], color="#1f4e5f")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Median primary candidate-background delta")
    ax.set_title("Base module deltas")
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    fig.suptitle("Configuration Driver Audit", fontsize=15, y=1.03)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(contrasts: pd.DataFrame, summary: pd.DataFrame, gene_priorities: pd.DataFrame) -> None:
    all_row = summary.loc[summary["summary_level"].eq("all")].iloc[0]
    primary_row = summary.loc[
        summary["summary_level"].eq("by_audit_scope")
        & summary["audit_scope"].eq("primary_core_candidate_background")
    ].iloc[0]
    local_row = summary.loc[
        summary["summary_level"].eq("by_audit_scope") & summary["audit_scope"].eq("local_named_comparator")
    ].iloc[0]
    primary_modules = summary.loc[summary["summary_level"].eq("primary_module_delta")].copy()
    top_downstream = gene_priorities.loc[
        gene_priorities["module_family"].eq("downstream_circuit_morphology")
    ].head(12)

    lines = [
        "# Configuration Driver Audit",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This audit decomposes the transcriptomic configuration score into downstream construction balance and regional fate polarity. It asks whether the assembly-plan signal is morphology-weighted, identity-weighted, or both.",
        "",
        "## Driver Classes",
        "",
        f"- All contrasts: {int(all_row['n_configuration_positive'])}/{int(all_row['n_contrasts'])} configuration-positive.",
        f"- Both components positive: {int(all_row['n_both_components_positive'])}.",
        f"- Fate-driven positive: {int(all_row['n_fate_driven_positive'])}.",
        f"- Construction-driven positive: {int(all_row['n_construction_driven_positive'])}.",
        f"- Configuration not positive: {int(all_row['n_configuration_not_positive'])}.",
        "",
        "Primary-core candidate-background layer:",
        f"- Configuration-positive: {int(primary_row['n_configuration_positive'])}/{int(primary_row['n_contrasts'])}.",
        f"- Both components positive: {int(primary_row['n_both_components_positive'])}.",
        f"- Fate-driven positive: {int(primary_row['n_fate_driven_positive'])}.",
        f"- Construction-driven positive: {int(primary_row['n_construction_driven_positive'])}.",
        f"- Median construction delta: {primary_row['median_delta_downstream_construction_balance']:.3f}.",
        f"- Median fate-polarity delta: {primary_row['median_delta_regional_fate_balance']:.3f}.",
        f"- Median combined configuration delta: {primary_row['median_delta_configuration_score']:.3f}.",
        "",
        "Local named-comparator layer:",
        f"- Configuration-positive: {int(local_row['n_configuration_positive'])}/{int(local_row['n_contrasts'])}.",
        f"- Both components positive: {int(local_row['n_both_components_positive'])}.",
        f"- Fate-driven positive: {int(local_row['n_fate_driven_positive'])}.",
        f"- Construction-driven positive: {int(local_row['n_construction_driven_positive'])}.",
        "",
        "## Base Module Delta Summary",
        "",
    ]
    for _, row in primary_modules.iterrows():
        lines.append(
            f"- {row['module_label']}: median candidate-background delta "
            f"{row['median_delta_candidate_minus_background']:.3f}; "
            f"{int(row['n_positive_delta'])}/{int(row['n_contrasts'])} positive."
        )

    lines.extend(
        [
            "",
            "## Top Downstream Configuration Genes",
            "",
        ]
    )
    for _, row in top_downstream.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['module_label']}): convergence {row['overall_convergence_delta']:.3f}, "
            f"branch bias {row['mean_branch_bias']:.3f}, pattern {row['branch_pattern']}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The broad primary-core configuration signal is real, but it is strongly identity-coupled: regional fate polarity is positive in most primary contrasts, while downstream construction balance is more selective.",
            "- This supports the phrase 'transcriptomic assembly configuration' more than a pure morphology-only transcriptomic code.",
            "- The most defensible claim is that distinct regional fate programs place granule cells into a permissive context, while shared downstream neurite/synaptic machinery contributes the morphology implementation layer.",
            "- A stronger morphology-specific test should add explicit pyramidal/Purkinje comparator labels in more datasets, morphology-linked datasets, or spatial/proteomic localization.",
            "",
            "## Outputs",
            "",
            f"- Driver contrast table: `{rel(OUT_DRIVER_CONTRASTS)}`",
            f"- Module delta table: `{rel(OUT_MODULE_DELTAS)}`",
            f"- Summary table: `{rel(OUT_SUMMARY)}`",
            f"- Gene priority table: `{rel(OUT_GENE_PRIORITIES)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    contrasts = read_and_classify_contrasts()
    module_deltas = build_primary_module_deltas()
    summary = summarize_drivers(contrasts, module_deltas)
    gene_priorities = build_gene_priorities()
    plot_audit(contrasts, module_deltas)
    write_report(contrasts, summary, gene_priorities)
    print(f"Wrote {rel(OUT_MD)}")
    print(
        summary.loc[
            summary["summary_level"].isin(["all", "by_audit_scope"]),
            [
                "summary_level",
                "audit_scope",
                "n_contrasts",
                "n_configuration_positive",
                "n_both_components_positive",
                "n_fate_driven_positive",
                "n_construction_driven_positive",
                "n_configuration_not_positive",
                "median_delta_downstream_construction_balance",
                "median_delta_regional_fate_balance",
                "median_delta_configuration_score",
            ],
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
