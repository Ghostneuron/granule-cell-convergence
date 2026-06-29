#!/usr/bin/env python3
"""Test transcriptomic configuration rather than single-module uniqueness.

The previous named-comparator analyses showed that individual downstream
neurite/synaptic modules are not uniquely granule-specific. This script asks a
slightly different question: do granule cells show a distinctive configuration
of common neuronal construction modules?

Configuration here means two interpretable balances:

- downstream construction balance:
  mean(neurite/morphology, synaptic/excitability) - neurogenic niche/progenitor
- regional fate balance:
  branch-matched fate module - branch-opposed fate module

The combined score is the sum of those two balances. It is not a causal model;
it is a transcriptomic "assembly-plan" score that tests whether module ratios
separate granule cells from named pyramidal/Purkinje comparators better than
single pathway membership does.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

IN_UNITS = RESULTS / "primary_core_niche_circuit_module_named_comparator_units.tsv"

OUT_UNITS = RESULTS / "primary_core_transcriptomic_configuration_units.tsv"
OUT_ROLE_SUMMARY = RESULTS / "primary_core_transcriptomic_configuration_role_summary.tsv"
OUT_CONTRASTS = RESULTS / "primary_core_transcriptomic_configuration_contrasts.tsv"
OUT_PLOT = RESULTS / "primary_core_transcriptomic_configuration_model.png"
OUT_MD = RESULTS / "primary_core_transcriptomic_configuration_model.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE_COLS = {
    "cerebellar_fate_rhombic_lip_shh": "cerebellar_fate_rank",
    "dentate_fate_wnt_prox1": "dentate_fate_rank",
    "shared_neurogenic_niche_state": "niche_rank",
    "downstream_neurite_morphology": "neurite_morphology_rank",
    "downstream_synaptic_excitability": "synaptic_excitability_rank",
}

ROLE_ORDER = [
    "dentate_granule",
    "pyramidal_comparator",
    "cerebellar_granule",
    "purkinje_comparator",
    "other_local_cell_type",
]

PLOT_ROLE_ORDER = [
    "dentate_granule",
    "pyramidal_comparator",
    "cerebellar_granule",
    "purkinje_comparator",
]

ROLE_LABELS = {
    "dentate_granule": "Dentate granule",
    "pyramidal_comparator": "Pyramidal comparator",
    "cerebellar_granule": "Cerebellar granule",
    "purkinje_comparator": "Purkinje comparator",
    "other_local_cell_type": "Other local cells",
}

SCORE_LABELS = {
    "downstream_construction_balance": "Construction over niche",
    "regional_fate_balance": "Regional fate polarity",
    "configuration_score": "Combined configuration",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def finite_median(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
    if values.size == 0:
        return np.nan
    return float(np.median(values))


def branch_for_dataset(dataset: str) -> str:
    if dataset == "GSE104323":
        return "dentate"
    if dataset == "GSE122357":
        return "cerebellar"
    return "unknown"


def build_configuration_units() -> pd.DataFrame:
    units = pd.read_csv(IN_UNITS, sep="\t")
    pivot = units.pivot_table(
        index=["dataset", "sample", "source_group", "specificity_role"],
        columns="module_id",
        values="within_sample_module_rank",
        aggfunc="median",
    ).reset_index()
    pivot.columns.name = None
    pivot = pivot.rename(columns=MODULE_COLS)
    for col in MODULE_COLS.values():
        if col not in pivot.columns:
            pivot[col] = np.nan

    pivot["dataset_branch"] = pivot["dataset"].map(branch_for_dataset)
    pivot["downstream_mean_rank"] = pivot[["neurite_morphology_rank", "synaptic_excitability_rank"]].mean(axis=1)
    pivot["downstream_construction_balance"] = pivot["downstream_mean_rank"] - pivot["niche_rank"]

    pivot["branch_matched_fate_rank"] = np.where(
        pivot["dataset_branch"].eq("dentate"),
        pivot["dentate_fate_rank"],
        np.where(pivot["dataset_branch"].eq("cerebellar"), pivot["cerebellar_fate_rank"], np.nan),
    )
    pivot["branch_opposed_fate_rank"] = np.where(
        pivot["dataset_branch"].eq("dentate"),
        pivot["cerebellar_fate_rank"],
        np.where(pivot["dataset_branch"].eq("cerebellar"), pivot["dentate_fate_rank"], np.nan),
    )
    pivot["regional_fate_balance"] = pivot["branch_matched_fate_rank"] - pivot["branch_opposed_fate_rank"]
    pivot["configuration_score"] = pivot["downstream_construction_balance"] + pivot["regional_fate_balance"]
    pivot["granule_role"] = pivot["specificity_role"].isin(["dentate_granule", "cerebellar_granule"])
    pivot["named_comparator_role"] = pivot["specificity_role"].isin(["pyramidal_comparator", "purkinje_comparator"])
    pivot = pivot.sort_values(["dataset", "sample", "specificity_role", "source_group"])
    pivot.to_csv(OUT_UNITS, sep="\t", index=False)
    return pivot


def summarize_roles(config: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for role, sub in config.groupby("specificity_role", sort=False):
        row = {
            "specificity_role": role,
            "role_label": ROLE_LABELS.get(role, role),
            "role_order": ROLE_ORDER.index(role) if role in ROLE_ORDER else 99,
            "n_group_units": int(len(sub)),
            "n_datasets": int(sub["dataset"].nunique()),
            "n_samples": int(sub[["dataset", "sample"]].drop_duplicates().shape[0]),
        }
        for col in [
            "cerebellar_fate_rank",
            "dentate_fate_rank",
            "niche_rank",
            "neurite_morphology_rank",
            "synaptic_excitability_rank",
            "downstream_mean_rank",
            "downstream_construction_balance",
            "regional_fate_balance",
            "configuration_score",
        ]:
            row[f"median_{col}"] = finite_median(sub[col])
        rows.append(row)
    out = pd.DataFrame(rows).sort_values("role_order")
    out.to_csv(OUT_ROLE_SUMMARY, sep="\t", index=False)
    return out


def summarize_contrasts(config: pd.DataFrame) -> pd.DataFrame:
    contrast_specs = [
        {
            "dataset": "GSE104323",
            "branch": "dentate",
            "granule_role": "dentate_granule",
            "comparator_role": "pyramidal_comparator",
            "contrast_label": "Dentate granule vs pyramidal comparator",
        },
        {
            "dataset": "GSE122357",
            "branch": "cerebellar",
            "granule_role": "cerebellar_granule",
            "comparator_role": "purkinje_comparator",
            "contrast_label": "Cerebellar granule vs Purkinje comparator",
        },
    ]
    score_cols = [
        "downstream_construction_balance",
        "regional_fate_balance",
        "configuration_score",
    ]
    rows: list[dict[str, object]] = []
    for spec in contrast_specs:
        dsub = config.loc[config["dataset"].eq(spec["dataset"])]
        for sample, sub in dsub.groupby("sample", sort=False):
            granule = sub.loc[sub["specificity_role"].eq(spec["granule_role"])]
            comparator = sub.loc[sub["specificity_role"].eq(spec["comparator_role"])]
            if granule.empty or comparator.empty:
                continue
            row: dict[str, object] = {
                "dataset": spec["dataset"],
                "branch": spec["branch"],
                "sample": sample,
                "contrast_label": spec["contrast_label"],
                "granule_role": spec["granule_role"],
                "comparator_role": spec["comparator_role"],
                "n_granule_groups": int(len(granule)),
                "n_comparator_groups": int(len(comparator)),
            }
            for col in score_cols:
                granule_median = finite_median(granule[col])
                comparator_median = finite_median(comparator[col])
                row[f"granule_median_{col}"] = granule_median
                row[f"comparator_median_{col}"] = comparator_median
                row[f"delta_{col}"] = granule_median - comparator_median
            rows.append(row)
    out = pd.DataFrame(rows)
    if not out.empty:
        out["granule_configuration_positive"] = out["delta_configuration_score"].gt(0)
    out.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    return out


def pvalue_for_deltas(contrasts: pd.DataFrame, col: str) -> float:
    values = contrasts[col].dropna().to_numpy(dtype=float)
    if values.size == 0:
        return np.nan
    return float(stats.wilcoxon(values, alternative="greater").pvalue) if values.size > 1 else np.nan


def plot_results(config: pd.DataFrame, role_summary: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    plot_roles = [role for role in PLOT_ROLE_ORDER if role in set(role_summary["specificity_role"])]
    module_cols = [
        "median_dentate_fate_rank",
        "median_cerebellar_fate_rank",
        "median_niche_rank",
        "median_neurite_morphology_rank",
        "median_synaptic_excitability_rank",
    ]
    module_labels = [
        "Dentate fate",
        "Cerebellar fate",
        "Niche/progenitor",
        "Neurite/morphology",
        "Synaptic/excitability",
    ]
    heat = (
        role_summary.set_index("specificity_role")
        .loc[plot_roles, module_cols]
        .to_numpy(dtype=float)
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.8), constrained_layout=True)
    ax = axes[0]
    im = ax.imshow(heat, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(module_labels)))
    ax.set_xticklabels(module_labels, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(plot_roles)))
    ax.set_yticklabels([ROLE_LABELS.get(role, role) for role in plot_roles])
    ax.set_title("Module-rank configuration")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            value = heat[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value < 0.55 else "black", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Median within-sample rank")

    ax = axes[1]
    score_cols = [
        "delta_downstream_construction_balance",
        "delta_regional_fate_balance",
        "delta_configuration_score",
    ]
    score_labels = [
        "Construction over niche",
        "Regional fate polarity",
        "Combined configuration",
    ]
    plot_contrasts = contrasts.copy()
    plot_contrasts["short_label"] = np.where(
        plot_contrasts["branch"].eq("dentate"),
        "DG vs Pyr",
        "CB vs Purk " + plot_contrasts["sample"].astype(str),
    )
    x = np.arange(len(plot_contrasts))
    width = 0.25
    colors = ["#2f7f8f", "#7b6d8d", "#1f4e5f"]
    for idx, col in enumerate(score_cols):
        ax.bar(x + (idx - 1) * width, plot_contrasts[col], width=width, label=score_labels[idx], color=colors[idx])
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_contrasts["short_label"], rotation=30, ha="right")
    ax.set_ylabel("Granule minus named-comparator delta")
    ax.set_title("Configuration contrasts")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    fig.suptitle("Transcriptomic Configuration Model", fontsize=15, y=1.03)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(config: pd.DataFrame, role_summary: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    n_positive = int(contrasts["granule_configuration_positive"].sum()) if "granule_configuration_positive" in contrasts else 0
    n_contrasts = int(len(contrasts))
    p_config = pvalue_for_deltas(contrasts, "delta_configuration_score")
    p_construct = pvalue_for_deltas(contrasts, "delta_downstream_construction_balance")
    p_fate = pvalue_for_deltas(contrasts, "delta_regional_fate_balance")

    role_view = role_summary.set_index("specificity_role")
    lines = [
        "# Transcriptomic Configuration Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis tests whether granule cells are distinguished by a configuration of shared neuronal construction programs rather than by unique single pathways.",
        "",
        "## Configuration Scores",
        "",
        "- Downstream construction balance = mean(neurite/morphology, synaptic/excitability) - neurogenic niche/progenitor rank.",
        "- Regional fate balance = branch-matched fate rank - branch-opposed fate rank.",
        "- Combined configuration score = downstream construction balance + regional fate balance.",
        "",
        "## Main Result",
        "",
        f"- Named granule-versus-comparator contrasts tested: {n_contrasts}.",
        f"- Contrasts with positive combined configuration score: {n_positive}/{n_contrasts}.",
        f"- Wilcoxon p, combined configuration greater than comparator: {p_config:.3g}.",
        f"- Wilcoxon p, construction-over-niche balance greater than comparator: {p_construct:.3g}.",
        f"- Wilcoxon p, regional fate balance greater than comparator: {p_fate:.3g}.",
        "",
        "Median role-level combined configuration scores:",
    ]
    for role in PLOT_ROLE_ORDER:
        if role in role_view.index:
            row = role_view.loc[role]
            lines.append(f"- {ROLE_LABELS[role]}: {row['median_configuration_score']:.3f}.")

    lines.extend(["", "Sample-level contrasts:"])
    for _, row in contrasts.iterrows():
        lines.append(
            f"- {row['contrast_label']} ({row['sample']}): construction delta "
            f"{row['delta_downstream_construction_balance']:.3f}, fate-polarity delta "
            f"{row['delta_regional_fate_balance']:.3f}, combined configuration delta "
            f"{row['delta_configuration_score']:.3f}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This supports the user's idea that morphology may be encoded in the transcriptome as a configuration of shared neuronal construction programs, not as a unique granule-only gene list.",
            "- The configuration signal is stronger than single downstream-module specificity because it accounts for module balance: construction over progenitor/niche state plus correct regional fate polarity.",
            "- This remains transcriptomic inference. Protein localization, local translation, post-translational control, activity, and physical circuit constraints still require external validation.",
            "",
            "## Outputs",
            "",
            f"- Configuration units: `{rel(OUT_UNITS)}`",
            f"- Role summary: `{rel(OUT_ROLE_SUMMARY)}`",
            f"- Contrast summary: `{rel(OUT_CONTRASTS)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    config = build_configuration_units()
    role_summary = summarize_roles(config)
    contrasts = summarize_contrasts(config)
    plot_results(config, role_summary, contrasts)
    write_report(config, role_summary, contrasts)
    print(f"Wrote {rel(OUT_MD)}")
    print(role_summary[["specificity_role", "median_downstream_construction_balance", "median_regional_fate_balance", "median_configuration_score"]].to_string(index=False))
    print(contrasts[["dataset", "sample", "delta_downstream_construction_balance", "delta_regional_fate_balance", "delta_configuration_score"]].to_string(index=False))


if __name__ == "__main__":
    main()
