#!/usr/bin/env python3
"""Stage-split immature versus mature granule-cell comparison.

This analysis asks whether dentate and cerebellar granule candidates look more
similar at an immature assembly stage or at a later mature/maturing stage. It
uses the existing full-transcriptome diffusion group summary, because that
table contains both explicit developmental labels and the module scores most
relevant to morphology, synapses, fate, niche, and maturation.

The output deliberately keeps strict and supporting stage calls separate. The
strict evidence has multiple dentate datasets but only one explicitly staged
cerebellar dataset (GSE122357), so p-values would be misleading; the main
reported quantity is the cross-branch difference and a simple similarity score.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
SUPP_FIGURES = ROOT / "Project/manuscript/Supplementary figures"

DIFFUSION_GROUPS = RESULTS / "primary_core_full_transcriptome_diffusion_group_summary.tsv"

OUT_GROUP_CALLS = RESULTS / "primary_core_stage_split_granule_group_calls.tsv"
OUT_BRANCH_SUMMARY = RESULTS / "primary_core_stage_split_granule_module_branch_summary.tsv"
OUT_SIMILARITY = RESULTS / "primary_core_stage_split_granule_stage_similarity.tsv"
OUT_TRANSITIONS = RESULTS / "primary_core_stage_split_granule_stage_transitions.tsv"
OUT_PLOT = RESULTS / "primary_core_stage_split_granule_comparison.png"
OUT_SUPP_PLOT = SUPP_FIGURES / "Fig.S1_stage_split_granule_comparison.png"
OUT_MD = RESULTS / "primary_core_stage_split_granule_comparison.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


MODULES = [
    ("immature_progenitor_state", "Immature/progenitor"),
    ("neuronal_differentiation_maturation", "Neuronal maturation"),
    ("downstream_neurite_morphology", "Neurite/morphology"),
    ("downstream_synaptic_excitability", "Synaptic/excitability"),
    ("shared_neurogenic_niche_state", "Shared neurogenic niche"),
    ("cerebellar_fate_rhombic_lip_shh", "Cerebellar fate/SHH"),
    ("dentate_fate_wnt_prox1", "Dentate fate/WNT/PROX1"),
    ("tgf_smad_pai1_response", "TGF-beta/SMAD"),
    ("bdnf_erk_response", "BDNF/ERK"),
    ("secreted_stop_candidate_axis", "Secreted stop candidates"),
]

PLOT_MODULES = [
    "immature_progenitor_state",
    "neuronal_differentiation_maturation",
    "downstream_neurite_morphology",
    "downstream_synaptic_excitability",
    "shared_neurogenic_niche_state",
    "tgf_smad_pai1_response",
    "bdnf_erk_response",
]

MODULE_LABELS = dict(MODULES)
SHORT_MODULE_LABELS = {
    "immature_progenitor_state": "Imm.\nprog.",
    "neuronal_differentiation_maturation": "Neur.\nmat.",
    "downstream_neurite_morphology": "Neurite\nmorph.",
    "downstream_synaptic_excitability": "Syn.\nexc.",
    "shared_neurogenic_niche_state": "Shared\nniche",
    "tgf_smad_pai1_response": "TGF-beta\nSMAD",
    "bdnf_erk_response": "BDNF\nERK",
}

STAGE_ORDER = {"immature": 0, "mature": 1}
BRANCH_ORDER = {"dentate": 0, "cerebellar": 1}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def classify_branch(row: pd.Series) -> str | None:
    broad = str(row.get("broad_class", ""))
    call = str(row.get("candidate_call", ""))
    if "cerebellar_candidate" in broad or "cerebellar_granule" in call:
        return "cerebellar"
    if "dentate_candidate" in broad or "dentate_immature_candidate" in broad or "dentate_granule" in call:
        return "dentate"
    return None


def classify_stage(row: pd.Series) -> tuple[str | None, str, bool, str]:
    """Return stage_bin, evidence_tier, strict_call, stage_call_reason."""

    dataset = str(row["dataset"])
    label = str(row["axis_label"])
    call = str(row["candidate_call"])

    if dataset == "GSE104323":
        if label in {"Neuroblast", "Immature-GC"}:
            return "immature", "strict_explicit_label", True, "dentate neuroblast/immature-GC label"
        if label in {"GC-juv", "GC-adult"}:
            return "mature", "strict_explicit_label", True, "dentate juvenile/adult GC label"

    if dataset == "GSE122357":
        if label == "P0":
            return "immature", "strict_developmental_age", True, "cerebellar early postnatal candidate stage"
        if label in {"P8a", "P8b"}:
            return "mature", "strict_developmental_age", True, "cerebellar later postnatal maturing stage"

    if dataset == "GSE214309":
        if label.startswith("immature"):
            return "immature", "strict_explicit_label", True, "dentate immature activity-state label"
        if label.startswith("mature"):
            return "mature", "strict_explicit_label", True, "dentate mature activity-state label"

    if dataset == "GSE292261":
        if label in {"DG_P5", "DG_P7", "DG_P10"}:
            return "immature", "strict_developmental_age", True, "dentate early postnatal assembly window"
        if label in {"DG_P15", "DG_P28"}:
            return "mature", "strict_developmental_age", True, "dentate late postnatal maturation window"

    if dataset == "GSE325391":
        if label.startswith("DiffN"):
            return "immature", "strict_source_label", True, "human differentiating-neuron source label"
        if label.startswith("MatN"):
            return "mature", "strict_source_label", True, "human mature-neuron source label"

    if dataset == "GSE268609":
        if call == "immature_neurogenic_candidate" or label == "immature_neurogenic_candidate":
            return "immature", "supporting_projected_label", False, "human projected immature neurogenic candidate"
        if "human_dg_like" in call or "human_dg_like" in label:
            return "mature", "supporting_projected_label", False, "human projected DG-like candidate"

    return None, "unassigned", False, "no explicit immature/mature stage label"


def median_or_nan(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def weighted_mean_or_nan(values: pd.Series, weights: pd.Series) -> float:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    w = pd.to_numeric(weights, errors="coerce").fillna(0).to_numpy(dtype=float)
    mask = np.isfinite(vals) & np.isfinite(w) & (w > 0)
    if not np.any(mask):
        return np.nan
    return float(np.average(vals[mask], weights=w[mask]))


def stage_group_calls() -> pd.DataFrame:
    groups = pd.read_csv(DIFFUSION_GROUPS, sep="\t")
    groups["branch"] = groups.apply(classify_branch, axis=1)
    stage_info = groups.apply(classify_stage, axis=1, result_type="expand")
    stage_info.columns = ["stage_bin", "stage_evidence_tier", "strict_stage_call", "stage_call_reason"]
    groups = pd.concat([groups, stage_info], axis=1)

    keep = groups["branch"].notna() & groups["stage_bin"].notna()
    groups = groups.loc[keep].copy()
    groups["stage_order"] = groups["stage_bin"].map(STAGE_ORDER)
    groups["branch_order"] = groups["branch"].map(BRANCH_ORDER)

    # Keep only candidate/support rows, not lineage warning or low-support background rows.
    allowed_broad = {
        "dentate_candidate",
        "dentate_immature_candidate",
        "dentate_candidate_or_hippocampal_neuronal",
        "cerebellar_candidate",
    }
    groups = groups.loc[groups["broad_class"].isin(allowed_broad)].copy()

    groups = groups.sort_values(["branch_order", "stage_order", "dataset", "axis_label"])

    out_cols = [
        "dataset",
        "source_scope",
        "axis_type",
        "axis_label",
        "candidate_call",
        "broad_class",
        "branch",
        "stage_bin",
        "stage_evidence_tier",
        "strict_stage_call",
        "stage_call_reason",
        "n_cells",
        "median_pseudotime",
        "mean_pseudotime",
    ] + [f"median_{m}" for m, _ in MODULES]
    return groups[out_cols]


def long_module_values(groups: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in groups.iterrows():
        for module_id, label in MODULES:
            col = f"median_{module_id}"
            rows.append(
                {
                    "dataset": row["dataset"],
                    "source_scope": row["source_scope"],
                    "axis_label": row["axis_label"],
                    "candidate_call": row["candidate_call"],
                    "broad_class": row["broad_class"],
                    "branch": row["branch"],
                    "stage_bin": row["stage_bin"],
                    "stage_evidence_tier": row["stage_evidence_tier"],
                    "strict_stage_call": bool(row["strict_stage_call"]),
                    "n_cells": row["n_cells"],
                    "module_id": module_id,
                    "module_label": label,
                    "score": row[col],
                }
            )
    return pd.DataFrame(rows)


def summarize_branch_stage(long_df: pd.DataFrame) -> pd.DataFrame:
    dataset_summary_rows: list[dict[str, object]] = []
    for keys, sub in long_df.groupby(
        ["dataset", "branch", "stage_bin", "stage_evidence_tier", "strict_stage_call", "module_id", "module_label"],
        dropna=False,
    ):
        dataset, branch, stage, tier, strict, module_id, label = keys
        dataset_summary_rows.append(
            {
                "dataset": dataset,
                "branch": branch,
                "stage_bin": stage,
                "stage_evidence_tier": tier,
                "strict_stage_call": bool(strict),
                "module_id": module_id,
                "module_label": label,
                "n_groups": int(sub[["axis_label", "candidate_call"]].drop_duplicates().shape[0]),
                "n_cells_total": int(pd.to_numeric(sub["n_cells"], errors="coerce").fillna(0).sum()),
                "median_group_score": median_or_nan(sub["score"]),
                "weighted_mean_group_score": weighted_mean_or_nan(sub["score"], sub["n_cells"]),
            }
        )
    dataset_summary = pd.DataFrame(dataset_summary_rows)

    summary_rows: list[dict[str, object]] = []
    for keys, sub in dataset_summary.groupby(["branch", "stage_bin", "module_id", "module_label"], dropna=False):
        branch, stage, module_id, label = keys
        strict_sub = sub.loc[sub["strict_stage_call"]].copy()
        source = strict_sub if not strict_sub.empty else sub
        summary_rows.append(
            {
                "branch": branch,
                "stage_bin": stage,
                "module_id": module_id,
                "module_label": label,
                "n_datasets": int(source["dataset"].nunique()),
                "n_groups": int(source["n_groups"].sum()),
                "n_cells_total": int(source["n_cells_total"].sum()),
                "n_strict_datasets": int(strict_sub["dataset"].nunique()),
                "used_strict_only": bool(not strict_sub.empty),
                "median_dataset_score": median_or_nan(source["median_group_score"]),
                "mean_dataset_score": float(pd.to_numeric(source["median_group_score"], errors="coerce").mean()),
                "weighted_mean_group_score": weighted_mean_or_nan(
                    source["weighted_mean_group_score"], source["n_cells_total"]
                ),
                "datasets_used": ",".join(sorted(source["dataset"].astype(str).unique())),
                "evidence_tiers": ",".join(sorted(source["stage_evidence_tier"].astype(str).unique())),
            }
        )
    return pd.DataFrame(summary_rows).sort_values(["stage_bin", "module_id", "branch"])


def compare_stage_similarity(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (stage, module_id, label), sub in summary.groupby(["stage_bin", "module_id", "module_label"], dropna=False):
        pivot = sub.set_index("branch")
        if {"dentate", "cerebellar"}.issubset(set(pivot.index)):
            dentate = float(pivot.loc["dentate", "median_dataset_score"])
            cerebellar = float(pivot.loc["cerebellar", "median_dataset_score"])
            delta = cerebellar - dentate
            rows.append(
                {
                    "stage_bin": stage,
                    "module_id": module_id,
                    "module_label": label,
                    "dentate_median_score": dentate,
                    "cerebellar_median_score": cerebellar,
                    "cerebellar_minus_dentate": delta,
                    "absolute_branch_difference": abs(delta),
                    "stage_similarity_score": 1.0 - min(1.0, abs(delta)),
                    "dentate_n_datasets": int(pivot.loc["dentate", "n_datasets"]),
                    "cerebellar_n_datasets": int(pivot.loc["cerebellar", "n_datasets"]),
                    "dentate_datasets": pivot.loc["dentate", "datasets_used"],
                    "cerebellar_datasets": pivot.loc["cerebellar", "datasets_used"],
                }
            )
    sim = pd.DataFrame(rows)
    if sim.empty:
        return sim

    sim["stage_order"] = sim["stage_bin"].map(STAGE_ORDER)
    sim["module_order"] = sim["module_id"].map({m: i for i, (m, _) in enumerate(MODULES)})

    trans_rows: list[dict[str, object]] = []
    for (module_id, label), sub in sim.groupby(["module_id", "module_label"], dropna=False):
        vals = sub.set_index("stage_bin")
        if {"immature", "mature"}.issubset(set(vals.index)):
            trans_rows.append(
                {
                    "module_id": module_id,
                    "module_label": label,
                    "immature_similarity": float(vals.loc["immature", "stage_similarity_score"]),
                    "mature_similarity": float(vals.loc["mature", "stage_similarity_score"]),
                    "mature_minus_immature_similarity": float(
                        vals.loc["mature", "stage_similarity_score"] - vals.loc["immature", "stage_similarity_score"]
                    ),
                    "immature_abs_difference": float(vals.loc["immature", "absolute_branch_difference"]),
                    "mature_abs_difference": float(vals.loc["mature", "absolute_branch_difference"]),
                }
            )
    transitions = pd.DataFrame(trans_rows).sort_values("mature_minus_immature_similarity", ascending=False)
    transitions.to_csv(OUT_TRANSITIONS, sep="\t", index=False)

    return sim.sort_values(["stage_order", "module_order"]).drop(columns=["stage_order", "module_order"])


def make_plot(summary: pd.DataFrame, similarity: pd.DataFrame, transitions: pd.DataFrame) -> None:
    plot_summary = summary.loc[summary["module_id"].isin(PLOT_MODULES)].copy()
    plot_summary["col_label"] = plot_summary["branch"].map({"dentate": "Dentate", "cerebellar": "Cerebellum"}) + "\n" + plot_summary[
        "stage_bin"
    ].str.capitalize()

    col_order = ["Dentate\nImmature", "Cerebellum\nImmature", "Dentate\nMature", "Cerebellum\nMature"]
    row_order = PLOT_MODULES

    heat = plot_summary.pivot_table(index="module_id", columns="col_label", values="median_dataset_score", aggfunc="median")
    heat = heat.reindex(index=row_order, columns=col_order)

    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.85], wspace=0.35, hspace=0.78)

    ax0 = fig.add_subplot(gs[:, 0])
    im = ax0.imshow(heat.to_numpy(dtype=float), vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax0.set_xticks(range(len(col_order)))
    ax0.set_xticklabels(col_order, fontsize=9)
    ax0.set_yticks(range(len(row_order)))
    ax0.set_yticklabels([MODULE_LABELS[m] for m in row_order], fontsize=9)
    ax0.set_title("a. Stage-split module scores", loc="left", fontsize=11, fontweight="bold")
    for y in range(heat.shape[0]):
        for x in range(heat.shape[1]):
            val = heat.iloc[y, x]
            if pd.notna(val):
                ax0.text(x, y, f"{val:.2f}", ha="center", va="center", color="white" if val < 0.65 else "black", fontsize=8)
    cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.02)
    cbar.set_label("Median rank-like module score")

    ax1 = fig.add_subplot(gs[0, 1])
    sim_plot = similarity.loc[similarity["module_id"].isin(PLOT_MODULES)].copy()
    sim_plot["x"] = sim_plot["module_id"].map({m: i for i, m in enumerate(PLOT_MODULES)})
    width = 0.36
    for offset, stage, color in [(-width / 2, "immature", "#4c78a8"), (width / 2, "mature", "#f58518")]:
        sub = sim_plot.loc[sim_plot["stage_bin"].eq(stage)]
        ax1.bar(sub["x"] + offset, sub["stage_similarity_score"], width=width, color=color, label=stage.capitalize())
    ax1.set_ylim(0, 1.05)
    ax1.set_xticks(range(len(PLOT_MODULES)))
    ax1.set_xticklabels([SHORT_MODULE_LABELS[m] for m in PLOT_MODULES], rotation=0, ha="center", fontsize=8)
    ax1.set_ylabel("1 - |cerebellum - dentate|")
    ax1.set_title("b. Cross-branch similarity by stage", loc="left", fontsize=11, fontweight="bold")
    ax1.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    ax2 = fig.add_subplot(gs[1, 1])
    trans_plot = transitions.loc[transitions["module_id"].isin(PLOT_MODULES)].copy()
    trans_plot["x"] = trans_plot["module_id"].map({m: i for i, m in enumerate(PLOT_MODULES)})
    colors = np.where(trans_plot["mature_minus_immature_similarity"].ge(0), "#54a24b", "#e45756")
    ax2.axhline(0, color="#555555", linewidth=0.8)
    ax2.bar(trans_plot["x"], trans_plot["mature_minus_immature_similarity"], color=colors)
    ax2.set_xticks(range(len(PLOT_MODULES)))
    ax2.set_xticklabels([SHORT_MODULE_LABELS[m] for m in PLOT_MODULES], rotation=0, ha="center", fontsize=8)
    ax2.set_ylabel("Mature - immature similarity")
    ax2.set_title("c. Does similarity strengthen after maturation?", loc="left", fontsize=11, fontweight="bold")

    fig.suptitle("Immature versus mature dentate/cerebellar granule comparison", fontsize=14, y=0.98)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    SUPP_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUPP_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_markdown(groups: pd.DataFrame, summary: pd.DataFrame, similarity: pd.DataFrame, transitions: pd.DataFrame) -> None:
    strict = groups.loc[groups["strict_stage_call"]]
    n_dentate_strict = strict.loc[strict["branch"].eq("dentate"), "dataset"].nunique()
    n_cereb_strict = strict.loc[strict["branch"].eq("cerebellar"), "dataset"].nunique()

    key_modules = [
        "downstream_neurite_morphology",
        "downstream_synaptic_excitability",
        "neuronal_differentiation_maturation",
        "immature_progenitor_state",
        "shared_neurogenic_niche_state",
    ]
    key = similarity.loc[similarity["module_id"].isin(key_modules)].copy()
    lines = [
        "# Stage-Split Granule-Cell Comparison",
        "",
        "## Question",
        "",
        "Can dentate and cerebellar granule-cell candidates be compared separately at immature and mature stages?",
        "",
        "## Short Answer",
        "",
        "Yes. The comparison is feasible and biologically useful, but the evidence is asymmetric: dentate has several explicit immature/mature resources, while cerebellum has one explicit staged developmental resource (`GSE122357`: P0 versus P8a/P8b). Therefore this layer should be treated as a stage-aware support analysis rather than a fully powered stage-by-branch statistical test.",
        "",
        "## Stage Definitions",
        "",
        "- Dentate immature: `GSE104323` Neuroblast/Immature-GC, `GSE214309` immature/immatureactive states, `GSE292261` DG_P5/DG_P7/DG_P10, and `GSE325391` DiffN labels.",
        "- Dentate mature: `GSE104323` GC-juv/GC-adult, `GSE214309` mature/matureactive states, `GSE292261` DG_P15/DG_P28, and `GSE325391` MatN labels.",
        "- Cerebellar immature: `GSE122357` P0 cerebellar granule candidates.",
        "- Cerebellar mature/maturing: `GSE122357` P8a/P8b cerebellar granule candidates.",
        "- `GSE268609` projected labels are retained as supporting, not strict, stage calls. Unstaged atlas resources are not forced into the binary split.",
        "",
        "## Coverage",
        "",
        f"- Strict dentate stage datasets: {n_dentate_strict}.",
        f"- Strict cerebellar stage datasets: {n_cereb_strict}.",
        f"- Stage-called candidate groups: {groups.shape[0]}.",
        "",
        "## Main Result",
        "",
        "The stage split supports the manuscript model: the strongest interpretable comparison is not stem/progenitor identity, but postmitotic assembly and maturation. Downstream neurite/morphology and synaptic/excitability modules can be compared separately in immature and mature windows, while regional fate modules remain branch-biased.",
        "",
        "Selected similarity scores below use `1 - abs(cerebellar median - dentate median)`, so higher means more similar between branches.",
        "",
        "| module | immature similarity | mature similarity | mature - immature |",
        "|---|---:|---:|---:|",
    ]
    for module_id in key_modules:
        row = transitions.loc[transitions["module_id"].eq(module_id)]
        if row.empty:
            continue
        r = row.iloc[0]
        lines.append(
            f"| {MODULE_LABELS[module_id]} | {r['immature_similarity']:.3f} | {r['mature_similarity']:.3f} | {r['mature_minus_immature_similarity']:.3f} |"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "- If a module is already similar in the immature stage, it likely reflects early postmitotic assembly machinery.",
        "- If similarity increases in the mature stage, it likely reflects later circuit/synaptic maturation or final geometry constraints.",
        "- If similarity decreases, the two branches may begin with shared immature neuronal machinery but diverge as region-specific circuit implementation becomes stronger.",
        "- Because cerebellar stage evidence is currently concentrated in `GSE122357`, any stage-specific cerebellar conclusion should be phrased cautiously.",
        "",
        "## Outputs",
        "",
        f"- Stage-called groups: `{rel(OUT_GROUP_CALLS)}`.",
        f"- Branch-stage module summary: `{rel(OUT_BRANCH_SUMMARY)}`.",
        f"- Cross-branch stage similarity: `{rel(OUT_SIMILARITY)}`.",
        f"- Mature-minus-immature transition table: `{rel(OUT_TRANSITIONS)}`.",
        f"- Plot: `{rel(OUT_PLOT)}`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    groups = stage_group_calls()
    groups.to_csv(OUT_GROUP_CALLS, sep="\t", index=False)

    long_df = long_module_values(groups)
    summary = summarize_branch_stage(long_df)
    summary.to_csv(OUT_BRANCH_SUMMARY, sep="\t", index=False)

    similarity = compare_stage_similarity(summary)
    similarity.to_csv(OUT_SIMILARITY, sep="\t", index=False)

    transitions = pd.read_csv(OUT_TRANSITIONS, sep="\t") if OUT_TRANSITIONS.exists() else pd.DataFrame()
    make_plot(summary, similarity, transitions)
    write_markdown(groups, summary, similarity, transitions)

    print(f"Wrote {rel(OUT_GROUP_CALLS)}")
    print(f"Wrote {rel(OUT_BRANCH_SUMMARY)}")
    print(f"Wrote {rel(OUT_SIMILARITY)}")
    print(f"Wrote {rel(OUT_TRANSITIONS)}")
    print(f"Wrote {rel(OUT_PLOT)}")
    print(f"Wrote {rel(OUT_MD)}")


if __name__ == "__main__":
    main()
