#!/usr/bin/env python3
"""Summarize manuscript candidate tiers into mechanism axes."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

IN_TIER = RESULTS / "primary_core_manuscript_candidate_tiers.tsv"

OUT_GENE_AXIS = RESULTS / "primary_core_mechanism_axis_gene_table.tsv"
OUT_AXIS_SUMMARY = RESULTS / "primary_core_mechanism_axis_summary.tsv"
OUT_BRANCH_SUMMARY = RESULTS / "primary_core_mechanism_axis_branch_summary.tsv"
OUT_PLOT = RESULTS / "primary_core_mechanism_axis_model.png"
OUT_MD = RESULTS / "primary_core_mechanism_axis_model.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BRANCHES = [
    ("selected_dentate", "median_delta_selected_dentate"),
    ("selected_cerebellar", "median_delta_selected_cerebellar"),
    ("full_dentate", "median_delta_full_matrix_dentate"),
    ("full_cerebellar", "median_delta_full_matrix_cerebellar"),
]

AXIS_ORDER = [
    "developmental_regulatory_control",
    "neurite_cytoskeleton_morphogenesis",
    "axon_guidance_adhesion",
    "synaptic_excitability_maturation",
    "exploratory_ortholog_completeness",
    "other_context",
]

AXIS_LABELS = {
    "developmental_regulatory_control": "Developmental regulatory control",
    "neurite_cytoskeleton_morphogenesis": "Neurite/cytoskeleton morphogenesis",
    "axon_guidance_adhesion": "Axon guidance and adhesion",
    "synaptic_excitability_maturation": "Synaptic/excitability maturation",
    "exploratory_ortholog_completeness": "Exploratory ortholog completeness",
    "other_context": "Other shared context",
}

GENE_AXIS_OVERRIDES = {
    "NFIA": "developmental_regulatory_control",
    "NFIB": "developmental_regulatory_control",
    "RFX3": "developmental_regulatory_control",
    "RFX7": "developmental_regulatory_control",
    "TCF4": "developmental_regulatory_control",
    "FOXN2": "developmental_regulatory_control",
    "BCL7A": "developmental_regulatory_control",
    "GPM6A": "neurite_cytoskeleton_morphogenesis",
    "DYNLL1": "neurite_cytoskeleton_morphogenesis",
    "BASP1": "neurite_cytoskeleton_morphogenesis",
    "MAPKAP1": "neurite_cytoskeleton_morphogenesis",
    "TUBA1A": "neurite_cytoskeleton_morphogenesis",
    "ACTB": "neurite_cytoskeleton_morphogenesis",
    "RTN3": "neurite_cytoskeleton_morphogenesis",
    "ACTG1": "neurite_cytoskeleton_morphogenesis",
    "TUBA1B": "neurite_cytoskeleton_morphogenesis",
    "MAP1B": "neurite_cytoskeleton_morphogenesis",
    "STMN2": "neurite_cytoskeleton_morphogenesis",
    "MAP3K13": "neurite_cytoskeleton_morphogenesis",
    "MAPK14": "neurite_cytoskeleton_morphogenesis",
    "TUBB": "neurite_cytoskeleton_morphogenesis",
    "ROBO2": "axon_guidance_adhesion",
    "CADM3": "axon_guidance_adhesion",
    "DCC": "axon_guidance_adhesion",
    "SEMA7A": "axon_guidance_adhesion",
    "KCNK1": "synaptic_excitability_maturation",
    "GABRA2": "synaptic_excitability_maturation",
    "GABRB3": "synaptic_excitability_maturation",
    "KCND2": "synaptic_excitability_maturation",
    "PPP3CA": "synaptic_excitability_maturation",
    "CACNA2D1": "synaptic_excitability_maturation",
    "KCNJ6": "synaptic_excitability_maturation",
    "GRIN2B": "synaptic_excitability_maturation",
    "KCNJ3": "synaptic_excitability_maturation",
    "STXBP5L": "synaptic_excitability_maturation",
    "KCNMB4": "synaptic_excitability_maturation",
    "SLC17A6": "synaptic_excitability_maturation",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def as_bool(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def assign_axis(row: pd.Series) -> str:
    gene = str(row["gene"])
    if int(row["tier_rank"]) == 5:
        return "exploratory_ortholog_completeness"
    if gene in GENE_AXIS_OVERRIDES:
        return GENE_AXIS_OVERRIDES[gene]

    text = " ".join(
        str(row.get(col, "")).lower()
        for col in ["mechanism_class", "suggested_role", "manuscript_tier", "mechanism_hit_tier"]
    )
    if "regulatory" in text or "transcription" in text or "ciliogenesis" in text:
        return "developmental_regulatory_control"
    if "axon" in text or "adhesion" in text or "guidance" in text:
        return "axon_guidance_adhesion"
    if (
        "synaptic" in text
        or "gaba" in text
        or "glutamat" in text
        or "potassium" in text
        or "channel" in text
        or "excitability" in text
        or "calcium" in text
    ):
        return "synaptic_excitability_maturation"
    if "cytoskeleton" in text or "neurite" in text or "structural" in text or "outgrowth" in text:
        return "neurite_cytoskeleton_morphogenesis"
    return "other_context"


def model_role(row: pd.Series) -> str:
    tier_rank = int(row["tier_rank"])
    axis = row["mechanism_axis"]
    if tier_rank == 1:
        if axis == "developmental_regulatory_control":
            return "Core upstream/regulatory seed"
        if axis == "neurite_cytoskeleton_morphogenesis":
            return "Core structural executor seed"
        if axis == "synaptic_excitability_maturation":
            return "Core excitability/synaptic seed"
        return "Core convergence seed"
    if tier_rank == 2:
        return "High-confidence pathway support"
    if tier_rank in {3, 4}:
        return "Supportive mechanism context"
    return "Exploratory ortholog-completeness candidate"


def support_for_branch(row: pd.Series, branch: str) -> bool:
    values = str(row.get("nominal_branch_support", "")).split(",")
    return branch in {value.strip() for value in values if value.strip()}


def load_gene_axis_table() -> pd.DataFrame:
    df = pd.read_csv(IN_TIER, sep="\t")
    df["mechanism_axis"] = df.apply(assign_axis, axis=1)
    df["mechanism_axis_label"] = df["mechanism_axis"].map(AXIS_LABELS).fillna(df["mechanism_axis"])
    df["model_role"] = df.apply(model_role, axis=1)
    df["axis_order"] = df["mechanism_axis"].map({axis: i for i, axis in enumerate(AXIS_ORDER)}).fillna(99).astype(int)
    for branch, delta_col in BRANCHES:
        df[f"{branch}_supported"] = df.apply(lambda row, b=branch: support_for_branch(row, b), axis=1)
        df[delta_col] = pd.to_numeric(df[delta_col], errors="coerce")
    return df.sort_values(
        ["axis_order", "tier_rank", "formal_n_nominal_branches", "formal_rank_priority_score"],
        ascending=[True, True, False, False],
    )


def summarize_axes(df: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for axis, sub in df.groupby("mechanism_axis", sort=False):
        tiers = {f"tier{tier}_genes": ",".join(sub.loc[sub["tier_rank"].eq(tier), "gene"].astype(str)) for tier in range(1, 6)}
        record = {
            "mechanism_axis": axis,
            "mechanism_axis_label": AXIS_LABELS.get(axis, axis),
            "n_genes": int(len(sub)),
            "n_tier1": int(sub["tier_rank"].eq(1).sum()),
            "n_tier2": int(sub["tier_rank"].eq(2).sum()),
            "n_tier3_4": int(sub["tier_rank"].isin([3, 4]).sum()),
            "n_tier5": int(sub["tier_rank"].eq(5).sum()),
            "median_formal_nominal_branches": float(pd.to_numeric(sub["formal_n_nominal_branches"], errors="coerce").median()),
            "median_formal_priority_score": float(pd.to_numeric(sub["formal_rank_priority_score"], errors="coerce").median()),
            "genes": ",".join(sub["gene"].astype(str)),
        }
        record.update(tiers)
        records.append(record)
    return pd.DataFrame(records)


def summarize_axis_branches(df: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for axis, sub in df.groupby("mechanism_axis", sort=False):
        for branch, delta_col in BRANCHES:
            supported = sub.loc[sub[f"{branch}_supported"]].copy()
            records.append(
                {
                    "mechanism_axis": axis,
                    "mechanism_axis_label": AXIS_LABELS.get(axis, axis),
                    "branch": branch,
                    "n_genes": int(len(sub)),
                    "n_supported_genes": int(len(supported)),
                    "supported_gene_fraction": float(len(supported) / len(sub)) if len(sub) else np.nan,
                    "median_delta_all_axis_genes": float(sub[delta_col].median()) if sub[delta_col].notna().any() else np.nan,
                    "median_delta_supported_genes": float(supported[delta_col].median())
                    if not supported.empty and supported[delta_col].notna().any()
                    else np.nan,
                    "supported_genes": ",".join(supported["gene"].astype(str)),
                }
            )
    return pd.DataFrame(records)


def plot_axis_model(df: pd.DataFrame) -> None:
    plot_df = df.loc[df["tier_rank"].le(4)].copy()
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values(
        ["axis_order", "tier_rank", "formal_n_nominal_branches", "formal_rank_priority_score"],
        ascending=[True, True, False, False],
    )
    values = plot_df[[col for _, col in BRANCHES]].to_numpy(dtype=float)
    labels = [f"T{int(row.tier_rank)} {row.gene}" for row in plot_df.itertuples()]

    fig_h = max(7.0, 0.31 * len(plot_df) + 2.2)
    fig, ax = plt.subplots(figsize=(9.0, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="PiYG", vmin=-0.5, vmax=0.5)
    ax.set_xticks(np.arange(len(BRANCHES)))
    ax.set_xticklabels([label.replace("_", " ") for label, _ in BRANCHES], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title("Mechanism axes across formal rank-meta branches")
    ax.set_xlabel("Formal screen/branch")
    ax.set_ylabel("Candidate genes grouped by mechanism axis")

    for i, (_, row) in enumerate(plot_df.iterrows()):
        for j, (branch, delta_col) in enumerate(BRANCHES):
            value = row[delta_col]
            support = as_bool(row.get(f"{branch}_supported", False))
            text = "NA" if pd.isna(value) else f"{float(value):.2f}{'*' if support else ''}"
            color = "#f6f6f6" if not pd.isna(value) and abs(float(value)) > 0.28 else "#202020"
            ax.text(j, i, text, ha="center", va="center", fontsize=7, color=color)

    last_axis = None
    for i, row in enumerate(plot_df.itertuples()):
        if last_axis is not None and row.mechanism_axis != last_axis:
            ax.axhline(i - 0.5, color="#202020", linewidth=0.8)
        last_axis = row.mechanism_axis

    axis_positions = (
        plot_df.reset_index(drop=True)
        .groupby("mechanism_axis", sort=False)
        .apply(lambda sub: (sub.index.min() + sub.index.max()) / 2, include_groups=False)
    )
    for axis, y in axis_positions.items():
        ax.text(
            len(BRANCHES) + 0.1,
            y,
            AXIS_LABELS.get(axis, axis),
            va="center",
            ha="left",
            fontsize=8,
        )

    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.14)
    cbar.set_label("Median dataset rank delta")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(df: pd.DataFrame, axis_summary: pd.DataFrame, branch_summary: pd.DataFrame) -> None:
    tier1 = df.loc[df["tier_rank"].eq(1)]
    tier2 = df.loc[df["tier_rank"].eq(2)]
    central = df.loc[df["tier_rank"].le(2)]

    lines = [
        "# Mechanism Axis Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This report translates the formal manuscript candidate tiers into biological mechanism axes. The goal is to support a paper argument that dentate and cerebellar granule cells do not share regional identity, but do converge on a downstream toolkit for compact neuronal morphology, neurite patterning, and excitability maturation.",
        "",
        "## Axis Definitions",
        "",
        "- Developmental regulatory control: transcriptional or regulatory candidates that may help coordinate granule-cell maturation programs.",
        "- Neurite/cytoskeleton morphogenesis: membrane, cytoskeletal, neurite-growth, and structural-executor candidates.",
        "- Axon guidance and adhesion: genes plausibly linking compact soma/neurite morphology to wiring and local circuit integration.",
        "- Synaptic/excitability maturation: receptor, ion-channel, calcium/signaling, and synaptic-release candidates.",
        "- Exploratory ortholog completeness: non-identical mouse/human ortholog hits retained for follow-up, not central claims.",
        "",
        "## Central Candidate Model",
        "",
    ]
    for axis in AXIS_ORDER:
        sub = central.loc[central["mechanism_axis"].eq(axis)]
        if sub.empty:
            continue
        genes = ", ".join(f"`{gene}`" for gene in sub["gene"].astype(str))
        lines.append(f"- {AXIS_LABELS.get(axis, axis)}: {genes}.")

    lines.extend(["", "## Tier 1 Seed Interpretation", ""])
    for _, row in tier1.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['mechanism_axis_label']}): {row['suggested_role']}; "
            f"{int(row['formal_n_nominal_branches'])}/{int(row['formal_n_available_branches'])} formal branches."
        )

    lines.extend(["", "## Tier 2 Pathway Support", ""])
    for _, row in tier2.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['mechanism_axis_label']}): {row['suggested_role']}; "
            f"{int(row['formal_n_nominal_branches'])}/{int(row['formal_n_available_branches'])} formal branches."
        )

    lines.extend(["", "## Axis Summary", ""])
    for _, row in axis_summary.iterrows():
        lines.append(
            f"- {row['mechanism_axis_label']}: {int(row['n_genes'])} genes "
            f"(Tier 1={int(row['n_tier1'])}, Tier 2={int(row['n_tier2'])}, "
            f"Tier 3/4={int(row['n_tier3_4'])}, Tier 5={int(row['n_tier5'])})."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The central model is not a single pathway. It is a layered toolkit: regulatory coordination, neurite/cytoskeleton execution, axon-guidance/adhesion, and synaptic/excitability maturation.",
            "- Tier 1 should anchor the manuscript model because it is both compact and robust across all available formal branches.",
            "- Tier 2 provides the strongest pathway-level biological breadth, especially synaptic/excitability maturation and ROBO2-mediated guidance.",
            "- Tier 5 should stay exploratory because non-identical ortholog recovery improves completeness but does not yet have the same manuscript-readiness as the curated mechanism tiers.",
            "",
            "## Outputs",
            "",
            f"- Gene-axis table: `{rel(OUT_GENE_AXIS)}`",
            f"- Axis summary: `{rel(OUT_AXIS_SUMMARY)}`",
            f"- Branch summary: `{rel(OUT_BRANCH_SUMMARY)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    gene_axis = load_gene_axis_table()
    axis_summary = summarize_axes(gene_axis)
    branch_summary = summarize_axis_branches(gene_axis)

    gene_axis.to_csv(OUT_GENE_AXIS, sep="\t", index=False)
    axis_summary.to_csv(OUT_AXIS_SUMMARY, sep="\t", index=False)
    branch_summary.to_csv(OUT_BRANCH_SUMMARY, sep="\t", index=False)
    plot_axis_model(gene_axis)
    write_report(gene_axis, axis_summary, branch_summary)

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Gene-axis rows: {len(gene_axis):,}")
    print(axis_summary[["mechanism_axis", "n_genes", "n_tier1", "n_tier2", "n_tier3_4", "n_tier5"]].to_string(index=False))


if __name__ == "__main__":
    main()
