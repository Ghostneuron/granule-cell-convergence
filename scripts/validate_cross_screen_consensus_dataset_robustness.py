#!/usr/bin/env python3
"""Dataset-aware validation for cross-screen consensus mechanism candidates."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

CONSENSUS = RESULTS / "primary_core_cross_screen_mechanism_consensus.tsv"
SELECTED_EXPR = RESULTS / "primary_core_expanded_gene_pseudobulk_expression.tsv.gz"
GENOME_EXPR = RESULTS / "primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz"

OUT_DELTAS = RESULTS / "primary_core_consensus_candidate_dataset_deltas.tsv"
OUT_SUMMARY = RESULTS / "primary_core_consensus_candidate_dataset_validation.tsv"
OUT_HEATMAP = RESULTS / "primary_core_consensus_candidate_dataset_validation_heatmap.png"
OUT_MD = RESULTS / "primary_core_consensus_candidate_dataset_validation.md"

MIN_CLASS_CELLS = 20

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DENTATE_BACKGROUNDS = {
    "non_dentate_background",
    "dentate_low_support",
    "other_or_ambiguous",
    "broad_neuronal_structural_warning",
}
CEREBELLAR_BACKGROUNDS = {"other_or_ambiguous", "broad_neuronal_structural_warning"}


def load_candidates() -> pd.DataFrame:
    consensus = pd.read_csv(CONSENSUS, sep="\t")
    candidates = consensus.loc[consensus["consensus_tier"].eq("consensus_figure_candidate")].copy()
    candidates = candidates[
        [
            "gene",
            "canonical_gene",
            "selected_mechanism_class",
            "genome_mechanism_class",
            "combined_priority_score",
        ]
    ].drop_duplicates("canonical_gene")
    return candidates


def load_expression(path: Path, screen: str, genes: set[str]) -> pd.DataFrame:
    usecols = [
        "dataset",
        "core_branch",
        "sample",
        "broad_class",
        "n_cells",
        "gene",
        "canonical_gene",
        "detection_fraction",
        "mean_log1p_expression",
        "eligible_class",
        "mean_log1p_rank_within_sample_gene",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    df = df.loc[df["canonical_gene"].isin(genes)].copy()
    df["screen"] = screen
    df["eligible_class"] = df["eligible_class"].astype(str).str.lower().isin(["true", "1", "yes"])
    return df


def summarize_one_unit(sub: pd.DataFrame, branch: str) -> dict[str, object] | None:
    if branch == "dentate":
        target_class = "dentate_candidate"
        backgrounds = DENTATE_BACKGROUNDS
    elif branch == "cerebellar":
        target_class = "cerebellar_candidate"
        backgrounds = CEREBELLAR_BACKGROUNDS
    else:
        raise ValueError(branch)

    eligible = sub.loc[sub["eligible_class"] & sub["n_cells"].ge(MIN_CLASS_CELLS)].copy()
    target = eligible.loc[eligible["broad_class"].eq(target_class)]
    background = eligible.loc[eligible["broad_class"].isin(backgrounds)]
    if target.empty or background.empty:
        return None
    target_rank = target["mean_log1p_rank_within_sample_gene"].median()
    background_rank = background["mean_log1p_rank_within_sample_gene"].median()
    return {
        "candidate_rank": target_rank,
        "background_rank": background_rank,
        "rank_delta": target_rank - background_rank,
        "candidate_detection": target["detection_fraction"].median(),
        "background_detection": background["detection_fraction"].median(),
        "detection_delta": target["detection_fraction"].median() - background["detection_fraction"].median(),
        "candidate_n_cells": int(target["n_cells"].sum()),
        "background_n_cells": int(background["n_cells"].sum()),
        "n_background_classes": int(background["broad_class"].nunique()),
    }


def build_deltas(expr: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    meta = candidates.set_index("canonical_gene").to_dict("index")
    rows: list[dict[str, object]] = []
    for (screen, dataset, sample, gene), sub in expr.groupby(["screen", "dataset", "sample", "canonical_gene"], sort=False):
        core_branch = sub["core_branch"].iloc[0]
        branches: list[str] = []
        if core_branch in {"mouse_dentate", "human_dentate_hippocampus"}:
            branches.append("dentate")
        if core_branch == "cerebellum":
            branches.append("cerebellar")
        for branch in branches:
            out = summarize_one_unit(sub, branch)
            if out is None:
                continue
            m = meta.get(gene, {})
            rows.append(
                {
                    "screen": screen,
                    "dataset": dataset,
                    "sample": sample,
                    "core_branch": core_branch,
                    "branch_tested": branch,
                    "gene": m.get("gene", sub["gene"].iloc[0]),
                    "canonical_gene": gene,
                    "selected_mechanism_class": m.get("selected_mechanism_class", ""),
                    "genome_mechanism_class": m.get("genome_mechanism_class", ""),
                    "combined_priority_score": m.get("combined_priority_score", np.nan),
                    **out,
                    "positive_delta": out["rank_delta"] > 0,
                }
            )
    return pd.DataFrame(rows)


def sign_p(n_positive: int, n_total: int) -> float:
    if n_total <= 0:
        return np.nan
    return float(stats.binomtest(n_positive, n_total, 0.5, alternative="greater").pvalue)


def summarize_deltas(deltas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (gene, screen, branch), sub in deltas.groupby(["canonical_gene", "screen", "branch_tested"], sort=False):
        n_units = len(sub)
        n_pos = int(sub["positive_delta"].sum())
        rows.append(
            {
                "gene": sub["gene"].iloc[0],
                "canonical_gene": gene,
                "screen": screen,
                "branch_tested": branch,
                "mechanism_class": sub["selected_mechanism_class"].iloc[0]
                if screen == "selected"
                else sub["genome_mechanism_class"].iloc[0],
                "n_units": n_units,
                "n_datasets": sub["dataset"].nunique(),
                "n_positive_units": n_pos,
                "positive_unit_fraction": n_pos / n_units if n_units else np.nan,
                "median_rank_delta": sub["rank_delta"].median(),
                "min_rank_delta": sub["rank_delta"].min(),
                "median_detection_delta": sub["detection_delta"].median(),
                "sign_test_p_greater": sign_p(n_pos, n_units),
                "datasets_positive": ",".join(sorted(sub.loc[sub["positive_delta"], "dataset"].unique())),
                "datasets_nonpositive": ",".join(sorted(sub.loc[~sub["positive_delta"], "dataset"].unique())),
            }
        )
    summary = pd.DataFrame(rows)
    summary["robust_branch_support"] = (
        summary["n_units"].ge(2)
        & summary["positive_unit_fraction"].ge(0.75)
        & summary["median_rank_delta"].gt(0)
    )
    wide = (
        summary.pivot_table(
            index=["gene", "canonical_gene"],
            columns=["screen", "branch_tested"],
            values=["positive_unit_fraction", "median_rank_delta", "robust_branch_support"],
            aggfunc="first",
        )
        .reset_index()
    )
    wide.columns = [
        "_".join(str(part) for part in col if str(part) != "") if isinstance(col, tuple) else col for col in wide.columns
    ]
    branch_cols = [col for col in wide.columns if col.startswith("robust_branch_support_")]
    if branch_cols:
        wide["robust_all_available_branches"] = wide[branch_cols].fillna(False).all(axis=1)
        wide["n_robust_screen_branches"] = wide[branch_cols].fillna(False).sum(axis=1)
    else:
        wide["robust_all_available_branches"] = False
        wide["n_robust_screen_branches"] = 0
    summary = summary.merge(
        wide[["canonical_gene", "robust_all_available_branches", "n_robust_screen_branches"]],
        on="canonical_gene",
        how="left",
    )
    summary = summary.sort_values(
        ["robust_all_available_branches", "n_robust_screen_branches", "median_rank_delta"],
        ascending=[False, False, False],
    )
    return summary


def plot_heatmap(summary: pd.DataFrame, candidates: pd.DataFrame) -> None:
    ordered_genes = candidates["canonical_gene"].tolist()
    matrix = summary.pivot_table(
        index="canonical_gene",
        columns=["screen", "branch_tested"],
        values="median_rank_delta",
        aggfunc="first",
    ).reindex(ordered_genes)
    if matrix.empty:
        return
    labels = [f"{a}:{b}" for a, b in matrix.columns]
    fig, ax = plt.subplots(figsize=(7.8, max(5.2, 0.28 * len(matrix))))
    data = matrix.to_numpy(dtype=float)
    im = ax.imshow(data, aspect="auto", cmap="PRGn", vmin=-0.5, vmax=0.5)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(matrix.index)))
    gene_labels = candidates.set_index("canonical_gene").loc[matrix.index, "gene"].tolist()
    ax.set_yticklabels(gene_labels)
    ax.set_title("Consensus candidate dataset-aware median rank deltas")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("candidate rank - local background rank")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_HEATMAP, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(deltas: pd.DataFrame, summary: pd.DataFrame, candidates: pd.DataFrame) -> None:
    robust_genes = (
        summary.loc[summary["robust_all_available_branches"], ["gene", "canonical_gene"]]
        .drop_duplicates("canonical_gene")
        .sort_values("gene")
    )
    top = (
        summary.groupby(["gene", "canonical_gene"], as_index=False)
        .agg(
            n_screen_branches=("branch_tested", "size"),
            n_robust_screen_branches=("robust_branch_support", "sum"),
            median_rank_delta=("median_rank_delta", "median"),
            min_positive_fraction=("positive_unit_fraction", "min"),
        )
        .sort_values(["n_robust_screen_branches", "median_rank_delta"], ascending=[False, False])
    )
    lines = [
        "# Consensus Candidate Dataset-Aware Validation",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This validation asks whether the 24 cross-screen consensus candidates remain positive across individual datasets/samples, rather than only in pooled pseudobulk summaries.",
        "",
        "## Summary",
        "",
        f"- Consensus candidates tested: {candidates['canonical_gene'].nunique()}.",
        f"- Dataset/sample/gene branch-delta rows: {len(deltas):,}.",
        f"- Genes robust across all available screen/branch tests: {len(robust_genes)}.",
        "",
        "## Robust Genes",
        "",
    ]
    if robust_genes.empty:
        lines.append("- No gene met the all-available-branch robustness rule.")
    else:
        lines.append("- " + ", ".join(f"`{gene}`" for gene in robust_genes["gene"]))
    lines.extend(["", "## Top Validation Summary", ""])
    for _, row in top.head(24).iterrows():
        lines.append(
            f"- `{row['gene']}`: {int(row['n_robust_screen_branches'])}/{int(row['n_screen_branches'])} "
            f"screen-branch tests robust; median delta {row['median_rank_delta']:.3f}; "
            f"minimum positive-unit fraction {row['min_positive_fraction']:.2f}."
        )
    lines.extend(
        [
            "",
            "## Robustness Rule",
            "",
            f"A screen/branch is robust if it has at least 2 dataset/sample units, >=75% positive deltas, and median rank delta > 0.",
            "",
            "## Outputs",
            "",
            f"- Dataset deltas: `{OUT_DELTAS.relative_to(ROOT)}`",
            f"- Validation summary: `{OUT_SUMMARY.relative_to(ROOT)}`",
            f"- Heatmap: `{OUT_HEATMAP.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    candidates = load_candidates()
    genes = set(candidates["canonical_gene"])
    selected = load_expression(SELECTED_EXPR, "selected", genes)
    genome = load_expression(GENOME_EXPR, "genome", genes)
    expr = pd.concat([selected, genome], ignore_index=True, sort=False)
    deltas = build_deltas(expr, candidates)
    summary = summarize_deltas(deltas)
    plot_heatmap(summary, candidates)
    deltas.to_csv(OUT_DELTAS, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    write_report(deltas, summary, candidates)
    print(f"Wrote {len(deltas):,} dataset-aware delta rows")
    print(f"Wrote {len(summary):,} validation summary rows")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
