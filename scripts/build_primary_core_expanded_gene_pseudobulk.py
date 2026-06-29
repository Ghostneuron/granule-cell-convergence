#!/usr/bin/env python3
"""Expanded primary-core pseudobulk screen over the 2,169 human-core selected genes.

This script reuses the primary-core matrix readers from the 67-gene candidate
screen, but swaps in the full selected-gene universe from the constructed human
core object. It is the bridge between the focused candidate-gene result and a
future whole-transcriptome, ortholog-curated mixed-effect DE analysis.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import build_primary_core_candidate_gene_pseudobulk as base


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"

SELECTED_VAR = PROCESSED / "human_core_normalized_reduced_object/var.tsv"
CANDIDATE_PACKET = RESULTS / "human_bridge_candidate_gene_packet.tsv"

OUT_EXPR = RESULTS / "primary_core_expanded_gene_pseudobulk_expression.tsv.gz"
OUT_COVERAGE = RESULTS / "primary_core_expanded_gene_pseudobulk_coverage.tsv"
OUT_STATS = RESULTS / "primary_core_expanded_gene_pseudobulk_statistics.tsv"
OUT_HITS = RESULTS / "primary_core_expanded_gene_pseudobulk_shared_hits.tsv"
OUT_BRANCH = RESULTS / "primary_core_expanded_gene_pseudobulk_branch_specific.tsv"
OUT_PLOT = RESULTS / "primary_core_expanded_gene_pseudobulk_shared_hits.png"
OUT_MD = RESULTS / "primary_core_expanded_gene_pseudobulk_analysis.md"

MIN_CLASS_CELLS = base.MIN_CLASS_CELLS
MIN_DETECTION = 0.05

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def canon_gene(gene: object) -> str:
    return base.canon_gene(gene)


def load_expanded_targets() -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    selected = pd.read_csv(SELECTED_VAR, sep="\t")
    selected["canonical_gene"] = selected["gene"].map(canon_gene)
    selected = selected.drop_duplicates("canonical_gene").copy()

    packet = pd.read_csv(CANDIDATE_PACKET, sep="\t")
    packet["canonical_gene"] = packet["gene"].map(canon_gene)
    packet_meta = packet.set_index("canonical_gene").to_dict("index")

    rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        canonical = row["canonical_gene"]
        packet_row = packet_meta.get(canonical, {})
        gene = packet_row.get("gene", row["gene"])
        selection_reason = str(row.get("selection_reason", "selected_gene"))
        rows.append(
            {
                "gene": gene,
                "canonical_gene": canonical,
                "human_symbol": canonical,
                "mouse_symbol": canonical[:1] + canonical[1:].lower(),
                "panel": packet_row.get("panel", f"expanded_{selection_reason}"),
                "candidate_role": packet_row.get("candidate_role", "expanded_selected_gene"),
                "support_tier": packet_row.get("support_tier", selection_reason),
                "selected_gene_reason": selection_reason,
                "is_original_candidate_gene": canonical in packet_meta,
                "is_marker_panel_gene": bool(row.get("is_marker_panel_gene", False)),
                "max_feature_score": row.get("max_feature_score", np.nan),
            }
        )
    targets = pd.DataFrame(rows)
    metadata = targets.set_index("canonical_gene").to_dict("index")
    return targets, metadata


def install_expanded_targets() -> pd.DataFrame:
    targets, metadata = load_expanded_targets()
    base.TARGETS = targets
    base.TARGET_META = metadata
    base.TARGET_GENES = set(targets["canonical_gene"])
    return targets


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    if not valid.any():
        return out
    idx = p_values.index[valid].to_numpy()
    order = np.argsort(p_values.loc[idx].to_numpy(dtype=float))
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    m = len(p_sorted)
    adjusted = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def class_delta(sub: pd.DataFrame, target_class: str, background_classes: set[str]) -> dict[str, object]:
    target = sub.loc[sub["broad_class"].eq(target_class), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    background = sub.loc[sub["broad_class"].isin(background_classes), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    if len(target) == 0 or len(background) == 0:
        return {
            "n_target_units": len(target),
            "n_background_units": len(background),
            "target_median_rank": np.nan,
            "background_median_rank": np.nan,
            "rank_delta_vs_background": np.nan,
            "p_greater": np.nan,
        }
    p_value = stats.mannwhitneyu(target, background, alternative="greater").pvalue
    return {
        "n_target_units": len(target),
        "n_background_units": len(background),
        "target_median_rank": float(np.median(target)),
        "background_median_rank": float(np.median(background)),
        "rank_delta_vs_background": float(np.median(target) - np.median(background)),
        "p_greater": float(p_value),
    }


def compute_expanded_stats(expr: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    eligible = expr.loc[expr["eligible_class"]].copy()
    target_meta = targets.set_index("canonical_gene")
    dentate_background = {"non_dentate_background", "other_or_ambiguous", "broad_neuronal_structural_warning"}
    cerebellar_background = {"other_or_ambiguous", "broad_neuronal_structural_warning"}

    rows: list[dict[str, object]] = []
    for canonical_gene, sub in eligible.groupby("canonical_gene", sort=False):
        if canonical_gene not in target_meta.index:
            continue
        meta = target_meta.loc[canonical_gene].to_dict()
        dentate_sub = sub.loc[sub["core_branch"].isin(["mouse_dentate", "human_dentate_hippocampus"])]
        cerebellar_sub = sub.loc[sub["core_branch"].eq("cerebellum")]
        d = class_delta(dentate_sub, "dentate_candidate", dentate_background)
        c = class_delta(cerebellar_sub, "cerebellar_candidate", cerebellar_background)
        dentate_candidates = dentate_sub.loc[dentate_sub["broad_class"].eq("dentate_candidate")]
        cerebellar_candidates = cerebellar_sub.loc[cerebellar_sub["broad_class"].eq("cerebellar_candidate")]
        rows.append(
            {
                "gene": meta["gene"],
                "canonical_gene": canonical_gene,
                "panel": meta["panel"],
                "candidate_role": meta["candidate_role"],
                "support_tier": meta["support_tier"],
                "selected_gene_reason": meta["selected_gene_reason"],
                "is_original_candidate_gene": bool(meta["is_original_candidate_gene"]),
                "is_marker_panel_gene": bool(meta["is_marker_panel_gene"]),
                "dentate_candidate_units": d["n_target_units"],
                "dentate_background_units": d["n_background_units"],
                "dentate_candidate_median_rank": d["target_median_rank"],
                "dentate_background_median_rank": d["background_median_rank"],
                "dentate_rank_delta_vs_background": d["rank_delta_vs_background"],
                "dentate_rank_p_greater": d["p_greater"],
                "cerebellar_candidate_units": c["n_target_units"],
                "cerebellar_background_units": c["n_background_units"],
                "cerebellar_candidate_median_rank": c["target_median_rank"],
                "cerebellar_background_median_rank": c["background_median_rank"],
                "cerebellar_rank_delta_vs_background": c["rank_delta_vs_background"],
                "cerebellar_rank_p_greater": c["p_greater"],
                "dentate_candidate_median_detection": dentate_candidates["detection_fraction"].median(),
                "cerebellar_candidate_median_detection": cerebellar_candidates["detection_fraction"].median(),
                "n_dentate_datasets_detected_5pct": int(
                    dentate_candidates.loc[dentate_candidates["detection_fraction"].ge(MIN_DETECTION), "dataset"].nunique()
                ),
                "n_cerebellar_datasets_detected_5pct": int(
                    cerebellar_candidates.loc[cerebellar_candidates["detection_fraction"].ge(MIN_DETECTION), "dataset"].nunique()
                ),
            }
        )

    stats_df = pd.DataFrame(rows)
    stats_df["dentate_rank_p_adj_bh"] = bh_adjust(stats_df["dentate_rank_p_greater"])
    stats_df["cerebellar_rank_p_adj_bh"] = bh_adjust(stats_df["cerebellar_rank_p_greater"])
    stats_df["shared_positive_rank_delta"] = (
        stats_df["dentate_rank_delta_vs_background"].gt(0)
        & stats_df["cerebellar_rank_delta_vs_background"].gt(0)
        & stats_df["dentate_candidate_units"].ge(3)
        & stats_df["cerebellar_candidate_units"].ge(3)
    )
    stats_df["shared_strict_bh_0_10"] = (
        stats_df["shared_positive_rank_delta"]
        & stats_df["dentate_rank_p_adj_bh"].lt(0.10)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.10)
    )
    stats_df["shared_strict_bh_0_20"] = (
        stats_df["shared_positive_rank_delta"]
        & stats_df["dentate_rank_p_adj_bh"].lt(0.20)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.20)
    )
    stats_df["combined_rank_delta"] = stats_df["dentate_rank_delta_vs_background"].fillna(0) + stats_df[
        "cerebellar_rank_delta_vs_background"
    ].fillna(0)
    stats_df["minimum_branch_detection"] = stats_df[
        ["dentate_candidate_median_detection", "cerebellar_candidate_median_detection"]
    ].min(axis=1)

    hits = stats_df.loc[stats_df["shared_positive_rank_delta"]].sort_values(
        [
            "shared_strict_bh_0_10",
            "shared_strict_bh_0_20",
            "combined_rank_delta",
            "minimum_branch_detection",
        ],
        ascending=[False, False, False, False],
    )

    dentate_specific = stats_df.loc[
        stats_df["dentate_rank_delta_vs_background"].gt(0)
        & stats_df["dentate_rank_p_adj_bh"].lt(0.10)
        & (
            stats_df["cerebellar_rank_delta_vs_background"].le(0)
            | stats_df["cerebellar_rank_delta_vs_background"].isna()
            | stats_df["cerebellar_rank_p_adj_bh"].ge(0.20)
        )
    ].copy()
    dentate_specific["branch_specificity"] = "dentate_biased"

    cerebellar_specific = stats_df.loc[
        stats_df["cerebellar_rank_delta_vs_background"].gt(0)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.10)
        & (
            stats_df["dentate_rank_delta_vs_background"].le(0)
            | stats_df["dentate_rank_delta_vs_background"].isna()
            | stats_df["dentate_rank_p_adj_bh"].ge(0.20)
        )
    ].copy()
    cerebellar_specific["branch_specificity"] = "cerebellar_biased"

    branch = pd.concat([dentate_specific, cerebellar_specific], ignore_index=True, sort=False)
    branch = branch.sort_values(["branch_specificity", "combined_rank_delta"], ascending=[True, False])
    return stats_df, hits, branch


def plot_hits(hits: pd.DataFrame) -> None:
    plot_df = hits.head(30).copy()
    if plot_df.empty:
        return
    colors = np.where(plot_df["is_original_candidate_gene"], "#168b7a", "#4f7fbf")
    fig, ax = plt.subplots(figsize=(10.2, 7.2))
    y = np.arange(len(plot_df))
    ax.barh(y - 0.18, plot_df["dentate_rank_delta_vs_background"], height=0.34, color="#168b7a", label="dentate candidate")
    ax.barh(y + 0.18, plot_df["cerebellar_rank_delta_vs_background"], height=0.34, color="#6d3bbd", label="cerebellar candidate")
    for yi, color in zip(y, colors):
        ax.scatter(plot_df["combined_rank_delta"].iloc[yi] / 2, yi, s=16, color=color, zorder=3)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Median within-sample gene-rank delta versus local background")
    ax.set_title("Expanded selected-gene shared pseudobulk hits")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(targets: pd.DataFrame, expr: pd.DataFrame, coverage: pd.DataFrame, stats_df: pd.DataFrame, hits: pd.DataFrame, branch: pd.DataFrame) -> None:
    strict_010 = int(stats_df["shared_strict_bh_0_10"].sum())
    strict_020 = int(stats_df["shared_strict_bh_0_20"].sum())
    shared_positive = int(stats_df["shared_positive_rank_delta"].sum())
    original_recovered = int(hits["is_original_candidate_gene"].sum())
    new_hits = hits.loc[~hits["is_original_candidate_gene"]].head(20)
    original_hits = hits.loc[hits["is_original_candidate_gene"]].head(20)

    lines = [
        "# Expanded Primary-Core Gene Pseudobulk Screen",
        "",
        "Date built: 2026-06-22",
        "",
        "## Scope",
        "",
        "This analysis expands the 67-gene candidate pseudobulk screen to the 2,169-gene selected human-core universe. It is not a final whole-transcriptome DE model, but it is the first broad discovery layer that keeps all 10 primary datasets in the same target-gene frame.",
        "",
        f"Minimum broad-class size for rank statistics: {MIN_CLASS_CELLS} cells/nuclei.",
        "",
        "## Coverage",
        "",
        f"- Selected-gene universe: {targets['canonical_gene'].nunique():,} genes.",
        f"- Pseudobulk expression rows: {len(expr):,}.",
        f"- Primary datasets represented: {coverage['dataset'].nunique()}/10.",
        f"- Genes tested in statistics: {stats_df['canonical_gene'].nunique():,}.",
        f"- Shared-positive rank genes: {shared_positive:,}.",
        f"- Shared-positive genes passing BH<0.10 in both branches: {strict_010:,}.",
        f"- Shared-positive genes passing BH<0.20 in both branches: {strict_020:,}.",
        f"- Original 67-gene packet genes recovered among shared-positive hits: {original_recovered}.",
        "",
    ]
    for _, row in coverage.sort_values(["dataset", "sample"]).iterrows():
        lines.append(
            f"- `{row['dataset']}` / `{row['sample']}`: "
            f"{int(row['target_genes_present'])}/{int(row['target_genes_total'])} selected genes, "
            f"{int(row['n_labeled_observations'])}/{int(row['n_matrix_observations'])} labeled observations "
            f"(`{row['source_layer']}`)."
        )

    lines.extend(["", "## Original Candidate Recovery", ""])
    if original_hits.empty:
        lines.append("- No original candidate genes met the shared-positive criterion in this expanded pass.")
    else:
        for _, row in original_hits.iterrows():
            lines.append(
                f"- `{row['gene']}` ({row['candidate_role']}): dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
                f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}."
            )

    lines.extend(["", "## New Shared-Positive Candidates", ""])
    if new_hits.empty:
        lines.append("- No non-packet selected genes met the shared-positive criterion.")
    else:
        for _, row in new_hits.iterrows():
            lines.append(
                f"- `{row['gene']}` ({row['selected_gene_reason']}): dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
                f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}, "
                f"BH<0.10 both branches={bool(row['shared_strict_bh_0_10'])}."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This is a discovery screen over the selected-gene bridge universe, not a final genome-wide model.",
            "- Genes that are shared-positive here should be divided into broad neuronal state, structural/morphogenesis executors, metabolic/supporting genes, and regional identity leakage before entering a manuscript mechanism figure.",
            "- The most important use is prioritization: it identifies new candidate genes to inspect in the future whole-transcriptome ortholog-aware DE model and tests whether the 67-gene packet is recovered in a broader feature space.",
            "- Branch-specific hits are useful too: they help separate shared morphology executors from region-specific wiring or maturation programs.",
            "",
            "## Outputs",
            "",
            f"- Expression table: `{OUT_EXPR.relative_to(ROOT)}`",
            f"- Coverage table: `{OUT_COVERAGE.relative_to(ROOT)}`",
            f"- Statistics table: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Shared hits: `{OUT_HITS.relative_to(ROOT)}`",
            f"- Branch-specific hits: `{OUT_BRANCH.relative_to(ROOT)}`",
            f"- Shared-hit plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    targets = install_expanded_targets()
    expr, coverage = base.collect_expression()
    expr = base.add_within_sample_gene_ranks(expr)
    stats_df, hits, branch = compute_expanded_stats(expr, targets)
    plot_hits(hits)
    expr.to_csv(OUT_EXPR, sep="\t", index=False, compression="gzip")
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    stats_df.to_csv(OUT_STATS, sep="\t", index=False)
    hits.to_csv(OUT_HITS, sep="\t", index=False)
    branch.to_csv(OUT_BRANCH, sep="\t", index=False)
    write_report(targets, expr, coverage, stats_df, hits, branch)
    print(f"Wrote {len(expr):,} expanded pseudobulk expression rows")
    print(f"Wrote {len(stats_df):,} expanded gene statistics")
    print(f"Wrote {len(hits):,} shared-positive hits")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
