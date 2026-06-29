#!/usr/bin/env python3
"""Run strict 10-dataset primary-core module analysis.

This analysis is ortholog-aware at the marker-panel level: each refined panel
gene is represented by a canonical root symbol, human symbol, and mouse symbol.
It tests module-level identity separation and structural convergence on the
rank-unit layer, not yet genome-wide differential expression.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "Project/config"
RESULTS = ROOT / "Project/results"

CORE_DATASETS = RESULTS / "integrated_primary_core_datasets.tsv"
PANELS = CONFIG / "granule_marker_panels_refined.tsv"
INTEGRATED_UNITS = RESULTS / "human_bridge_backbone_rank_units.tsv"

OUT_ORTHOLOG = RESULTS / "primary_core_marker_panel_ortholog_map.tsv"
OUT_UNITS = RESULTS / "primary_core_integrated_module_units.tsv"
OUT_STATS = RESULTS / "primary_core_ortholog_module_statistics.tsv"
OUT_LOO = RESULTS / "primary_core_ortholog_module_leave_one_dataset_out.tsv"
OUT_PLOT = RESULTS / "primary_core_identity_structural_modules.png"
OUT_MD = RESULTS / "primary_core_ortholog_module_analysis.md"

MIN_UNIT_CELLS = 20

MODULES = [
    "dentate_identity",
    "cerebellar_identity",
    "shared_granule_neuronal",
    "morphogenesis_cytoskeleton",
    "axon_guidance_synapse",
    "structural_program_mean",
    "identity_rank_contrast",
    "dentate_rank",
    "cerebellar_rank",
    "shared_rank",
    "structural_rank",
]


def human_symbol(gene: str) -> str:
    return str(gene).strip().upper()


def mouse_symbol(gene: str) -> str:
    g = str(gene).strip()
    return g[:1].upper() + g[1:].lower() if g else g


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    idx = p_values.index[valid].to_numpy()
    if len(idx) == 0:
        return out
    order = p_values.loc[idx].to_numpy(dtype=float).argsort()
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    m = len(p_sorted)
    adjusted = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def load_core() -> pd.DataFrame:
    return pd.read_csv(CORE_DATASETS, sep="\t")


def build_ortholog_map(core: pd.DataFrame) -> pd.DataFrame:
    panels = pd.read_csv(PANELS, sep="\t")
    out = panels.copy()
    out["canonical_symbol"] = out["gene"].map(human_symbol)
    out["human_symbol"] = out["gene"].map(human_symbol)
    out["mouse_symbol"] = out["gene"].map(mouse_symbol)
    out["ortholog_scope"] = "refined_marker_panel"
    out["ortholog_basis"] = "same-root mouse-human marker symbol curated for panel-level module scoring"
    out["used_in_primary_core_species"] = ";".join(sorted(core["species"].str.split(";").explode().str.strip().unique()))
    return out[
        [
            "panel",
            "canonical_symbol",
            "human_symbol",
            "mouse_symbol",
            "gene",
            "notes",
            "ortholog_scope",
            "ortholog_basis",
            "used_in_primary_core_species",
        ]
    ].sort_values(["panel", "canonical_symbol"])


def load_primary_units(core: pd.DataFrame) -> pd.DataFrame:
    units = pd.read_csv(INTEGRATED_UNITS, sep="\t")
    units = units.loc[units["dataset"].isin(core["dataset"])].copy()
    units = units.merge(
        core[
            [
                "dataset",
                "core_branch",
                "core_role",
                "primary_reason",
                "current_use",
                "caveat",
            ]
        ],
        on="dataset",
        how="left",
        validate="many_to_one",
    )
    units["eligible_for_primary_stats"] = units["n_cells_or_spots"] >= MIN_UNIT_CELLS
    return units.sort_values(["core_branch", "dataset", "source_layer", "sample", "broad_class", "label"])


def values(units: pd.DataFrame, broad_class: str, metric: str, branch: str | None = None) -> np.ndarray:
    col = f"{metric}_median"
    mask = units["eligible_for_primary_stats"] & units["broad_class"].eq(broad_class)
    if branch is not None:
        mask &= units["core_branch"].eq(branch)
    return pd.to_numeric(units.loc[mask, col], errors="coerce").dropna().to_numpy(dtype=float)


def compare(
    units: pd.DataFrame,
    name: str,
    class_a: str,
    class_b: str,
    metric: str,
    note: str,
    branch_a: str | None = None,
    branch_b: str | None = None,
) -> dict[str, object]:
    a = values(units, class_a, metric, branch_a)
    b = values(units, class_b, metric, branch_b)
    p_value = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) and len(b) else np.nan
    return {
        "test": name,
        "metric": metric,
        "class_a": class_a,
        "class_b": class_b,
        "branch_a": branch_a or "all",
        "branch_b": branch_b or "all",
        "n_units_a": len(a),
        "n_units_b": len(b),
        "median_a": float(np.median(a)) if len(a) else np.nan,
        "median_b": float(np.median(b)) if len(b) else np.nan,
        "delta_a_minus_b": float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan,
        "p_value": p_value,
        "note": note,
    }


def sign_test(
    units: pd.DataFrame,
    name: str,
    broad_class: str,
    metric: str,
    threshold: float,
    direction: str,
    note: str,
    branch: str | None = None,
) -> dict[str, object]:
    vals = values(units, broad_class, metric, branch)
    if direction == "greater":
        successes = int(np.sum(vals > threshold))
    elif direction == "less":
        successes = int(np.sum(vals < threshold))
    else:
        raise ValueError(direction)
    p_value = stats.binomtest(successes, len(vals), 0.5, alternative="greater").pvalue if len(vals) else np.nan
    return {
        "test": name,
        "metric": metric,
        "class_a": broad_class,
        "class_b": f"{direction}_{threshold}",
        "branch_a": branch or "all",
        "branch_b": "threshold",
        "n_units_a": len(vals),
        "n_units_b": "",
        "median_a": float(np.median(vals)) if len(vals) else np.nan,
        "median_b": threshold,
        "delta_a_minus_b": float(np.median(vals) - threshold) if len(vals) else np.nan,
        "p_value": p_value,
        "note": f"{note}; successes={successes}/{len(vals)}",
    }


def compute_stats(units: pd.DataFrame) -> pd.DataFrame:
    rows = [
        compare(
            units,
            "primary_core_dentate_vs_cerebellar_identity_separation",
            "dentate_candidate",
            "cerebellar_candidate",
            "identity_rank_contrast",
            "Primary test: dentate candidates should be dentate-rank high and cerebellar candidates cerebellar-rank high.",
        ),
        compare(
            units,
            "primary_core_dentate_vs_background_identity_separation",
            "dentate_candidate",
            "non_dentate_background",
            "identity_rank_contrast",
            "Dentate candidates should separate from non-dentate/background units.",
        ),
        compare(
            units,
            "primary_core_dentate_vs_cerebellar_structural_convergence",
            "dentate_candidate",
            "cerebellar_candidate",
            "structural_rank",
            "Structural rank tests magnitude difference after establishing that both groups are above the within-sample structural median.",
        ),
        sign_test(
            units,
            "primary_core_dentate_identity_contrast_above_zero",
            "dentate_candidate",
            "identity_rank_contrast",
            0.0,
            "greater",
            "Dentate candidates should have positive dentate-minus-cerebellar rank.",
        ),
        sign_test(
            units,
            "primary_core_cerebellar_identity_contrast_below_zero",
            "cerebellar_candidate",
            "identity_rank_contrast",
            0.0,
            "less",
            "Cerebellar candidates should have negative dentate-minus-cerebellar rank.",
        ),
        sign_test(
            units,
            "primary_core_dentate_structural_rank_above_median",
            "dentate_candidate",
            "structural_rank",
            0.5,
            "greater",
            "Dentate candidates should be above the within-sample structural median.",
        ),
        sign_test(
            units,
            "primary_core_cerebellar_structural_rank_above_median",
            "cerebellar_candidate",
            "structural_rank",
            0.5,
            "greater",
            "Cerebellar candidates should be above the within-sample structural median.",
        ),
        compare(
            units,
            "mouse_dentate_vs_human_dentate_structural_rank",
            "dentate_candidate",
            "dentate_candidate",
            "structural_rank",
            "Checks whether human dentate bridge and mouse dentate units differ strongly on structural rank.",
            branch_a="mouse_dentate",
            branch_b="human_dentate_hippocampus",
        ),
        compare(
            units,
            "mouse_dentate_vs_human_dentate_identity_rank_contrast",
            "dentate_candidate",
            "dentate_candidate",
            "identity_rank_contrast",
            "DG-enriched human sources can compress identity contrast; interpret source-aware.",
            branch_a="mouse_dentate",
            branch_b="human_dentate_hippocampus",
        ),
        sign_test(
            units,
            "human_dentate_branch_structural_rank_above_median",
            "dentate_candidate",
            "structural_rank",
            0.5,
            "greater",
            "Human dentate/hippocampal candidate units should retain structural-program support.",
            branch="human_dentate_hippocampus",
        ),
    ]
    out = pd.DataFrame(rows)
    out["p_adj_bh"] = bh_adjust(out["p_value"])
    out["min_unit_cells"] = MIN_UNIT_CELLS
    return out.sort_values(["p_adj_bh", "test"]).reset_index(drop=True)


def compute_leave_one_dataset_out(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    datasets = sorted(units["dataset"].dropna().unique())
    for held_out in ["none"] + datasets:
        sub = units.copy() if held_out == "none" else units.loc[~units["dataset"].eq(held_out)].copy()
        stats_df = compute_stats(sub)
        keep = stats_df.loc[
            stats_df["test"].isin(
                [
                    "primary_core_dentate_vs_cerebellar_identity_separation",
                    "primary_core_dentate_structural_rank_above_median",
                    "primary_core_cerebellar_structural_rank_above_median",
                    "human_dentate_branch_structural_rank_above_median",
                ]
            )
        ].copy()
        keep.insert(0, "held_out_dataset", held_out)
        keep.insert(1, "n_eligible_units", int(sub["eligible_for_primary_stats"].sum()))
        rows.append(keep)
    return pd.concat(rows, ignore_index=True, sort=False)


def plot_units(units: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_df = units.loc[units["eligible_for_primary_stats"]].copy()
    plot_df = plot_df.loc[
        plot_df["broad_class"].isin(
            [
                "dentate_candidate",
                "cerebellar_candidate",
                "non_dentate_background",
                "broad_neuronal_structural_warning",
                "other_or_ambiguous",
                "dentate_low_support",
            ]
        )
    ]
    colors = {
        "dentate_candidate": "#168b7a",
        "cerebellar_candidate": "#6d3bbd",
        "non_dentate_background": "#7f8790",
        "broad_neuronal_structural_warning": "#c97938",
        "other_or_ambiguous": "#c9c9c9",
        "dentate_low_support": "#71b7ad",
    }
    markers = {
        "mouse_dentate": "o",
        "cerebellum": "^",
        "human_dentate_hippocampus": "s",
    }
    fig, ax = plt.subplots(figsize=(9.0, 6.2))
    for (branch, broad), sub in plot_df.groupby(["core_branch", "broad_class"], dropna=False):
        ax.scatter(
            sub["identity_rank_contrast_median"],
            sub["structural_rank_median"],
            s=18 + np.clip(np.sqrt(sub["n_cells_or_spots"].astype(float)), 2, 72),
            c=colors.get(broad, "#999999"),
            marker=markers.get(branch, "o"),
            alpha=0.78,
            edgecolors="white",
            linewidths=0.5,
            label=f"{branch} | {broad}",
        )
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.axhline(0.5, color="#777777", linewidth=0.9, linestyle=":")
    ax.set_xlabel("Identity rank contrast (dentate - cerebellar)")
    ax.set_ylabel("Structural-program rank")
    ax.set_title("Strict 10-dataset primary core module units")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), frameon=False, fontsize=6.6, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.grid(True, linewidth=0.4, color="#d9d9d9")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 0.74, 1])
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_md(core: pd.DataFrame, units: pd.DataFrame, stats_df: pd.DataFrame, loo: pd.DataFrame) -> None:
    eligible = units.loc[units["eligible_for_primary_stats"]].copy()
    counts = (
        eligible.groupby(["core_branch", "broad_class"], dropna=False)
        .agg(n_units=("unit_id", "size"), n_cells_or_spots=("n_cells_or_spots", "sum"))
        .reset_index()
        .sort_values(["core_branch", "broad_class"])
    )
    top_tests = stats_df.loc[
        stats_df["test"].isin(
            [
                "primary_core_dentate_vs_cerebellar_identity_separation",
                "primary_core_dentate_vs_cerebellar_structural_convergence",
                "primary_core_dentate_structural_rank_above_median",
                "primary_core_cerebellar_structural_rank_above_median",
                "human_dentate_branch_structural_rank_above_median",
            ]
        )
    ]
    lines = [
        "# Primary Core Ortholog-Aware Module Analysis",
        "",
        "Date built: 2026-06-22",
        "",
        "## Scope",
        "",
        "This analysis freezes the strict 10-dataset primary core and tests the current module-level model: regional identity separation with structural-program convergence.",
        "",
        "The ortholog-aware layer is currently marker-panel level, using curated same-root mouse-human symbols from the refined marker panel. It is not yet genome-wide ortholog-aware differential expression.",
        "",
        "## Primary Core Composition",
        "",
    ]
    branch_counts = core.groupby("core_branch")["dataset"].nunique().reset_index()
    for _, row in branch_counts.iterrows():
        lines.append(f"- `{row['core_branch']}`: {int(row['dataset'])} datasets.")
    lines.extend(["", "## Eligible Module Units", ""])
    for _, row in counts.iterrows():
        lines.append(
            f"- `{row['core_branch']}` / `{row['broad_class']}`: "
            f"{int(row['n_units'])} units, {int(row['n_cells_or_spots'])} cells/spots."
        )
    lines.extend(["", "## Main Tests", ""])
    for _, row in top_tests.sort_values("test").iterrows():
        lines.append(
            f"- `{row['test']}` / `{row['metric']}`: median delta {row['delta_a_minus_b']:.4f}; "
            f"n={row['n_units_a']} vs {row['n_units_b']}; BH-adjusted p={row['p_adj_bh']:.3g}."
        )
    identity_loo = loo.loc[loo["test"].eq("primary_core_dentate_vs_cerebellar_identity_separation")].copy()
    stable = identity_loo["p_adj_bh"].lt(0.05).all() and (identity_loo["delta_a_minus_b"] > 0).all()
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Leave-one-dataset-out identity separation is {'stable' if stable else 'not uniformly stable'} at BH-adjusted p<0.05 with positive dentate-minus-cerebellar delta.",
            "- The strongest model remains: dentate and cerebellar granule candidates are identity-distinct, while structural/morphogenesis modules provide a shared elevated executor axis.",
            "- The structural axis should not be framed as equal magnitude: in this strict-core rank layer, cerebellar candidates are higher than dentate candidates, while both candidate classes remain above the within-sample structural median.",
            "- Human DG-enriched sources are most useful for anchoring human dentate state and structural programs; they can compress within-sample dentate-versus-cerebellar contrast.",
            "- The next stricter step is genome-wide ortholog-aware pseudobulk or mixed-effect differential expression within the 10-dataset core.",
            "",
            "## Outputs",
            "",
            f"- Ortholog marker map: `{OUT_ORTHOLOG.relative_to(ROOT)}`",
            f"- Primary core module units: `{OUT_UNITS.relative_to(ROOT)}`",
            f"- Module statistics: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Leave-one-dataset-out checks: `{OUT_LOO.relative_to(ROOT)}`",
            f"- Primary core module plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))
    core = load_core()
    ortholog = build_ortholog_map(core)
    units = load_primary_units(core)
    stats_df = compute_stats(units)
    loo = compute_leave_one_dataset_out(units)

    ortholog.to_csv(OUT_ORTHOLOG, sep="\t", index=False)
    units.to_csv(OUT_UNITS, sep="\t", index=False, float_format="%.6g")
    stats_df.to_csv(OUT_STATS, sep="\t", index=False, float_format="%.6g")
    loo.to_csv(OUT_LOO, sep="\t", index=False, float_format="%.6g")
    plot_units(units)
    write_md(core, units, stats_df, loo)

    print(f"Wrote {OUT_ORTHOLOG}")
    print(f"Wrote {OUT_UNITS}")
    print(f"Wrote {OUT_STATS}")
    print(f"Wrote {OUT_LOO}")
    print(f"Wrote {OUT_PLOT}")
    print(f"Wrote {OUT_MD}")
    print(f"datasets={core['dataset'].nunique()}; units={len(units)}; eligible={int(units['eligible_for_primary_stats'].sum())}")


if __name__ == "__main__":
    main()
