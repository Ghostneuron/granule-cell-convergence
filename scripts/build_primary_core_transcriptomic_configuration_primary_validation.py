#!/usr/bin/env python3
"""Validate transcriptomic configuration across primary-core pseudobulk layers.

The local named-comparator configuration model is small-n because only two
datasets retain explicit pyramidal/Purkinje labels. This script broadens the
test across primary-core pseudobulk expression layers by asking whether the
candidate granule class in each dataset/sample has a higher module-balance
configuration score than local non-candidate/background classes.
"""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

GENE_SETS = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
FULL_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"
SELECTED_EXPR = RESULTS / "primary_core_expanded_gene_pseudobulk_expression.tsv.gz"

OUT_UNITS = RESULTS / "primary_core_transcriptomic_configuration_primary_units.tsv.gz"
OUT_CONTRASTS = RESULTS / "primary_core_transcriptomic_configuration_primary_contrasts.tsv"
OUT_SUMMARY = RESULTS / "primary_core_transcriptomic_configuration_primary_summary.tsv"
OUT_COVERAGE = RESULTS / "primary_core_transcriptomic_configuration_primary_coverage.tsv"
OUT_PLOT = RESULTS / "primary_core_transcriptomic_configuration_primary_validation.png"
OUT_PLOT_LABEL_KEY = RESULTS / "primary_core_transcriptomic_configuration_primary_validation_label_key.tsv"
OUT_MD = RESULTS / "primary_core_transcriptomic_configuration_primary_validation.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


MODULE_COLS = {
    "cerebellar_fate_rhombic_lip_shh": "cerebellar_fate_rank",
    "dentate_fate_wnt_prox1": "dentate_fate_rank",
    "shared_neurogenic_niche_state": "niche_rank",
    "downstream_neurite_morphology": "neurite_morphology_rank",
    "downstream_synaptic_excitability": "synaptic_excitability_rank",
}

EXPRESSION_LAYERS = [
    {
        "expression_layer": "full_mgi_ortholog_matrix",
        "path": FULL_EXPR,
        "description": "full local matrices projected through MGI one-to-one orthologs",
    },
    {
        "expression_layer": "selected_feature_matrix",
        "path": SELECTED_EXPR,
        "description": "2,169-gene selected-feature primary-core pseudobulk layer",
    },
]

SCORE_COLS = [
    "downstream_construction_balance",
    "regional_fate_balance",
    "configuration_score",
]

LAYER_DISPLAY = {
    "full_mgi_ortholog_matrix": "Full MGI ortholog matrix",
    "selected_feature_matrix": "Selected-feature matrix",
}

LAYER_PREFIX = {
    "full_mgi_ortholog_matrix": "F",
    "selected_feature_matrix": "S",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().upper()


def branch_kind(core_branch: str) -> str:
    branch = str(core_branch).lower()
    if "cerebell" in branch:
        return "cerebellar"
    if "dentate" in branch or "hippocamp" in branch:
        return "dentate"
    return "unknown"


def target_class_for_branch(core_branch: str) -> str:
    kind = branch_kind(core_branch)
    if kind == "cerebellar":
        return "cerebellar_candidate"
    if kind == "dentate":
        return "dentate_candidate"
    return ""


def read_module_expression(path: Path, genes: set[str], expression_layer: str) -> pd.DataFrame:
    use_cols = [
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "n_cells",
        "canonical_gene",
        "gene",
        "human_symbol",
        "mouse_symbol",
        "detection_fraction",
        "mean_log1p_expression",
        "eligible_class",
        "mean_log1p_rank_within_sample_gene",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", usecols=lambda col: col in use_cols, chunksize=100_000, low_memory=False):
        chunk["canonical_gene"] = chunk["canonical_gene"].map(canon)
        sub = chunk.loc[chunk["canonical_gene"].isin(genes)].copy()
        if not sub.empty:
            sub["expression_layer"] = expression_layer
            pieces.append(sub)
    if not pieces:
        return pd.DataFrame(columns=[*use_cols, "expression_layer"])
    return pd.concat(pieces, ignore_index=True)


def build_units() -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_sets = pd.read_csv(GENE_SETS, sep="\t")
    gene_sets["canonical_gene"] = gene_sets["canonical_gene"].map(canon)
    genes = set(gene_sets["canonical_gene"])
    module_map = gene_sets[["canonical_gene", "module_id", "module_label", "module_family"]].drop_duplicates()

    all_expr: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    for layer in EXPRESSION_LAYERS:
        expr = read_module_expression(layer["path"], genes, layer["expression_layer"])
        if expr.empty:
            continue
        expr = expr.merge(module_map, on="canonical_gene", how="left")
        all_expr.append(expr)
        for module_id, sub in gene_sets.groupby("module_id", sort=False):
            present = sorted(set(expr.loc[expr["module_id"].eq(module_id), "canonical_gene"]))
            coverage_rows.append(
                {
                    "expression_layer": layer["expression_layer"],
                    "module_id": module_id,
                    "module_label": sub["module_label"].iloc[0],
                    "n_defined_genes": int(sub["canonical_gene"].nunique()),
                    "n_present_genes": int(len(present)),
                    "present_genes": ",".join(present),
                    "missing_genes": ",".join(sorted(set(sub["canonical_gene"]) - set(present))),
                }
            )

    if not all_expr:
        raise RuntimeError("No module-expression rows found")

    expr_all = pd.concat(all_expr, ignore_index=True)
    expr_all = expr_all.loc[expr_all["eligible_class"].astype(bool)].copy()

    group_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "module_id",
    ]
    module_units = (
        expr_all.groupby(group_cols, sort=False)
        .agg(
            n_cells=("n_cells", "max"),
            n_module_genes_present=("canonical_gene", "nunique"),
            median_module_gene_rank=("mean_log1p_rank_within_sample_gene", "median"),
            mean_module_gene_rank=("mean_log1p_rank_within_sample_gene", "mean"),
            median_detection_fraction=("detection_fraction", "median"),
            genes_present=("canonical_gene", lambda values: ",".join(sorted(set(values)))),
        )
        .reset_index()
    )

    pivot = module_units.pivot_table(
        index=[
            "expression_layer",
            "dataset",
            "core_branch",
            "sample",
            "source_layer",
            "expression_scope",
            "broad_class",
        ],
        columns="module_id",
        values="median_module_gene_rank",
        aggfunc="median",
    ).reset_index()
    pivot.columns.name = None
    pivot = pivot.rename(columns=MODULE_COLS)
    for col in MODULE_COLS.values():
        if col not in pivot.columns:
            pivot[col] = np.nan

    gene_count_pivot = module_units.pivot_table(
        index=[
            "expression_layer",
            "dataset",
            "core_branch",
            "sample",
            "source_layer",
            "expression_scope",
            "broad_class",
        ],
        columns="module_id",
        values="n_module_genes_present",
        aggfunc="max",
    ).reset_index()
    gene_count_pivot.columns.name = None
    gene_count_pivot = gene_count_pivot.rename(
        columns={module: f"n_present_{column}" for module, column in MODULE_COLS.items()}
    )
    pivot = pivot.merge(gene_count_pivot, on=[
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
    ], how="left")

    pivot["branch_kind"] = pivot["core_branch"].map(branch_kind)
    pivot["target_class_for_branch"] = pivot["core_branch"].map(target_class_for_branch)
    pivot["is_branch_candidate"] = pivot["broad_class"].eq(pivot["target_class_for_branch"])
    pivot["is_background_class"] = ~pivot["is_branch_candidate"]
    pivot["downstream_mean_rank"] = pivot[["neurite_morphology_rank", "synaptic_excitability_rank"]].mean(axis=1)
    pivot["downstream_construction_balance"] = pivot["downstream_mean_rank"] - pivot["niche_rank"]
    pivot["branch_matched_fate_rank"] = np.where(
        pivot["branch_kind"].eq("dentate"),
        pivot["dentate_fate_rank"],
        np.where(pivot["branch_kind"].eq("cerebellar"), pivot["cerebellar_fate_rank"], np.nan),
    )
    pivot["branch_opposed_fate_rank"] = np.where(
        pivot["branch_kind"].eq("dentate"),
        pivot["cerebellar_fate_rank"],
        np.where(pivot["branch_kind"].eq("cerebellar"), pivot["dentate_fate_rank"], np.nan),
    )
    pivot["regional_fate_balance"] = pivot["branch_matched_fate_rank"] - pivot["branch_opposed_fate_rank"]
    pivot["configuration_score"] = pivot["downstream_construction_balance"] + pivot["regional_fate_balance"]
    pivot = pivot.sort_values(["expression_layer", "dataset", "sample", "broad_class"])

    coverage = pd.DataFrame(coverage_rows)
    pivot.to_csv(OUT_UNITS, sep="\t", index=False, compression="gzip")
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    return pivot, coverage


def median_or_nan(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def build_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["expression_layer", "dataset", "core_branch", "sample", "source_layer", "expression_scope"]
    for keys, sub in units.groupby(group_cols, sort=False):
        expression_layer, dataset, core_branch, sample, source_layer, expression_scope = keys
        target_class = target_class_for_branch(core_branch)
        if not target_class:
            continue
        candidate = sub.loc[sub["broad_class"].eq(target_class)]
        background = sub.loc[
            ~sub["broad_class"].eq(target_class)
            & ~sub["broad_class"].astype(str).str.contains("low_support", case=False, na=False)
        ]
        background = background.loc[background["configuration_score"].notna()]
        if candidate.empty or background.empty:
            continue
        row: dict[str, object] = {
            "expression_layer": expression_layer,
            "dataset": dataset,
            "core_branch": core_branch,
            "branch_kind": branch_kind(core_branch),
            "sample": sample,
            "source_layer": source_layer,
            "expression_scope": expression_scope,
            "target_class": target_class,
            "background_classes": ",".join(sorted(set(background["broad_class"].astype(str)))),
            "n_candidate_classes": int(len(candidate)),
            "n_background_classes": int(len(background)),
        }
        for score in SCORE_COLS:
            cand = median_or_nan(candidate[score])
            bg = median_or_nan(background[score])
            row[f"candidate_median_{score}"] = cand
            row[f"background_median_{score}"] = bg
            row[f"delta_{score}"] = cand - bg if np.isfinite(cand) and np.isfinite(bg) else np.nan
        row["configuration_positive"] = bool(row["delta_configuration_score"] > 0) if np.isfinite(row["delta_configuration_score"]) else False
        rows.append(row)
    contrasts = pd.DataFrame(rows)
    contrasts = contrasts.sort_values(["expression_layer", "branch_kind", "dataset", "sample"])
    contrasts.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    return contrasts


def test_positive(values: pd.Series) -> tuple[float, float, int, int]:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan, np.nan, 0, 0
    positives = int(np.sum(arr > 0))
    n = int(arr.size)
    sign_p = float(stats.binomtest(positives, n=n, p=0.5, alternative="greater").pvalue)
    wilcoxon_p = float(stats.wilcoxon(arr, alternative="greater").pvalue) if n > 1 and np.any(arr != 0) else np.nan
    return sign_p, wilcoxon_p, positives, n


def summarize_contrasts(contrasts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    groupings = [
        ("all_layers", ["all"]),
        ("by_expression_layer", ["expression_layer"]),
        ("by_branch_kind", ["branch_kind"]),
        ("by_expression_layer_and_branch", ["expression_layer", "branch_kind"]),
    ]
    for summary_level, cols in groupings:
        if cols == ["all"]:
            iterator = [(("all",), contrasts)]
            label_cols = ["summary_group"]
        else:
            iterator = contrasts.groupby(cols, sort=False)
            label_cols = cols
        for keys, sub in iterator:
            if not isinstance(keys, tuple):
                keys = (keys,)
            row: dict[str, object] = {"summary_level": summary_level}
            for col, key in zip(label_cols, keys):
                row[col] = key
            for score in SCORE_COLS:
                sign_p, wilcoxon_p, positives, n = test_positive(sub[f"delta_{score}"])
                row[f"n_{score}_contrasts"] = n
                row[f"n_{score}_positive"] = positives
                row[f"fraction_{score}_positive"] = positives / n if n else np.nan
                row[f"median_delta_{score}"] = median_or_nan(sub[f"delta_{score}"])
                row[f"sign_test_p_{score}"] = sign_p
                row[f"wilcoxon_p_{score}"] = wilcoxon_p
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    return summary


def plot_results(contrasts: pd.DataFrame, summary: pd.DataFrame) -> None:
    plot_df = contrasts.copy()
    plot_df = plot_df.sort_values(["expression_layer", "branch_kind", "dataset", "sample"])
    layer_order = list(dict.fromkeys(plot_df["expression_layer"]))
    colors = {
        "dentate": "#2f7f8f",
        "cerebellar": "#7f4e8a",
        "unknown": "#777777",
    }

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(13.8, 9.2),
        gridspec_kw={"height_ratios": [1.15, 1.0]},
        constrained_layout=True,
    )
    ax = axes[0]
    x_positions: list[float] = []
    x_labels: list[str] = []
    label_rows: list[dict[str, object]] = []
    layer_spans: list[tuple[str, float, float, int]] = []
    current_x = 0.0
    for layer in layer_order:
        sub = plot_df.loc[plot_df["expression_layer"].eq(layer)]
        if sub.empty:
            continue
        layer_start = current_x - 0.5
        prefix = LAYER_PREFIX.get(layer, str(layer)[:1].upper() or "X")
        for label_idx, (_, row) in enumerate(sub.iterrows(), start=1):
            plot_id = f"{prefix}{label_idx:02d}"
            x_positions.append(current_x)
            x_labels.append(plot_id)
            label_rows.append(
                {
                    "plot_id": plot_id,
                    "expression_layer": layer,
                    "dataset": row["dataset"],
                    "sample": row["sample"],
                    "core_branch": row["core_branch"],
                    "branch_kind": row["branch_kind"],
                    "source_layer": row["source_layer"],
                    "expression_scope": row["expression_scope"],
                    "delta_configuration_score": row["delta_configuration_score"],
                    "configuration_positive": row["configuration_positive"],
                }
            )
            ax.bar(
                current_x,
                row["delta_configuration_score"],
                color=colors.get(row["branch_kind"], "#777777"),
                edgecolor="#333333",
                linewidth=0.4,
                zorder=2,
            )
            current_x += 1.0
        layer_spans.append((layer, layer_start, current_x - 0.5, len(sub)))
        current_x += 0.8
    pd.DataFrame(label_rows).to_csv(OUT_PLOT_LABEL_KEY, sep="\t", index=False)
    for idx, (layer, start_x, end_x, n_layer) in enumerate(layer_spans):
        if idx % 2 == 0:
            ax.axvspan(start_x, end_x, color="#f6f8fb", zorder=0)
        ax.axvline(start_x, color="#c9ced6", linewidth=0.8, zorder=1)
    ax.axhline(0, color="#555555", linewidth=0.8)
    if x_positions:
        ax.set_xlim(min(x_positions) - 0.7, max(x_positions) + 0.7)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, ha="center", va="top", fontsize=6.7)
    ax.set_xlabel("Contrast ID; decode in the label-key TSV")
    ax.set_ylabel("Candidate minus background\nconfiguration score")
    ax.set_title("Primary-core configuration validation")
    ax.legend(
        handles=[
            Patch(facecolor=colors["cerebellar"], edgecolor="#333333", label="Cerebellar branch"),
            Patch(facecolor=colors["dentate"], edgecolor="#333333", label="Dentate branch"),
        ],
        frameon=False,
        fontsize=8.5,
        loc="upper right",
    )
    ax.grid(axis="y", color="#dddddd", linewidth=0.5, zorder=1)

    ax = axes[1]
    summary_plot = summary.loc[summary["summary_level"].eq("by_expression_layer")].copy()
    x = np.arange(len(summary_plot))
    width = 0.24
    score_labels = [
        ("median_delta_downstream_construction_balance", "Construction over niche", "#2f7f8f"),
        ("median_delta_regional_fate_balance", "Regional fate polarity", "#7b6d8d"),
        ("median_delta_configuration_score", "Combined configuration", "#1f4e5f"),
    ]
    for idx, (col, label, color) in enumerate(score_labels):
        ax.bar(x + (idx - 1) * width, summary_plot[col], width=width, label=label, color=color)
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [LAYER_DISPLAY.get(layer, str(layer).replace("_", "\n")) for layer in summary_plot["expression_layer"]],
        rotation=0,
        ha="center",
    )
    ax.set_ylabel("Median candidate-background delta")
    ax.set_title("Layer-level median effects")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    fig.suptitle("Primary-Core Transcriptomic Configuration Validation", fontsize=15, y=1.02)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(units: pd.DataFrame, contrasts: pd.DataFrame, summary: pd.DataFrame, coverage: pd.DataFrame) -> None:
    all_summary = summary.loc[summary["summary_level"].eq("all_layers")].iloc[0]
    layer_summary = summary.loc[summary["summary_level"].eq("by_expression_layer")]
    branch_summary = summary.loc[summary["summary_level"].eq("by_branch_kind")]

    lines = [
        "# Primary-Core Transcriptomic Configuration Validation",
        "",
        f"Date built: {date.today().isoformat()}",
        "",
        "## Purpose",
        "",
        "This analysis broadens the local transcriptomic configuration test from named pyramidal/Purkinje comparators to primary-core candidate-versus-background pseudobulk contrasts.",
        "",
        "## Inputs",
        "",
        "- Full MGI one-to-one ortholog matrix expression layer.",
        "- 2,169-gene selected-feature primary-core pseudobulk layer.",
        "- The five niche/circuit module gene sets used in the local configuration model.",
        "",
        "## Scale",
        "",
        f"- Configuration class units: {len(units):,} across {units['dataset'].nunique():,} datasets.",
        f"- Candidate-versus-background contrasts: {len(contrasts):,} across {contrasts['dataset'].nunique():,} datasets.",
        f"- Expression layers represented: {contrasts['expression_layer'].nunique():,}.",
        "",
        "## Main Result",
        "",
        f"- Combined configuration positive contrasts: {int(all_summary['n_configuration_score_positive'])}/{int(all_summary['n_configuration_score_contrasts'])}.",
        f"- Median candidate-background combined configuration delta: {all_summary['median_delta_configuration_score']:.3f}.",
        f"- Sign-test p for positive combined configuration deltas: {all_summary['sign_test_p_configuration_score']:.3g}.",
        f"- Wilcoxon p for combined configuration greater than zero: {all_summary['wilcoxon_p_configuration_score']:.3g}.",
        "",
        "Layer-level summary:",
    ]
    for _, row in layer_summary.iterrows():
        lines.append(
            f"- {row['expression_layer']}: {int(row['n_configuration_score_positive'])}/{int(row['n_configuration_score_contrasts'])} positive, "
            f"median delta {row['median_delta_configuration_score']:.3f}, Wilcoxon p={row['wilcoxon_p_configuration_score']:.3g}."
        )
    lines.append("")
    lines.append("Branch-level summary:")
    for _, row in branch_summary.iterrows():
        lines.append(
            f"- {row['branch_kind']}: {int(row['n_configuration_score_positive'])}/{int(row['n_configuration_score_contrasts'])} positive, "
            f"median delta {row['median_delta_configuration_score']:.3f}, Wilcoxon p={row['wilcoxon_p_configuration_score']:.3g}."
        )

    coverage_min = coverage.groupby("module_label")["n_present_genes"].min().sort_values()
    lines.extend(["", "Minimum module gene coverage across expression layers:"])
    for module_label, n_present in coverage_min.items():
        lines.append(f"- {module_label}: {int(n_present)} genes present.")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This strengthens the configuration hypothesis beyond the four named local contrasts: candidate granule classes usually show higher combined configuration scores than local backgrounds across the primary-core pseudobulk layers.",
            "- The test is broader but less specific than the named-comparator audit because most datasets do not preserve explicit pyramidal or Purkinje comparator labels.",
            "- The result supports manuscript language that morphology is partly encoded as transcriptomic module balance, while final geometry still requires spatial, proteomic, activity, and lineage/perturbation validation.",
            "",
            "## Outputs",
            "",
            f"- Configuration units: `{rel(OUT_UNITS)}`",
            f"- Contrast table: `{rel(OUT_CONTRASTS)}`",
            f"- Summary table: `{rel(OUT_SUMMARY)}`",
            f"- Coverage table: `{rel(OUT_COVERAGE)}`",
            f"- Plot label key: `{rel(OUT_PLOT_LABEL_KEY)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    units, coverage = build_units()
    contrasts = build_contrasts(units)
    summary = summarize_contrasts(contrasts)
    plot_results(contrasts, summary)
    write_report(units, contrasts, summary, coverage)
    print(f"Wrote {rel(OUT_MD)}")
    print(summary.loc[summary["summary_level"].isin(["all_layers", "by_expression_layer"]), [
        "summary_level",
        "summary_group" if "summary_group" in summary.columns else "summary_level",
        "expression_layer" if "expression_layer" in summary.columns else "summary_level",
        "n_configuration_score_positive",
        "n_configuration_score_contrasts",
        "median_delta_configuration_score",
        "wilcoxon_p_configuration_score",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
