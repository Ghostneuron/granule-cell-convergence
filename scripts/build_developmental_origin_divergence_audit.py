#!/usr/bin/env python3
"""Audit deep neural origin, regional divergence, and later granule convergence.

This is a deliberately conservative support analysis. It does not attempt to
prove clonal lineage history. Instead, it asks whether the locally available
dentate and cerebellar full matrices are consistent with a model in which the
two granule-cell lineages share a deep neural-progenitor origin, diverge into
distinct anterior/telencephalic versus hindbrain/rhombic-lip territories, and
later reuse shared postmitotic/construction programs.
"""

from __future__ import annotations

import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

import build_primary_core_granule_specificity_named_comparators as base


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
SUPP_FIGURES = ROOT / "Project/manuscript/Supplementary figures"

OUT_GENE_SETS = RESULTS / "developmental_origin_divergence_audit_gene_sets.tsv"
OUT_UNITS = RESULTS / "developmental_origin_divergence_audit_units.tsv"
OUT_STATE_SUMMARY = RESULTS / "developmental_origin_divergence_audit_state_summary.tsv"
OUT_METRICS = RESULTS / "developmental_origin_divergence_audit_branch_metrics.tsv"
OUT_PLOT = RESULTS / "developmental_origin_divergence_audit.png"
OUT_SUPP_PLOT = SUPP_FIGURES / "Fig.S4_developmental_origin_divergence_audit.png"
OUT_MD = RESULTS / "developmental_origin_divergence_audit.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


MODULES = [
    {
        "module_id": "deep_neural_progenitor",
        "module_label": "Deep neural progenitor",
        "module_family": "deep_origin",
        "interpretation_role": "shared neural ectoderm/neuroepithelium competence, not a direct clonal marker",
        "genes": [
            "SOX1",
            "SOX2",
            "SOX3",
            "NES",
            "VIM",
            "HES1",
            "HES5",
            "PAX6",
            "HMGA2",
            "MSI1",
            "PROM1",
            "NOTCH1",
            "ASCL1",
        ],
    },
    {
        "module_id": "anterior_telencephalic_pattern",
        "module_label": "Anterior/telencephalic pattern",
        "module_family": "regional_origin",
        "interpretation_role": "forebrain/telencephalic regional identity upstream of hippocampal/dentate lineage",
        "genes": [
            "FOXG1",
            "EMX1",
            "EMX2",
            "LHX2",
            "SIX3",
            "OTX2",
            "PAX6",
            "EOMES",
            "TBR1",
            "ZBTB20",
            "TCF7L2",
            "LEF1",
            "WNT3A",
            "WNT7A",
        ],
    },
    {
        "module_id": "medial_pallium_dentate_lineage",
        "module_label": "Medial pallium/dentate lineage",
        "module_family": "regional_origin",
        "interpretation_role": "dentate/hippocampal granule developmental trajectory and WNT/PROX1-associated fate",
        "genes": [
            "PROX1",
            "NEUROD1",
            "NEUROD2",
            "EOMES",
            "TBR1",
            "LEF1",
            "LHX2",
            "EMX2",
            "ZBTB20",
            "TCF7L2",
            "WNT3A",
            "WNT7A",
            "DKK3",
            "BCL11B",
            "CALB1",
            "C1QL3",
            "GLIS3",
        ],
    },
    {
        "module_id": "hindbrain_rhombic_lip_pattern",
        "module_label": "Hindbrain/rhombic-lip pattern",
        "module_family": "regional_origin",
        "interpretation_role": "posterior neural-tube, isthmic, rhombic-lip, and SHH-associated cerebellar granule lineage",
        "genes": [
            "GBX2",
            "EN1",
            "EN2",
            "WNT1",
            "FGF8",
            "PAX2",
            "PAX5",
            "PAX8",
            "ATOH1",
            "BARHL1",
            "ZIC1",
            "ZIC2",
            "ZIC3",
            "MEIS1",
            "MYCN",
            "PTCH1",
            "GLI1",
            "GLI2",
        ],
    },
    {
        "module_id": "shared_postmitotic_granule_maturation",
        "module_label": "Shared postmitotic/granule maturation",
        "module_family": "later_convergence",
        "interpretation_role": "reused postmitotic neuronal maturation and granule-lineage toolkit",
        "genes": [
            "NFIA",
            "NEUROD1",
            "RBFOX3",
            "HMGN2",
            "DCX",
            "TUBB3",
            "STMN2",
            "STMN3",
            "GPM6A",
            "GAP43",
            "CALB2",
        ],
    },
    {
        "module_id": "downstream_neurite_synapse_construction",
        "module_label": "Downstream neurite/synapse construction",
        "module_family": "later_convergence",
        "interpretation_role": "neurite, adhesion, synaptic, and excitability implementation machinery",
        "genes": [
            "GPM6A",
            "ROBO2",
            "DCC",
            "CADM3",
            "STMN2",
            "STMN3",
            "GAP43",
            "DPYSL2",
            "DPYSL3",
            "MAP1B",
            "BASP1",
            "NCAM1",
            "L1CAM",
            "KCNK1",
            "GABRA2",
            "GABRB3",
            "GRIN2B",
            "SLC17A6",
            "SLC17A7",
            "SNAP25",
            "SYT1",
        ],
    },
]

MODULE_ORDER = [module["module_id"] for module in MODULES]
MODULE_LABELS = {module["module_id"]: module["module_label"] for module in MODULES}
SHORT_MODULE_LABELS = {
    "deep_neural_progenitor": "Deep neural\nprogenitor",
    "anterior_telencephalic_pattern": "Anterior /\ntelencephalic",
    "medial_pallium_dentate_lineage": "Medial pallium /\ndentate",
    "hindbrain_rhombic_lip_pattern": "Hindbrain /\nrhombic lip",
    "shared_postmitotic_granule_maturation": "Postmitotic\nmaturation",
    "downstream_neurite_synapse_construction": "Neurite /\nsynapse",
}

DENTATE_FOCUS_ORDER = {
    "RGL_young": (0.0, "early_progenitor", "RGL_young"),
    "RGL": (1.0, "early_progenitor", "RGL"),
    "nIPC": (2.0, "intermediate_progenitor", "nIPC"),
    "nIPC-perin": (3.0, "intermediate_progenitor", "nIPC-perin"),
    "Neuroblast": (4.0, "postmitotic_immature", "Neuroblast"),
    "Immature-GC": (5.0, "postmitotic_immature", "Immature-GC"),
    "GC-juv": (6.0, "mature_granule", "GC-juv"),
    "GC-adult": (7.0, "mature_granule", "GC-adult"),
}

CEREBELLAR_SAMPLE_ORDER = {"P0": 0.0, "P8a": 8.1, "P8b": 8.2}
CEREBELLAR_FOCUS_OFFSET = {"Granule precursor": 0.0, "Granule cells": 0.25}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().upper()


def mouse_case(symbol: str) -> str:
    if not symbol:
        return symbol
    return symbol[:1].upper() + symbol[1:].lower()


def write_gene_sets() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module in MODULES:
        for order, gene in enumerate(module["genes"], start=1):
            rows.append(
                {
                    "module_id": module["module_id"],
                    "module_label": module["module_label"],
                    "module_family": module["module_family"],
                    "interpretation_role": module["interpretation_role"],
                    "gene": gene,
                    "canonical_gene": canon(gene),
                    "default_mouse_symbol": mouse_case(gene),
                    "gene_order": order,
                }
            )
    gene_sets = pd.DataFrame(rows)
    gene_sets.to_csv(OUT_GENE_SETS, sep="\t", index=False)
    return gene_sets


def compute_module_units(
    *,
    dataset: str,
    sample: str,
    branch: str,
    expression: pd.DataFrame,
    cell_groups: pd.Series,
    gene_sets: pd.DataFrame,
) -> list[dict[str, object]]:
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    index_by_canon = {canon(gene): gene for gene in expression.index}
    common_cells = [cell for cell in expression.columns if cell in set(cell_groups.index)]
    expression = expression[common_cells]
    cell_groups = cell_groups.loc[common_cells]

    records: list[dict[str, object]] = []
    for module_id, sub in gene_sets.groupby("module_id", sort=False):
        meta = sub.iloc[0]
        genes = sub["gene"].astype(str).tolist()
        present = [index_by_canon[canon(gene)] for gene in genes if canon(gene) in index_by_canon]
        if not present:
            continue
        cell_scores = np.log1p(expression.loc[present].to_numpy(dtype=float)).mean(axis=0)
        score_df = pd.DataFrame(
            {
                "cell_id": common_cells,
                "source_group": cell_groups.to_numpy(),
                "module_score": cell_scores,
            }
        )
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("module_score", "size"),
                median_module_score=("module_score", "median"),
                mean_module_score=("module_score", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_module_rank"] = group_summary["median_module_score"].rank(pct=True, method="average")
        for _, row in group_summary.iterrows():
            records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "branch": branch,
                    "source_group": str(row["source_group"]),
                    "module_id": module_id,
                    "module_label": meta["module_label"],
                    "module_family": meta["module_family"],
                    "n_defined_genes": int(len(genes)),
                    "n_present_genes": int(len(present)),
                    "present_genes": ",".join(present),
                    "missing_genes": ",".join(gene for gene in genes if canon(gene) not in index_by_canon),
                    "n_cells": int(row["n_cells"]),
                    "median_module_score": float(row["median_module_score"]),
                    "mean_module_score": float(row["mean_module_score"]),
                    "within_sample_module_rank": float(row["within_sample_module_rank"]),
                }
            )
    return records


def load_dentate(gene_sets: pd.DataFrame, wanted_genes: set[str]) -> list[dict[str, object]]:
    meta = pd.read_csv(base.GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    expression = base.read_selected_rows(base.GSE104323_EXPR, "\t", wanted_genes)
    return compute_module_units(
        dataset="GSE104323",
        sample="10X_all_cells",
        branch="dentate",
        expression=expression,
        cell_groups=cell_groups,
        gene_sets=gene_sets,
    )


def load_cerebellum(gene_sets: pd.DataFrame, wanted_genes: set[str]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in base.GSE122357_FILES.items():
            labels = base.load_gse122357_label_map(prefix)
            expression = base.read_selected_rows_from_tar(base.GSE122357_TAR, member, wanted_genes, tmpdir)
            records.extend(
                compute_module_units(
                    dataset="GSE122357",
                    sample=sample,
                    branch="cerebellar",
                    expression=expression,
                    cell_groups=labels,
                    gene_sets=gene_sets,
                )
            )
    return records


def annotate_focus(units: pd.DataFrame) -> pd.DataFrame:
    def annotate(row: pd.Series) -> pd.Series:
        branch = str(row["branch"])
        group = str(row["source_group"])
        sample = str(row["sample"])
        if branch == "dentate" and group in DENTATE_FOCUS_ORDER:
            order, bin_id, label = DENTATE_FOCUS_ORDER[group]
            row["focus_state"] = True
            row["state_order"] = order
            row["state_bin"] = bin_id
            row["state_label"] = label
        elif branch == "cerebellar" and group in CEREBELLAR_FOCUS_OFFSET:
            row["focus_state"] = True
            row["state_order"] = CEREBELLAR_SAMPLE_ORDER[sample] + CEREBELLAR_FOCUS_OFFSET[group]
            row["state_bin"] = "early_precursor" if sample == "P0" or group == "Granule precursor" else "maturing_granule"
            row["state_label"] = f"{sample} {'GP' if group == 'Granule precursor' else 'GC'}"
        else:
            row["focus_state"] = False
            row["state_order"] = np.nan
            row["state_bin"] = "other"
            row["state_label"] = group
        return row

    return units.apply(annotate, axis=1)


def build_state_summary(units: pd.DataFrame) -> pd.DataFrame:
    focus = units.loc[units["focus_state"]].copy()
    pivot = focus.pivot_table(
        index=["branch", "dataset", "sample", "source_group", "state_label", "state_order", "state_bin", "n_cells"],
        columns="module_id",
        values="within_sample_module_rank",
        aggfunc="median",
    ).reset_index()
    for module_id in MODULE_ORDER:
        if module_id not in pivot:
            pivot[module_id] = np.nan

    pivot["branch_origin_rank"] = np.where(
        pivot["branch"].eq("dentate"),
        pivot[["anterior_telencephalic_pattern", "medial_pallium_dentate_lineage"]].mean(axis=1),
        pivot["hindbrain_rhombic_lip_pattern"],
    )
    pivot["opposed_origin_rank"] = np.where(
        pivot["branch"].eq("dentate"),
        pivot["hindbrain_rhombic_lip_pattern"],
        pivot[["anterior_telencephalic_pattern", "medial_pallium_dentate_lineage"]].mean(axis=1),
    )
    pivot["regional_origin_polarity"] = pivot["branch_origin_rank"] - pivot["opposed_origin_rank"]
    pivot["shared_convergence_rank"] = pivot[
        ["shared_postmitotic_granule_maturation", "downstream_neurite_synapse_construction"]
    ].mean(axis=1)
    pivot["origin_plus_convergence_index"] = pivot["regional_origin_polarity"] + pivot["shared_convergence_rank"]
    return pivot.sort_values(["branch", "state_order"])


def median_for(summary: pd.DataFrame, branch: str, bins: set[str], column: str) -> float:
    values = summary.loc[summary["branch"].eq(branch) & summary["state_bin"].isin(bins), column].dropna()
    return float(values.median()) if not values.empty else np.nan


def build_metrics(summary: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("dentate", "early_to_postmitotic", {"early_progenitor", "intermediate_progenitor"}, {"postmitotic_immature", "mature_granule"}),
        ("dentate", "early_to_mature", {"early_progenitor", "intermediate_progenitor"}, {"mature_granule"}),
        ("cerebellar", "P0_or_precursor_to_maturing", {"early_precursor"}, {"maturing_granule"}),
    ]
    columns = [
        "deep_neural_progenitor",
        "branch_origin_rank",
        "opposed_origin_rank",
        "regional_origin_polarity",
        "shared_postmitotic_granule_maturation",
        "downstream_neurite_synapse_construction",
        "shared_convergence_rank",
        "origin_plus_convergence_index",
    ]
    rows: list[dict[str, object]] = []
    for branch, comparison, early_bins, late_bins in comparisons:
        for column in columns:
            early = median_for(summary, branch, early_bins, column)
            late = median_for(summary, branch, late_bins, column)
            rows.append(
                {
                    "branch": branch,
                    "comparison": comparison,
                    "metric": column,
                    "early_median": early,
                    "late_median": late,
                    "late_minus_early": late - early if pd.notna(early) and pd.notna(late) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def add_box(ax: plt.Axes, xy: tuple[float, float], w: float, h: float, text: str, fc: str, ec: str) -> None:
    patch = Rectangle(xy, w, h, facecolor=fc, edgecolor=ec, linewidth=1.1)
    ax.add_patch(patch)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=7.0, color="#202020", linespacing=1.0)


def add_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#4a4a4a") -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=1.0,
            color=color,
            shrinkA=3,
            shrinkB=3,
            connectionstyle="arc3,rad=0",
        )
    )


def add_orthogonal_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    corner: tuple[float, float],
    end: tuple[float, float],
    color: str = "#4a4a4a",
) -> None:
    ax.plot([start[0], corner[0]], [start[1], corner[1]], color=color, linewidth=1.0)
    add_arrow(ax, corner, end, color=color)


def draw_route_model(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("a. Developmental route model", loc="left", fontsize=11, fontweight="bold")

    left_x, right_x, box_w = 0.025, 0.595, 0.38
    add_box(ax, (0.30, 0.84), 0.40, 0.10, "Neural ectoderm /\nneuroepithelium", "#eef4ff", "#5975a4")
    add_box(ax, (left_x, 0.64), box_w, 0.12, "Anterior\nneural tube /\ntelencephalon", "#eaf5ee", "#498c61")
    add_box(ax, (right_x, 0.64), box_w, 0.12, "Posterior\nneural tube /\nhindbrain", "#fff1e5", "#c87830")
    add_box(ax, (left_x, 0.43), box_w, 0.12, "Medial pallium /\ndentate stream", "#eaf5ee", "#498c61")
    add_box(ax, (right_x, 0.43), box_w, 0.12, "Rhombic lip /\nexternal\ngranule layer", "#fff1e5", "#c87830")
    add_box(ax, (left_x, 0.22), box_w, 0.12, "Dentate granule\nlineage", "#d8eadf", "#498c61")
    add_box(ax, (right_x, 0.22), box_w, 0.12, "Cerebellar\ngranule lineage", "#ffe1c2", "#c87830")
    add_box(ax, (0.30, 0.04), 0.40, 0.10, "Reused maturation +\nconstruction toolkit", "#f3f0fb", "#7460a8")

    left_center, right_center = left_x + box_w / 2, right_x + box_w / 2
    ax.plot([0.50, 0.50], [0.84, 0.80], color="#4a4a4a", linewidth=1.0)
    add_orthogonal_arrow(ax, (0.50, 0.80), (left_center, 0.80), (left_center, 0.76))
    add_orthogonal_arrow(ax, (0.50, 0.80), (right_center, 0.80), (right_center, 0.76))
    add_arrow(ax, (left_center, 0.64), (left_center, 0.55))
    add_arrow(ax, (right_center, 0.64), (right_center, 0.55))
    add_arrow(ax, (left_center, 0.43), (left_center, 0.34))
    add_arrow(ax, (right_center, 0.43), (right_center, 0.34))
    add_orthogonal_arrow(ax, (left_center, 0.22), (left_center, 0.17), (0.42, 0.14), "#7460a8")
    add_orthogonal_arrow(ax, (right_center, 0.22), (right_center, 0.17), (0.58, 0.14), "#7460a8")


def plot_heatmap(ax: plt.Axes, summary: pd.DataFrame, branch: str, title: str) -> None:
    sub = summary.loc[summary["branch"].eq(branch)].sort_values("state_order")
    matrix = sub[MODULE_ORDER].to_numpy(dtype=float).T
    im = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xticks(np.arange(len(sub)))
    ax.set_xticklabels(sub["state_label"].astype(str), rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(np.arange(len(MODULE_ORDER)))
    ax.set_yticklabels([SHORT_MODULE_LABELS.get(module_id, MODULE_LABELS[module_id]) for module_id in MODULE_ORDER], fontsize=7.3)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text_color = "#f8f8f8" if np.isfinite(value) and value < 0.35 else "#202020"
            ax.text(j, i, "NA" if pd.isna(value) else f"{value:.2f}", ha="center", va="center", fontsize=6.5, color=text_color)
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    return im


def plot_index_lines(ax: plt.Axes, summary: pd.DataFrame) -> None:
    colors = {"dentate": "#498c61", "cerebellar": "#c87830"}
    metrics = [
        ("regional_origin_polarity", "regional-origin polarity", "-", "o"),
        ("shared_convergence_rank", "shared-convergence rank", (0, (5, 2)), "^"),
    ]
    for branch in ["dentate", "cerebellar"]:
        sub = summary.loc[summary["branch"].eq(branch)].sort_values("state_order").reset_index(drop=True)
        x = np.arange(len(sub))
        for metric, _label, linestyle, marker in metrics:
            ax.plot(
                x,
                sub[metric],
                marker=marker,
                markersize=5.5,
                linewidth=1.9,
                linestyle=linestyle,
                color=colors[branch],
            )
    ax.axhline(0, color="#555555", linewidth=0.7)
    ax.set_title("d. Indices", loc="left", fontsize=11, fontweight="bold")
    ax.set_ylabel("Within-sample rank / rank delta")
    ax.set_xticks([])
    ax.set_ylim(-0.55, 1.08)
    branch_handles = [
        Line2D([0], [0], color=colors["dentate"], marker="o", linewidth=1.8, label="dentate"),
        Line2D([0], [0], color=colors["cerebellar"], marker="o", linewidth=1.8, label="cerebellar"),
    ]
    metric_handles = [
        Line2D([0], [0], color="#333333", marker="o", markersize=5.5, linewidth=2.0, linestyle="-", label="origin polarity"),
        Line2D([0], [0], color="#333333", marker="^", markersize=6.0, linewidth=2.0, linestyle=(0, (5, 2)), label="shared convergence"),
    ]
    branch_legend = ax.legend(
        handles=branch_handles,
        title="Branch",
        loc="lower left",
        bbox_to_anchor=(0.00, 0.02),
        fontsize=7.0,
        title_fontsize=7.2,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.6,
    )
    ax.add_artist(branch_legend)
    ax.legend(
        handles=metric_handles,
        title="Line meaning",
        loc="lower left",
        bbox_to_anchor=(0.42, 0.02),
        fontsize=7.0,
        title_fontsize=7.2,
        frameon=False,
        borderaxespad=0.0,
        handlelength=3.2,
    )
    ax.spines[["top", "right"]].set_visible(False)


def plot_summary(summary: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(18.0, 9.4))
    grid = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.65, 1.12], height_ratios=[1.0, 1.0], wspace=0.52, hspace=0.58)
    ax_model = fig.add_subplot(grid[:, 0])
    draw_route_model(ax_model)
    ax_dentate = fig.add_subplot(grid[0, 1:])
    im = plot_heatmap(ax_dentate, summary, "dentate", "b. Dentate ordered states")
    ax_cereb = fig.add_subplot(grid[1, 1])
    plot_heatmap(ax_cereb, summary, "cerebellar", "c. Cerebellar ordered states")
    ax_lines = fig.add_subplot(grid[1, 2])
    plot_index_lines(ax_lines, summary)
    pos = ax_lines.get_position()
    ax_lines.set_position([pos.x0 - 0.025, pos.y0, pos.width, pos.height])
    cbar = fig.colorbar(im, ax=[ax_dentate, ax_cereb], fraction=0.024, pad=0.01)
    cbar.set_label("Within-sample module rank")
    fig.suptitle("Developmental-origin control: deep origin, regional split, later toolkit reuse", fontsize=14.5, fontweight="bold")
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    SUPP_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUPP_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def metric_value(metrics: pd.DataFrame, branch: str, comparison: str, metric: str, column: str = "late_minus_early") -> float:
    sub = metrics.loc[
        metrics["branch"].eq(branch) & metrics["comparison"].eq(comparison) & metrics["metric"].eq(metric),
        column,
    ]
    return float(sub.iloc[0]) if not sub.empty else np.nan


def write_report(summary: pd.DataFrame, metrics: pd.DataFrame) -> None:
    dentate_origin_early = median_for(
        summary, "dentate", {"early_progenitor", "intermediate_progenitor"}, "regional_origin_polarity"
    )
    dentate_origin_late = median_for(
        summary, "dentate", {"postmitotic_immature", "mature_granule"}, "regional_origin_polarity"
    )
    cereb_origin_early = median_for(summary, "cerebellar", {"early_precursor"}, "regional_origin_polarity")
    cereb_origin_late = median_for(summary, "cerebellar", {"maturing_granule"}, "regional_origin_polarity")
    dentate_shared_shift = metric_value(metrics, "dentate", "early_to_postmitotic", "shared_convergence_rank")
    cereb_shared_shift = metric_value(metrics, "cerebellar", "P0_or_precursor_to_maturing", "shared_convergence_rank")
    dentate_deep_shift = metric_value(metrics, "dentate", "early_to_postmitotic", "deep_neural_progenitor")
    cereb_deep_shift = metric_value(metrics, "cerebellar", "P0_or_precursor_to_maturing", "deep_neural_progenitor")

    lines = [
        "# Developmental-Origin Divergence Audit",
        "",
        "## Purpose",
        "",
        "This support analysis asks whether the available local full-expression matrices are consistent with a deep-origin/regional-divergence/later-convergence model. It is not a clonal lineage-tracing analysis and cannot prove that dentate and cerebellar granule cells descend from the same individual embryonic progenitor.",
        "",
        "## Marker Modules",
        "",
        "- Deep neural progenitor: broad neural ectoderm/neuroepithelium and progenitor-competence markers.",
        "- Anterior/telencephalic and medial pallium/dentate modules: forebrain and hippocampal/dentate lineage markers.",
        "- Hindbrain/rhombic-lip module: posterior neural-tube, isthmic, rhombic-lip, and cerebellar granule-lineage markers.",
        "- Shared postmitotic/construction modules: reused maturation, neurite, synaptic, and excitability machinery.",
        "",
        "## Main Findings",
        "",
        f"- Dentate ordered states show branch-matched regional-origin polarity that shifts from median {dentate_origin_early:.3f} in early/intermediate progenitor states to {dentate_origin_late:.3f} in postmitotic/mature granule states.",
        f"- Cerebellar ordered states retain cerebellar/rhombic-lip polarity, with median {cereb_origin_early:.3f} in P0 or precursor states and {cereb_origin_late:.3f} in maturing granule-cell states.",
        f"- Shared postmitotic/construction convergence rank changes by {dentate_shared_shift:.3f} from dentate progenitor to postmitotic/mature states and by {cereb_shared_shift:.3f} from cerebellar P0/precursor to maturing states.",
        f"- Deep neural-progenitor rank changes by {dentate_deep_shift:.3f} in the dentate trajectory and by {cereb_deep_shift:.3f} in the cerebellar trajectory; this module is a competence/origin control, not evidence of a shared recent clone.",
        "",
        "## Manuscript-Safe Interpretation",
        "",
        "- The first common ancestor of the two lineages is best placed at the deep neural ectoderm/neuroepithelium level, before anterior-posterior regional patterning.",
        "- The available transcriptomic data support distinct regional paths: telencephalic/medial-pallial/dentate versus hindbrain/rhombic-lip/cerebellar.",
        "- Similar granule-cell morphology is therefore better framed as later reuse of postmitotic maturation and construction programs within distinct regional identities, not as migration of one recent granule progenitor to both sites.",
        "- True migration/clonal tracking would require early embryonic lineage tracing, barcode lineage reconstruction, or spatial transcriptomic atlases that include neural plate/tube patterning through both branch origins.",
        "",
        "## Outputs",
        "",
        f"- Gene-set table: `{rel(OUT_GENE_SETS)}`",
        f"- Module unit table: `{rel(OUT_UNITS)}`",
        f"- State summary: `{rel(OUT_STATE_SUMMARY)}`",
        f"- Branch metrics: `{rel(OUT_METRICS)}`",
        f"- Figure: `{rel(OUT_PLOT)}`",
        f"- Supplementary figure copy: `{rel(OUT_SUPP_PLOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    gene_sets = write_gene_sets()
    wanted_genes = set(gene_sets["gene"].astype(str)) | set(gene_sets["default_mouse_symbol"].astype(str))
    records = []
    records.extend(load_dentate(gene_sets, wanted_genes))
    records.extend(load_cerebellum(gene_sets, wanted_genes))

    units = annotate_focus(pd.DataFrame(records))
    summary = build_state_summary(units)
    metrics = build_metrics(summary)

    units.to_csv(OUT_UNITS, sep="\t", index=False)
    summary.to_csv(OUT_STATE_SUMMARY, sep="\t", index=False)
    metrics.to_csv(OUT_METRICS, sep="\t", index=False)
    plot_summary(summary)
    write_report(summary, metrics)

    print(f"Wrote {rel(OUT_MD)}")
    print(summary[["branch", "state_label", "deep_neural_progenitor", "regional_origin_polarity", "shared_convergence_rank"]].to_string(index=False))
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
