#!/usr/bin/env python3
"""Build streamlined main figures for the comparator-relative analysis."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
OUT = ROOT / "Project/manuscript/main_figures"
OUT.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import seaborn as sns


TIER1 = "Tier 1 core convergent program"
TIER2 = "Tier 2 high-confidence wiring/synaptic executor"
GREEN = "#1B9E77"
PURPLE = "#6C63B5"
ORANGE = "#D95F02"
BLUE = "#3B6FB6"
GRAY = "#64748B"
LIGHT = "#F5F7FA"


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.06, 1.03, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")


def box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, color: str) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.2,
        edgecolor=color,
        facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=10)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = GRAY) -> None:
    ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": color})


def build_figure1() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    ax = axes[0, 0]
    ax.set_axis_off()
    panel_label(ax, "a")
    events = [
        ("Histological naming", "small, densely packed neurons"),
        ("Altman and Das, 1966", "postnatal microneuron production"),
        ("RU49/Zipro1, 1996", "a proposed shared developmental element"),
        ("Lu et al., 2005", "cerebellar medium alters dentate proliferation"),
        ("Public single-cell era", "test identity, recurrence and specificity"),
    ]
    xs = np.linspace(0.08, 0.92, len(events))
    ax.plot([xs[0], xs[-1]], [0.52, 0.52], color=GRAY, lw=1.5)
    for index, ((title, subtitle), x) in enumerate(zip(events, xs, strict=True)):
        ax.scatter(x, 0.52, s=70, color=GREEN if index in {2, 4} else BLUE, zorder=3)
        y = 0.72 if index % 2 == 0 else 0.29
        ax.text(x, y, title, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(x, y - 0.085, subtitle, ha="center", va="center", fontsize=8.5, color="#334155", wrap=True)
        ax.plot([x, x], [0.55 if y > 0.5 else 0.49, y - 0.04 if y > 0.5 else y + 0.04], color="#94A3B8", lw=0.9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax = axes[0, 1]
    ax.set_axis_off()
    panel_label(ax, "b")
    box(ax, (0.05, 0.72), 0.37, 0.15, "Telencephalic hippocampal lineage\nDentate granule cell", BLUE)
    box(ax, (0.58, 0.72), 0.37, 0.15, "Hindbrain rhombic-lip lineage\nCerebellar granule cell", ORANGE)
    arrow(ax, (0.235, 0.695), (0.235, 0.565), BLUE)
    arrow(ax, (0.765, 0.695), (0.765, 0.565), ORANGE)
    box(ax, (0.13, 0.38), 0.74, 0.16, "Do independently specified neurons reuse a limited\ncomparator-relative differentiation program?", GREEN)
    arrow(ax, (0.5, 0.355), (0.5, 0.235), GRAY)
    box(ax, (0.21, 0.08), 0.58, 0.12, "A recent shared migratory progenitor is not assumed", GRAY)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax = axes[1, 0]
    panel_label(ax, "c")
    labels = ["Mouse dentate", "Mouse cerebellum", "Human dentate", "Human cerebellum"]
    counts = [4, 2, 3, 1]
    colors = [BLUE, ORANGE, BLUE, ORANGE]
    bars = ax.barh(labels[::-1], counts[::-1], color=colors[::-1])
    for bar, value in zip(bars, counts[::-1], strict=True):
        ax.text(value + 0.08, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=11)
    ax.set_xlim(0, 4.6)
    ax.set_xlabel("Strict primary-core datasets")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 1]
    ax.set_axis_off()
    panel_label(ax, "d")
    steps = [
        (0.025, "Mouse-first\ndiscovery", BLUE),
        (0.275, "Dataset-level\nrobustness", GREEN),
        (0.525, "Allen\ncommon-matrix\nexternal test", ORANGE),
        (0.775, "Narrowed\nbiological\ninterpretation", PURPLE),
    ]
    for x, text, color in steps:
        box(ax, (x, 0.50), 0.18, 0.20, text, color)
    for left, right in zip(steps[:-1], steps[1:], strict=True):
        arrow(ax, (left[0] + 0.195, 0.60), (right[0] - 0.015, 0.60))
    ax.text(
        0.5,
        0.25,
        "Direct adult module convergence and developmental causality remain unproven",
        ha="center",
        va="center",
        fontsize=10,
        color="#7F1D1D",
        bbox={"facecolor": "#FEF2F2", "edgecolor": "#FCA5A5", "boxstyle": "round,pad=0.5"},
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    fig.savefig(OUT / "Figure_1.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure2() -> None:
    shared = pd.read_csv(RESULTS / "dgd_species_shared_support.tsv", sep="\t")
    lodo = pd.read_csv(RESULTS / "dgd_candidate_lodo_summary.tsv", sep="\t")
    null = pd.read_csv(RESULTS / "dgd_matched_null_summary.tsv", sep="\t")
    tiers = pd.read_csv(RESULTS / "primary_core_manuscript_candidate_tiers.tsv", sep="\t")
    tiers = tiers[tiers["manuscript_tier"].isin([TIER1, TIER2])][
        ["gene", "manuscript_tier", "mechanism_class"]
    ].copy()

    route_map = {
        ("mouse", "selected"): "Mouse selected",
        ("mouse", "full_matrix"): "Mouse full matrix",
        ("human", "selected"): "Human selected bridge",
    }
    shared["route"] = [route_map.get((row.species, row.screen)) for row in shared.itertuples()]
    heat = shared.dropna(subset=["route"]).pivot(
        index="canonical_gene", columns="route", values="shared_minimum_median_delta"
    )
    gene_order = tiers["gene"].tolist()
    heat = heat.reindex(index=gene_order, columns=list(route_map.values()))

    sns.set_theme(style="whitegrid", context="notebook")
    fig = plt.figure(figsize=(16.5, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.25, 1.25, 1.1])

    ax = fig.add_subplot(grid[:, 0])
    panel_label(ax, "a")
    sns.heatmap(
        heat,
        ax=ax,
        cmap="BrBG",
        center=0,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        cbar_kws={"label": "Minimum dentate/cerebellar median rank delta", "shrink": 0.55},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=30)

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "b")
    plot = lodo.merge(tiers, left_on=["canonical_gene", "manuscript_tier"], right_on=["gene", "manuscript_tier"])
    plot = plot.sort_values(["manuscript_tier", "minimum_lodo_median_delta"])
    y = np.arange(len(plot))
    colors = plot["manuscript_tier"].map({TIER1: GREEN, TIER2: PURPLE})
    ax.scatter(plot["minimum_lodo_median_delta"], y, c=colors, s=55)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y, plot["canonical_gene"])
    ax.set_xlabel("Minimum median delta across leave-one-dataset tests")
    ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(grid[0, 2])
    panel_label(ax, "c")
    metric = null[null["metric"].eq("mean_minimum_branch_median_delta")].copy()
    xpos = np.arange(len(metric))
    ax.errorbar(
        xpos,
        metric["null_median"],
        yerr=[metric["null_median"] - metric["null_q025"], metric["null_q975"] - metric["null_median"]],
        fmt="o",
        color=GRAY,
        capsize=5,
        label="Matched null median and 95% interval",
    )
    ax.scatter(xpos, metric["observed_value"], marker="D", s=70, color=ORANGE, label="Observed candidates", zorder=3)
    ax.set_xticks(xpos, ["Tier 1", "Tier 1+2"])
    ax.set_ylabel("Mean minimum branch median delta")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(grid[1, 1:])
    panel_label(ax, "d")
    class_labels = {
        "regulatory_morphogenesis_candidate": "Regulatory/morphogenesis",
        "synaptic_wiring": "Synaptic/wiring",
        "cytoskeleton_morphogenesis": "Cytoskeleton/morphogenesis",
        "curated_shared_structural_executor": "Curated structural executor",
    }
    tiers["mechanism_label"] = tiers["mechanism_class"].map(class_labels).fillna(tiers["mechanism_class"])
    class_counts = tiers.groupby(["mechanism_label", "manuscript_tier"]).size().unstack(fill_value=0)
    class_counts.columns = ["Tier 1" if col == TIER1 else "Tier 2" for col in class_counts.columns]
    class_counts = class_counts.sort_values(class_counts.columns.tolist(), ascending=False)
    class_counts.plot(kind="barh", stacked=True, ax=ax, color=[GREEN, PURPLE])
    ax.invert_yaxis()
    ax.set_xlabel("Candidate genes")
    ax.set_ylabel("")
    ax.legend(frameon=False, title="")
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(OUT / "Figure_2.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure3() -> None:
    datasets = pd.read_csv(RESULTS / "dgd_dataset_level_configuration.tsv", sep="\t")
    summaries = pd.read_csv(RESULTS / "dgd_dataset_level_configuration_summary.tsv", sep="\t")
    modules = pd.read_csv(RESULTS / "dgd_module_level_inference.tsv", sep="\t")
    module_rows = modules[modules["interpretation"].eq("module_summary")].copy()

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), constrained_layout=True)

    ax = axes[0]
    panel_label(ax, "a")
    data = datasets.sort_values("median_delta_configuration_score")
    y = np.arange(len(data))
    colors = data["region"].map({"dentate": BLUE, "cerebellum": ORANGE})
    ax.scatter(data["median_delta_configuration_score"], y, c=colors, s=65)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y, data["dataset"])
    ax.set_xlabel("Dataset median configuration delta")
    ax.text(0.02, 0.98, "7/7 positive\nexact sign p=0.0078", transform=ax.transAxes, va="top", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    panel_label(ax, "b")
    order = ["all", "region:dentate", "region:cerebellum", "species:mouse", "species:human"]
    labels = ["All", "Dentate", "Cerebellum", "Mouse", "Human"]
    s = summaries.set_index("summary_group").reindex(order)
    x = np.arange(len(s))
    yerr = [s["median_dataset_delta"] - s["bootstrap_95ci_low"], s["bootstrap_95ci_high"] - s["median_dataset_delta"]]
    ax.errorbar(x, s["median_dataset_delta"], yerr=yerr, fmt="o", color=GREEN, capsize=5)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, labels, rotation=25)
    ax.set_ylabel("Median dataset delta, bootstrap 95% CI")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[2]
    panel_label(ax, "c")
    module_rows = module_rows.sort_values("median_overall_convergence_delta")
    colors = module_rows["inference_group"].map({"downstream": GREEN, "upstream_or_niche": GRAY})
    ax.barh(module_rows["module_label"], module_rows["median_overall_convergence_delta"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Median convergence delta")
    ax.text(
        0.02,
        0.03,
        "Downstream > upstream: exact p=0.10\nDirectional, not conventionally significant",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        bbox={"facecolor": "white", "edgecolor": "#CBD5E1", "boxstyle": "round,pad=0.4"},
    )
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(OUT / "Figure_3.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    build_figure1()
    build_figure2()
    build_figure3()
    print(f"Wrote Figures 1-3 to {OUT}")


if __name__ == "__main__":
    main()
