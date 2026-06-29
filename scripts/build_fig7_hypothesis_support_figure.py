#!/usr/bin/env python3
"""Build final Fig. 7: hypothesis support matrix and integrated model.

The figure is generated from the hypothesis-support outputs built by
build_hypothesis_support_score_matrix.py. It is intended as a manuscript-facing
summary of evidence alignment, not as a causal model-selection result.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
FINAL_FIGURES = ROOT / "Project/manuscript/Final figures"
PANELS = ROOT / "Project/manuscript/Figure panels"

SCORE_SCRIPT = ROOT / "Project/scripts/build_hypothesis_support_score_matrix.py"
MATRIX = RESULTS / "hypothesis_support_score_matrix.tsv"
SCORES = RESULTS / "hypothesis_support_scores.tsv"

OUT_PNG = RESULTS / "manuscript_figures/figure7_hypothesis_support_model.png"
OUT_PDF = RESULTS / "manuscript_figures/figure7_hypothesis_support_model.pdf"
OUT_JPG = FINAL_FIGURES / "Fig.7.jpg"
OUT_PANEL = PANELS / "figure7_hypothesis_support_model_panel.png"
OUT_PANEL_A = PANELS / "figure7a_hypothesis_relationship_panel.png"


TERM_LABELS = {
    "Y": "granule\nconfiguration",
    "F": "regional\nfate polarity",
    "C": "construction\nbalance",
    "I": "fate x\nconstruction",
    "T": "stage /\npseudotime",
    "N": "niche\nsignal",
    "E": "regulatory\ncompatibility",
    "M": "morphology\nsampling",
    "A": "activity\nsparsity",
    "R": "resource\nconstraint",
}

DOMAIN_COLORS = {
    "direct_transcriptomic_configuration": "#4C78A8",
    "stage_niche_regulatory": "#59A14F",
    "external_morphology_activity_circuit": "#D88435",
}

HYPOTHESIS_ORDER = [
    "H1_shared_granule_fate",
    "H2_identity_coupled_assembly",
    "H3_niche_circuit_constraint",
    "H2_H3_integrated_model",
]

HYPOTHESIS_SHORT = {
    "H1_shared_granule_fate": "H1 fate",
    "H2_identity_coupled_assembly": "H2 assembly",
    "H3_niche_circuit_constraint": "H3 niche/circuit",
    "H2_H3_integrated_model": "H2+H3 integrated",
}


def wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.045,
        1.045,
        label,
        transform=ax.transAxes,
        fontsize=28,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def add_box(ax, xy, width, height, text, fc, ec="#222222", fontsize=12, lw=1.5):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.02",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        linespacing=1.12,
    )
    return box


def add_arrow(ax, start, end, color="#333333", lw=1.6, rad=0.0):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=6,
        shrinkB=6,
    )
    ax.add_patch(arrow)


def add_orthogonal_arrow(ax, points, color="#333333", lw=1.6):
    """Draw an arrow composed only of horizontal and vertical segments."""
    for start, end in zip(points[:-2], points[1:-1]):
        ax.plot([start[0], end[0]], [start[1], end[1]], color=color, linewidth=lw)
    arrow = FancyArrowPatch(
        points[-2],
        points[-1],
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=lw,
        color=color,
        connectionstyle="arc3,rad=0",
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(arrow)


def draw_panel_a(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    panel_label(ax, "a")
    ax.text(
        0.02,
        0.97,
        "Three hypotheses are related but not equivalent",
        ha="left",
        va="top",
        fontsize=20,
        fontweight="bold",
    )

    add_box(
        ax,
        (0.05, 0.65),
        0.32,
        0.20,
        "H1\nhidden shared\nfate identity",
        "#F3D3D3",
        "#9A4B4B",
        fontsize=15,
    )
    add_box(
        ax,
        (0.52, 0.65),
        0.36,
        0.20,
        "H2\nidentity-coupled\nassembly convergence",
        "#DDEAF7",
        "#2F5F91",
        fontsize=15,
    )
    ax.text(
        0.46,
        0.89,
        "transcriptomic alternatives",
        ha="center",
        va="center",
        fontsize=12,
        color="#333333",
    )
    add_arrow(ax, (0.38, 0.74), (0.51, 0.74), "#555555", lw=1.4)
    add_arrow(ax, (0.51, 0.70), (0.38, 0.70), "#555555", lw=1.4)

    add_box(
        ax,
        (0.24, 0.31),
        0.50,
        0.18,
        "H3\nstage, niche and\ncircuit constraints",
        "#F7E1C7",
        "#A15D1D",
        fontsize=15,
    )
    ax.text(
        0.86,
        0.54,
        "higher-level\nexplanatory layer",
        ha="center",
        va="center",
        fontsize=12,
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.92),
    )
    add_orthogonal_arrow(ax, [(0.21, 0.615), (0.21, 0.56), (0.34, 0.56), (0.34, 0.535)], "#555555", lw=1.5)
    add_orthogonal_arrow(ax, [(0.70, 0.615), (0.70, 0.56), (0.62, 0.56), (0.62, 0.535)], "#555555", lw=1.5)

    add_box(
        ax,
        (0.14, 0.04),
        0.72,
        0.18,
        "Final model\nnot H1 alone: distinct lineages + shared assembly state\nfiltered by niche and sparse-expansion constraints",
        "#E8F2E4",
        "#3F7F3A",
        fontsize=12.5,
    )
    add_orthogonal_arrow(ax, [(0.49, 0.295), (0.49, 0.255)], "#555555", lw=1.5)


def draw_panel_b(ax, matrix):
    panel_label(ax, "b")
    terms = matrix["term_symbol"].tolist()
    observed = matrix["observed_term_score"].astype(float).to_numpy()
    colors = [DOMAIN_COLORS[d] for d in matrix["interpretive_domain"]]
    x = np.arange(len(terms))
    ax.bar(x, observed, color=colors, width=0.72)
    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.set_ylim(-0.15, 0.95)
    ax.set_ylabel("Observed evidence score", fontsize=15)
    ax.set_xticks(x)
    ax.set_xticklabels([TERM_LABELS[t] for t in terms], fontsize=12, rotation=28, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="y", labelsize=13)
    ax.set_title("Observed evidence terms from hierarchical synthesis", loc="left", fontsize=21, fontweight="bold")
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.8, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for xi, yi in zip(x, observed):
        ax.text(xi, yi + 0.025, f"{yi:.2f}", ha="center", va="bottom", fontsize=12)
    legend_items = [
        mpl.patches.Patch(color=DOMAIN_COLORS["direct_transcriptomic_configuration"], label="direct transcriptomic"),
        mpl.patches.Patch(color=DOMAIN_COLORS["stage_niche_regulatory"], label="stage/niche/regulatory"),
        mpl.patches.Patch(color=DOMAIN_COLORS["external_morphology_activity_circuit"], label="external/circuit"),
    ]
    ax.legend(handles=legend_items, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.02), ncol=3, fontsize=12)


def draw_panel_c(ax, matrix):
    panel_label(ax, "c")
    terms = matrix["term_symbol"].tolist()
    coeff = np.array(
        [
            matrix[f"{hypothesis_id}_coefficient"].astype(float).to_numpy()
            for hypothesis_id in HYPOTHESIS_ORDER
        ]
    )
    im = ax.imshow(coeff, aspect=1.38, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_anchor("N")
    ax.set_xticks(range(len(terms)))
    ax.set_xticklabels(terms, fontsize=11)
    ax.set_yticks(range(len(HYPOTHESIS_ORDER)))
    ax.set_yticklabels([HYPOTHESIS_SHORT[h] for h in HYPOTHESIS_ORDER], fontsize=11)
    ax.set_title("Prediction coefficients used for hypothesis scoring", loc="left", fontsize=18, fontweight="bold")
    for i in range(coeff.shape[0]):
        for j in range(coeff.shape[1]):
            text_color = "white" if abs(coeff[i, j]) >= 0.68 else "#111111"
            ax.text(j, i, f"{coeff[i, j]:.1f}", ha="center", va="center", fontsize=11, color=text_color, fontweight="bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = plt.colorbar(im, ax=ax, fraction=0.032, pad=0.014, shrink=0.72)
    cbar.set_label("coefficient", fontsize=11)
    cbar.ax.tick_params(labelsize=10.5)


def draw_panel_d(ax, scores):
    panel_label(ax, "d")
    score_map = scores.set_index("hypothesis_id")["support_index_0_to100"].to_dict()
    labels = [HYPOTHESIS_SHORT[h] for h in HYPOTHESIS_ORDER]
    values = [score_map[h] for h in HYPOTHESIS_ORDER]
    colors = ["#B95F5F", "#4C78A8", "#D88435", "#59A14F"]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, height=0.62)
    ax.axvline(50, color="#555555", linestyle="--", linewidth=1.3)
    ax.set_xlim(0, 100)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=14)
    ax.set_xlabel("Support index (0-100; 50 = neutral)", fontsize=15)
    ax.tick_params(axis="x", labelsize=13)
    ax.set_title("Support index favors convergence over H1 alone", loc="left", fontsize=21, fontweight="bold")
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.8, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for yi, value in zip(y, values):
        ax.text(value + 1.3, yi, f"{value:.1f}", va="center", fontsize=14)
    ax.invert_yaxis()
    note = (
        "Evidence-index result:\n"
        "H1 is disfavored.\n"
        "H2 explains the RNA configuration.\n"
        "H3 explains why the design is favored.\n"
        "The conclusion is integrated H2+H3."
    )
    ax.text(
        0.985,
        0.97,
        note,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F7F7F7", edgecolor="#BBBBBB"),
        linespacing=1.25,
    )


def main() -> int:
    import subprocess
    import sys

    subprocess.run([sys.executable, str(SCORE_SCRIPT)], check=True)
    matrix = pd.read_csv(MATRIX, sep="\t")
    scores = pd.read_csv(SCORES, sep="\t")

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(18, 14), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.05], width_ratios=[0.98, 1.08])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b, matrix)
    draw_panel_c(ax_c, matrix)
    draw_panel_d(ax_d, scores)

    fig.subplots_adjust(left=0.075, right=0.975, top=0.96, bottom=0.07, wspace=0.30, hspace=0.34)

    for out in [OUT_PNG, OUT_PDF, OUT_JPG, OUT_PANEL]:
        out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_JPG, dpi=300, pil_kwargs={"quality": 95})
    fig.savefig(OUT_PANEL, dpi=300)
    plt.close(fig)

    panel_fig, panel_ax = plt.subplots(figsize=(8.5, 6.6))
    draw_panel_a(panel_ax)
    panel_fig.subplots_adjust(left=0.03, right=0.98, top=0.97, bottom=0.04)
    panel_fig.savefig(OUT_PANEL_A, dpi=300)
    plt.close(panel_fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
