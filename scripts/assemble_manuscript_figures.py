#!/usr/bin/env python3
"""Assemble draft manuscript figures for the granule-cell convergence project."""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
FIG_DIR = RESULTS / "manuscript_figures"

OUT_MANIFEST = FIG_DIR / "manuscript_figure_manifest.tsv"
OUT_MD = FIG_DIR / "manuscript_figure_assembly.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import patches


COLORS = {
    "dentate": "#2f7f8f",
    "cerebellar": "#7f4e8a",
    "shared": "#526f4f",
    "neutral": "#4b5563",
    "light": "#f3f4f6",
    "dark": "#111827",
    "warning": "#a16207",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def wrap(text: str, width: int = 42) -> str:
    wrapped_lines: list[str] = []
    for line in str(text).splitlines():
        if line.strip():
            wrapped_lines.extend(textwrap.wrap(line, width=width, break_long_words=False))
        else:
            wrapped_lines.append("")
    return "\n".join(wrapped_lines)


def no_axes(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.03,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold",
        color=COLORS["dark"],
    )


def add_box(
    ax,
    xy,
    width,
    height,
    text,
    fc="#ffffff",
    ec="#333333",
    fontsize=10,
    lw=1.2,
    wrap_width: int | None = None,
    pad: float = 0.025,
):
    box = patches.FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad={pad},rounding_size=0.025",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        wrap(text, wrap_width) if wrap_width else text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=COLORS["dark"],
    )
    return box


def add_arrow(ax, start, end, color="#333333", lw=1.2):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="->", color=color, linewidth=lw, shrinkA=4, shrinkB=4),
    )


def add_ortho_arrow(ax, points, color="#333333", lw=1.2, linestyle="-", alpha=1.0):
    """Draw an arrow using only horizontal and vertical segments."""
    if len(points) < 2:
        return
    clean_points = [points[0]]
    for point in points[1:]:
        if point != clean_points[-1]:
            clean_points.append(point)
    if len(clean_points) < 2:
        return

    for start, end in zip(clean_points[:-2], clean_points[1:-1]):
        if start[0] != end[0] and start[1] != end[1]:
            raise ValueError(f"Non-orthogonal connector segment: {start} -> {end}")
        ax.plot(
            [start[0], end[0]],
            [start[1], end[1]],
            color=color,
            linewidth=lw,
            linestyle=linestyle,
            alpha=alpha,
            solid_capstyle="round",
        )

    start = clean_points[-2]
    end = clean_points[-1]
    if start[0] != end[0] and start[1] != end[1]:
        raise ValueError(f"Non-orthogonal connector segment: {start} -> {end}")
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            linewidth=lw,
            linestyle=linestyle,
            alpha=alpha,
            shrinkA=4,
            shrinkB=4,
        ),
    )


def add_image(ax, path: Path, title: str | None = None) -> None:
    image = mpimg.imread(path)
    ax.imshow(image)
    no_axes(ax)
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", pad=5)


def save_figure(fig, name: str, title: str, sources: list[str], purpose: str, rows: list[dict[str, str]]) -> None:
    png = FIG_DIR / f"{name}.png"
    pdf = FIG_DIR / f"{name}.pdf"
    fig.savefig(png, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    rows.append(
        {
            "figure_id": name,
            "title": title,
            "png": rel(png),
            "pdf": rel(pdf),
            "purpose": purpose,
            "source_outputs": ";".join(sources),
        }
    )


def draw_mini_neuron(ax, center, color, label):
    cx, cy = center
    soma = patches.Circle((cx, cy), 0.045, facecolor=color, edgecolor="#222222", linewidth=1)
    ax.add_patch(soma)
    for angle in np.linspace(40, 140, 4):
        rad = np.deg2rad(angle)
        x2 = cx + 0.16 * np.cos(rad)
        y2 = cy + 0.16 * np.sin(rad)
        ax.plot([cx, x2], [cy, y2], color=color, linewidth=2)
        ax.plot([x2, x2 + 0.04 * np.cos(rad + 0.7)], [y2, y2 + 0.04 * np.sin(rad + 0.7)], color=color, linewidth=1.4)
    ax.plot([cx, cx + 0.22], [cy, cy - 0.16], color=color, linewidth=2)
    ax.text(cx, cy + 0.25, label, ha="center", va="bottom", fontsize=10, fontweight="bold")


def draw_figure1a_biological_puzzle(ax, label: str | None = "A") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Biological Puzzle", fontsize=13, fontweight="bold", pad=8)
    draw_mini_neuron(ax, (0.30, 0.52), COLORS["cerebellar"], "Cerebellar\ngranule")
    draw_mini_neuron(ax, (0.70, 0.52), COLORS["dentate"], "Dentate\ngranule")
    ax.text(
        0.50,
        0.13,
        wrap("Similar compact excitatory-neuron design, but distinct origins, anatomy, and circuit roles.", 42),
        ha="center",
        va="center",
        fontsize=10.4,
    )


def draw_figure1b_working_model(ax, label: str | None = "B") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Working Model", fontsize=13, fontweight="bold", pad=8)
    add_box(
        ax,
        (0.06, 0.73),
        0.34,
        0.14,
        "Cerebellar fate\nrhombic lip / SHH",
        fc="#f3e8ff",
        ec=COLORS["cerebellar"],
        fontsize=9.8,
        wrap_width=28,
    )
    add_box(
        ax,
        (0.60, 0.73),
        0.34,
        0.14,
        "Dentate fate\nWNT / PROX1",
        fc="#e0f2fe",
        ec=COLORS["dentate"],
        fontsize=9.8,
        wrap_width=28,
    )
    add_box(
        ax,
        (0.25, 0.46),
        0.50,
        0.14,
        "Downstream assembly\nneurites + synapses + excitability",
        fc="#ecfdf5",
        ec=COLORS["shared"],
        fontsize=9.8,
        wrap_width=42,
    )
    add_box(
        ax,
        (0.25, 0.17),
        0.50,
        0.13,
        "Compact input-expansion\nneuron design",
        fc="#fff7ed",
        ec=COLORS["warning"],
        fontsize=9.8,
        wrap_width=38,
    )
    add_ortho_arrow(ax, [(0.23, 0.73), (0.23, 0.66), (0.42, 0.66), (0.42, 0.60)], COLORS["cerebellar"])
    add_ortho_arrow(ax, [(0.77, 0.73), (0.77, 0.66), (0.58, 0.66), (0.58, 0.60)], COLORS["dentate"])
    add_ortho_arrow(ax, [(0.50, 0.46), (0.50, 0.30)], COLORS["shared"])


def draw_figure2a_workflow(ax, label: str | None = "A") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Ortholog-Aware Rank-Meta Workflow", fontsize=13, fontweight="bold", pad=8)
    workflow = [
        ("10 primary datasets", 0.82),
        ("Pseudobulk ranks\nwithin sample/gene", 0.66),
        ("MGI one-to-one\northolog frame", 0.50),
        ("Formal branch tests\nand shared hits", 0.34),
        ("Manuscript tiers", 0.18),
    ]
    for text, y in workflow:
        add_box(
            ax,
            (0.14, y - 0.055),
            0.72,
            0.11,
            text,
            fc="#f9fafb",
            ec=COLORS["neutral"],
            fontsize=9.1,
            wrap_width=36,
        )
    for (_, y1), (_, y2) in zip(workflow[:-1], workflow[1:]):
        add_ortho_arrow(ax, [(0.50, y1 - 0.055), (0.50, y2 + 0.055)], COLORS["neutral"])
    ax.text(
        0.5,
        0.045,
        wrap("Formal MGI model: 17,611 one-to-one genes; 1,370 shared hits; 36 mechanism-prioritized genes.", 66),
        ha="center",
        va="center",
        fontsize=8.6,
        color=COLORS["neutral"],
    )


def draw_figure4a_configuration(ax, label: str | None = "A") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Transcriptomic Configuration Score", fontsize=13, fontweight="bold", pad=8)
    add_box(
        ax,
        (0.08, 0.68),
        0.84,
        0.14,
        "Construction balance\n= downstream neurite/synaptic - niche",
        fc="#ecfdf5",
        ec=COLORS["shared"],
        fontsize=9.6,
        wrap_width=42,
    )
    add_box(
        ax,
        (0.08, 0.43),
        0.84,
        0.14,
        "Regional fate polarity\n= branch-matched fate - opposed fate",
        fc="#f5f3ff",
        ec="#6d28d9",
        fontsize=9.6,
        wrap_width=42,
    )
    add_box(
        ax,
        (0.08, 0.18),
        0.84,
        0.14,
        "Configuration score\n= construction balance + fate polarity",
        fc="#fff7ed",
        ec=COLORS["warning"],
        fontsize=9.6,
        wrap_width=42,
    )
    add_ortho_arrow(ax, [(0.50, 0.68), (0.50, 0.57)], COLORS["neutral"])
    add_ortho_arrow(ax, [(0.50, 0.43), (0.50, 0.32)], COLORS["neutral"])


def draw_figure5d_final_model(ax, label: str | None = "D") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Revised Final Model", fontsize=13, fontweight="bold", pad=8)
    add_box(ax, (0.05, 0.80), 0.36, 0.10, "Cerebellar lineage\nSHH / rhombic lip", fc="#f3e8ff", ec=COLORS["cerebellar"], fontsize=9.2, wrap_width=28)
    add_box(ax, (0.59, 0.80), 0.36, 0.10, "Dentate lineage\nWNT / PROX1", fc="#e0f2fe", ec=COLORS["dentate"], fontsize=9.2, wrap_width=28)
    add_box(ax, (0.26, 0.62), 0.48, 0.10, "Pseudotime / stage window\nmaturation readiness", fc="#fff7ed", ec=COLORS["warning"], fontsize=9.2, wrap_width=38)
    add_box(ax, (0.23, 0.43), 0.54, 0.10, "Identity-coupled\nassembly configuration", fc="#ecfdf5", ec=COLORS["shared"], fontsize=9.5, wrap_width=38)
    add_box(ax, (0.05, 0.24), 0.39, 0.10, "TGF/BDNF/secreted cues\nbranch-specific timing", fc="#fef3c7", ec="#b45309", fontsize=8.4, wrap_width=34)
    add_box(ax, (0.56, 0.24), 0.39, 0.10, "Resource-constrained\nsparse expansion", fc="#fefce8", ec="#854d0e", fontsize=8.4, wrap_width=34)
    add_box(ax, (0.25, 0.07), 0.50, 0.10, "Convergent compact\ngranule-cell design", fc="#ffffff", ec=COLORS["dark"], fontsize=9.0, wrap_width=36)
    add_ortho_arrow(ax, [(0.23, 0.80), (0.23, 0.74), (0.42, 0.74), (0.42, 0.72)], COLORS["cerebellar"])
    add_ortho_arrow(ax, [(0.77, 0.80), (0.77, 0.74), (0.58, 0.74), (0.58, 0.72)], COLORS["dentate"])
    add_ortho_arrow(ax, [(0.50, 0.62), (0.50, 0.53)], COLORS["warning"])
    add_ortho_arrow(ax, [(0.41, 0.43), (0.41, 0.38), (0.245, 0.38), (0.245, 0.34)], "#b45309")
    add_ortho_arrow(ax, [(0.59, 0.43), (0.59, 0.38), (0.755, 0.38), (0.755, 0.34)], "#854d0e")
    add_ortho_arrow(ax, [(0.245, 0.24), (0.245, 0.19), (0.40, 0.19), (0.40, 0.17)], "#b45309")
    add_ortho_arrow(ax, [(0.755, 0.24), (0.755, 0.19), (0.60, 0.19), (0.60, 0.17)], "#854d0e")


def draw_figure6_working_model(ax, label: str | None = "A") -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Integrated Working Model", fontsize=16, fontweight="bold", pad=12)

    add_box(ax, (0.04, 0.80), 0.22, 0.11, "Primary-core evidence\n10 datasets\nmouse + human", fc="#f9fafb", ec=COLORS["neutral"], fontsize=9.0, wrap_width=30)
    add_box(ax, (0.04, 0.62), 0.22, 0.11, "Ortholog rank-meta\ncandidate tiers\nshared toolkit", fc="#eef2ff", ec="#4338ca", fontsize=9.0, wrap_width=30)
    add_box(ax, (0.04, 0.44), 0.22, 0.11, "Trajectory layer\npseudotime / stage\nwindowed signals", fc="#fff7ed", ec=COLORS["warning"], fontsize=9.0, wrap_width=30)
    add_box(ax, (0.04, 0.26), 0.22, 0.10, "Computation layer\nsparse expansion\ncircuit constraint", fc="#fefce8", ec="#854d0e", fontsize=8.9, wrap_width=30)
    add_box(ax, (0.04, 0.09), 0.22, 0.10, "External validation\nNeuroMorpho morphology\nDANDI 26 units", fc="#f0fdf4", ec=COLORS["shared"], fontsize=8.7, wrap_width=30)

    # Top-to-bottom evidence workflow on the left.
    for y1, y2 in [(0.80, 0.73), (0.62, 0.55), (0.44, 0.36), (0.26, 0.19)]:
        add_ortho_arrow(ax, [(0.15, y1), (0.15, y2)], COLORS["neutral"], lw=0.9)

    add_box(ax, (0.34, 0.78), 0.25, 0.12, "Cerebellar granule lineage\nrhombic lip / SHH\npostnatal expansion", fc="#f3e8ff", ec=COLORS["cerebellar"], fontsize=9.0, wrap_width=34)
    add_box(ax, (0.66, 0.78), 0.25, 0.12, "Dentate granule lineage\nWNT / PROX1\nlifelong neurogenesis", fc="#e0f2fe", ec=COLORS["dentate"], fontsize=9.0, wrap_width=34)
    add_box(ax, (0.39, 0.59), 0.49, 0.11, "Stage-windowed maturation readiness\nTGF/BDNF/SMAD/ERK and secreted cues are overlays", fc="#fff7ed", ec=COLORS["warning"], fontsize=9.0, wrap_width=58)
    add_box(ax, (0.38, 0.40), 0.50, 0.11, "Identity-coupled transcriptomic assembly configuration\nregional fate polarity + downstream construction balance", fc="#ecfdf5", ec=COLORS["shared"], fontsize=9.0, wrap_width=58)
    add_box(ax, (0.33, 0.23), 0.27, 0.10, "Shared construction toolkit\nneurites, synapses,\nexcitability, guidance", fc="#ffffff", ec=COLORS["shared"], fontsize=8.8, wrap_width=34, pad=0.012)
    add_box(ax, (0.66, 0.23), 0.27, 0.10, "Circuit-level filter\nsparse expansion favors\npattern separation", fc="#ffffff", ec="#854d0e", fontsize=8.8, wrap_width=34, pad=0.012)
    add_box(ax, (0.42, 0.07), 0.45, 0.10, "Convergent compact excitatory granule-cell morphology\nsimilar design, different origins", fc="#f8fafc", ec=COLORS["dark"], fontsize=9.0, wrap_width=54, pad=0.012)

    feed_style = "--"
    feed_alpha = 0.72
    add_ortho_arrow(ax, [(0.275, 0.855), (0.31, 0.855), (0.31, 0.84), (0.335, 0.84)], COLORS["neutral"], lw=0.80, linestyle=feed_style, alpha=feed_alpha)
    add_ortho_arrow(ax, [(0.275, 0.675), (0.31, 0.675), (0.31, 0.455), (0.375, 0.455)], "#4338ca", lw=0.80, linestyle=feed_style, alpha=feed_alpha)
    add_ortho_arrow(ax, [(0.275, 0.495), (0.32, 0.495), (0.32, 0.645), (0.385, 0.645)], COLORS["warning"], lw=0.80, linestyle=feed_style, alpha=feed_alpha)
    add_ortho_arrow(ax, [(0.275, 0.31), (0.30, 0.31), (0.30, 0.195), (0.64, 0.195), (0.64, 0.28), (0.655, 0.28)], "#854d0e", lw=0.80, linestyle=feed_style, alpha=feed_alpha)
    add_ortho_arrow(ax, [(0.275, 0.14), (0.34, 0.14), (0.34, 0.12), (0.415, 0.12)], COLORS["shared"], lw=0.80, linestyle=feed_style, alpha=feed_alpha)

    add_ortho_arrow(ax, [(0.465, 0.765), (0.465, 0.725), (0.51, 0.725), (0.51, 0.705)], COLORS["cerebellar"], lw=1.05)
    add_ortho_arrow(ax, [(0.785, 0.765), (0.785, 0.725), (0.76, 0.725), (0.76, 0.705)], COLORS["dentate"], lw=1.05)
    add_ortho_arrow(ax, [(0.635, 0.575), (0.635, 0.525)], COLORS["warning"], lw=1.05)
    add_ortho_arrow(ax, [(0.50, 0.385), (0.50, 0.355), (0.465, 0.355), (0.465, 0.345)], COLORS["shared"], lw=1.05)
    add_ortho_arrow(ax, [(0.76, 0.385), (0.76, 0.355), (0.795, 0.355), (0.795, 0.345)], "#854d0e", lw=1.05)
    add_ortho_arrow(ax, [(0.465, 0.215), (0.465, 0.185)], COLORS["shared"], lw=1.05)
    add_ortho_arrow(ax, [(0.795, 0.215), (0.795, 0.185)], "#854d0e", lw=1.05)


def compact_core_role(role: str) -> str:
    role_map = {
        "primary mouse dentate reference": "Mouse DG reference",
        "mouse dentate maturation validation": "Dentate maturation",
        "mouse postnatal dentate developmental validation": "Postnatal DG development",
        "mouse adult/activity-state dentate validation": "Adult/activity DG",
        "primary mouse cerebellar developmental comparison": "Mouse cerebellar development",
        "primary human cerebellar validation": "Human cerebellar validation",
        "mouse cerebellar perturbation validation": "Cerebellar perturbation",
        "human DG taxonomy anchor": "Human DG taxonomy",
        "primary adult human dentate anchor": "Adult human DG",
        "broader human hippocampal aging/AD RNA expansion": "Human hippocampal aging/AD",
    }
    return role_map.get(str(role), str(role))


def draw_figure1d_dataset_map(
    ax,
    core: pd.DataFrame,
    species_counts: pd.Series,
    label: str | None = "D",
) -> None:
    no_axes(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if label:
        add_panel_label(ax, label)
    ax.set_title("Primary Dataset Map", fontsize=13, fontweight="bold", pad=8)
    branch_order = [
        ("mouse_dentate", "Mouse dentate", COLORS["dentate"]),
        ("cerebellum", "Cerebellum", COLORS["cerebellar"]),
        ("human_dentate_hippocampus", "Human dentate / hippocampus", "#0891b2"),
    ]
    x_starts = [0.055, 0.370, 0.685]
    col_w = 0.260
    header_y = 0.80
    item_top = 0.66
    item_h = 0.108
    gap = 0.022
    for (branch, branch_label, color), x0 in zip(branch_order, x_starts):
        add_box(ax, (x0, header_y), col_w, 0.105, branch_label, fc="#ffffff", ec=color, fontsize=9.4, wrap_width=26, pad=0.012)
        sub = core.loc[core["core_branch"].eq(branch)].copy()
        y = item_top
        for _, row in sub.iterrows():
            text = f"{row['dataset']}\n{compact_core_role(row['core_role'])}"
            add_box(
                ax,
                (x0, y),
                col_w,
                item_h,
                text,
                fc="#f9fafb",
                ec=color,
                fontsize=7.4,
                lw=0.85,
                wrap_width=26,
                pad=0.012,
            )
            y -= item_h + gap
    ax.text(
        0.5,
        0.030,
        f"Species represented: {', '.join(f'{k}={v}' for k, v in species_counts.items())}. "
        "Scaffolds/supporting resources tracked separately.",
        ha="center",
        va="bottom",
        fontsize=8.8,
        color=COLORS["neutral"],
    )


def draw_figure1c_primary_core_counts(ax, core: pd.DataFrame, label: str | None = "C") -> None:
    branch_counts = core.groupby("core_branch")["dataset"].nunique().reindex(
        ["mouse_dentate", "cerebellum", "human_dentate_hippocampus"]
    )
    if label:
        add_panel_label(ax, label)
    ax.barh(range(len(branch_counts)), branch_counts.values, color=[COLORS["dentate"], COLORS["cerebellar"], "#0891b2"])
    ax.set_yticks(range(len(branch_counts)))
    ax.set_yticklabels(["Mouse DG", "Cerebellum", "Human DG/HC"], fontsize=9.5)
    ax.set_xlabel("Datasets")
    ax.set_title("Strict 10-Dataset Primary Core", fontsize=13, fontweight="bold")
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)
    for i, v in enumerate(branch_counts.values):
        ax.text(v + 0.05, i, str(int(v)), va="center", fontsize=11)
    ax.set_xlim(0, max(branch_counts.values) + 1)


def load_figure2_data() -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    tiers = pd.read_csv(RESULTS / "primary_core_manuscript_candidate_tiers.tsv", sep="\t")
    axis_summary = pd.read_csv(RESULTS / "primary_core_mechanism_axis_summary.tsv", sep="\t")
    axis_summary = axis_summary.loc[~axis_summary["mechanism_axis"].eq("exploratory_ortholog_completeness")].copy()
    axis_summary = axis_summary.sort_values(["n_tier1", "n_tier2", "n_genes"], ascending=False)
    top = tiers.loc[
        tiers["manuscript_tier"].isin(
            ["Tier 1 core convergent program", "Tier 2 high-confidence wiring/synaptic executor"]
        )
    ].copy()
    top = top.sort_values(["tier_rank", "gene"]).head(15)
    heat_cols = [
        "median_delta_selected_dentate",
        "median_delta_selected_cerebellar",
        "median_delta_full_matrix_dentate",
        "median_delta_full_matrix_cerebellar",
    ]
    return top, heat_cols, axis_summary


def draw_figure2b_candidate_deltas(
    ax,
    top: pd.DataFrame,
    heat_cols: list[str],
    label: str | None = "B",
):
    heat = top[heat_cols].to_numpy(dtype=float)
    if label:
        add_panel_label(ax, label)
    im = ax.imshow(heat, vmin=-0.5, vmax=0.5, cmap="RdBu_r", aspect="auto")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{g} ({'T1' if r == 1 else 'T2'})" for g, r in zip(top["gene"], top["tier_rank"])], fontsize=9)
    ax.set_xticks(range(len(heat_cols)))
    ax.set_xticklabels(["Sel dentate", "Sel cereb", "Full dentate", "Full cereb"], rotation=25, ha="right", fontsize=9)
    ax.set_title("Tier 1/2 Candidate Branch Deltas", fontsize=13, fontweight="bold")
    return im


def draw_figure2c_mechanism_axes(ax, axis_summary: pd.DataFrame, label: str | None = "C") -> None:
    if label:
        add_panel_label(ax, label)
    y_pos = np.arange(len(axis_summary))
    left = np.zeros(len(axis_summary))
    tier_cols = [
        ("n_tier1", "Tier 1", COLORS["dentate"]),
        ("n_tier2", "Tier 2", COLORS["cerebellar"]),
        ("n_tier3_4", "Tier 3/4", COLORS["shared"]),
    ]
    for col, tier_label, color in tier_cols:
        vals = axis_summary[col].to_numpy(dtype=float)
        ax.barh(y_pos, vals, left=left, color=color, edgecolor="white", height=0.62, label=tier_label)
        for i, val in enumerate(vals):
            if val > 0:
                ax.text(left[i] + val / 2, i, str(int(val)), ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        left += vals
    ax.set_yticks(y_pos)
    ax.set_yticklabels(axis_summary["mechanism_axis_label"].tolist(), fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Mechanism-prioritized genes")
    ax.set_title("Mechanism-axis Organization", fontsize=13, fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.6)
    ax.legend(frameon=False, loc="lower right", ncol=3, fontsize=9)
    ax.set_xlim(0, max(left) + 5.0)
    for i, row in enumerate(axis_summary.itertuples(index=False)):
        genes = ", ".join(
            ", ".join(g.split(","))
            for g in [row.tier1_genes, row.tier2_genes]
            if isinstance(g, str) and g
        )
        if genes:
            ax.text(left[i] + 0.25, i, wrap(genes, 36), ha="left", va="center", fontsize=8.4, color=COLORS["neutral"])


def figure1(rows: list[dict[str, str]]) -> None:
    core = pd.read_csv(RESULTS / "integrated_primary_core_datasets.tsv", sep="\t")
    species_counts = core.groupby("species")["dataset"].nunique()

    fig = plt.figure(figsize=(15.8, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], width_ratios=[1.1, 1.2, 1.1], hspace=0.28, wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    draw_figure1a_biological_puzzle(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    draw_figure1b_working_model(ax, "B")

    ax = fig.add_subplot(gs[0, 2])
    draw_figure1c_primary_core_counts(ax, core, "C")

    ax = fig.add_subplot(gs[1, :])
    draw_figure1d_dataset_map(ax, core, species_counts, "D")
    fig.suptitle("Figure 1. Granule-cell convergence question and primary-core design", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure1_primary_core_concept",
        "Granule-cell convergence question and primary-core design",
        [rel(RESULTS / "integrated_primary_core_datasets.tsv"), rel(RESULTS / "primary_core_niche_circuit_module_model.md")],
        "Introduces the biological puzzle, working model, and strict 10-dataset core.",
        rows,
    )


def figure2(rows: list[dict[str, str]]) -> None:
    top, heat_cols, axis_summary = load_figure2_data()

    fig = plt.figure(figsize=(16.0, 10.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.2], height_ratios=[0.9, 1.2], hspace=0.28, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0])
    draw_figure2a_workflow(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    im = draw_figure2b_candidate_deltas(ax, top, heat_cols, "B")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="Median candidate-background rank delta")

    ax = fig.add_subplot(gs[1, :])
    draw_figure2c_mechanism_axes(ax, axis_summary, "C")

    fig.suptitle("Figure 2. Ortholog-aware candidate tiers and mechanism axes", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure2_ortholog_candidate_tiers",
        "Ortholog-aware candidate tiers and mechanism axes",
        [
            rel(RESULTS / "primary_core_manuscript_candidate_tiers.tsv"),
            rel(RESULTS / "primary_core_mechanism_axis_summary.tsv"),
        ],
        "Shows the analysis workflow, Tier 1/2 candidates, and biological mechanism axes.",
        rows,
    )


def figure3(rows: list[dict[str, str]]) -> None:
    fig = plt.figure(figsize=(16.0, 11.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.72, 1.0], hspace=0.20)
    ax = fig.add_subplot(gs[0, 0])
    add_panel_label(ax, "A")
    add_image(ax, RESULTS / "primary_core_granule_specificity_named_comparators.png", "Named-comparator specificity audit")
    ax = fig.add_subplot(gs[1, 0])
    add_panel_label(ax, "B")
    add_image(ax, RESULTS / "primary_core_niche_circuit_module_model.png", "Niche/fate versus circuit/morphology model")
    fig.suptitle("Figure 3. Specificity constraints and downstream convergence", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure3_specificity_niche_circuit",
        "Specificity constraints and downstream convergence",
        [
            rel(RESULTS / "primary_core_granule_specificity_named_comparators.png"),
            rel(RESULTS / "primary_core_niche_circuit_module_model.png"),
        ],
        "Places the named-comparator caveat beside the upstream-versus-downstream convergence result.",
        rows,
    )


def figure4(rows: list[dict[str, str]]) -> None:
    fig = plt.figure(figsize=(16.5, 12.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.2], width_ratios=[0.8, 1.2], hspace=0.22, wspace=0.16)

    ax = fig.add_subplot(gs[0, 0])
    draw_figure4a_configuration(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    add_panel_label(ax, "B")
    add_image(ax, RESULTS / "primary_core_transcriptomic_configuration_model.png", "Local named-comparator configuration")

    ax = fig.add_subplot(gs[1, 0])
    add_panel_label(ax, "C")
    add_image(ax, RESULTS / "primary_core_transcriptomic_configuration_primary_validation.png", "Primary-core validation")

    ax = fig.add_subplot(gs[1, 1])
    add_panel_label(ax, "D")
    add_image(ax, RESULTS / "primary_core_configuration_driver_audit.png", "Driver decomposition")

    fig.suptitle("Figure 4. Identity-coupled transcriptomic assembly configuration", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure4_transcriptomic_configuration",
        "Identity-coupled transcriptomic assembly configuration",
        [
            rel(RESULTS / "primary_core_transcriptomic_configuration_model.png"),
            rel(RESULTS / "primary_core_transcriptomic_configuration_primary_validation.png"),
            rel(RESULTS / "primary_core_configuration_driver_audit.png"),
        ],
        "Defines, validates, and decomposes the identity-coupled configuration score; interpret with the stage-window trajectory layer in Figure 5.",
        rows,
    )


def figure5(rows: list[dict[str, str]]) -> None:
    fig = plt.figure(figsize=(17.2, 13.2))
    gs = fig.add_gridspec(2, 2, hspace=0.24, wspace=0.16)

    ax = fig.add_subplot(gs[0, 0])
    add_panel_label(ax, "A")
    add_image(ax, RESULTS / "primary_core_aim2_niche_pathway_model.png", "Aim 2: niche/pathway readiness")

    ax = fig.add_subplot(gs[0, 1])
    add_panel_label(ax, "B")
    add_image(
        ax,
        RESULTS / "aim2_stage_window_model.png",
        "Fitted stage-window model: TGF/BDNF timing differs by branch",
    )

    ax = fig.add_subplot(gs[1, 0])
    add_panel_label(ax, "C")
    add_image(ax, RESULTS / "aim3_empirical_calibration.png", "Aim 3: empirical sparse-coding calibration")

    ax = fig.add_subplot(gs[1, 1])
    draw_figure5d_final_model(ax, "D")

    fig.suptitle("Figure 5. Stage-windowed niche, computation, and final convergence model", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure5_stage_window_sparse_coding_final_model",
        "Stage-windowed niche, computation, and final convergence model",
        [
            rel(RESULTS / "primary_core_aim2_niche_pathway_model.png"),
            rel(RESULTS / "aim2_stage_window_model.png"),
            rel(RESULTS / "aim3_empirical_calibration.png"),
            rel(RESULTS / "primary_core_full_transcriptome_diffusion.md"),
            rel(RESULTS / "dandi_000003_multisession_spatial_extension.md"),
            rel(RESULTS / "specific_aims_completion_audit.md"),
        ],
        "Integrates pathway-readiness scoring, fitted stage-window timing, empirical sparse-coding calibration, and the revised resource-constrained final model.",
        rows,
    )


def figure6(rows: list[dict[str, str]]) -> None:
    fig = plt.figure(figsize=(17.2, 10.6))
    ax = fig.add_subplot(111)
    draw_figure6_working_model(ax, "A")

    fig.suptitle("Graphical abstract. Integrated working model for granule-cell convergence", fontsize=16, fontweight="bold")
    save_figure(
        fig,
        "figure6_integrated_working_model",
        "Graphical abstract source: integrated working model for granule-cell convergence",
        [
            rel(RESULTS / "integrated_primary_core_datasets.tsv"),
            rel(RESULTS / "primary_core_manuscript_candidate_tiers.tsv"),
            rel(RESULTS / "primary_core_transcriptomic_configuration_model.png"),
            rel(RESULTS / "primary_core_full_transcriptome_diffusion.md"),
            rel(RESULTS / "aim2_stage_window_model.md"),
            rel(RESULTS / "primary_core_aim3_sparse_coding_model.md"),
            rel(RESULTS / "aim3_empirical_calibration.md"),
            rel(RESULTS / "neuromorpho_granule_morphometry_validation.md"),
            rel(RESULTS / "dandi_000003_multisession_spatial_extension.md"),
        ],
        "Integrates the Fig1 primary-core frame, Fig2 candidate tiers, Fig4 configuration model, Fig5 fitted stage-window/computational model, and morphology/activity validation into one working hypothesis.",
        rows,
    )


def standalone_flow_panels(rows: list[dict[str, str]]) -> None:
    core = pd.read_csv(RESULTS / "integrated_primary_core_datasets.tsv", sep="\t")
    species_counts = core.groupby("species")["dataset"].nunique()
    top, heat_cols, axis_summary = load_figure2_data()
    panel_specs = [
        (
            "figure1a_biological_puzzle_panel",
            "Figure 1A standalone biological-puzzle panel",
            (5.4, 4.2),
            draw_figure1a_biological_puzzle,
            [rel(RESULTS / "primary_core_niche_circuit_module_model.md")],
            "Standalone conceptual export for Figure 1A.",
        ),
        (
            "figure1b_working_model_panel",
            "Figure 1B standalone working-model panel",
            (5.6, 4.2),
            draw_figure1b_working_model,
            [rel(RESULTS / "primary_core_niche_circuit_module_model.md")],
            "Standalone orthogonal flow-chart export for Figure 1B.",
        ),
        (
            "figure1c_primary_core_counts_panel",
            "Figure 1C standalone primary-core counts panel",
            (5.2, 4.0),
            lambda ax, label: draw_figure1c_primary_core_counts(ax, core, label),
            [rel(RESULTS / "integrated_primary_core_datasets.tsv")],
            "Standalone bar-chart export for Figure 1C.",
        ),
        (
            "figure1d_primary_dataset_map_panel",
            "Figure 1D standalone primary-dataset map panel",
            (9.4, 5.4),
            lambda ax, label: draw_figure1d_dataset_map(ax, core, species_counts, label),
            [rel(RESULTS / "integrated_primary_core_datasets.tsv")],
            "Standalone cleaned dataset-map export for Figure 1D.",
        ),
        (
            "figure2a_ortholog_workflow_panel",
            "Figure 2A standalone ortholog-aware workflow panel",
            (5.4, 5.0),
            draw_figure2a_workflow,
            [rel(RESULTS / "primary_core_mgi_ortholog_formal_rank_model.md")],
            "Standalone orthogonal flow-chart export for Figure 2A.",
        ),
        (
            "figure2b_candidate_branch_deltas_panel",
            "Figure 2B standalone candidate branch-deltas panel",
            (6.4, 5.4),
            lambda ax, label: draw_figure2b_candidate_deltas(ax, top, heat_cols, label),
            [rel(RESULTS / "primary_core_manuscript_candidate_tiers.tsv")],
            "Standalone heatmap export for Figure 2B.",
        ),
        (
            "figure2c_mechanism_axis_panel",
            "Figure 2C standalone mechanism-axis panel",
            (10.8, 5.8),
            lambda ax, label: draw_figure2c_mechanism_axes(ax, axis_summary, label),
            [rel(RESULTS / "primary_core_mechanism_axis_summary.tsv")],
            "Standalone mechanism-axis export for Figure 2C.",
        ),
        (
            "figure4a_configuration_score_panel",
            "Figure 4A standalone transcriptomic-configuration panel",
            (5.6, 4.4),
            draw_figure4a_configuration,
            [rel(RESULTS / "primary_core_transcriptomic_configuration_model.md")],
            "Standalone orthogonal flow-chart export for Figure 4A.",
        ),
        (
            "figure5d_revised_final_model_panel",
            "Figure 5D standalone revised final-model panel",
            (6.1, 5.1),
            draw_figure5d_final_model,
            [
                rel(RESULTS / "aim2_stage_window_model.md"),
                rel(RESULTS / "aim3_empirical_calibration.md"),
            ],
            "Standalone orthogonal flow-chart export for Figure 5D.",
        ),
        (
            "figure6_integrated_working_model_panel",
            "Graphical abstract standalone integrated working-model panel",
            (12.8, 7.6),
            draw_figure6_working_model,
            [
                rel(RESULTS / "integrated_primary_core_datasets.tsv"),
                rel(RESULTS / "primary_core_manuscript_candidate_tiers.tsv"),
                rel(RESULTS / "primary_core_full_transcriptome_diffusion.md"),
                rel(RESULTS / "aim2_stage_window_model.md"),
                rel(RESULTS / "aim3_empirical_calibration.md"),
                rel(RESULTS / "neuromorpho_granule_morphometry_validation.md"),
                rel(RESULTS / "dandi_000003_multisession_spatial_extension.md"),
            ],
            "Standalone orthogonal flow-chart export for the graphical abstract/source working model.",
        ),
    ]
    for name, title, figsize, drawer, sources, purpose in panel_specs:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        result = drawer(ax, None)
        if name == "figure2b_candidate_branch_deltas_panel":
            fig.colorbar(result, ax=ax, fraction=0.046, pad=0.035, label="Median candidate-background rank delta")
        save_figure(fig, name, title, sources, purpose, rows)


def write_report(manifest: pd.DataFrame) -> None:
    lines = [
        "# Manuscript Figure Assembly",
        "",
        "Date updated: 2026-06-24",
        "",
        "## Purpose",
        "",
        "This packet assembles draft manuscript-scale composite figures and source panels, including primary-core diffusion/pseudotime refinement, the fitted stage-window model, empirical sparse-coding calibration, and morphology/activity validation.",
        "",
        "## Figure Set",
        "",
    ]
    for _, row in manifest.iterrows():
        lines.append(f"- `{row['figure_id']}`: {row['title']}.")
        lines.append(f"  - PNG: `{row['png']}`")
        lines.append(f"  - PDF: `{row['pdf']}`")
        lines.append(f"  - Purpose: {row['purpose']}")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- These are draft assembly figures, not final journal-polished layouts.",
            "- Panels that reuse existing analysis plots preserve the current quantitative outputs.",
            "- Conceptual panels were redrawn from project conclusions rather than copied from source publications.",
            "- Figure 5 now uses the fitted stage-window model and the empirical sparse-coding calibration.",
            "- The integrated working-model panel is retained as the graphical-abstract/source panel; final Figure 6 is the focused sender-receiver ligand-receptor figure in the manuscript figure folder.",
            "- Standalone cleaned panel exports are included for Figure 1A, Figure 1B, Figure 1C, Figure 1D, Figure 2A, Figure 2B, Figure 2C, Figure 4A, Figure 5D, and the graphical abstract.",
            "",
            "## Manifest",
            "",
            f"- `{rel(OUT_MANIFEST)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    figure1(rows)
    figure2(rows)
    figure3(rows)
    figure4(rows)
    figure5(rows)
    figure6(rows)
    standalone_flow_panels(rows)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)
    write_report(manifest)
    print(f"Wrote {rel(OUT_MD)}")
    print(manifest[["figure_id", "png"]].to_string(index=False))


if __name__ == "__main__":
    main()
