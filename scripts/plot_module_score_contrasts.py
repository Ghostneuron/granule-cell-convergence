#!/usr/bin/env python3
"""Plot identity contrast versus structural module score."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CONTRAST_IN = RESULTS / "module_score_identity_structural_contrasts.tsv"
PLOT_OUT = RESULTS / "module_score_identity_structural_scatter.png"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {
    "curated_dentate_granule_lineage": "#0077b6",
    "dentate_granule_enriched": "#2a9d8f",
    "cerebellum_mixed_or_spatial": "#7b2cbf",
    "organoid_mixed": "#b56576",
    "mixed_or_reference": "#8d99ae",
}

LABEL_ROWS = {
    ("GSE104323", "10X_all_cells", "GC-adult"),
    ("GSE104323", "10X_all_cells", "Immature-GC"),
    ("GSE104323", "10X_all_cells", "Neuroblast"),
    ("GSE122357", "GSM3464549_P0", "all"),
    ("GSE165657", "Cerebellum_aggr", "all"),
    ("GSE214905", "patch_RNA_QC_counts", "all"),
    ("GSE292261", "SS2_filtered_counts", "DG_P5"),
    ("GSE292261", "SS2_filtered_counts", "DG_P28"),
    ("GSE312658", "Ctrl", "all"),
}


def short_label(row: dict[str, str]) -> str:
    dataset = row["dataset"]
    sample = row["sample"]
    group = row["group"]
    if dataset == "GSE104323":
        return group
    if dataset == "GSE292261":
        return group
    if dataset == "GSE122357":
        return f"{dataset} P0"
    if dataset == "GSE312658":
        return "GSE312658 Ctrl"
    return dataset


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    with CONTRAST_IN.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    for context, color in COLORS.items():
        subset = [row for row in rows if row["context_note"] == context]
        if not subset:
            continue
        xs = [float(row["identity_contrast_dentate_minus_cerebellar"]) for row in subset]
        ys = [float(row["structural_program_mean_log1p"]) for row in subset]
        sizes = [22 + min(float(row["n_cells_or_spots"]), 5000) ** 0.5 * 2.0 for row in subset]
        ax.scatter(xs, ys, s=sizes, c=color, alpha=0.78, edgecolors="white", linewidths=0.6, label=context.replace("_", " "))

    for row in rows:
        key = (row["dataset"], row["sample"], row["group"])
        if key not in LABEL_ROWS:
            continue
        x = float(row["identity_contrast_dentate_minus_cerebellar"])
        y = float(row["structural_program_mean_log1p"])
        offset = (-8, 4) if x > 1.8 else (5, 4)
        ha = "right" if x > 1.8 else "left"
        ax.annotate(short_label(row), (x, y), xytext=offset, textcoords="offset points", fontsize=8, ha=ha)

    ax.axvline(0, color="#555555", linewidth=1, linestyle="--")
    ax.set_xlabel("Dentate identity minus cerebellar identity, mean log1p score")
    ax.set_ylabel("Structural program mean log1p score")
    ax.set_title("Granule-cell identity and shared structural-program modules")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, linewidth=0.4, color="#d9d9d9", alpha=0.7)
    fig.tight_layout(rect=[0, 0, 0.8, 1])
    fig.savefig(PLOT_OUT, dpi=180)
    print(f"Wrote {PLOT_OUT}")


if __name__ == "__main__":
    main()
