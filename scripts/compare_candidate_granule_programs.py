#!/usr/bin/env python3
"""Compare candidate dentate and cerebellar granule-cell programs."""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CALLS_IN = Path(os.environ.get("CANDIDATE_CALLS_INPUT", RESULTS / "candidate_granule_cell_calls.tsv.gz"))
OUTPUT_PREFIX = os.environ.get("PROGRAM_COMPARE_PREFIX", "candidate_granule_program")
CLASS_SUMMARY_OUT = RESULTS / f"{OUTPUT_PREFIX}_class_summary.tsv"
STATS_OUT = RESULTS / f"{OUTPUT_PREFIX}_statistics.tsv"
FIG_OUT = RESULTS / f"{OUTPUT_PREFIX}_comparison.png"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


NUMERIC_COLS = [
    "dentate_identity",
    "cerebellar_identity",
    "identity_contrast",
    "shared_granule_neuronal",
    "morphogenesis_cytoskeleton",
    "axon_guidance_synapse",
    "structural_program_mean",
    "dentate_rank",
    "cerebellar_rank",
    "shared_rank",
    "structural_rank",
]


def analysis_class(row: pd.Series) -> str:
    call = row["candidate_call"]
    if call == "candidate_dentate_granule":
        return "dentate_candidate"
    if call == "candidate_cerebellar_granule":
        return "cerebellar_candidate"
    if call == "known_non_dentate_reference":
        return "known_non_dentate_reference"
    if call == "cerebellum_dentate_panel_high_warning":
        return "cerebellum_warning"
    if call == "dentate_like_low_support":
        return "dentate_low_support"
    if call == "organoid_granule_like_candidate":
        return "organoid_granule_like"
    return "other_or_ambiguous"


def cliff_delta(x: np.ndarray, y: np.ndarray, max_n: int = 5000) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    if len(x) > max_n:
        x = np.random.default_rng(17).choice(x, max_n, replace=False)
    if len(y) > max_n:
        y = np.random.default_rng(23).choice(y, max_n, replace=False)
    diffs = x[:, None] - y[None, :]
    return float((np.sum(diffs > 0) - np.sum(diffs < 0)) / diffs.size)


def summarize_group(df: pd.DataFrame, label: str) -> dict[str, str]:
    out = {"analysis_class": label, "n_cells_or_spots": str(len(df))}
    for col in NUMERIC_COLS:
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        out[f"{col}_mean"] = f"{values.mean():.6g}" if len(values) else "nan"
        out[f"{col}_median"] = f"{values.median():.6g}" if len(values) else "nan"
        out[f"{col}_iqr"] = f"{(values.quantile(0.75) - values.quantile(0.25)):.6g}" if len(values) else "nan"
    datasets = df[["dataset", "sample"]].drop_duplicates()
    out["n_dataset_samples"] = str(len(datasets))
    out["dataset_samples"] = ";".join(f"{row.dataset}:{row.sample}" for row in datasets.itertuples(index=False))
    return out


def mann_whitney_row(df: pd.DataFrame, class_a: str, class_b: str, metric: str, interpretation: str) -> dict[str, str]:
    a = pd.to_numeric(df.loc[df["analysis_class"] == class_a, metric], errors="coerce").dropna().to_numpy()
    b = pd.to_numeric(df.loc[df["analysis_class"] == class_b, metric], errors="coerce").dropna().to_numpy()
    if len(a) == 0 or len(b) == 0:
        return {
            "test": "mann_whitney_u",
            "class_a": class_a,
            "class_b": class_b,
            "metric": metric,
            "n_a": str(len(a)),
            "n_b": str(len(b)),
            "median_a": "nan",
            "median_b": "nan",
            "effect_delta_a_minus_b": "nan",
            "p_value": "nan",
            "interpretation": interpretation,
        }
    test = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "test": "mann_whitney_u",
        "class_a": class_a,
        "class_b": class_b,
        "metric": metric,
        "n_a": str(len(a)),
        "n_b": str(len(b)),
        "median_a": f"{np.median(a):.6g}",
        "median_b": f"{np.median(b):.6g}",
        "effect_delta_a_minus_b": f"{cliff_delta(a, b):.6g}",
        "p_value": f"{test.pvalue:.6g}",
        "interpretation": interpretation,
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_comparison(df: pd.DataFrame) -> None:
    plot_classes = [
        "dentate_candidate",
        "cerebellar_candidate",
        "known_non_dentate_reference",
        "other_or_ambiguous",
        "cerebellum_warning",
    ]
    plot_df = df[df["analysis_class"].isin(plot_classes)].copy()
    colors = {
        "dentate_candidate": "#2a9d8f",
        "cerebellar_candidate": "#7b2cbf",
        "known_non_dentate_reference": "#8d99ae",
        "other_or_ambiguous": "#c7c7c7",
        "cerebellum_warning": "#d08c60",
    }

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    metrics = [
        ("identity_contrast", "Dentate minus cerebellar identity"),
        ("structural_rank", "Within-sample structural-program rank"),
    ]
    for ax, (metric, title) in zip(axes, metrics):
        data = [plot_df.loc[plot_df["analysis_class"] == cls, metric].dropna().to_numpy() for cls in plot_classes]
        parts = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
        for body, cls in zip(parts["bodies"], plot_classes):
            body.set_facecolor(colors[cls])
            body.set_edgecolor("white")
            body.set_alpha(0.78)
        parts["cmedians"].set_color("#222222")
        ax.axhline(0, color="#555555", linewidth=0.8, linestyle="--") if metric == "identity_contrast" else None
        ax.set_xticks(range(1, len(plot_classes) + 1))
        ax.set_xticklabels([cls.replace("_", "\n") for cls in plot_classes], fontsize=8)
        ax.set_title(title)
        ax.grid(True, axis="y", linewidth=0.4, alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Candidate granule-cell program comparison", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=180, bbox_inches="tight")


def main() -> None:
    with gzip.open(CALLS_IN, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["analysis_class"] = df.apply(analysis_class, axis=1)

    class_order = [
        "dentate_candidate",
        "cerebellar_candidate",
        "known_non_dentate_reference",
        "dentate_low_support",
        "cerebellum_warning",
        "organoid_granule_like",
        "other_or_ambiguous",
    ]
    summary_rows = [summarize_group(df[df["analysis_class"] == cls], cls) for cls in class_order if (df["analysis_class"] == cls).any()]
    write_tsv(CLASS_SUMMARY_OUT, summary_rows)

    stat_rows = [
        mann_whitney_row(
            df,
            "dentate_candidate",
            "cerebellar_candidate",
            "identity_contrast",
            "Positive effect means dentate candidates are dentate-identity shifted relative to cerebellar candidates.",
        ),
        mann_whitney_row(
            df,
            "dentate_candidate",
            "known_non_dentate_reference",
            "identity_contrast",
            "Positive effect means dentate candidates separate from curated non-dentate reference cells.",
        ),
        mann_whitney_row(
            df,
            "cerebellar_candidate",
            "other_or_ambiguous",
            "cerebellar_identity",
            "Positive effect means cerebellar candidates are cerebellar-identity enriched relative to ambiguous cells.",
        ),
        mann_whitney_row(
            df,
            "dentate_candidate",
            "known_non_dentate_reference",
            "structural_rank",
            "Positive effect means dentate candidates are structurally enriched relative to curated non-dentate references.",
        ),
        mann_whitney_row(
            df,
            "cerebellar_candidate",
            "other_or_ambiguous",
            "structural_rank",
            "Positive effect means cerebellar candidates are structurally enriched relative to ambiguous cells/spots.",
        ),
        mann_whitney_row(
            df,
            "dentate_candidate",
            "cerebellar_candidate",
            "structural_rank",
            "Small effect with both medians high would support convergence of structural-program rank.",
        ),
    ]
    write_tsv(STATS_OUT, stat_rows)
    plot_comparison(df)
    print(f"Wrote {CLASS_SUMMARY_OUT}")
    print(f"Wrote {STATS_OUT}")
    print(f"Wrote {FIG_OUT}")


if __name__ == "__main__":
    main()
