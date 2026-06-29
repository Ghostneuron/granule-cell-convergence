#!/usr/bin/env python3
"""Dataset-aware analysis of refined candidate granule-cell programs."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CALLS_IN = Path(os.environ.get("CANDIDATE_CALLS_INPUT", RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"))
OUTPUT_PREFIX = os.environ.get("DATASET_LEVEL_PREFIX", "refined_dataset_level_granule_program")
UNITS_OUT = RESULTS / f"{OUTPUT_PREFIX}_units.tsv"
STATS_OUT = RESULTS / f"{OUTPUT_PREFIX}_statistics.tsv"
LEAVE_ONE_OUT = RESULTS / f"{OUTPUT_PREFIX}_leave_one_dataset_out.tsv"
FIG_OUT = RESULTS / f"{OUTPUT_PREFIX}_identity_structural_units.png"

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


PLOT_CLASSES = [
    "dentate_candidate",
    "cerebellar_candidate",
    "known_non_dentate_reference",
    "cerebellum_warning",
    "other_or_ambiguous",
]

COLORS = {
    "dentate_candidate": "#2a9d8f",
    "cerebellar_candidate": "#7b2cbf",
    "known_non_dentate_reference": "#8d99ae",
    "cerebellum_warning": "#d08c60",
    "other_or_ambiguous": "#c7c7c7",
}


def analysis_class(call: str) -> str:
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


def load_calls() -> pd.DataFrame:
    with gzip.open(CALLS_IN, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["analysis_class"] = df["candidate_call"].map(analysis_class)
    return df


def summarize_units(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["dataset", "sample", "group", "region", "platform", "species", "analysis_class"]
    rows = []
    for key, sub in df.groupby(group_cols, dropna=False):
        row = dict(zip(group_cols, key))
        row["n_cells_or_spots"] = len(sub)
        row["n_high_confidence"] = int((sub["call_confidence"] == "high").sum())
        row["n_medium_confidence"] = int((sub["call_confidence"] == "medium").sum())
        row["n_reference_confidence"] = int((sub["call_confidence"] == "reference").sum())
        for col in NUMERIC_COLS:
            vals = sub[col].dropna()
            row[f"{col}_mean"] = vals.mean() if len(vals) else np.nan
            row[f"{col}_median"] = vals.median() if len(vals) else np.nan
            row[f"{col}_q25"] = vals.quantile(0.25) if len(vals) else np.nan
            row[f"{col}_q75"] = vals.quantile(0.75) if len(vals) else np.nan
        rows.append(row)
    units = pd.DataFrame(rows)
    units["unit_id"] = (
        units["dataset"].astype(str)
        + "|"
        + units["sample"].astype(str)
        + "|"
        + units["group"].astype(str)
        + "|"
        + units["analysis_class"].astype(str)
    )
    return units


def compare_units(units: pd.DataFrame, class_a: str, class_b: str, metric: str, note: str) -> dict[str, str]:
    col = f"{metric}_median"
    a = units.loc[units["analysis_class"] == class_a, col].dropna().to_numpy()
    b = units.loc[units["analysis_class"] == class_b, col].dropna().to_numpy()
    if len(a) == 0 or len(b) == 0:
        return {
            "comparison": f"{class_a}_vs_{class_b}",
            "metric": metric,
            "n_units_a": str(len(a)),
            "n_units_b": str(len(b)),
            "median_of_unit_medians_a": "nan",
            "median_of_unit_medians_b": "nan",
            "difference_a_minus_b": "nan",
            "mann_whitney_p": "nan",
            "note": note,
        }
    pval = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return {
        "comparison": f"{class_a}_vs_{class_b}",
        "metric": metric,
        "n_units_a": str(len(a)),
        "n_units_b": str(len(b)),
        "median_of_unit_medians_a": f"{np.median(a):.6g}",
        "median_of_unit_medians_b": f"{np.median(b):.6g}",
        "difference_a_minus_b": f"{np.median(a) - np.median(b):.6g}",
        "mann_whitney_p": f"{pval:.6g}",
        "note": note,
    }


def sign_test_row(units: pd.DataFrame, cls: str, metric: str, direction: str, threshold: float, note: str) -> dict[str, str]:
    vals = units.loc[units["analysis_class"] == cls, f"{metric}_median"].dropna().to_numpy()
    if direction == "greater":
        successes = int(np.sum(vals > threshold))
    elif direction == "less":
        successes = int(np.sum(vals < threshold))
    else:
        raise ValueError(direction)
    n = len(vals)
    pval = stats.binomtest(successes, n, 0.5, alternative="greater").pvalue if n else np.nan
    return {
        "comparison": f"{cls}_sign_test_{direction}_{threshold}",
        "metric": metric,
        "n_units_a": str(n),
        "n_units_b": "",
        "median_of_unit_medians_a": f"{np.median(vals):.6g}" if n else "nan",
        "median_of_unit_medians_b": str(threshold),
        "difference_a_minus_b": f"{np.median(vals) - threshold:.6g}" if n else "nan",
        "mann_whitney_p": f"{pval:.6g}" if n else "nan",
        "note": note + f"; successes={successes}/{n}",
    }


def leave_one_dataset_out(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for excluded in ["none"] + sorted(units["dataset"].dropna().unique()):
        sub = units if excluded == "none" else units[units["dataset"] != excluded]
        dentate = sub.loc[sub["analysis_class"] == "dentate_candidate", "identity_contrast_median"].dropna()
        cerebellar = sub.loc[sub["analysis_class"] == "cerebellar_candidate", "identity_contrast_median"].dropna()
        dentate_struct = sub.loc[sub["analysis_class"] == "dentate_candidate", "structural_rank_median"].dropna()
        cerebellar_struct = sub.loc[sub["analysis_class"] == "cerebellar_candidate", "structural_rank_median"].dropna()
        rows.append(
            {
                "excluded_dataset": excluded,
                "n_dentate_units": len(dentate),
                "n_cerebellar_units": len(cerebellar),
                "dentate_identity_contrast_median": f"{dentate.median():.6g}" if len(dentate) else "nan",
                "cerebellar_identity_contrast_median": f"{cerebellar.median():.6g}" if len(cerebellar) else "nan",
                "identity_median_difference": f"{dentate.median() - cerebellar.median():.6g}" if len(dentate) and len(cerebellar) else "nan",
                "dentate_structural_rank_median": f"{dentate_struct.median():.6g}" if len(dentate_struct) else "nan",
                "cerebellar_structural_rank_median": f"{cerebellar_struct.median():.6g}" if len(cerebellar_struct) else "nan",
                "both_structural_medians_above_0_5": str(
                    bool(len(dentate_struct) and len(cerebellar_struct) and dentate_struct.median() > 0.5 and cerebellar_struct.median() > 0.5)
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_units(units: pd.DataFrame) -> None:
    plot_df = units[units["analysis_class"].isin(PLOT_CLASSES)].copy()
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    for cls in PLOT_CLASSES:
        sub = plot_df[plot_df["analysis_class"] == cls]
        if sub.empty:
            continue
        sizes = 20 + np.clip(np.sqrt(sub["n_cells_or_spots"].astype(float)), 2, 85)
        ax.scatter(
            sub["identity_contrast_median"],
            sub["structural_rank_median"],
            s=sizes,
            c=COLORS[cls],
            label=cls.replace("_", " "),
            alpha=0.82,
            edgecolors="white",
            linewidths=0.6,
        )
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.axhline(0.5, color="#777777", linewidth=0.9, linestyle=":")
    ax.set_xlabel("Unit median: dentate identity minus cerebellar identity")
    ax.set_ylabel("Unit median: within-sample structural-program rank")
    ax.set_title("Dataset-level candidate granule-cell program units")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    ax.grid(True, linewidth=0.4, color="#d9d9d9", alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(FIG_OUT, dpi=180)


def main() -> None:
    df = load_calls()
    units = summarize_units(df)
    units.to_csv(UNITS_OUT, sep="\t", index=False, float_format="%.6g")

    stat_rows = [
        compare_units(
            units,
            "dentate_candidate",
            "cerebellar_candidate",
            "identity_contrast",
            "Dataset-level identity separation; positive difference supports distinct regional identity.",
        ),
        compare_units(
            units,
            "dentate_candidate",
            "known_non_dentate_reference",
            "identity_contrast",
            "Dataset-level dentate candidate separation from curated non-dentate references.",
        ),
        compare_units(
            units,
            "cerebellar_candidate",
            "other_or_ambiguous",
            "cerebellar_identity",
            "Dataset-level cerebellar enrichment relative to ambiguous cells/spots.",
        ),
        compare_units(
            units,
            "dentate_candidate",
            "cerebellar_candidate",
            "structural_rank",
            "Structural rank difference is not expected to separate identity; both groups being above 0.5 supports partial convergence.",
        ),
        sign_test_row(
            units,
            "dentate_candidate",
            "identity_contrast",
            "greater",
            0.0,
            "Dentate candidate units should have positive identity contrast.",
        ),
        sign_test_row(
            units,
            "cerebellar_candidate",
            "identity_contrast",
            "less",
            0.0,
            "Cerebellar candidate units should have negative identity contrast.",
        ),
        sign_test_row(
            units,
            "dentate_candidate",
            "structural_rank",
            "greater",
            0.5,
            "Dentate candidate units above median within-sample structural rank support structural program activity.",
        ),
        sign_test_row(
            units,
            "cerebellar_candidate",
            "structural_rank",
            "greater",
            0.5,
            "Cerebellar candidate units above median within-sample structural rank support structural program activity.",
        ),
    ]
    pd.DataFrame(stat_rows).to_csv(STATS_OUT, sep="\t", index=False)

    loo = leave_one_dataset_out(units)
    loo.to_csv(LEAVE_ONE_OUT, sep="\t", index=False)
    plot_units(units)
    print(f"Wrote {UNITS_OUT}")
    print(f"Wrote {STATS_OUT}")
    print(f"Wrote {LEAVE_ONE_OUT}")
    print(f"Wrote {FIG_OUT}")


if __name__ == "__main__":
    main()
