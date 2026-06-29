#!/usr/bin/env python3
"""Compare broad and cluster-supported GSE322785 epigenomic contrasts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

BROAD = RESULTS / "gse322785_human_h5_epigenomic_marker_group_contrasts.tsv"
STRICT = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_contrasts.tsv"

OUT_COMPARISON = RESULTS / "gse322785_epigenomic_broad_vs_cluster_supported_contrasts.tsv"
OUT_ROBUST = RESULTS / "gse322785_epigenomic_robust_positive_contrasts.tsv"
OUT_SUMMARY = RESULTS / "gse322785_epigenomic_robust_contrast_summary.tsv"
OUT_MD = RESULTS / "gse322785_epigenomic_robust_contrast_summary.md"

KEYS = ["comparator", "feature_type", "model_term_supported", "target_set", "peak_category"]

MIN_POSITIVE_DELTA = 0.10
MIN_STRONG_DELTA = 0.25


def pooled(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    p = df.loc[df["donor_id"].eq("pooled")].copy()
    keep = KEYS + [
        "n_donors",
        "mean_delta_score",
        "median_delta_score",
        "mean_delta_detection",
        "median_delta_detection",
    ]
    p = p[keep].copy()
    return p.rename(
        columns={
            "n_donors": f"{prefix}_n_donors",
            "mean_delta_score": f"{prefix}_mean_delta_score",
            "median_delta_score": f"{prefix}_median_delta_score",
            "mean_delta_detection": f"{prefix}_mean_delta_detection",
            "median_delta_detection": f"{prefix}_median_delta_detection",
        }
    )


def classify(row: pd.Series) -> str:
    broad = row["broad_mean_delta_score"]
    strict = row["cluster_supported_mean_delta_score"]
    if pd.isna(broad) or pd.isna(strict):
        return "missing_layer"
    if broad >= MIN_POSITIVE_DELTA and strict >= MIN_POSITIVE_DELTA:
        if broad >= MIN_STRONG_DELTA and strict >= MIN_STRONG_DELTA:
            return "robust_positive_strong"
        return "robust_positive"
    if broad >= MIN_POSITIVE_DELTA and strict < -MIN_POSITIVE_DELTA:
        return "discordant_broad_positive_strict_negative"
    if strict >= MIN_POSITIVE_DELTA and broad < -MIN_POSITIVE_DELTA:
        return "discordant_strict_positive_broad_negative"
    if broad >= MIN_POSITIVE_DELTA:
        return "broad_only_positive"
    if strict >= MIN_POSITIVE_DELTA:
        return "cluster_supported_only_positive"
    if broad <= -MIN_POSITIVE_DELTA and strict <= -MIN_POSITIVE_DELTA:
        return "robust_negative"
    return "weak_or_neutral"


def main() -> None:
    broad = pooled(BROAD, "broad")
    strict = pooled(STRICT, "cluster_supported")
    merged = broad.merge(strict, on=KEYS, how="outer")
    merged["delta_shift_strict_minus_broad"] = (
        merged["cluster_supported_mean_delta_score"] - merged["broad_mean_delta_score"]
    )
    merged["detection_shift_strict_minus_broad"] = (
        merged["cluster_supported_mean_delta_detection"] - merged["broad_mean_delta_detection"]
    )
    merged["concordance_class"] = merged.apply(classify, axis=1)
    merged["robust_rank_score"] = np.minimum(
        merged["broad_mean_delta_score"].fillna(-np.inf),
        merged["cluster_supported_mean_delta_score"].fillna(-np.inf),
    )
    merged["robust_detection_score"] = np.minimum(
        merged["broad_mean_delta_detection"].fillna(-np.inf),
        merged["cluster_supported_mean_delta_detection"].fillna(-np.inf),
    )
    merged = merged.sort_values(
        ["concordance_class", "robust_rank_score", "robust_detection_score"],
        ascending=[True, False, False],
    )
    merged.to_csv(OUT_COMPARISON, sep="\t", index=False)

    robust = merged.loc[merged["concordance_class"].isin(["robust_positive", "robust_positive_strong"])].copy()
    robust = robust.sort_values(["robust_rank_score", "robust_detection_score"], ascending=False)
    robust.to_csv(OUT_ROBUST, sep="\t", index=False)

    summary = (
        merged.groupby(["concordance_class", "comparator", "feature_type", "model_term_supported"], dropna=False)
        .agg(
            n_contrasts=("target_set", "size"),
            mean_broad_delta=("broad_mean_delta_score", "mean"),
            mean_cluster_supported_delta=("cluster_supported_mean_delta_score", "mean"),
            mean_shift_strict_minus_broad=("delta_shift_strict_minus_broad", "mean"),
            max_robust_rank_score=("robust_rank_score", "max"),
        )
        .reset_index()
        .sort_values(["concordance_class", "n_contrasts", "max_robust_rank_score"], ascending=[True, False, False])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    class_counts = merged["concordance_class"].value_counts().to_dict()
    robust_by_comp = robust.groupby("comparator").size().sort_values(ascending=False) if not robust.empty else pd.Series(dtype=int)
    robust_by_term = (
        robust.groupby(["feature_type", "model_term_supported"]).size().sort_values(ascending=False).head(10)
        if not robust.empty
        else pd.Series(dtype=int)
    )
    top = robust.head(8)
    top_lines = [
        "- `{}` / `{}` / `{}` / `{}` / `{}`: broad {:.3g}, strict {:.3g}.".format(
            row.comparator,
            row.feature_type,
            row.model_term_supported,
            row.target_set,
            row.peak_category,
            float(row.broad_mean_delta_score),
            float(row.cluster_supported_mean_delta_score),
        )
        for row in top.itertuples()
    ]
    comp_lines = [f"- `{idx}`: {int(value)} robust-positive contrasts." for idx, value in robust_by_comp.items()]
    term_lines = [
        f"- `{idx[0]}` / `{idx[1]}`: {int(value)} robust-positive contrasts."
        for idx, value in robust_by_term.items()
    ]
    class_lines = [f"- `{key}`: {int(value)} contrasts." for key, value in sorted(class_counts.items())]

    lines = [
        "# GSE322785 Epigenomic Robust Contrast Summary",
        "",
        "Date built: 2026-06-26",
        "",
        "## Rule",
        "",
        f"Robust-positive contrasts require both broad and cluster-supported pooled granule-minus-comparator mean delta scores >= {MIN_POSITIVE_DELTA}. Strong robust-positive contrasts require both layers >= {MIN_STRONG_DELTA}.",
        "",
        "## Concordance Classes",
        "",
        *class_lines,
        "",
        "## Robust-Positive Counts by Comparator",
        "",
        *(comp_lines or ["- No robust-positive contrasts."]),
        "",
        "## Robust-Positive Counts by Term",
        "",
        *(term_lines or ["- No robust-positive contrasts."]),
        "",
        "## Top Robust-Positive Contrasts",
        "",
        *(top_lines or ["- No robust-positive contrasts."]),
        "",
        "## Interpretation",
        "",
        "This table identifies GSE322785 epigenomic signals that are stable to the stricter cluster-supported sensitivity filter. These robust contrasts are stronger candidates for discussion than broad-only or strict-only effects, while still remaining provisional until source taxonomy or full multimodal clustering is available.",
        "",
        "## Outputs",
        "",
        f"- Full comparison: `{OUT_COMPARISON.relative_to(ROOT)}`",
        f"- Robust-positive contrasts: `{OUT_ROBUST.relative_to(ROOT)}`",
        f"- Summary: `{OUT_SUMMARY.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
