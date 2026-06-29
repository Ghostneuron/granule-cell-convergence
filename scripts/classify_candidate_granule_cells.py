#!/usr/bin/env python3
"""Classify candidate granule-like cells from per-cell marker module scores.

The calls are intentionally conservative and should be treated as candidate
labels for triage, not final cell-type annotations.
"""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
SCORES_IN = Path(os.environ.get("CANDIDATE_SCORE_INPUT", RESULTS / "per_cell_marker_module_scores.tsv.gz"))
OUTPUT_PREFIX = os.environ.get("CANDIDATE_CALL_PREFIX", "candidate_granule_cell")
CALLS_OUT = RESULTS / f"{OUTPUT_PREFIX}_calls.tsv.gz"
SUMMARY_OUT = RESULTS / f"{OUTPUT_PREFIX}_call_summary.tsv"

KNOWN_DG_GROUPS = {"GC-adult", "GC-juv", "Immature-GC", "Neuroblast"}

SCORE_PANELS = [
    "dentate_identity",
    "cerebellar_identity",
    "shared_granule_neuronal",
    "morphogenesis_cytoskeleton",
    "axon_guidance_synapse",
]


def percentile_rank(series: pd.Series) -> pd.Series:
    if len(series) <= 1:
        return pd.Series([1.0] * len(series), index=series.index)
    return series.rank(method="average", pct=True)


def call_candidate(row: pd.Series) -> tuple[str, str, str]:
    dataset = row["dataset"]
    region = row["region"]
    group = row["group"]
    dentate = float(row["dentate_identity"])
    cerebellar = float(row["cerebellar_identity"])
    shared_rank = float(row["shared_rank"])
    structural_rank = float(row["structural_rank"])
    dentate_rank = float(row["dentate_rank"])
    cerebellar_rank = float(row["cerebellar_rank"])
    contrast = float(row["identity_contrast"])

    if dataset == "GSE104323" and group in KNOWN_DG_GROUPS:
        return ("candidate_dentate_granule", "reference", "curated GSE104323 dentate granule lineage group")
    if dataset == "GSE104323":
        return ("known_non_dentate_reference", "reference", "curated GSE104323 non-dentate-granule metadata group")

    if region == "dentate_gyrus":
        if contrast > 0 and dentate_rank >= 0.55 and shared_rank >= 0.35:
            confidence = "high" if structural_rank >= 0.60 else "medium"
            return ("candidate_dentate_granule", confidence, "dentate identity higher with neuronal/structural support")
        if dentate > cerebellar:
            return ("dentate_like_low_support", "low", "dentate identity higher but below rank thresholds")
        return ("non_granule_or_ambiguous", "low", "dentate sample without dentate identity dominance")

    if region in {"cerebellum", "cerebellum_spatial"}:
        if contrast < 0 and cerebellar_rank >= 0.55 and shared_rank >= 0.35:
            confidence = "high" if structural_rank >= 0.60 else "medium"
            return ("candidate_cerebellar_granule", confidence, "cerebellar identity higher with neuronal/structural support")
        if contrast > 0 and dentate_rank >= 0.70:
            return ("cerebellum_dentate_panel_high_warning", "low", "cerebellar source has higher dentate-panel score; panel needs refinement")
        return ("non_granule_or_ambiguous", "low", "cerebellar source without strong cerebellar granule candidate pattern")

    if region == "organoid":
        if max(dentate_rank, cerebellar_rank) >= 0.70 and shared_rank >= 0.40:
            return ("organoid_granule_like_candidate", "low", "organoid score resembles granule-like program but metadata are incomplete")
        return ("non_granule_or_ambiguous", "low", "organoid source without strong candidate pattern")

    return ("non_granule_or_ambiguous", "low", "unsupported region for granule candidate call")


def main() -> None:
    usecols = [
        "dataset",
        "sample",
        "cell_id",
        "group",
        "species",
        "region",
        "platform",
        "panel",
        "genes_detected_in_cell",
        "detection_fraction_panel",
        "mean_log1p_expression_panel",
        "source_path",
    ]
    long = pd.read_csv(SCORES_IN, sep="\t", usecols=usecols)
    long["mean_log1p_expression_panel"] = pd.to_numeric(long["mean_log1p_expression_panel"], errors="coerce").fillna(0.0)
    long["genes_detected_in_cell"] = pd.to_numeric(long["genes_detected_in_cell"], errors="coerce").fillna(0).astype(int)
    long["detection_fraction_panel"] = pd.to_numeric(long["detection_fraction_panel"], errors="coerce").fillna(0.0)

    index_cols = ["dataset", "sample", "cell_id", "group", "species", "region", "platform", "source_path"]
    score_wide = long.pivot_table(
        index=index_cols,
        columns="panel",
        values="mean_log1p_expression_panel",
        aggfunc="first",
        fill_value=0.0,
    ).reset_index()
    detected_wide = long.pivot_table(
        index=index_cols,
        columns="panel",
        values="genes_detected_in_cell",
        aggfunc="first",
        fill_value=0,
    ).reset_index()

    for panel in SCORE_PANELS:
        if panel not in score_wide:
            score_wide[panel] = 0.0
        if panel not in detected_wide:
            detected_wide[panel] = 0

    calls = score_wide.copy()
    calls["identity_contrast"] = calls["dentate_identity"] - calls["cerebellar_identity"]
    calls["structural_program_mean"] = calls[["shared_granule_neuronal", "morphogenesis_cytoskeleton", "axon_guidance_synapse"]].mean(axis=1)
    calls["structural_detected_genes"] = (
        detected_wide["shared_granule_neuronal"]
        + detected_wide["morphogenesis_cytoskeleton"]
        + detected_wide["axon_guidance_synapse"]
    )

    group_cols = ["dataset", "sample"]
    calls["dentate_rank"] = calls.groupby(group_cols)["dentate_identity"].transform(percentile_rank)
    calls["cerebellar_rank"] = calls.groupby(group_cols)["cerebellar_identity"].transform(percentile_rank)
    calls["shared_rank"] = calls.groupby(group_cols)["shared_granule_neuronal"].transform(percentile_rank)
    calls["structural_rank"] = calls.groupby(group_cols)["structural_program_mean"].transform(percentile_rank)

    call_tuples = calls.apply(call_candidate, axis=1)
    calls["candidate_call"] = [item[0] for item in call_tuples]
    calls["call_confidence"] = [item[1] for item in call_tuples]
    calls["call_reason"] = [item[2] for item in call_tuples]

    out_cols = [
        "dataset",
        "sample",
        "cell_id",
        "group",
        "species",
        "region",
        "platform",
        "candidate_call",
        "call_confidence",
        "call_reason",
        "dentate_identity",
        "cerebellar_identity",
        "identity_contrast",
        "shared_granule_neuronal",
        "morphogenesis_cytoskeleton",
        "axon_guidance_synapse",
        "structural_program_mean",
        "structural_detected_genes",
        "dentate_rank",
        "cerebellar_rank",
        "shared_rank",
        "structural_rank",
        "source_path",
    ]
    with gzip.open(CALLS_OUT, "wt", newline="") as fh:
        calls.to_csv(fh, sep="\t", index=False, columns=out_cols, float_format="%.6g")

    summary = (
        calls.groupby(["dataset", "sample", "group", "region", "platform", "candidate_call", "call_confidence"], dropna=False)
        .agg(
            n_cells_or_spots=("cell_id", "size"),
            mean_dentate_identity=("dentate_identity", "mean"),
            mean_cerebellar_identity=("cerebellar_identity", "mean"),
            mean_identity_contrast=("identity_contrast", "mean"),
            mean_structural_program=("structural_program_mean", "mean"),
            median_structural_rank=("structural_rank", "median"),
        )
        .reset_index()
    )
    total_by_group = summary.groupby(["dataset", "sample", "group"])["n_cells_or_spots"].transform("sum")
    summary["fraction_of_group"] = summary["n_cells_or_spots"] / total_by_group
    summary = summary.sort_values(["dataset", "sample", "group", "candidate_call", "call_confidence"])
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)

    print(f"Wrote {len(calls)} candidate call rows to {CALLS_OUT}")
    print(f"Wrote {len(summary)} summary rows to {SUMMARY_OUT}")


if __name__ == "__main__":
    main()
