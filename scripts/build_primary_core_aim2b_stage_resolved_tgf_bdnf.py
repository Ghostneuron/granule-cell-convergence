#!/usr/bin/env python3
"""Stage-resolved TGF-beta/BDNF audit for Specific Aim 2.

This add-on asks whether the historical TGF-beta/BDNF maturation/stop
mechanism is stage-dependent in dentate and cerebellar granule-cell systems.
It deliberately scores within each dataset/axis using rank-normalized
pseudobulk expression, because raw expression scales are not comparable across
10x, Smart-seq2, bulked single-cell count tables, and published summaries.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_primary_core_aim2_niche_pathway_model import PATHWAY_MODULES, SIGNATURES  # noqa: E402


GSE104323_EXPR = EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz"
GSE104323_META = EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz"
GSE292261_EXPR = EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz"
GSE292261_META = EXTERNAL / "GEO/GSE292261/GSE292261_sample_data_SS2_filtered.csv.gz"
GSE214309_EXPR = EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz"
REFINED_CALLS = RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"
FULL_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"

OUT_GENE_UNITS = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_gene_units.tsv"
OUT_PATHWAY_UNITS = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_pathway_units.tsv"
OUT_SIGNATURE_UNITS = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_signature_units.tsv"
OUT_TRANSITIONS = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_transitions.tsv"
OUT_SUMMARY = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_summary.tsv"
OUT_PLOT = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_plot.png"
OUT_MD = RESULTS / "primary_core_aim2b_stage_tgf_bdnf.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


csv.field_size_limit(sys.maxsize)


DENTATE_LINEAGE_ORDER = {
    "RGL_young": 0,
    "RGL": 1,
    "nIPC": 2,
    "nIPC-perin": 3,
    "Neuroblast": 4,
    "Immature-GC": 5,
    "GC-juv": 6,
    "GC-adult": 7,
}

GSE214309_ORDER = {
    "immature_1hr": 0,
    "immatureactive_1hr": 1,
    "mature_1hr": 2,
    "matureactive_1hr": 3,
    "immature_4hr": 4,
    "immatureactive_4hr": 5,
    "mature_4hr": 6,
    "matureactive_4hr": 7,
}

SIGNATURE_IDS = [
    "tgf_bdnf_2005_index",
    "differentiation_stop_index",
    "neurogenic_permissive_index",
    "stop_minus_permissive_index",
]

SIGNATURE_COLORS = {
    "tgf_bdnf_2005_index": "#31588c",
    "differentiation_stop_index": "#b35806",
    "neurogenic_permissive_index": "#1b7837",
    "stop_minus_permissive_index": "#6a3d9a",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').strip("'").upper()


def finite_mean(values: list[float] | pd.Series) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.mean(arr))


def finite_median(values: list[float] | pd.Series) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def pathway_gene_sets() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module in PATHWAY_MODULES:
        for order, (gene, role) in enumerate(module["genes"].items(), start=1):
            rows.append(
                {
                    "pathway_id": module["pathway_id"],
                    "pathway_label": module["pathway_label"],
                    "pathway_family": module["pathway_family"],
                    "hypothesis_role": module["hypothesis_role"],
                    "canonical_gene": canon(gene),
                    "gene": gene,
                    "gene_role": role,
                    "gene_order": order,
                }
            )
    return pd.DataFrame(rows)


def target_gene_set() -> set[str]:
    return set(pathway_gene_sets()["canonical_gene"])


def summarize_group_values(
    *,
    dataset: str,
    species: str,
    region: str,
    axis_type: str,
    axis_label: str,
    axis_order: float,
    comparison_group: str,
    n_cells: int,
    source_path: str,
    source_gene_symbol: str,
    canonical_gene: str,
    values: list[float],
) -> dict[str, object]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        mean_expr = np.nan
        mean_log = np.nan
        detection = np.nan
        nonzero = 0
    else:
        mean_expr = float(np.mean(arr))
        mean_log = float(np.mean(np.log1p(arr)))
        nonzero = int(np.sum(arr > 0))
        detection = float(nonzero / arr.size)
    return {
        "dataset": dataset,
        "species": species,
        "region": region,
        "axis_type": axis_type,
        "axis_label": axis_label,
        "axis_order": axis_order,
        "comparison_group": comparison_group,
        "n_cells": int(n_cells),
        "canonical_gene": canonical_gene,
        "source_gene_symbol": source_gene_symbol,
        "nonzero_cells": nonzero,
        "detection_fraction": detection,
        "mean_expression": mean_expr,
        "mean_log1p_expression": mean_log,
        "source_path": source_path,
    }


def add_within_axis_gene_ranks(gene_units: pd.DataFrame) -> pd.DataFrame:
    ranked = gene_units.copy()
    ranked["gene_rank_within_axis"] = np.nan
    group_cols = ["dataset", "axis_type", "comparison_group", "canonical_gene"]
    for _, idx in ranked.groupby(group_cols, sort=False).groups.items():
        sub = ranked.loc[idx, "mean_log1p_expression"]
        ranked.loc[idx, "gene_rank_within_axis"] = sub.rank(method="average", pct=True)
    return ranked


def load_gse104323_lineage_units(target_genes: set[str]) -> pd.DataFrame:
    meta = pd.read_csv(GSE104323_META, sep="\t", dtype=str)
    cell_col = "Sample name (24185 single cells)"
    cluster_col = "characteristics: cell cluster"
    meta = meta.loc[meta[cluster_col].isin(DENTATE_LINEAGE_ORDER)].copy()
    cell_to_group = dict(zip(meta[cell_col].astype(str), meta[cluster_col].astype(str)))

    rows: list[dict[str, object]] = []
    with open_text(GSE104323_EXPR) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cell_ids = [str(cell).strip() for cell in header[1:]]
        selected = [(idx + 1, cell, cell_to_group[cell]) for idx, cell in enumerate(cell_ids) if cell in cell_to_group]
        group_cells = Counter(group for _, _, group in selected)
        selected_by_group: dict[str, list[int]] = defaultdict(list)
        for idx, _, group in selected:
            selected_by_group[group].append(idx)

        for row in reader:
            if not row:
                continue
            gene = str(row[0]).strip()
            cg = canon(gene)
            if cg not in target_genes:
                continue
            for group, col_indices in selected_by_group.items():
                vals: list[float] = []
                for col_idx in col_indices:
                    try:
                        vals.append(float(row[col_idx]))
                    except (IndexError, ValueError):
                        vals.append(0.0)
                rows.append(
                    summarize_group_values(
                        dataset="GSE104323",
                        species="Mus musculus",
                        region="dentate_gyrus",
                        axis_type="adult_dentate_lineage_state",
                        axis_label=group,
                        axis_order=DENTATE_LINEAGE_ORDER[group],
                        comparison_group="curated_dentate_granule_lineage",
                        n_cells=group_cells[group],
                        source_path="External_Data/GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
                        source_gene_symbol=gene,
                        canonical_gene=cg,
                        values=vals,
                    )
                )
    return pd.DataFrame(rows)


def gse292261_candidate_cells() -> set[str]:
    cols = ["dataset", "sample", "cell_id", "candidate_call"]
    calls = pd.read_csv(REFINED_CALLS, sep="\t", usecols=cols, dtype=str, low_memory=False)
    sub = calls.loc[
        calls["dataset"].eq("GSE292261") & calls["candidate_call"].eq("candidate_dentate_granule")
    ].copy()
    return set(sub["cell_id"].astype(str))


def load_gse292261_stage_units(target_genes: set[str]) -> pd.DataFrame:
    meta = pd.read_csv(GSE292261_META, dtype=str)
    cell_col = meta.columns[0]
    meta = meta.rename(columns={cell_col: "cell_id"})
    meta["stage"] = meta["Sample"].astype(str)
    stage_order = {f"DG_P{age}": float(age) for age in [5, 7, 10, 15, 28]}
    meta = meta.loc[meta["stage"].isin(stage_order)].copy()
    candidate_cells = gse292261_candidate_cells()

    header = pd.read_csv(GSE292261_EXPR, nrows=0).columns.tolist()
    gene_cols = [col for col in header[1:] if canon(col) in target_genes]
    df = pd.read_csv(GSE292261_EXPR, usecols=[header[0], *gene_cols], low_memory=False)
    df = df.rename(columns={header[0]: "cell_id"})
    df["cell_id"] = df["cell_id"].astype(str)
    df = df.merge(meta[["cell_id", "stage"]], on="cell_id", how="inner")
    df["is_candidate"] = df["cell_id"].isin(candidate_cells)

    rows: list[dict[str, object]] = []
    for comparison_group, sub_df in [
        ("all_DG_cells", df),
        ("candidate_dentate_granule_only", df.loc[df["is_candidate"]].copy()),
    ]:
        for stage, stage_df in sub_df.groupby("stage", sort=False):
            for gene_col in gene_cols:
                vals = pd.to_numeric(stage_df[gene_col], errors="coerce").fillna(0.0).astype(float).tolist()
                rows.append(
                    summarize_group_values(
                        dataset="GSE292261",
                        species="Mus musculus",
                        region="dentate_gyrus",
                        axis_type="postnatal_dentate_age",
                        axis_label=stage.replace("DG_", ""),
                        axis_order=stage_order[stage],
                        comparison_group=comparison_group,
                        n_cells=len(stage_df),
                        source_path="External_Data/GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
                        source_gene_symbol=gene_col,
                        canonical_gene=canon(gene_col),
                        values=vals,
                    )
                )
    return pd.DataFrame(rows)


def load_gse214309_state_units(target_genes: set[str]) -> pd.DataFrame:
    cols = ["dataset", "sample", "cell_id", "group"]
    calls = pd.read_csv(REFINED_CALLS, sep="\t", usecols=cols, dtype=str, low_memory=False)
    sub = calls.loc[calls["dataset"].eq("GSE214309") & calls["group"].isin(GSE214309_ORDER)].copy()
    cell_to_group = dict(zip(sub["cell_id"].astype(str), sub["group"].astype(str)))

    rows: list[dict[str, object]] = []
    with open_text(GSE214309_EXPR) as fh:
        reader = csv.reader(fh, delimiter=",")
        header = next(reader)
        cell_ids = [str(cell).strip() for cell in header]
        selected = [(idx + 1, cell, cell_to_group[cell]) for idx, cell in enumerate(cell_ids) if cell in cell_to_group]
        selected_by_group: dict[str, list[int]] = defaultdict(list)
        group_cells = Counter(group for _, _, group in selected)
        for idx, _, group in selected:
            selected_by_group[group].append(idx)

        for row in reader:
            if not row:
                continue
            gene = str(row[0]).strip()
            cg = canon(gene)
            if cg not in target_genes:
                continue
            for group, col_indices in selected_by_group.items():
                vals: list[float] = []
                for col_idx in col_indices:
                    try:
                        vals.append(float(row[col_idx]))
                    except (IndexError, ValueError):
                        vals.append(0.0)
                rows.append(
                    summarize_group_values(
                        dataset="GSE214309",
                        species="Mus musculus",
                        region="dentate_gyrus",
                        axis_type="adult_dentate_activity_maturation_state",
                        axis_label=group,
                        axis_order=GSE214309_ORDER[group],
                        comparison_group="adult_DGC_maturation_activity_state",
                        n_cells=group_cells[group],
                        source_path="External_Data/GEO/GSE214309/GSE214309_counts.txt.gz",
                        source_gene_symbol=gene,
                        canonical_gene=cg,
                        values=vals,
                    )
                )
    return pd.DataFrame(rows)


def load_gse122357_cerebellar_units(target_genes: set[str]) -> pd.DataFrame:
    use_cols = [
        "dataset",
        "sample",
        "broad_class",
        "n_cells",
        "canonical_gene",
        "source_gene_symbol",
        "detection_fraction",
        "mean_expression",
        "mean_log1p_expression",
        "eligible_class",
        "source_path",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(FULL_EXPR, sep="\t", usecols=use_cols, chunksize=100_000, low_memory=False):
        chunk["canonical_gene"] = chunk["canonical_gene"].map(canon)
        sub = chunk.loc[
            chunk["dataset"].eq("GSE122357")
            & chunk["broad_class"].eq("cerebellar_candidate")
            & chunk["canonical_gene"].isin(target_genes)
            & chunk["eligible_class"].astype(str).str.lower().isin({"true", "1", "yes"})
        ].copy()
        if not sub.empty:
            pieces.append(sub)
    if not pieces:
        return pd.DataFrame()
    expr = pd.concat(pieces, ignore_index=True)
    sample_map = {
        "GSM3464549_P0": ("P0", 0.0),
        "GSM3464550_P8a": ("P8a", 8.1),
        "GSM3464551_P8b": ("P8b", 8.2),
    }
    expr = expr.loc[expr["sample"].isin(sample_map)].copy()
    expr["axis_label"] = expr["sample"].map(lambda value: sample_map[value][0])
    expr["axis_order"] = expr["sample"].map(lambda value: sample_map[value][1])
    return pd.DataFrame(
        {
            "dataset": "GSE122357",
            "species": "Mus musculus",
            "region": "cerebellum",
            "axis_type": "postnatal_cerebellar_age",
            "axis_label": expr["axis_label"],
            "axis_order": expr["axis_order"],
            "comparison_group": "candidate_cerebellar_granule_only",
            "n_cells": expr["n_cells"].astype(int),
            "canonical_gene": expr["canonical_gene"],
            "source_gene_symbol": expr["source_gene_symbol"],
            "nonzero_cells": np.nan,
            "detection_fraction": pd.to_numeric(expr["detection_fraction"], errors="coerce"),
            "mean_expression": pd.to_numeric(expr["mean_expression"], errors="coerce"),
            "mean_log1p_expression": pd.to_numeric(expr["mean_log1p_expression"], errors="coerce"),
            "source_path": expr["source_path"],
        }
    )


def build_pathway_units(gene_units: pd.DataFrame, gene_sets: pd.DataFrame) -> pd.DataFrame:
    module_cols = [
        "canonical_gene",
        "pathway_id",
        "pathway_label",
        "pathway_family",
        "hypothesis_role",
        "gene_role",
    ]
    expr = gene_units.merge(gene_sets[module_cols].drop_duplicates(), on="canonical_gene", how="inner")
    group_cols = [
        "dataset",
        "species",
        "region",
        "axis_type",
        "axis_label",
        "axis_order",
        "comparison_group",
        "pathway_id",
        "pathway_label",
        "pathway_family",
        "hypothesis_role",
    ]
    units = (
        expr.groupby(group_cols, sort=False)
        .agg(
            n_cells=("n_cells", "max"),
            n_genes_present=("canonical_gene", "nunique"),
            genes_present=("canonical_gene", lambda values: ",".join(sorted(set(values)))),
            gene_roles_present=("gene_role", lambda values: ",".join(sorted(set(values)))),
            median_gene_rank=("gene_rank_within_axis", "median"),
            mean_gene_rank=("gene_rank_within_axis", "mean"),
            median_detection_fraction=("detection_fraction", "median"),
            median_mean_log1p_expression=("mean_log1p_expression", "median"),
        )
        .reset_index()
    )
    defined = gene_sets.groupby("pathway_id")["canonical_gene"].nunique().to_dict()
    units["n_genes_defined"] = units["pathway_id"].map(defined).astype(int)
    units["pathway_gene_coverage"] = units["n_genes_present"] / units["n_genes_defined"]
    units = units.sort_values(["dataset", "comparison_group", "axis_order", "pathway_id"])
    return units


def build_signature_units(pathway_units: pd.DataFrame) -> pd.DataFrame:
    idx_cols = [
        "dataset",
        "species",
        "region",
        "axis_type",
        "axis_label",
        "axis_order",
        "comparison_group",
    ]
    pivot = pathway_units.pivot_table(
        index=idx_cols,
        columns="pathway_id",
        values="median_gene_rank",
        aggfunc="median",
    ).reset_index()
    pivot.columns.name = None
    rows: list[dict[str, object]] = []
    for _, row in pivot.iterrows():
        base = {col: row[col] for col in idx_cols}
        computed: dict[str, float] = {}
        for sig in SIGNATURES:
            sig_id = sig["signature_id"]
            if "pathways" in sig:
                values = [row[p] for p in sig["pathways"] if p in row.index and pd.notna(row[p])]
                score = finite_mean(values)
                n_pathways = len(values)
            else:
                score = computed.get("differentiation_stop_index", np.nan) - computed.get(
                    "neurogenic_permissive_index", np.nan
                )
                n_pathways = 2
            computed[sig_id] = score
            rows.append(
                {
                    **base,
                    "signature_id": sig_id,
                    "signature_label": sig["signature_label"],
                    "signature_score": score,
                    "n_pathways_present": n_pathways,
                }
            )
    out = pd.DataFrame(rows)
    out = out.sort_values(["dataset", "comparison_group", "axis_order", "signature_id"])
    return out


def add_transition(
    rows: list[dict[str, object]],
    table: pd.DataFrame,
    *,
    metric_type: str,
    metric_id: str,
    label_col: str,
    score_col: str,
    dataset: str,
    comparison_group: str,
    start_label: str,
    end_label: str,
    transition_id: str,
    transition_label: str,
) -> None:
    sub = table.loc[
        table["dataset"].eq(dataset)
        & table["comparison_group"].eq(comparison_group)
        & table[label_col].eq(metric_id)
        & table["axis_label"].isin([start_label, end_label])
    ].copy()
    if set(sub["axis_label"]) != {start_label, end_label}:
        return
    start = float(sub.loc[sub["axis_label"].eq(start_label), score_col].iloc[0])
    end = float(sub.loc[sub["axis_label"].eq(end_label), score_col].iloc[0])
    meta = sub.iloc[0]
    rows.append(
        {
            "dataset": dataset,
            "region": meta["region"],
            "axis_type": meta["axis_type"],
            "comparison_group": comparison_group,
            "transition_id": transition_id,
            "transition_label": transition_label,
            "metric_type": metric_type,
            "metric_id": metric_id,
            "metric_label": meta.get("signature_label", meta.get("pathway_label", metric_id)),
            "start_label": start_label,
            "end_label": end_label,
            "start_score": start,
            "end_score": end,
            "delta_end_minus_start": end - start,
            "direction": "increased" if end > start else "decreased" if end < start else "unchanged",
        }
    )


def build_transitions(pathway_units: pd.DataFrame, signature_units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    comparisons = [
        ("GSE104323", "curated_dentate_granule_lineage", "RGL_young", "GC-adult", "rgl_to_gc_adult", "RGL_young to GC-adult"),
        ("GSE104323", "curated_dentate_granule_lineage", "Neuroblast", "Immature-GC", "neuroblast_to_immature_gc", "Neuroblast to Immature-GC"),
        ("GSE104323", "curated_dentate_granule_lineage", "Immature-GC", "GC-adult", "immature_gc_to_adult_gc", "Immature-GC to GC-adult"),
        ("GSE292261", "all_DG_cells", "P5", "P28", "p5_to_p28_all", "P5 to P28 all DG cells"),
        ("GSE292261", "candidate_dentate_granule_only", "P5", "P15", "p5_to_p15_candidate", "P5 to P15 candidate dentate GCs"),
        ("GSE292261", "candidate_dentate_granule_only", "P15", "P28", "p15_to_p28_candidate", "P15 to P28 candidate dentate GCs"),
        ("GSE122357", "candidate_cerebellar_granule_only", "P0", "P8a", "p0_to_p8a_cerebellum", "P0 to P8a cerebellar candidates"),
        ("GSE122357", "candidate_cerebellar_granule_only", "P0", "P8b", "p0_to_p8b_cerebellum", "P0 to P8b cerebellar candidates"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "immature_1hr", "mature_1hr", "immature_to_mature_1hr", "Immature to mature, 1 hr"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "immature_4hr", "mature_4hr", "immature_to_mature_4hr", "Immature to mature, 4 hr"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "immature_1hr", "immatureactive_1hr", "immature_activity_1hr", "Immature activity, 1 hr"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "mature_1hr", "matureactive_1hr", "mature_activity_1hr", "Mature activity, 1 hr"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "immature_4hr", "immatureactive_4hr", "immature_activity_4hr", "Immature activity, 4 hr"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "mature_4hr", "matureactive_4hr", "mature_activity_4hr", "Mature activity, 4 hr"),
    ]

    key_pathways = ["tgf_beta_smad", "bdnf_trkb_mapk", "shh_granule_expansion", "notch_hes"]
    for comp in comparisons:
        dataset, group, start, end, transition_id, transition_label = comp
        for sig_id in SIGNATURE_IDS:
            add_transition(
                rows,
                signature_units,
                metric_type="signature",
                metric_id=sig_id,
                label_col="signature_id",
                score_col="signature_score",
                dataset=dataset,
                comparison_group=group,
                start_label=start,
                end_label=end,
                transition_id=transition_id,
                transition_label=transition_label,
            )
        for pathway_id in key_pathways:
            add_transition(
                rows,
                pathway_units,
                metric_type="pathway",
                metric_id=pathway_id,
                label_col="pathway_id",
                score_col="median_gene_rank",
                dataset=dataset,
                comparison_group=group,
                start_label=start,
                end_label=end,
                transition_id=transition_id,
                transition_label=transition_label,
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["dataset", "comparison_group", "transition_id", "metric_type", "metric_id"])
    return out


def build_summary(signature_units: pd.DataFrame, transitions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (dataset, axis_type, comparison_group, signature_id), sub in signature_units.groupby(
        ["dataset", "axis_type", "comparison_group", "signature_id"], sort=False
    ):
        sub = sub.sort_values("axis_order")
        first = sub.iloc[0]
        last = sub.iloc[-1]
        values = pd.to_numeric(sub["signature_score"], errors="coerce")
        rows.append(
            {
                "dataset": dataset,
                "region": first["region"],
                "axis_type": axis_type,
                "comparison_group": comparison_group,
                "signature_id": signature_id,
                "signature_label": first["signature_label"],
                "n_axis_units": int(len(sub)),
                "start_label": first["axis_label"],
                "end_label": last["axis_label"],
                "start_score": float(first["signature_score"]),
                "end_score": float(last["signature_score"]),
                "delta_end_minus_start": float(last["signature_score"] - first["signature_score"]),
                "min_score": float(values.min()),
                "max_score": float(values.max()),
                "peak_label": sub.iloc[int(values.to_numpy().argmax())]["axis_label"] if values.notna().any() else "",
                "mean_score": float(values.mean()),
            }
        )
    summary = pd.DataFrame(rows)

    if not transitions.empty:
        tgf = transitions.loc[
            transitions["metric_type"].eq("signature") & transitions["metric_id"].eq("tgf_bdnf_2005_index")
        ]
        for (dataset, comparison_group), sub in tgf.groupby(["dataset", "comparison_group"], sort=False):
            rows.append(
                {
                    "dataset": dataset,
                    "region": sub["region"].iloc[0],
                    "axis_type": sub["axis_type"].iloc[0],
                    "comparison_group": comparison_group,
                    "signature_id": "tgf_bdnf_2005_index_transition_set",
                    "signature_label": "TGF-beta/BDNF 2005 mechanism transition set",
                    "n_axis_units": int(len(sub)),
                    "start_label": ";".join(sub["start_label"].astype(str)),
                    "end_label": ";".join(sub["end_label"].astype(str)),
                    "start_score": np.nan,
                    "end_score": np.nan,
                    "delta_end_minus_start": finite_median(sub["delta_end_minus_start"]),
                    "min_score": float(pd.to_numeric(sub["delta_end_minus_start"], errors="coerce").min()),
                    "max_score": float(pd.to_numeric(sub["delta_end_minus_start"], errors="coerce").max()),
                    "peak_label": "",
                    "mean_score": float(pd.to_numeric(sub["delta_end_minus_start"], errors="coerce").mean()),
                }
            )
        summary = pd.DataFrame(rows)
    return summary.sort_values(["dataset", "comparison_group", "signature_id"])


def plot_signature_trajectories(signature_units: pd.DataFrame) -> None:
    panels = [
        ("GSE104323", "curated_dentate_granule_lineage", "Adult dentate lineage"),
        ("GSE292261", "candidate_dentate_granule_only", "Postnatal dentate candidates"),
        ("GSE122357", "candidate_cerebellar_granule_only", "Postnatal cerebellar candidates"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "Adult dentate maturation/activity"),
    ]
    label_map = {
        "immature_1hr": "immature\n1 hr",
        "immatureactive_1hr": "imm active\n1 hr",
        "mature_1hr": "mature\n1 hr",
        "matureactive_1hr": "mat active\n1 hr",
        "immature_4hr": "immature\n4 hr",
        "immatureactive_4hr": "imm active\n4 hr",
        "mature_4hr": "mature\n4 hr",
        "matureactive_4hr": "mat active\n4 hr",
        "RGL_young": "RGL\nyoung",
        "Immature-GC": "Immature\nGC",
        "GC-adult": "GC\nadult",
        "GC-juv": "GC\njuv",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for ax, (dataset, group, title) in zip(axes, panels):
        sub = signature_units.loc[
            signature_units["dataset"].eq(dataset)
            & signature_units["comparison_group"].eq(group)
            & signature_units["signature_id"].isin(SIGNATURE_IDS)
        ].copy()
        if sub.empty:
            ax.axis("off")
            ax.set_title(title)
            continue
        sub = sub.sort_values("axis_order")
        for sig_id in SIGNATURE_IDS:
            sig = sub.loc[sub["signature_id"].eq(sig_id)].copy()
            x_labels = [label_map.get(str(label), str(label)) for label in sig["axis_label"]]
            ax.plot(
                x_labels,
                sig["signature_score"],
                marker="o",
                linewidth=2,
                markersize=4,
                color=SIGNATURE_COLORS.get(sig_id, "#333333"),
                label=sig["signature_label"].iloc[0],
            )
        ax.axhline(0.5, color="#777777", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.set_title(title)
        ax.set_ylabel("Within-axis rank score")
        ax.set_ylim(-0.15, 1.05)
        ax.tick_params(axis="x", rotation=0, labelsize=9)
        ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.025))
    fig.suptitle("Aim 2b stage-resolved TGF-beta/BDNF and niche-pathway readiness", fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.16, hspace=0.40, wspace=0.22)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def fmt(value: object, digits: int = 3) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(val):
        return "NA"
    return f"{val:.{digits}f}"


def signature_lookup(signature_units: pd.DataFrame, dataset: str, group: str, sig_id: str) -> pd.DataFrame:
    return signature_units.loc[
        signature_units["dataset"].eq(dataset)
        & signature_units["comparison_group"].eq(group)
        & signature_units["signature_id"].eq(sig_id)
    ].sort_values("axis_order")


def transition_delta(transitions: pd.DataFrame, dataset: str, transition_id: str, metric_id: str) -> str:
    sub = transitions.loc[
        transitions["dataset"].eq(dataset)
        & transitions["transition_id"].eq(transition_id)
        & transitions["metric_id"].eq(metric_id)
        & transitions["metric_type"].eq("signature")
    ]
    if sub.empty:
        return "NA"
    return fmt(sub["delta_end_minus_start"].iloc[0])


def write_markdown(
    gene_units: pd.DataFrame,
    pathway_units: pd.DataFrame,
    signature_units: pd.DataFrame,
    transitions: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    g104_tgf = signature_lookup(signature_units, "GSE104323", "curated_dentate_granule_lineage", "tgf_bdnf_2005_index")
    g292_tgf = signature_lookup(signature_units, "GSE292261", "candidate_dentate_granule_only", "tgf_bdnf_2005_index")
    g122_tgf = signature_lookup(signature_units, "GSE122357", "candidate_cerebellar_granule_only", "tgf_bdnf_2005_index")
    g214_tgf = signature_lookup(signature_units, "GSE214309", "adult_DGC_maturation_activity_state", "tgf_bdnf_2005_index")

    def peak_text(df: pd.DataFrame) -> str:
        if df.empty:
            return "NA"
        row = df.loc[df["signature_score"].idxmax()]
        return f"{row['axis_label']} ({fmt(row['signature_score'])})"

    p28_cells = gene_units.loc[
        gene_units["dataset"].eq("GSE292261")
        & gene_units["comparison_group"].eq("candidate_dentate_granule_only")
        & gene_units["axis_label"].eq("P28")
    ]["n_cells"]
    p28_note = int(p28_cells.iloc[0]) if not p28_cells.empty else 0

    lines = [
        "# Aim 2b Stage-Resolved TGF-beta/BDNF Audit",
        "",
        "## Question",
        "",
        "Does the historical TGF-beta/BDNF maturation or stop mechanism behave differently across granule-cell developmental stage, especially given that cerebellar granule cells develop mostly postnatally while dentate granule cells retain adult neurogenesis?",
        "",
        "## Scope",
        "",
        "This is a pathway-readiness audit, not a direct niche sender-receiver assay. Scores are percentile ranks within each dataset and axis for each gene, then summarized at pathway and signature level. Therefore values test whether a module is relatively higher at a stage/state inside the same dataset, not whether raw expression is larger across platforms.",
        "",
        "Datasets included:",
        "",
        "- `GSE104323`: adult mouse dentate lineage states from RGL/RGL_young through neuroblast, immature GC, juvenile GC, and adult GC.",
        "- `GSE292261`: mouse postnatal dentate stages P5, P7, P10, P15, and P28, scored both as all DG cells and candidate dentate granule cells.",
        "- `GSE122357`: mouse cerebellar candidate granule cells at P0 and P8 replicates.",
        "- `GSE214309`: adult mouse dentate granule-cell maturation/activity states, including immature, mature, active immature, and active mature cells at 1 hr and 4 hr.",
        "",
        "## Main Findings",
        "",
        f"- Dentate lineage (`GSE104323`) TGF-beta/BDNF peaks at {peak_text(g104_tgf)}; RGL_young to GC-adult delta is {transition_delta(transitions, 'GSE104323', 'rgl_to_gc_adult', 'tgf_bdnf_2005_index')}, while neuroblast to immature-GC delta is {transition_delta(transitions, 'GSE104323', 'neuroblast_to_immature_gc', 'tgf_bdnf_2005_index')}.",
        f"- Postnatal dentate candidates (`GSE292261`) TGF-beta/BDNF peaks at {peak_text(g292_tgf)}. The P15 to P28 candidate-only delta is {transition_delta(transitions, 'GSE292261', 'p15_to_p28_candidate', 'tgf_bdnf_2005_index')}, but P28 candidate-cell count is only {p28_note}, so P28 must be treated as a low-support endpoint.",
        f"- Cerebellar candidates (`GSE122357`) TGF-beta/BDNF peaks at {peak_text(g122_tgf)}; P0 to P8a delta is {transition_delta(transitions, 'GSE122357', 'p0_to_p8a_cerebellum', 'tgf_bdnf_2005_index')} and P0 to P8b delta is {transition_delta(transitions, 'GSE122357', 'p0_to_p8b_cerebellum', 'tgf_bdnf_2005_index')}.",
        f"- Adult dentate activity/maturation (`GSE214309`) TGF-beta/BDNF peaks at {peak_text(g214_tgf)}. Immature-to-mature deltas are {transition_delta(transitions, 'GSE214309', 'immature_to_mature_1hr', 'tgf_bdnf_2005_index')} at 1 hr and {transition_delta(transitions, 'GSE214309', 'immature_to_mature_4hr', 'tgf_bdnf_2005_index')} at 4 hr.",
        "",
        "## Interpretation",
        "",
        "The 2005 TGF-beta/BDNF mechanism should not be framed as a simple cerebellum-versus-dentate regional effect. In the current primary-core data, it is better framed as a stage- and state-sensitive maturation/readiness axis. Dentate evidence is strongest because adult lineage and postnatal stages can be directly ordered, while cerebellar evidence is limited to P0 versus P8 candidate granule cells in this first pass.",
        "",
        "This supports the project hypothesis in a more precise form: cerebellar and dentate granule cells can converge on similar morphology through shared downstream maturation and wiring modules, but the timing and upstream niche logic are region-specific. Dentate retains a lifelong progenitor-to-granule continuum, whereas cerebellum has a more developmentally bounded expansion-and-differentiation program.",
        "",
        "## Manuscript Use",
        "",
        "- Add Aim 2b as a stage-aware refinement under the niche pathway aim.",
        "- Phrase the TGF-beta/BDNF result as `stage-dependent maturation/readiness`, not as a universal granule-cell stop switch.",
        "- Use `GSE104323`, `GSE292261`, and `GSE214309` as dentate neurogenesis anchors; use `GSE122357` as the cerebellar postnatal comparator.",
        "- Keep direct ligand-source claims for a future spatial or ligand-receptor sender analysis.",
        "",
        "## Outputs",
        "",
        f"- Gene units: `{OUT_GENE_UNITS.relative_to(ROOT)}` ({len(gene_units):,} rows).",
        f"- Pathway units: `{OUT_PATHWAY_UNITS.relative_to(ROOT)}` ({len(pathway_units):,} rows).",
        f"- Signature units: `{OUT_SIGNATURE_UNITS.relative_to(ROOT)}` ({len(signature_units):,} rows).",
        f"- Transitions: `{OUT_TRANSITIONS.relative_to(ROOT)}` ({len(transitions):,} rows).",
        f"- Summary: `{OUT_SUMMARY.relative_to(ROOT)}` ({len(summary):,} rows).",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    gene_sets = pathway_gene_sets()
    targets = target_gene_set()

    pieces = [
        load_gse104323_lineage_units(targets),
        load_gse292261_stage_units(targets),
        load_gse122357_cerebellar_units(targets),
        load_gse214309_state_units(targets),
    ]
    gene_units = pd.concat([piece for piece in pieces if not piece.empty], ignore_index=True)
    gene_units = add_within_axis_gene_ranks(gene_units)
    gene_units = gene_units.sort_values(
        ["dataset", "comparison_group", "axis_order", "canonical_gene"]
    ).reset_index(drop=True)

    pathway_units = build_pathway_units(gene_units, gene_sets)
    signature_units = build_signature_units(pathway_units)
    transitions = build_transitions(pathway_units, signature_units)
    summary = build_summary(signature_units, transitions)

    gene_units.to_csv(OUT_GENE_UNITS, sep="\t", index=False)
    pathway_units.to_csv(OUT_PATHWAY_UNITS, sep="\t", index=False)
    signature_units.to_csv(OUT_SIGNATURE_UNITS, sep="\t", index=False)
    transitions.to_csv(OUT_TRANSITIONS, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    plot_signature_trajectories(signature_units)
    write_markdown(gene_units, pathway_units, signature_units, transitions, summary)

    print(f"Wrote {OUT_MD}")
    print(f"Gene units: {len(gene_units):,}")
    print(f"Pathway units: {len(pathway_units):,}")
    print(f"Signature units: {len(signature_units):,}")
    print(f"Transitions: {len(transitions):,}")


if __name__ == "__main__":
    main()
