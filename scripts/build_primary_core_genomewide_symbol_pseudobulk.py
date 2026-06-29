#!/usr/bin/env python3
"""Genome-wide same-symbol pseudobulk screen for full-matrix primary-core data.

This is the first true full-matrix step after the selected-gene screens. It uses
the full human DG taxonomy gene list as the target symbol universe, then streams
local full matrices from primary-core datasets. Mouse genes are mapped by
same-root upper-case symbol only, so this remains "ortholog-ready" rather than a
curated ortholog model.

The contrast statistics require at least two eligible broad classes within a
dataset/sample/gene, preventing single-class human DG anchor samples from
inflating dentate candidate ranks.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import build_primary_core_candidate_gene_pseudobulk as base


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"

HUMAN_DG_FULL = PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates"
GENE_METADATA = HUMAN_DG_FULL / "gene_metadata.tsv.gz"

OUT_EXPR = RESULTS / "primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz"
OUT_COVERAGE = RESULTS / "primary_core_genomewide_symbol_pseudobulk_coverage.tsv"
OUT_STATS = RESULTS / "primary_core_genomewide_symbol_pseudobulk_statistics.tsv"
OUT_SHARED = RESULTS / "primary_core_genomewide_symbol_pseudobulk_shared_hits.tsv"
OUT_BRANCH = RESULTS / "primary_core_genomewide_symbol_pseudobulk_branch_specific.tsv"
OUT_PLOT = RESULTS / "primary_core_genomewide_symbol_pseudobulk_shared_hits.png"
OUT_MD = RESULTS / "primary_core_genomewide_symbol_pseudobulk_analysis.md"

MIN_CLASS_CELLS = 20
MIN_DATASETS_DETECTED = 2

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def canon_gene(gene: object) -> str:
    return base.canon_gene(gene)


def load_symbol_targets() -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    genes = pd.read_csv(GENE_METADATA, sep="\t")
    genes["canonical_gene"] = genes["gene"].map(canon_gene)
    genes = genes.loc[genes["canonical_gene"].ne("")].drop_duplicates("canonical_gene").copy()
    rows = []
    for _, row in genes.iterrows():
        canonical = row["canonical_gene"]
        rows.append(
            {
                "gene": row["gene"],
                "canonical_gene": canonical,
                "human_symbol": canonical,
                "mouse_symbol": canonical[:1] + canonical[1:].lower(),
                "panel": "genomewide_same_symbol",
                "candidate_role": "genomewide_same_symbol_gene",
                "support_tier": "full_matrix_symbol_universe",
                "human_dg_full_n_counts": row.get("n_counts", np.nan),
                "human_dg_full_n_cells": row.get("n_cells", np.nan),
            }
        )
    targets = pd.DataFrame(rows)
    return targets, targets.set_index("canonical_gene").to_dict("index")


def install_symbol_targets() -> pd.DataFrame:
    targets, metadata = load_symbol_targets()
    base.TARGETS = targets
    base.TARGET_META = metadata
    base.TARGET_GENES = set(targets["canonical_gene"])
    return targets


def rel(path: Path | str) -> str:
    return base.rel(path)


def open_text(path: Path):
    return base.open_text(path)


def finalize_rows(
    *,
    dataset: str,
    sample: str,
    source_layer: str,
    expression_scope: str,
    expression_scale: str,
    source_path: str,
    groups: Counter[str],
    totals: dict[str, dict[str, float]],
    log_totals: dict[str, dict[str, float]],
    nonzeros: dict[str, dict[str, int]],
    source_gene_symbol: dict[str, str],
    present_genes: set[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for canonical_gene in sorted(present_genes):
        meta = base.gene_meta(canonical_gene)
        for group, n_cells in sorted(groups.items()):
            total = float(totals.get(group, {}).get(canonical_gene, 0.0))
            log_total = float(log_totals.get(group, {}).get(canonical_gene, 0.0))
            nz = int(nonzeros.get(group, {}).get(canonical_gene, 0))
            rows.append(
                {
                    "dataset": dataset,
                    "core_branch": base.CORE_BRANCH.get(dataset, "unknown"),
                    "sample": sample,
                    "source_layer": source_layer,
                    "expression_scope": expression_scope,
                    "expression_scale": expression_scale,
                    "broad_class": group,
                    "n_cells": int(n_cells),
                    **meta,
                    "source_gene_symbol": source_gene_symbol.get(canonical_gene, meta["gene"]),
                    "nonzero_cells": nz,
                    "detection_fraction": (nz / n_cells) if n_cells else np.nan,
                    "total_expression": total,
                    "mean_expression": (total / n_cells) if n_cells else np.nan,
                    "mean_log1p_expression": (log_total / n_cells) if n_cells else np.nan,
                    "source_path": source_path,
                }
            )
    return rows


def extract_wide_by_cell_table_chunked(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_genomewide_same_symbol",
    expression_scale: str = "raw_counts_or_reported_counts",
    tar_member: str | None = None,
    chunksize: int = 64,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    if tar_member:
        return base.extract_wide_by_cell_table(
            dataset=dataset,
            sample=sample,
            path=path,
            delimiter=delimiter,
            source_layer=source_layer,
            expression_scope=expression_scope,
            expression_scale=expression_scale,
            tar_member=tar_member,
            source_path=f"{rel(path)}:{tar_member}",
        )

    label_map = base.BACKBONE_LABELS.get((dataset, sample), {})
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split(delimiter)
    cells = [cell.strip().strip('"') for cell in header[1:]]
    labels = np.array([label_map.get(cell, "unmapped") for cell in cells], dtype=object)
    groups = Counter(labels.tolist())
    group_masks = {group: np.flatnonzero(labels == group) for group in groups}

    totals: dict[str, dict[str, float]] = {group: defaultdict(float) for group in groups}
    log_totals: dict[str, dict[str, float]] = {group: defaultdict(float) for group in groups}
    nonzeros: dict[str, dict[str, int]] = {group: defaultdict(int) for group in groups}
    source_gene_symbol: dict[str, str] = {}
    present_genes: set[str] = set()
    rows_scanned = 0

    reader = pd.read_csv(path, sep=delimiter, compression="gzip" if path.suffix == ".gz" else None, chunksize=chunksize)
    gene_col = reader._engine.names[0]
    for chunk in reader:
        rows_scanned += len(chunk)
        gene_symbols = chunk[gene_col].astype(str)
        canonical = gene_symbols.map(canon_gene)
        keep = canonical.isin(base.TARGET_GENES)
        if not keep.any():
            continue
        kept = chunk.loc[keep].copy()
        kept_canonical = canonical.loc[keep].to_numpy()
        kept_symbols = gene_symbols.loc[keep].to_numpy()
        values = kept.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(dtype=np.float32)
        for row_idx, canonical_gene in enumerate(kept_canonical):
            present_genes.add(canonical_gene)
            source_gene_symbol.setdefault(canonical_gene, kept_symbols[row_idx])
        log_values = np.log1p(values)
        nz_values = values > 0
        for group, cols in group_masks.items():
            if len(cols) == 0:
                continue
            sums = values[:, cols].sum(axis=1, dtype=np.float64)
            log_sums = log_values[:, cols].sum(axis=1, dtype=np.float64)
            nz = nz_values[:, cols].sum(axis=1)
            for i, canonical_gene in enumerate(kept_canonical):
                totals[group][canonical_gene] += float(sums[i])
                log_totals[group][canonical_gene] += float(log_sums[i])
                nonzeros[group][canonical_gene] += int(nz[i])

    rows = finalize_rows(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=rel(path),
        groups=groups,
        totals=totals,
        log_totals=log_totals,
        nonzeros=nonzeros,
        source_gene_symbol=source_gene_symbol,
        present_genes=present_genes,
    )
    coverage = {
        "dataset": dataset,
        "sample": sample,
        "source_layer": source_layer,
        "expression_scope": expression_scope,
        "expression_scale": expression_scale,
        "source_path": rel(path),
        "n_matrix_observations": len(cells),
        "n_labeled_observations": sum(v for k, v in groups.items() if k != "unmapped"),
        "n_unmapped_observations": groups.get("unmapped", 0),
        "n_gene_rows_scanned": rows_scanned,
        "target_genes_present": len(present_genes),
        "target_genes_total": len(base.TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(base.TARGET_GENES),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
    }
    return rows, coverage


def gse186538_full_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    labels = pd.read_csv(HUMAN_DG_FULL / "cell_metadata.tsv.gz", sep="\t", low_memory=False)
    labels["sample_for_pb"] = labels["samplename"].fillna("GSE186538").astype(str)
    labels["broad_class"] = "dentate_candidate"
    return base.aggregate_selected_sparse(
        dataset="GSE186538",
        source_layer="full_sparse_subset",
        expression_scope="full_human_dg_taxonomy_genes",
        expression_scale="raw_counts_full_dg_gc_subset",
        matrix_path=HUMAN_DG_FULL / "matrix_cells_by_genes.npz",
        var_path=HUMAN_DG_FULL / "gene_metadata.tsv.gz",
        labels=labels,
        sample_col="sample_for_pb",
        class_col="broad_class",
        source_path=rel(HUMAN_DG_FULL / "matrix_cells_by_genes.npz"),
    )


def collect_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    jobs = [
        lambda: extract_wide_by_cell_table_chunked(
            dataset="GSE104323",
            sample="10X_all_cells",
            path=EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
            delimiter="\t",
        ),
        lambda: extract_wide_by_cell_table_chunked(
            dataset="GSE95752",
            sample="C1_all_cells",
            path=EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz",
            delimiter="\t",
        ),
        lambda: base.extract_obs_by_gene_table(
            dataset="GSE292261",
            sample="SS2_filtered_counts",
            path=EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
            delimiter=",",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_wide_by_cell_table(
            dataset="GSE214309",
            sample="snRNA_counts",
            path=EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
            delimiter=",",
            header_has_gene_col=False,
            expression_scope="full_matrix_genomewide_same_symbol",
            expression_scale="reported_counts_gene_symbols",
        ),
        lambda: base.extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464549_P0",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464549_P0.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464549_P0.csv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464550_P8a",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464550_P8a.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464550_P8a.csv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464551_P8b",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464551_P8b.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464551_P8b.csv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_10x_matrix(
            dataset="GSE165657",
            sample="Cerebellum_aggr",
            matrix=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_genes.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_10x_matrix(
            dataset="GSE312658",
            sample="Ctrl",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_barcodes.tsv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        lambda: base.extract_10x_matrix(
            dataset="GSE312658",
            sample="cKO",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_barcodes.tsv.gz",
            expression_scope="full_matrix_genomewide_same_symbol",
        ),
        gse186538_full_rows,
    ]

    for job in jobs:
        job_rows, job_coverage = job()
        rows.extend(job_rows)
        coverage_rows.append(job_coverage)
        print(
            f"{job_coverage['dataset']} {job_coverage['sample']}: "
            f"{job_coverage['target_genes_present']}/{job_coverage['target_genes_total']} full-symbol target genes; "
            f"{job_coverage['n_labeled_observations']} labeled observations",
            flush=True,
        )

    return pd.DataFrame(rows), pd.DataFrame(coverage_rows)


def percentile_rank(series: pd.Series) -> pd.Series:
    if series.notna().sum() <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def add_strict_ranks(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    excluded = {"unmapped", "excluded_low_qc"}
    prelim = expr["n_cells"].ge(MIN_CLASS_CELLS) & ~expr["broad_class"].isin(excluded)
    expr["eligible_class"] = False
    expr["sample_gene_n_eligible_classes"] = 0
    expr["mean_log1p_rank_within_sample_gene"] = np.nan
    group_cols = ["dataset", "sample", "canonical_gene"]
    class_counts = expr.loc[prelim].groupby(group_cols)["broad_class"].transform("nunique")
    expr.loc[prelim, "sample_gene_n_eligible_classes"] = class_counts.to_numpy()
    eligible_idx = class_counts.index[class_counts.ge(2)]
    expr.loc[eligible_idx, "eligible_class"] = True
    expr.loc[eligible_idx, "mean_log1p_rank_within_sample_gene"] = (
        expr.loc[eligible_idx].groupby(group_cols)["mean_log1p_expression"].rank(method="average", pct=True)
    )
    return expr


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    if not valid.any():
        return out
    idx = p_values.index[valid].to_numpy()
    order = np.argsort(p_values.loc[idx].to_numpy(dtype=float))
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    adjusted = np.minimum.accumulate((p_sorted * len(p_sorted) / np.arange(1, len(p_sorted) + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def class_delta(sub: pd.DataFrame, target_class: str, background_classes: set[str]) -> dict[str, object]:
    target = sub.loc[sub["broad_class"].eq(target_class), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    background = sub.loc[sub["broad_class"].isin(background_classes), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    if len(target) == 0 or len(background) == 0:
        return {
            "n_target_units": len(target),
            "n_background_units": len(background),
            "target_median_rank": np.nan,
            "background_median_rank": np.nan,
            "rank_delta_vs_background": np.nan,
            "p_greater": np.nan,
        }
    return {
        "n_target_units": len(target),
        "n_background_units": len(background),
        "target_median_rank": float(np.median(target)),
        "background_median_rank": float(np.median(background)),
        "rank_delta_vs_background": float(np.median(target) - np.median(background)),
        "p_greater": float(stats.mannwhitneyu(target, background, alternative="greater").pvalue),
    }


def compute_stats(expr: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    eligible = expr.loc[expr["eligible_class"]].copy()
    target_meta = targets.set_index("canonical_gene")
    dentate_detection = (
        expr.loc[
            expr["broad_class"].eq("dentate_candidate")
            & expr["core_branch"].isin(["mouse_dentate", "human_dentate_hippocampus"])
        ]
        .groupby("canonical_gene")
        .agg(
            dentate_candidate_median_detection=("detection_fraction", "median"),
            n_dentate_datasets_detected_5pct=("dataset", lambda s: 0),
        )
    )
    dentate_detected = expr.loc[
        expr["broad_class"].eq("dentate_candidate")
        & expr["core_branch"].isin(["mouse_dentate", "human_dentate_hippocampus"])
        & expr["detection_fraction"].ge(0.05)
    ].groupby("canonical_gene")["dataset"].nunique()
    dentate_detection["n_dentate_datasets_detected_5pct"] = dentate_detected.reindex(dentate_detection.index).fillna(0).astype(int)

    cerebellar_detection = (
        expr.loc[expr["broad_class"].eq("cerebellar_candidate") & expr["core_branch"].eq("cerebellum")]
        .groupby("canonical_gene")
        .agg(
            cerebellar_candidate_median_detection=("detection_fraction", "median"),
            n_cerebellar_datasets_detected_5pct=("dataset", lambda s: 0),
        )
    )
    cerebellar_detected = expr.loc[
        expr["broad_class"].eq("cerebellar_candidate")
        & expr["core_branch"].eq("cerebellum")
        & expr["detection_fraction"].ge(0.05)
    ].groupby("canonical_gene")["dataset"].nunique()
    cerebellar_detection["n_cerebellar_datasets_detected_5pct"] = (
        cerebellar_detected.reindex(cerebellar_detection.index).fillna(0).astype(int)
    )

    dentate_background = {
        "non_dentate_background",
        "dentate_low_support",
        "other_or_ambiguous",
        "broad_neuronal_structural_warning",
    }
    cerebellar_background = {"other_or_ambiguous", "broad_neuronal_structural_warning"}
    rows: list[dict[str, object]] = []
    for canonical_gene, sub in eligible.groupby("canonical_gene", sort=False):
        meta = target_meta.loc[canonical_gene].to_dict() if canonical_gene in target_meta.index else base.gene_meta(canonical_gene)
        dentate_sub = sub.loc[sub["core_branch"].eq("mouse_dentate")]
        cerebellar_sub = sub.loc[sub["core_branch"].eq("cerebellum")]
        d = class_delta(dentate_sub, "dentate_candidate", dentate_background)
        c = class_delta(cerebellar_sub, "cerebellar_candidate", cerebellar_background)
        d_det = dentate_detection.loc[canonical_gene] if canonical_gene in dentate_detection.index else None
        c_det = cerebellar_detection.loc[canonical_gene] if canonical_gene in cerebellar_detection.index else None
        rows.append(
            {
                "gene": meta.get("gene", canonical_gene),
                "canonical_gene": canonical_gene,
                "candidate_role": meta.get("candidate_role", "genomewide_same_symbol_gene"),
                "dentate_candidate_units": d["n_target_units"],
                "dentate_background_units": d["n_background_units"],
                "dentate_candidate_median_rank": d["target_median_rank"],
                "dentate_background_median_rank": d["background_median_rank"],
                "dentate_rank_delta_vs_background": d["rank_delta_vs_background"],
                "dentate_rank_p_greater": d["p_greater"],
                "cerebellar_candidate_units": c["n_target_units"],
                "cerebellar_background_units": c["n_background_units"],
                "cerebellar_candidate_median_rank": c["target_median_rank"],
                "cerebellar_background_median_rank": c["background_median_rank"],
                "cerebellar_rank_delta_vs_background": c["rank_delta_vs_background"],
                "cerebellar_rank_p_greater": c["p_greater"],
                "dentate_candidate_median_detection": np.nan
                if d_det is None
                else float(d_det["dentate_candidate_median_detection"]),
                "cerebellar_candidate_median_detection": np.nan
                if c_det is None
                else float(c_det["cerebellar_candidate_median_detection"]),
                "n_dentate_datasets_detected_5pct": 0
                if d_det is None
                else int(d_det["n_dentate_datasets_detected_5pct"]),
                "n_cerebellar_datasets_detected_5pct": 0
                if c_det is None
                else int(c_det["n_cerebellar_datasets_detected_5pct"]),
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["dentate_rank_p_adj_bh"] = bh_adjust(stats_df["dentate_rank_p_greater"])
    stats_df["cerebellar_rank_p_adj_bh"] = bh_adjust(stats_df["cerebellar_rank_p_greater"])
    stats_df["shared_positive_rank_delta"] = (
        stats_df["dentate_rank_delta_vs_background"].gt(0)
        & stats_df["cerebellar_rank_delta_vs_background"].gt(0)
        & stats_df["dentate_candidate_units"].ge(2)
        & stats_df["cerebellar_candidate_units"].ge(3)
        & stats_df["n_dentate_datasets_detected_5pct"].ge(MIN_DATASETS_DETECTED)
        & stats_df["n_cerebellar_datasets_detected_5pct"].ge(MIN_DATASETS_DETECTED)
    )
    stats_df["shared_strict_bh_0_10"] = (
        stats_df["shared_positive_rank_delta"]
        & stats_df["dentate_rank_p_adj_bh"].lt(0.10)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.10)
    )
    stats_df["combined_rank_delta"] = stats_df["dentate_rank_delta_vs_background"].fillna(0) + stats_df[
        "cerebellar_rank_delta_vs_background"
    ].fillna(0)
    stats_df["minimum_branch_detection"] = stats_df[
        ["dentate_candidate_median_detection", "cerebellar_candidate_median_detection"]
    ].min(axis=1)
    shared = stats_df.loc[stats_df["shared_positive_rank_delta"]].sort_values(
        ["shared_strict_bh_0_10", "combined_rank_delta", "minimum_branch_detection"],
        ascending=[False, False, False],
    )
    dentate = stats_df.loc[
        stats_df["dentate_rank_delta_vs_background"].gt(0)
        & stats_df["dentate_rank_p_adj_bh"].lt(0.10)
        & (
            stats_df["cerebellar_rank_delta_vs_background"].le(0)
            | stats_df["cerebellar_rank_delta_vs_background"].isna()
            | stats_df["cerebellar_rank_p_adj_bh"].ge(0.20)
        )
    ].copy()
    dentate["branch_specificity"] = "dentate_biased"
    cerebellar = stats_df.loc[
        stats_df["cerebellar_rank_delta_vs_background"].gt(0)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.10)
        & (
            stats_df["dentate_rank_delta_vs_background"].le(0)
            | stats_df["dentate_rank_delta_vs_background"].isna()
            | stats_df["dentate_rank_p_adj_bh"].ge(0.20)
        )
    ].copy()
    cerebellar["branch_specificity"] = "cerebellar_biased"
    branch = pd.concat([dentate, cerebellar], ignore_index=True, sort=False).sort_values(
        ["branch_specificity", "combined_rank_delta"], ascending=[True, False]
    )
    return stats_df, shared, branch


def plot_shared(shared: pd.DataFrame) -> None:
    plot_df = shared.head(30).copy()
    if plot_df.empty:
        return
    fig, ax = plt.subplots(figsize=(10.2, 7.2))
    y = np.arange(len(plot_df))
    ax.barh(y - 0.18, plot_df["dentate_rank_delta_vs_background"], height=0.34, color="#168b7a", label="dentate candidate")
    ax.barh(y + 0.18, plot_df["cerebellar_rank_delta_vs_background"], height=0.34, color="#6d3bbd", label="cerebellar candidate")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Median within-sample rank delta versus local background")
    ax.set_title("Full-matrix same-symbol shared pseudobulk hits")
    ax.legend(frameon=False)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(targets: pd.DataFrame, expr: pd.DataFrame, coverage: pd.DataFrame, stats_df: pd.DataFrame, shared: pd.DataFrame, branch: pd.DataFrame) -> None:
    contrast_datasets = sorted(expr.loc[expr["eligible_class"], "dataset"].unique())
    lines = [
        "# Genome-Wide Same-Symbol Full-Matrix Pseudobulk Screen",
        "",
        "Date built: 2026-06-22",
        "",
        "## Scope",
        "",
        "This analysis streams local full matrices for primary-core datasets and aggregates a full human-symbol target universe from the GSE186538 human DG taxonomy gene list. Mouse genes are mapped by same-root upper-case symbol, so this is an ortholog-ready screen rather than a curated ortholog model.",
        "",
        "Contrast statistics require at least two eligible broad classes within a dataset/sample/gene. Single-class human DG anchor samples are retained for expression/detection but excluded from candidate-versus-background rank tests.",
        "",
        "## Coverage",
        "",
        f"- Full-symbol target universe: {targets['canonical_gene'].nunique():,} genes.",
        f"- Pseudobulk expression rows: {len(expr):,}.",
        f"- Datasets with full-matrix expression represented: {coverage['dataset'].nunique()}/10 primary datasets.",
        f"- Datasets contributing to rank contrasts: {len(contrast_datasets)} ({', '.join(contrast_datasets)}).",
        f"- Genes tested in contrast statistics: {stats_df['canonical_gene'].nunique():,}.",
        f"- Shared-positive rank genes: {len(shared):,}.",
        f"- Shared-positive genes passing BH<0.10 in both branches: {int(shared['shared_strict_bh_0_10'].sum()):,}.",
        f"- Branch-specific genes: {len(branch):,}.",
        "",
    ]
    for _, row in coverage.sort_values(["dataset", "sample"]).iterrows():
        lines.append(
            f"- `{row['dataset']}` / `{row['sample']}`: "
            f"{int(row['target_genes_present'])}/{int(row['target_genes_total'])} target symbols, "
            f"{int(row['n_labeled_observations'])}/{int(row['n_matrix_observations'])} labeled observations "
            f"(`{row['source_layer']}`)."
        )
    lines.extend(["", "## Top Shared Hits", ""])
    for _, row in shared.head(30).iterrows():
        lines.append(
            f"- `{row['gene']}`: dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
            f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}, "
            f"BH<0.10 both branches={bool(row['shared_strict_bh_0_10'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This is the first full-matrix discovery pass, but it still uses same-symbol mapping rather than a curated mouse-human ortholog table.",
            "- The result should be used to prioritize candidates and to design the final mixed-effect/ortholog DE model, not as the final gene-level claim.",
            "- The current full-matrix contrast layer excludes `GSE325391` and `GSE268609` from genome-wide tests because they are represented locally by selected-gene bridge objects or very large source objects requiring dedicated export.",
            "- `GSE186538` contributes full human DG expression/detection but not rank contrast because the extracted local object is a DG GC subset without local non-DG background.",
            "",
            "## Outputs",
            "",
            f"- Expression table: `{OUT_EXPR.relative_to(ROOT)}`",
            f"- Coverage table: `{OUT_COVERAGE.relative_to(ROOT)}`",
            f"- Statistics table: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Shared hits: `{OUT_SHARED.relative_to(ROOT)}`",
            f"- Branch-specific hits: `{OUT_BRANCH.relative_to(ROOT)}`",
            f"- Shared-hit plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    targets = install_symbol_targets()
    expr, coverage = collect_expression()
    expr = add_strict_ranks(expr)
    stats_df, shared, branch = compute_stats(expr, targets)
    plot_shared(shared)
    expr.to_csv(OUT_EXPR, sep="\t", index=False, compression="gzip")
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    stats_df.to_csv(OUT_STATS, sep="\t", index=False)
    shared.to_csv(OUT_SHARED, sep="\t", index=False)
    branch.to_csv(OUT_BRANCH, sep="\t", index=False)
    write_report(targets, expr, coverage, stats_df, shared, branch)
    print(f"Wrote {len(expr):,} full-matrix pseudobulk expression rows")
    print(f"Wrote {len(stats_df):,} same-symbol gene statistics")
    print(f"Wrote {len(shared):,} shared-positive hits")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
