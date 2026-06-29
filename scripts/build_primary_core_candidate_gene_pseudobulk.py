#!/usr/bin/env python3
"""Build candidate-gene pseudobulk summaries across the 10 primary-core datasets.

This is the first gene-level layer after the marker-module analysis. It focuses
on the 67 candidate genes selected from the human bridge packet, then aggregates
expression by the refined candidate-cell calls in each primary dataset.

The output is intentionally explicit about matrix depth:
- full_raw_matrix: local full matrix is used, but only candidate genes are read.
- selected_gene_bridge: local selected-gene bridge matrix is used.
- selected_norm_bridge: selected log1p(CP10K) human-core bridge is used.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats


csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"

CANDIDATE_PACKET = RESULTS / "human_bridge_candidate_gene_packet.tsv"
CORE_DATASETS = RESULTS / "integrated_primary_core_datasets.tsv"
BACKBONE_CALLS = RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"
HUMAN_CORE_LABELS = RESULTS / "human_core_tuned_labels.tsv.gz"
GSE325391_LABELS = RESULTS / "gse325391_human_core_label_projection.tsv.gz"
GSE268609_LABELS = RESULTS / "gse268609_human_core_label_projection.tsv.gz"

OUT_EXPR = RESULTS / "primary_core_candidate_gene_pseudobulk_expression.tsv"
OUT_COVERAGE = RESULTS / "primary_core_candidate_gene_pseudobulk_coverage.tsv"
OUT_STATS = RESULTS / "primary_core_candidate_gene_pseudobulk_statistics.tsv"
OUT_HITS = RESULTS / "primary_core_candidate_gene_pseudobulk_hits.tsv"
OUT_PLOT = RESULTS / "primary_core_candidate_gene_pseudobulk_effects.png"
OUT_MD = RESULTS / "primary_core_candidate_gene_pseudobulk_analysis.md"

MIN_CLASS_CELLS = 20

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def rel(path: Path | str) -> str:
    path = Path(path)
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def clean_gene(gene: object) -> str:
    if pd.isna(gene):
        return ""
    return str(gene).strip().strip('"').strip("'")


def canon_gene(gene: object) -> str:
    return clean_gene(gene).upper()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def tar_member_text(tf: tarfile.TarFile, member_name: str):
    raw = tf.extractfile(member_name)
    if raw is None:
        raise FileNotFoundError(member_name)
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")
    return io.TextIOWrapper(raw, newline="")


def load_targets() -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    packet = pd.read_csv(CANDIDATE_PACKET, sep="\t")
    packet["canonical_gene"] = packet["gene"].map(canon_gene)
    packet["human_symbol"] = packet["canonical_gene"]
    packet["mouse_symbol"] = packet["canonical_gene"].str.slice(0, 1) + packet["canonical_gene"].str.slice(1).str.lower()
    keep_cols = [
        "gene",
        "canonical_gene",
        "human_symbol",
        "mouse_symbol",
        "panel",
        "candidate_role",
        "support_tier",
        "notes",
    ]
    targets = packet[keep_cols].drop_duplicates("canonical_gene").copy()
    metadata = targets.set_index("canonical_gene").to_dict("index")
    return targets, metadata


TARGETS, TARGET_META = load_targets()
TARGET_GENES = set(TARGETS["canonical_gene"])


def resolve_gene_from_parts(parts: list[str]) -> str | None:
    for part in parts[:3]:
        cg = canon_gene(part)
        if cg in TARGET_GENES:
            return cg
    return None


def broad_class_from_backbone_call(call: str) -> str:
    if call == "candidate_dentate_granule":
        return "dentate_candidate"
    if call == "candidate_cerebellar_granule":
        return "cerebellar_candidate"
    if call == "known_non_dentate_reference":
        return "non_dentate_background"
    if call == "cerebellum_dentate_panel_high_warning":
        return "broad_neuronal_structural_warning"
    if call == "dentate_like_low_support":
        return "dentate_low_support"
    return "other_or_ambiguous"


def classify_human_core(label: str) -> str:
    if label == "curated_human_dg_gc_anchor":
        return "dentate_candidate"
    if label == "human_dg_like_high_confidence":
        return "dentate_candidate"
    if label == "human_dg_like_candidate":
        return "dentate_low_support"
    if label == "immature_neurogenic_candidate":
        return "dentate_candidate"
    if label == "immature_neurogenic_candidate_low_support":
        return "dentate_low_support"
    if label == "non_neuronal_background":
        return "non_dentate_background"
    if label == "broad_neuronal_structural_warning":
        return "broad_neuronal_structural_warning"
    return "other_or_ambiguous"


def classify_gse325391(label: str) -> str:
    if label in {"adult_human_dg_mature_anchor", "adult_human_dg_differentiating_anchor"}:
        return "dentate_candidate"
    if label == "adult_dg_background_warning":
        return "non_dentate_background"
    return "other_or_ambiguous"


def classify_gse268609(label: str) -> str:
    if label in {"human_dg_like_high_confidence", "human_dg_like_candidate"}:
        return "dentate_candidate"
    if label == "immature_neurogenic_candidate":
        return "dentate_candidate"
    if label == "non_neuronal_background":
        return "non_dentate_background"
    if label == "broad_neuronal_structural_warning":
        return "broad_neuronal_structural_warning"
    return "other_or_ambiguous"


def load_core_branches() -> dict[str, str]:
    core = pd.read_csv(CORE_DATASETS, sep="\t")
    return dict(zip(core["dataset"], core["core_branch"]))


CORE_BRANCH = load_core_branches()


def load_backbone_labels() -> dict[tuple[str, str], dict[str, str]]:
    cols = ["dataset", "sample", "cell_id", "candidate_call"]
    calls = pd.read_csv(BACKBONE_CALLS, sep="\t", usecols=cols, low_memory=False)
    calls = calls.loc[calls["dataset"].isin(CORE_BRANCH)].copy()
    calls["broad_class"] = calls["candidate_call"].map(broad_class_from_backbone_call)
    labels: dict[tuple[str, str], dict[str, str]] = {}
    for (dataset, sample), sub in calls.groupby(["dataset", "sample"], sort=False):
        labels[(dataset, sample)] = dict(zip(sub["cell_id"].astype(str), sub["broad_class"].astype(str)))
    return labels


BACKBONE_LABELS = load_backbone_labels()


def gene_meta(canonical_gene: str) -> dict[str, object]:
    meta = TARGET_META.get(canonical_gene, {})
    return {
        "gene": meta.get("gene", canonical_gene),
        "canonical_gene": canonical_gene,
        "human_symbol": meta.get("human_symbol", canonical_gene),
        "mouse_symbol": meta.get("mouse_symbol", canonical_gene.capitalize()),
        "panel": meta.get("panel", ""),
        "candidate_role": meta.get("candidate_role", ""),
        "support_tier": meta.get("support_tier", ""),
    }


def empty_accumulator(groups: Counter[str]) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, int]]]:
    totals = {group: defaultdict(float) for group in groups}
    nonzeros = {group: defaultdict(int) for group in groups}
    return totals, nonzeros


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
        meta = gene_meta(canonical_gene)
        for group, n_cells in sorted(groups.items()):
            total = float(totals.get(group, {}).get(canonical_gene, 0.0))
            log_total = float(log_totals.get(group, {}).get(canonical_gene, 0.0))
            nz = int(nonzeros.get(group, {}).get(canonical_gene, 0))
            rows.append(
                {
                    "dataset": dataset,
                    "core_branch": CORE_BRANCH.get(dataset, "unknown"),
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


def extract_wide_by_cell_table(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_candidate_genes",
    expression_scale: str = "raw_counts_or_reported_counts",
    header_has_gene_col: bool = True,
    source_path: str | None = None,
    tar_member: str | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    label_map = BACKBONE_LABELS.get((dataset, sample), {})
    if tar_member:
        with tarfile.open(path) as tf:
            with tar_member_text(tf, tar_member) as fh:
                return _extract_wide_reader(
                    dataset=dataset,
                    sample=sample,
                    fh=fh,
                    delimiter=delimiter,
                    label_map=label_map,
                    source_layer=source_layer,
                    expression_scope=expression_scope,
                    expression_scale=expression_scale,
                    header_has_gene_col=header_has_gene_col,
                    source_path=source_path or f"{rel(path)}:{tar_member}",
                )
    with open_text(path) as fh:
        return _extract_wide_reader(
            dataset=dataset,
            sample=sample,
            fh=fh,
            delimiter=delimiter,
            label_map=label_map,
            source_layer=source_layer,
            expression_scope=expression_scope,
            expression_scale=expression_scale,
            header_has_gene_col=header_has_gene_col,
            source_path=source_path or rel(path),
        )


def _extract_wide_reader(
    *,
    dataset: str,
    sample: str,
    fh,
    delimiter: str,
    label_map: dict[str, str],
    source_layer: str,
    expression_scope: str,
    expression_scale: str,
    header_has_gene_col: bool,
    source_path: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    reader = csv.reader(fh, delimiter=delimiter)
    header = next(reader)
    cell_fields = header[1:] if header_has_gene_col else header
    cells = [str(cell).strip().strip('"') for cell in cell_fields]
    labels = [label_map.get(cell, "unmapped") for cell in cells]
    groups = Counter(labels)
    totals, nonzeros = empty_accumulator(groups)
    log_totals = {group: defaultdict(float) for group in groups}
    present_genes: set[str] = set()
    source_gene_symbol: dict[str, str] = {}
    rows_scanned = 0

    for row in reader:
        if not row:
            continue
        rows_scanned += 1
        canonical_gene = canon_gene(row[0])
        if canonical_gene not in TARGET_GENES:
            continue
        present_genes.add(canonical_gene)
        source_gene_symbol.setdefault(canonical_gene, clean_gene(row[0]))
        values = row[1:] if header_has_gene_col else row[1:]
        for idx, raw in enumerate(values):
            if idx >= len(labels):
                break
            try:
                value = float(raw)
            except ValueError:
                value = 0.0
            if value <= 0:
                continue
            group = labels[idx]
            totals[group][canonical_gene] += value
            log_totals[group][canonical_gene] += math.log1p(value)
            nonzeros[group][canonical_gene] += 1

    rows = finalize_rows(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
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
        "source_path": source_path,
        "n_matrix_observations": len(cells),
        "n_labeled_observations": sum(v for k, v in groups.items() if k != "unmapped"),
        "n_unmapped_observations": groups.get("unmapped", 0),
        "n_gene_rows_scanned": rows_scanned,
        "target_genes_present": len(present_genes),
        "target_genes_total": len(TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(TARGET_GENES),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
    }
    return rows, coverage


def extract_obs_by_gene_table(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_candidate_genes",
    expression_scale: str = "raw_counts_or_reported_counts",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    label_map = BACKBONE_LABELS.get((dataset, sample), {})
    source_path = rel(path)
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        target_cols: list[tuple[int, str, str]] = []
        for idx, gene in enumerate(header[1:], start=1):
            canonical_gene = canon_gene(gene)
            if canonical_gene in TARGET_GENES:
                target_cols.append((idx, canonical_gene, clean_gene(gene)))
        source_gene_symbol = {canonical_gene: raw for _, canonical_gene, raw in target_cols}
        present_genes = {canonical_gene for _, canonical_gene, _ in target_cols}
        groups: Counter[str] = Counter()
        totals: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        log_totals: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        nonzeros: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        n_cells = 0
        for row in reader:
            if not row:
                continue
            n_cells += 1
            cell = str(row[0]).strip().strip('"')
            group = label_map.get(cell, "unmapped")
            groups[group] += 1
            for col_idx, canonical_gene, _ in target_cols:
                if col_idx >= len(row):
                    continue
                try:
                    value = float(row[col_idx])
                except ValueError:
                    value = 0.0
                if value <= 0:
                    continue
                totals[group][canonical_gene] += value
                log_totals[group][canonical_gene] += math.log1p(value)
                nonzeros[group][canonical_gene] += 1

    rows = finalize_rows(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
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
        "source_path": source_path,
        "n_matrix_observations": n_cells,
        "n_labeled_observations": sum(v for k, v in groups.items() if k != "unmapped"),
        "n_unmapped_observations": groups.get("unmapped", 0),
        "n_gene_rows_scanned": len(header) - 1,
        "target_genes_present": len(present_genes),
        "target_genes_total": len(TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(TARGET_GENES),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
    }
    return rows, coverage


def read_barcodes(path: Path) -> list[str]:
    with open_text(path) as fh:
        return [line.rstrip("\n").split("\t")[0] for line in fh if line.strip()]


def read_10x_target_features(path: Path) -> tuple[dict[int, str], dict[str, str], int]:
    target_rows: dict[int, str] = {}
    source_gene_symbol: dict[str, str] = {}
    n_features = 0
    with open_text(path) as fh:
        for idx, line in enumerate(fh, start=1):
            n_features = idx
            parts = line.rstrip("\n").split("\t")
            canonical_gene = resolve_gene_from_parts(parts)
            if canonical_gene:
                target_rows[idx] = canonical_gene
                source_gene_symbol.setdefault(canonical_gene, parts[1] if len(parts) > 1 else parts[0])
    return target_rows, source_gene_symbol, n_features


def extract_10x_matrix(
    *,
    dataset: str,
    sample: str,
    matrix: Path,
    features: Path,
    barcodes: Path,
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_candidate_genes",
    expression_scale: str = "raw_counts",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    target_rows, source_gene_symbol, n_features = read_10x_target_features(features)
    barcode_list = read_barcodes(barcodes)
    label_map = BACKBONE_LABELS.get((dataset, sample), {})
    labels = [label_map.get(barcode, "unmapped") for barcode in barcode_list]
    groups = Counter(labels)
    totals, nonzeros = empty_accumulator(groups)
    log_totals = {group: defaultdict(float) for group in groups}
    source_path = rel(matrix)

    n_features_seen = n_obs_seen = n_nonzero_seen = None
    shape_seen = False
    with open_text(matrix) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if not shape_seen:
                n_features_seen, n_obs_seen, n_nonzero_seen = map(int, parts[:3])
                shape_seen = True
                continue
            feature_idx = int(parts[0])
            if feature_idx not in target_rows:
                continue
            cell_idx = int(parts[1]) - 1
            if cell_idx < 0 or cell_idx >= len(labels):
                continue
            value = float(parts[2])
            if value <= 0:
                continue
            group = labels[cell_idx]
            canonical_gene = target_rows[feature_idx]
            totals[group][canonical_gene] += value
            log_totals[group][canonical_gene] += math.log1p(value)
            nonzeros[group][canonical_gene] += 1

    present_genes = set(target_rows.values())
    rows = finalize_rows(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
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
        "source_path": source_path,
        "n_matrix_observations": int(n_obs_seen or len(barcode_list)),
        "n_labeled_observations": sum(v for k, v in groups.items() if k != "unmapped"),
        "n_unmapped_observations": groups.get("unmapped", 0),
        "n_gene_rows_scanned": int(n_features_seen or n_features),
        "target_genes_present": len(present_genes),
        "target_genes_total": len(TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(TARGET_GENES),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
        "matrix_nonzero_entries": int(n_nonzero_seen or 0),
    }
    return rows, coverage


def load_sparse_npz(path: Path):
    obj = np.load(path, allow_pickle=False)
    if {"data", "indices", "indptr", "shape"}.issubset(set(obj.files)):
        return sparse.csr_matrix((obj["data"], obj["indices"], obj["indptr"]), shape=tuple(obj["shape"]))
    return sparse.load_npz(path)


def aggregate_selected_sparse(
    *,
    dataset: str,
    source_layer: str,
    expression_scope: str,
    expression_scale: str,
    matrix_path: Path,
    var_path: Path,
    labels: pd.DataFrame,
    sample_col: str,
    class_col: str,
    source_path: str | None = None,
    dataset_filter: str | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    var = pd.read_csv(var_path, sep="\t")
    var["canonical_gene"] = var["gene"].map(canon_gene)
    target_cols = [(idx, row["canonical_gene"], row["gene"]) for idx, row in var.iterrows() if row["canonical_gene"] in TARGET_GENES]
    col_indices = [idx for idx, _, _ in target_cols]
    present_genes = {gene for _, gene, _ in target_cols}
    source_gene_symbol = {gene: raw for _, gene, raw in target_cols}
    X = load_sparse_npz(matrix_path).tocsr()
    labels = labels.copy()
    if dataset_filter and "dataset" in labels:
        keep = labels["dataset"].astype(str).eq(dataset_filter)
        labels = labels.loc[keep].copy()
        X = X[keep.to_numpy(), :]
    labels = labels.reset_index(drop=True)
    if len(labels) != X.shape[0]:
        raise ValueError(f"{dataset}: label rows ({len(labels)}) do not match matrix rows ({X.shape[0]})")

    if col_indices and col_indices == list(range(X.shape[1])):
        X_target = X
    else:
        X_target = X[:, col_indices] if col_indices else X[:, []]
    out_rows: list[dict[str, object]] = []
    groups = labels[[sample_col, class_col]].fillna("missing").astype(str)
    group_counts = groups.value_counts(sort=False).to_dict()
    for (sample, broad_class), n_cells in sorted(group_counts.items()):
        idx = np.flatnonzero((groups[sample_col].to_numpy() == sample) & (groups[class_col].to_numpy() == broad_class))
        if len(idx) == 0:
            continue
        Xg = X_target[idx, :]
        sums = np.asarray(Xg.sum(axis=0)).ravel()
        nonzero_counts = np.asarray((Xg > 0).sum(axis=0)).ravel()
        log_Xg = Xg.copy()
        log_Xg.data = np.log1p(log_Xg.data)
        log_sums = np.asarray(log_Xg.sum(axis=0)).ravel()
        for j, (_, canonical_gene, _) in enumerate(target_cols):
            meta = gene_meta(canonical_gene)
            out_rows.append(
                {
                    "dataset": dataset,
                    "core_branch": CORE_BRANCH.get(dataset, "unknown"),
                    "sample": sample,
                    "source_layer": source_layer,
                    "expression_scope": expression_scope,
                    "expression_scale": expression_scale,
                    "broad_class": broad_class,
                    "n_cells": int(n_cells),
                    **meta,
                    "source_gene_symbol": source_gene_symbol.get(canonical_gene, meta["gene"]),
                    "nonzero_cells": int(nonzero_counts[j]),
                    "detection_fraction": float(nonzero_counts[j] / n_cells) if n_cells else np.nan,
                    "total_expression": float(sums[j]),
                    "mean_expression": float(sums[j] / n_cells) if n_cells else np.nan,
                    "mean_log1p_expression": float(log_sums[j] / n_cells) if n_cells else np.nan,
                    "source_path": source_path or rel(matrix_path),
                }
            )

    broad_counts = Counter(groups[class_col])
    coverage = {
        "dataset": dataset,
        "sample": "selected_bridge_all_samples",
        "source_layer": source_layer,
        "expression_scope": expression_scope,
        "expression_scale": expression_scale,
        "source_path": source_path or rel(matrix_path),
        "n_matrix_observations": int(X.shape[0]),
        "n_labeled_observations": int(X.shape[0]),
        "n_unmapped_observations": 0,
        "n_gene_rows_scanned": int(var.shape[0]),
        "target_genes_present": len(present_genes),
        "target_genes_total": len(TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(TARGET_GENES),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(broad_counts.items())),
    }
    return out_rows, coverage


def human_core_selected_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    obs = pd.read_csv(PROCESSED / "human_core_normalized_reduced_object/obs.tsv.gz", sep="\t", low_memory=False)
    labels = pd.read_csv(HUMAN_CORE_LABELS, sep="\t", usecols=["cell_id", "dataset", "replicate_unit", "tuned_label"], low_memory=False)
    labels["broad_class"] = labels["tuned_label"].map(classify_human_core)
    merged = obs[["cell_id", "dataset"]].merge(labels, on=["cell_id", "dataset"], how="left")
    merged = merged.loc[merged["dataset"].eq("GSE186538")].copy()
    merged["sample"] = merged["replicate_unit"].fillna("GSE186538")
    return aggregate_selected_sparse(
        dataset="GSE186538",
        source_layer="selected_norm_bridge",
        expression_scope="selected_human_core_candidate_genes",
        expression_scale="log1p_cp10k_selected_genes",
        matrix_path=PROCESSED / "human_core_normalized_reduced_object/X_log1p_cp10k_selected_genes.npz",
        var_path=PROCESSED / "human_core_normalized_reduced_object/var.tsv",
        labels=obs[["cell_id", "dataset"]].merge(labels, on=["cell_id", "dataset"], how="left").assign(
            sample=lambda d: d["replicate_unit"].fillna("missing"),
            broad_class=lambda d: d["broad_class"].fillna("other_or_ambiguous"),
        ),
        sample_col="sample",
        class_col="broad_class",
        source_path=rel(PROCESSED / "human_core_normalized_reduced_object/X_log1p_cp10k_selected_genes.npz"),
        dataset_filter="GSE186538",
    )


def gse325391_selected_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    labels = pd.read_csv(GSE325391_LABELS, sep="\t", low_memory=False)
    keep = labels["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    labels["broad_class"] = labels["source_anchor_label"].map(classify_gse325391)
    labels.loc[~keep, "broad_class"] = "excluded_low_qc"
    labels["sample_for_pb"] = labels["sample"].astype(str)
    return aggregate_selected_sparse(
        dataset="GSE325391",
        source_layer="selected_gene_bridge",
        expression_scope="selected_bridge_candidate_genes",
        expression_scale="raw_counts_selected_genes",
        matrix_path=PROCESSED / "gse325391_adult_dg_selected/matrix_cells_by_selected_genes.npz",
        var_path=PROCESSED / "gse325391_adult_dg_selected/var_selected_features.tsv",
        labels=labels,
        sample_col="sample_for_pb",
        class_col="broad_class",
        source_path=rel(PROCESSED / "gse325391_adult_dg_selected/matrix_cells_by_selected_genes.npz"),
    )


def gse268609_selected_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    labels = pd.read_csv(GSE268609_LABELS, sep="\t", low_memory=False, dtype={"sample_id": str})
    keep = labels["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    labels["broad_class"] = labels["projected_label"].map(classify_gse268609)
    labels.loc[~keep, "broad_class"] = "excluded_low_qc"
    labels["sample_for_pb"] = labels["sample_id"].astype(str)
    return aggregate_selected_sparse(
        dataset="GSE268609",
        source_layer="selected_gene_bridge",
        expression_scope="selected_bridge_candidate_genes",
        expression_scale="raw_counts_selected_genes",
        matrix_path=PROCESSED / "gse268609_rna_selected/matrix_cells_by_selected_genes.npz",
        var_path=PROCESSED / "gse268609_rna_selected/var_selected_features.tsv",
        labels=labels,
        sample_col="sample_for_pb",
        class_col="broad_class",
        source_path=rel(PROCESSED / "gse268609_rna_selected/matrix_cells_by_selected_genes.npz"),
    )


def collect_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    jobs = [
        lambda: extract_wide_by_cell_table(
            dataset="GSE104323",
            sample="10X_all_cells",
            path=EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
            delimiter="\t",
        ),
        lambda: extract_wide_by_cell_table(
            dataset="GSE95752",
            sample="C1_all_cells",
            path=EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz",
            delimiter="\t",
        ),
        lambda: extract_obs_by_gene_table(
            dataset="GSE292261",
            sample="SS2_filtered_counts",
            path=EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
            delimiter=",",
        ),
        lambda: extract_wide_by_cell_table(
            dataset="GSE214309",
            sample="snRNA_counts",
            path=EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
            delimiter=",",
            header_has_gene_col=False,
            expression_scale="reported_counts_gene_symbols",
        ),
        lambda: extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464549_P0",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464549_P0.csv.gz",
            tar_member="GSM3464549_P0.csv.gz",
        ),
        lambda: extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464550_P8a",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464550_P8a.csv.gz",
            tar_member="GSM3464550_P8a.csv.gz",
        ),
        lambda: extract_wide_by_cell_table(
            dataset="GSE122357",
            sample="GSM3464551_P8b",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464551_P8b.csv.gz",
            tar_member="GSM3464551_P8b.csv.gz",
        ),
        lambda: extract_10x_matrix(
            dataset="GSE165657",
            sample="Cerebellum_aggr",
            matrix=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_genes.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
        ),
        lambda: extract_10x_matrix(
            dataset="GSE312658",
            sample="Ctrl",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_barcodes.tsv.gz",
        ),
        lambda: extract_10x_matrix(
            dataset="GSE312658",
            sample="cKO",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_barcodes.tsv.gz",
        ),
        human_core_selected_rows,
        gse325391_selected_rows,
        gse268609_selected_rows,
    ]

    for job in jobs:
        job_rows, job_coverage = job()
        rows.extend(job_rows)
        coverage_rows.append(job_coverage)
        print(
            f"{job_coverage['dataset']} {job_coverage['sample']}: "
            f"{job_coverage['target_genes_present']}/{job_coverage['target_genes_total']} target genes; "
            f"{job_coverage['n_labeled_observations']} labeled observations"
        )

    expr = pd.DataFrame(rows)
    coverage = pd.DataFrame(coverage_rows)
    return expr, coverage


def percentile_rank(series: pd.Series) -> pd.Series:
    if series.notna().sum() <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def add_within_sample_gene_ranks(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    excluded_classes = {"unmapped", "excluded_low_qc"}
    eligible = expr["n_cells"].ge(MIN_CLASS_CELLS) & ~expr["broad_class"].isin(excluded_classes)
    expr["eligible_class"] = eligible
    expr["mean_log1p_rank_within_sample_gene"] = np.nan
    for _, idx in expr.loc[eligible].groupby(["dataset", "sample", "canonical_gene"]).groups.items():
        values = expr.loc[idx, "mean_log1p_expression"]
        expr.loc[idx, "mean_log1p_rank_within_sample_gene"] = percentile_rank(values)
    return expr


def class_delta(sub: pd.DataFrame, target_class: str, background_classes: set[str]) -> tuple[float, int, int, float, float]:
    target = sub.loc[sub["broad_class"].eq(target_class), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    background = sub.loc[sub["broad_class"].isin(background_classes), "mean_log1p_rank_within_sample_gene"].dropna().to_numpy(dtype=float)
    if len(target) == 0 or len(background) == 0:
        return np.nan, len(target), len(background), np.nan, np.nan
    p = stats.mannwhitneyu(target, background, alternative="greater").pvalue
    return float(np.median(target) - np.median(background)), len(target), len(background), float(np.median(target)), float(p)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    if not valid.any():
        return out
    idx = p_values.index[valid].to_numpy()
    order = np.argsort(p_values.loc[idx].to_numpy(dtype=float))
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    m = len(p_sorted)
    adjusted = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def compute_gene_stats(expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = expr.loc[expr["eligible_class"]].copy()
    dentate_background = {"non_dentate_background", "other_or_ambiguous", "broad_neuronal_structural_warning"}
    cerebellar_background = {"other_or_ambiguous", "broad_neuronal_structural_warning"}
    rows: list[dict[str, object]] = []
    for canonical_gene, sub in eligible.groupby("canonical_gene"):
        meta = gene_meta(canonical_gene)
        dentate_sub = sub.loc[sub["core_branch"].isin(["mouse_dentate", "human_dentate_hippocampus"])]
        cerebellar_sub = sub.loc[sub["core_branch"].eq("cerebellum")]
        d_delta, d_n, d_bg_n, d_med, d_p = class_delta(dentate_sub, "dentate_candidate", dentate_background)
        c_delta, c_n, c_bg_n, c_med, c_p = class_delta(cerebellar_sub, "cerebellar_candidate", cerebellar_background)
        all_dentate = dentate_sub.loc[dentate_sub["broad_class"].eq("dentate_candidate")]
        all_cerebellar = cerebellar_sub.loc[cerebellar_sub["broad_class"].eq("cerebellar_candidate")]
        rows.append(
            {
                **meta,
                "dentate_candidate_units": d_n,
                "dentate_background_units": d_bg_n,
                "dentate_candidate_median_rank": d_med,
                "dentate_rank_delta_vs_background": d_delta,
                "dentate_rank_p_greater": d_p,
                "cerebellar_candidate_units": c_n,
                "cerebellar_background_units": c_bg_n,
                "cerebellar_candidate_median_rank": c_med,
                "cerebellar_rank_delta_vs_background": c_delta,
                "cerebellar_rank_p_greater": c_p,
                "dentate_candidate_median_detection": all_dentate["detection_fraction"].median(),
                "cerebellar_candidate_median_detection": all_cerebellar["detection_fraction"].median(),
                "n_datasets_detected_dentate_candidate_5pct": int(
                    all_dentate.loc[all_dentate["detection_fraction"].ge(0.05), "dataset"].nunique()
                ),
                "n_datasets_detected_cerebellar_candidate_5pct": int(
                    all_cerebellar.loc[all_cerebellar["detection_fraction"].ge(0.05), "dataset"].nunique()
                ),
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["dentate_rank_p_adj_bh"] = bh_adjust(stats_df["dentate_rank_p_greater"])
    stats_df["cerebellar_rank_p_adj_bh"] = bh_adjust(stats_df["cerebellar_rank_p_greater"])
    stats_df["shared_structural_support"] = (
        stats_df["candidate_role"].eq("shared_structural_executor")
        &
        stats_df["dentate_rank_delta_vs_background"].gt(0)
        & stats_df["cerebellar_rank_delta_vs_background"].gt(0)
        & stats_df["dentate_candidate_units"].ge(3)
        & stats_df["cerebellar_candidate_units"].ge(3)
    )
    stats_df["strict_shared_support"] = (
        stats_df["shared_structural_support"]
        & stats_df["dentate_rank_p_adj_bh"].lt(0.2)
        & stats_df["cerebellar_rank_p_adj_bh"].lt(0.2)
    )

    hit_cols = [
        "gene",
        "canonical_gene",
        "panel",
        "candidate_role",
        "support_tier",
        "dentate_rank_delta_vs_background",
        "cerebellar_rank_delta_vs_background",
        "dentate_candidate_median_detection",
        "cerebellar_candidate_median_detection",
        "dentate_rank_p_adj_bh",
        "cerebellar_rank_p_adj_bh",
        "shared_structural_support",
        "strict_shared_support",
    ]
    hits = stats_df.sort_values(
        [
            "strict_shared_support",
            "shared_structural_support",
            "dentate_rank_delta_vs_background",
            "cerebellar_rank_delta_vs_background",
        ],
        ascending=[False, False, False, False],
    )[hit_cols].copy()
    return stats_df, hits


def plot_gene_effects(stats_df: pd.DataFrame) -> None:
    plot_df = stats_df.loc[stats_df["candidate_role"].eq("shared_structural_executor")].copy()
    plot_df = plot_df.dropna(subset=["dentate_rank_delta_vs_background", "cerebellar_rank_delta_vs_background"])
    if plot_df.empty:
        return
    plot_df["score"] = plot_df["dentate_rank_delta_vs_background"] + plot_df["cerebellar_rank_delta_vs_background"]
    plot_df = plot_df.sort_values("score", ascending=False).head(24)
    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    y = np.arange(len(plot_df))
    ax.barh(y - 0.18, plot_df["dentate_rank_delta_vs_background"], height=0.34, color="#168b7a", label="dentate candidate")
    ax.barh(y + 0.18, plot_df["cerebellar_rank_delta_vs_background"], height=0.34, color="#6d3bbd", label="cerebellar candidate")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Median within-sample gene-rank delta versus local background")
    ax.set_title("Candidate structural-executor gene support across primary core")
    ax.legend(frameon=False)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(expr: pd.DataFrame, coverage: pd.DataFrame, stats_df: pd.DataFrame, hits: pd.DataFrame) -> None:
    datasets_with_expression = coverage.loc[coverage["target_genes_present"].gt(0), "dataset"].nunique()
    datasets_total = pd.read_csv(CORE_DATASETS, sep="\t")["dataset"].nunique()
    shared = int(stats_df["shared_structural_support"].sum())
    strict = int(stats_df["strict_shared_support"].sum())
    top_hits = hits.loc[
        hits["candidate_role"].eq("shared_structural_executor") & hits["shared_structural_support"]
    ].head(15)
    lines = [
        "# Primary-Core Candidate-Gene Pseudobulk Analysis",
        "",
        "Date built: 2026-06-22",
        "",
        "## Scope",
        "",
        "This layer aggregates the 67 candidate genes from the human bridge packet by refined candidate-cell class across the strict 10-dataset primary core.",
        "",
        "It is DE-adjacent but not yet a genome-wide mixed-effect differential-expression model. It tests whether the proposed structural-executor genes are elevated in dentate and cerebellar candidate granule populations relative to local backgrounds.",
        "",
        f"Minimum broad-class size used for rank statistics: {MIN_CLASS_CELLS} cells/nuclei.",
        "",
        "## Coverage",
        "",
        f"- Primary datasets represented by at least one candidate gene: {datasets_with_expression}/{datasets_total}.",
        f"- Pseudobulk expression rows: {len(expr):,}.",
        f"- Candidate genes tested: {stats_df['canonical_gene'].nunique():,}.",
        f"- Structural-executor genes with positive dentate and cerebellar candidate deltas: {shared}.",
        f"- Structural-executor genes passing the stricter exploratory BH<0.2 rule in both branches: {strict}.",
        "",
    ]
    for _, row in coverage.sort_values(["dataset", "sample"]).iterrows():
        note = ""
        if row["target_genes_present"] == 0:
            note = " Gene-symbol cleanup needed before this source can contribute gene-level evidence."
        lines.append(
            f"- `{row['dataset']}` / `{row['sample']}`: "
            f"{int(row['target_genes_present'])}/{int(row['target_genes_total'])} candidate genes, "
            f"{int(row['n_labeled_observations'])}/{int(row['n_matrix_observations'])} labeled observations "
            f"(`{row['source_layer']}`).{note}"
        )
    lines.extend(["", "## Top Shared-Executor Signals", ""])
    for _, row in top_hits.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['panel']}): dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
            f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}, "
            f"shared={bool(row['shared_structural_support'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The useful positive result is not that dentate and cerebellar granule cells are transcriptionally identical. The module analysis already argues against that. The stronger claim is that both lineages repeatedly use an elevated structural-executor gene axis on top of distinct regional identity programs.",
            "- This candidate-gene pass supports prioritizing shared structural genes whose rank is positive in both dentate-candidate and cerebellar-candidate pseudobulks.",
            "- `GSE214309` contributes gene-symbol-resolved candidate-gene pseudobulk evidence in this focused pass. Its broader whole-transcriptome cleanup should still be revisited before genome-wide DE.",
            "- The next stricter analysis should expand from these 67 genes to genome-wide ortholog-aware pseudobulk DE, with donor/sample/stage modeled explicitly.",
            "",
            "## Outputs",
            "",
            f"- Expression table: `{OUT_EXPR.relative_to(ROOT)}`",
            f"- Coverage table: `{OUT_COVERAGE.relative_to(ROOT)}`",
            f"- Gene statistics: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Ranked hits: `{OUT_HITS.relative_to(ROOT)}`",
            f"- Effect plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    expr, coverage = collect_expression()
    expr = add_within_sample_gene_ranks(expr)
    stats_df, hits = compute_gene_stats(expr)
    plot_gene_effects(stats_df)
    expr.to_csv(OUT_EXPR, sep="\t", index=False)
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    stats_df.to_csv(OUT_STATS, sep="\t", index=False)
    hits.to_csv(OUT_HITS, sep="\t", index=False)
    write_report(expr, coverage, stats_df, hits)
    print(f"Wrote {len(expr):,} pseudobulk expression rows")
    print(f"Wrote {len(stats_df):,} candidate-gene statistics")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
