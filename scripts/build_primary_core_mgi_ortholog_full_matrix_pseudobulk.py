#!/usr/bin/env python3
"""Full-matrix pseudobulk screen using MGI one-to-one ortholog symbols.

This expands the earlier same-symbol full-matrix screen by resolving mouse
matrix rows through the MGI mouse symbol side of one-to-one human-mouse
homology classes. Human matrices are resolved through the human symbol side.

The output canonical gene is always the human symbol. Mouse-only source symbols
that differ from their human ortholog are now retained instead of being missed
or incorrectly assigned by same-symbol matching.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import build_primary_core_candidate_gene_pseudobulk as base
import build_primary_core_genomewide_symbol_pseudobulk as same_symbol


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"

MGI_MAP = RESULTS / "primary_core_mgi_ortholog_meta_model_map.tsv"
HUMAN_DG_FULL = PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates"
GENE_METADATA = HUMAN_DG_FULL / "gene_metadata.tsv.gz"

OUT_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"
OUT_COVERAGE = RESULTS / "primary_core_mgi_ortholog_full_matrix_coverage.tsv"
OUT_STATS = RESULTS / "primary_core_mgi_ortholog_full_matrix_statistics.tsv"
OUT_SHARED = RESULTS / "primary_core_mgi_ortholog_full_matrix_shared_hits.tsv"
OUT_BRANCH = RESULTS / "primary_core_mgi_ortholog_full_matrix_branch_specific.tsv"
OUT_NONIDENTICAL = RESULTS / "primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv"
OUT_PLOT = RESULTS / "primary_core_mgi_ortholog_full_matrix_shared_hits.png"
OUT_MD = RESULTS / "primary_core_mgi_ortholog_full_matrix_analysis.md"

MIN_CLASS_CELLS = 20

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def canon_gene(gene: object) -> str:
    if pd.isna(gene):
        return ""
    return str(gene).strip().strip('"').strip("'").upper()


def rel(path: Path | str) -> str:
    return base.rel(path)


def open_text(path: Path):
    return base.open_text(path)


def tar_member_text(tf: tarfile.TarFile, member_name: str):
    raw = tf.extractfile(member_name)
    if raw is None:
        raise FileNotFoundError(member_name)
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")
    return io.TextIOWrapper(raw, newline="")


def load_targets() -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    targets = pd.read_csv(MGI_MAP, sep="\t", dtype={"db_class_key": str}).copy()
    targets = targets.loc[targets["mgi_one_to_one_human_mouse"].astype(str).str.lower().eq("true")].copy()
    targets["canonical_gene"] = targets["canonical_gene"].map(canon_gene)
    targets["canonical_human_symbol"] = targets["human_symbol"].map(canon_gene)
    targets["canonical_mouse_symbol"] = targets["mouse_symbol"].map(canon_gene)
    targets["ortholog_symbol_class"] = np.where(
        targets["canonical_human_symbol"].eq(targets["canonical_mouse_symbol"]),
        "same_symbol",
        "nonidentical_symbol",
    )
    targets["panel"] = "mgi_one_to_one_ortholog"
    targets["candidate_role"] = "mgi_one_to_one_ortholog_gene"
    targets["support_tier"] = np.where(
        targets["ortholog_symbol_class"].eq("same_symbol"),
        "mgi_one_to_one_same_symbol",
        "mgi_one_to_one_nonidentical_symbol",
    )
    human_dg = pd.read_csv(GENE_METADATA, sep="\t")
    human_dg["canonical_gene"] = human_dg["gene"].map(canon_gene)
    targets = targets.merge(
        human_dg[["canonical_gene", "n_counts", "n_cells"]].rename(
            columns={"n_counts": "human_dg_full_n_counts", "n_cells": "human_dg_full_n_cells"}
        ),
        on="canonical_gene",
        how="left",
    )
    targets["present_in_human_dg_full_gene_list"] = targets["human_dg_full_n_cells"].notna()

    human_map = dict(zip(targets["canonical_human_symbol"], targets["canonical_gene"]))
    mouse_map = dict(zip(targets["canonical_mouse_symbol"], targets["canonical_gene"]))

    target_cols = [
        "gene",
        "canonical_gene",
        "human_symbol",
        "mouse_symbol",
        "panel",
        "candidate_role",
        "support_tier",
        "db_class_key",
        "canonical_mouse_symbol",
        "ortholog_symbol_class",
        "human_entrez_id",
        "mouse_entrez_id",
        "human_hgnc_id",
        "mouse_mgi_id",
        "human_dg_full_n_counts",
        "human_dg_full_n_cells",
        "present_in_human_dg_full_gene_list",
    ]
    targets["gene"] = targets["human_symbol"]
    targets = targets[target_cols].drop_duplicates("canonical_gene").sort_values("canonical_gene")

    base.TARGETS = targets
    base.TARGET_META = targets.set_index("canonical_gene").to_dict("index")
    base.TARGET_GENES = set(targets["canonical_gene"])
    return targets, human_map, mouse_map


def resolve_source_gene(parts: list[str], source_species: str, human_map: dict[str, str], mouse_map: dict[str, str]) -> str | None:
    source_map = human_map if source_species == "human" else mouse_map
    for part in parts[:3]:
        canonical_source = canon_gene(part)
        if canonical_source in source_map:
            return source_map[canonical_source]
    return None


def finalize_rows(
    *,
    dataset: str,
    sample: str,
    source_layer: str,
    expression_scope: str,
    expression_scale: str,
    source_path: str,
    source_species: str,
    groups: Counter[str],
    totals: dict[str, dict[str, float]],
    log_totals: dict[str, dict[str, float]],
    nonzeros: dict[str, dict[str, int]],
    source_gene_symbol: dict[str, str],
    present_genes: set[str],
) -> list[dict[str, object]]:
    rows = same_symbol.finalize_rows(
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
    for row in rows:
        row["ortholog_source_species"] = source_species
    return rows


def coverage_row(
    *,
    dataset: str,
    sample: str,
    source_layer: str,
    expression_scope: str,
    expression_scale: str,
    source_path: str,
    source_species: str,
    n_matrix_observations: int,
    n_labeled_observations: int,
    n_unmapped_observations: int,
    n_gene_rows_scanned: int,
    present_genes: set[str],
    groups: Counter[str],
    matrix_nonzero_entries: int | None = None,
) -> dict[str, object]:
    targets = base.TARGETS
    nonidentical = set(targets.loc[targets["ortholog_symbol_class"].eq("nonidentical_symbol"), "canonical_gene"])
    present_nonidentical = present_genes & nonidentical
    row = {
        "dataset": dataset,
        "sample": sample,
        "source_layer": source_layer,
        "expression_scope": expression_scope,
        "expression_scale": expression_scale,
        "source_species": source_species,
        "source_path": source_path,
        "n_matrix_observations": int(n_matrix_observations),
        "n_labeled_observations": int(n_labeled_observations),
        "n_unmapped_observations": int(n_unmapped_observations),
        "n_gene_rows_scanned": int(n_gene_rows_scanned),
        "target_genes_present": len(present_genes),
        "target_genes_total": len(base.TARGET_GENES),
        "target_gene_coverage_fraction": len(present_genes) / len(base.TARGET_GENES),
        "nonidentical_target_genes_present": len(present_nonidentical),
        "nonidentical_target_genes_total": len(nonidentical),
        "nonidentical_target_gene_coverage_fraction": len(present_nonidentical) / len(nonidentical),
        "broad_class_counts": ";".join(f"{k}:{v}" for k, v in sorted(groups.items())),
    }
    if matrix_nonzero_entries is not None:
        row["matrix_nonzero_entries"] = int(matrix_nonzero_entries)
    return row


def extract_wide_by_cell_table_chunked_mapped(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_mgi_one_to_one_ortholog",
    expression_scale: str = "raw_counts_or_reported_counts",
    chunksize: int = 64,
) -> tuple[list[dict[str, object]], dict[str, object]]:
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
        canonical = gene_symbols.map(lambda gene: resolve_source_gene([gene], source_species, human_map, mouse_map))
        keep = canonical.notna()
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
        source_species=source_species,
        groups=groups,
        totals=totals,
        log_totals=log_totals,
        nonzeros=nonzeros,
        source_gene_symbol=source_gene_symbol,
        present_genes=present_genes,
    )
    coverage = coverage_row(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=rel(path),
        source_species=source_species,
        n_matrix_observations=len(cells),
        n_labeled_observations=sum(v for k, v in groups.items() if k != "unmapped"),
        n_unmapped_observations=groups.get("unmapped", 0),
        n_gene_rows_scanned=rows_scanned,
        present_genes=present_genes,
        groups=groups,
    )
    return rows, coverage


def extract_wide_by_cell_table_mapped(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_mgi_one_to_one_ortholog",
    expression_scale: str = "raw_counts_or_reported_counts",
    header_has_gene_col: bool = True,
    source_path: str | None = None,
    tar_member: str | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    label_map = base.BACKBONE_LABELS.get((dataset, sample), {})
    if tar_member:
        with tarfile.open(path) as tf:
            with tar_member_text(tf, tar_member) as fh:
                return _extract_wide_reader_mapped(
                    dataset=dataset,
                    sample=sample,
                    fh=fh,
                    delimiter=delimiter,
                    label_map=label_map,
                    source_species=source_species,
                    human_map=human_map,
                    mouse_map=mouse_map,
                    source_layer=source_layer,
                    expression_scope=expression_scope,
                    expression_scale=expression_scale,
                    header_has_gene_col=header_has_gene_col,
                    source_path=source_path or f"{rel(path)}:{tar_member}",
                )
    with open_text(path) as fh:
        return _extract_wide_reader_mapped(
            dataset=dataset,
            sample=sample,
            fh=fh,
            delimiter=delimiter,
            label_map=label_map,
            source_species=source_species,
            human_map=human_map,
            mouse_map=mouse_map,
            source_layer=source_layer,
            expression_scope=expression_scope,
            expression_scale=expression_scale,
            header_has_gene_col=header_has_gene_col,
            source_path=source_path or rel(path),
        )


def _extract_wide_reader_mapped(
    *,
    dataset: str,
    sample: str,
    fh,
    delimiter: str,
    label_map: dict[str, str],
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
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
    totals = {group: defaultdict(float) for group in groups}
    log_totals = {group: defaultdict(float) for group in groups}
    nonzeros = {group: defaultdict(int) for group in groups}
    present_genes: set[str] = set()
    source_gene_symbol: dict[str, str] = {}
    rows_scanned = 0

    for row in reader:
        if not row:
            continue
        rows_scanned += 1
        raw_gene = str(row[0]).strip().strip('"')
        canonical_gene = resolve_source_gene([raw_gene], source_species, human_map, mouse_map)
        if not canonical_gene:
            continue
        present_genes.add(canonical_gene)
        source_gene_symbol.setdefault(canonical_gene, raw_gene)
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
        source_species=source_species,
        groups=groups,
        totals=totals,
        log_totals=log_totals,
        nonzeros=nonzeros,
        source_gene_symbol=source_gene_symbol,
        present_genes=present_genes,
    )
    coverage = coverage_row(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
        source_species=source_species,
        n_matrix_observations=len(cells),
        n_labeled_observations=sum(v for k, v in groups.items() if k != "unmapped"),
        n_unmapped_observations=groups.get("unmapped", 0),
        n_gene_rows_scanned=rows_scanned,
        present_genes=present_genes,
        groups=groups,
    )
    return rows, coverage


def extract_obs_by_gene_table_mapped(
    *,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_mgi_one_to_one_ortholog",
    expression_scale: str = "raw_counts_or_reported_counts",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    label_map = base.BACKBONE_LABELS.get((dataset, sample), {})
    source_path = rel(path)
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        target_cols: list[tuple[int, str, str]] = []
        for idx, gene in enumerate(header[1:], start=1):
            canonical_gene = resolve_source_gene([gene], source_species, human_map, mouse_map)
            if canonical_gene:
                target_cols.append((idx, canonical_gene, str(gene).strip().strip('"')))
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
        source_species=source_species,
        groups=groups,
        totals=totals,
        log_totals=log_totals,
        nonzeros=nonzeros,
        source_gene_symbol=source_gene_symbol,
        present_genes=present_genes,
    )
    coverage = coverage_row(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
        source_species=source_species,
        n_matrix_observations=n_cells,
        n_labeled_observations=sum(v for k, v in groups.items() if k != "unmapped"),
        n_unmapped_observations=groups.get("unmapped", 0),
        n_gene_rows_scanned=len(header) - 1,
        present_genes=present_genes,
        groups=groups,
    )
    return rows, coverage


def read_barcodes(path: Path) -> list[str]:
    return base.read_barcodes(path)


def read_10x_target_features_mapped(
    path: Path,
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
) -> tuple[dict[int, str], dict[str, str], int]:
    target_rows: dict[int, str] = {}
    source_gene_symbol: dict[str, str] = {}
    n_features = 0
    with open_text(path) as fh:
        for idx, line in enumerate(fh, start=1):
            n_features = idx
            parts = line.rstrip("\n").split("\t")
            canonical_gene = resolve_source_gene(parts, source_species, human_map, mouse_map)
            if canonical_gene:
                target_rows[idx] = canonical_gene
                source_gene_symbol.setdefault(canonical_gene, parts[1] if len(parts) > 1 else parts[0])
    return target_rows, source_gene_symbol, n_features


def extract_10x_matrix_mapped(
    *,
    dataset: str,
    sample: str,
    matrix: Path,
    features: Path,
    barcodes: Path,
    source_species: str,
    human_map: dict[str, str],
    mouse_map: dict[str, str],
    source_layer: str = "full_raw_matrix",
    expression_scope: str = "full_matrix_mgi_one_to_one_ortholog",
    expression_scale: str = "raw_counts",
) -> tuple[list[dict[str, object]], dict[str, object]]:
    target_rows, source_gene_symbol, n_features = read_10x_target_features_mapped(
        features, source_species, human_map, mouse_map
    )
    barcode_list = read_barcodes(barcodes)
    label_map = base.BACKBONE_LABELS.get((dataset, sample), {})
    labels = [label_map.get(barcode, "unmapped") for barcode in barcode_list]
    groups = Counter(labels)
    totals = {group: defaultdict(float) for group in groups}
    log_totals = {group: defaultdict(float) for group in groups}
    nonzeros = {group: defaultdict(int) for group in groups}
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
        source_species=source_species,
        groups=groups,
        totals=totals,
        log_totals=log_totals,
        nonzeros=nonzeros,
        source_gene_symbol=source_gene_symbol,
        present_genes=present_genes,
    )
    coverage = coverage_row(
        dataset=dataset,
        sample=sample,
        source_layer=source_layer,
        expression_scope=expression_scope,
        expression_scale=expression_scale,
        source_path=source_path,
        source_species=source_species,
        n_matrix_observations=int(n_obs_seen or len(barcode_list)),
        n_labeled_observations=sum(v for k, v in groups.items() if k != "unmapped"),
        n_unmapped_observations=groups.get("unmapped", 0),
        n_gene_rows_scanned=int(n_features_seen or n_features),
        present_genes=present_genes,
        groups=groups,
        matrix_nonzero_entries=int(n_nonzero_seen or 0),
    )
    return rows, coverage


def gse186538_full_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    labels = pd.read_csv(HUMAN_DG_FULL / "cell_metadata.tsv.gz", sep="\t", low_memory=False)
    labels["sample_for_pb"] = labels["samplename"].fillna("GSE186538").astype(str)
    labels["broad_class"] = "dentate_candidate"
    rows, coverage = base.aggregate_selected_sparse(
        dataset="GSE186538",
        source_layer="full_sparse_subset",
        expression_scope="full_human_dg_mgi_one_to_one_orthologs",
        expression_scale="raw_counts_full_dg_gc_subset",
        matrix_path=HUMAN_DG_FULL / "matrix_cells_by_genes.npz",
        var_path=HUMAN_DG_FULL / "gene_metadata.tsv.gz",
        labels=labels,
        sample_col="sample_for_pb",
        class_col="broad_class",
        source_path=rel(HUMAN_DG_FULL / "matrix_cells_by_genes.npz"),
    )
    for row in rows:
        row["ortholog_source_species"] = "human"
    coverage["source_species"] = "human"
    targets = base.TARGETS
    nonidentical = set(targets.loc[targets["ortholog_symbol_class"].eq("nonidentical_symbol"), "canonical_gene"])
    present = set(row["canonical_gene"] for row in rows)
    coverage["nonidentical_target_genes_present"] = len(present & nonidentical)
    coverage["nonidentical_target_genes_total"] = len(nonidentical)
    coverage["nonidentical_target_gene_coverage_fraction"] = len(present & nonidentical) / len(nonidentical)
    return rows, coverage


def collect_expression(human_map: dict[str, str], mouse_map: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    jobs = [
        lambda: extract_wide_by_cell_table_chunked_mapped(
            dataset="GSE104323",
            sample="10X_all_cells",
            path=EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
            delimiter="\t",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_wide_by_cell_table_chunked_mapped(
            dataset="GSE95752",
            sample="C1_all_cells",
            path=EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz",
            delimiter="\t",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_obs_by_gene_table_mapped(
            dataset="GSE292261",
            sample="SS2_filtered_counts",
            path=EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
            delimiter=",",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_wide_by_cell_table_mapped(
            dataset="GSE214309",
            sample="snRNA_counts",
            path=EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
            delimiter=",",
            header_has_gene_col=False,
            expression_scale="reported_counts_symbol_rows_after_ensembl_preface",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_wide_by_cell_table_mapped(
            dataset="GSE122357",
            sample="GSM3464549_P0",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464549_P0.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464549_P0.csv.gz",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_wide_by_cell_table_mapped(
            dataset="GSE122357",
            sample="GSM3464550_P8a",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464550_P8a.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464550_P8a.csv.gz",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_wide_by_cell_table_mapped(
            dataset="GSE122357",
            sample="GSM3464551_P8b",
            path=EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar",
            delimiter=",",
            tar_member="GSM3464551_P8b.csv.gz",
            source_path="External_Data/GEO/GSE122357/GSE122357_RAW.tar:GSM3464551_P8b.csv.gz",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_10x_matrix_mapped(
            dataset="GSE165657",
            sample="Cerebellum_aggr",
            matrix=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_genes.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
            source_species="human",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_10x_matrix_mapped(
            dataset="GSE312658",
            sample="Ctrl",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350909_Ctrl_barcodes.tsv.gz",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        lambda: extract_10x_matrix_mapped(
            dataset="GSE312658",
            sample="cKO",
            matrix=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_matrix.mtx.gz",
            features=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_features.tsv.gz",
            barcodes=EXTERNAL / "GEO/GSE312658/GSM9350910_cKO_barcodes.tsv.gz",
            source_species="mouse",
            human_map=human_map,
            mouse_map=mouse_map,
        ),
        gse186538_full_rows,
    ]

    for job in jobs:
        job_rows, job_coverage = job()
        rows.extend(job_rows)
        coverage_rows.append(job_coverage)
        print(
            f"{job_coverage['dataset']} {job_coverage['sample']}: "
            f"{job_coverage['target_genes_present']}/{job_coverage['target_genes_total']} MGI targets; "
            f"{job_coverage.get('nonidentical_target_genes_present', 0)}/"
            f"{job_coverage.get('nonidentical_target_genes_total', 0)} nonidentical targets; "
            f"{job_coverage['n_labeled_observations']} labeled observations",
            flush=True,
        )
    return pd.DataFrame(rows), pd.DataFrame(coverage_rows)


def add_ortholog_columns(df: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    ortho_cols = [
        "canonical_gene",
        "db_class_key",
        "human_symbol",
        "mouse_symbol",
        "canonical_mouse_symbol",
        "ortholog_symbol_class",
        "human_entrez_id",
        "mouse_entrez_id",
        "human_hgnc_id",
        "mouse_mgi_id",
        "present_in_human_dg_full_gene_list",
    ]
    keep = targets[ortho_cols].drop_duplicates("canonical_gene")
    overlap = [col for col in keep.columns if col in df.columns and col != "canonical_gene"]
    df = df.drop(columns=overlap, errors="ignore")
    return df.merge(keep, on="canonical_gene", how="left")


def plot_shared(shared: pd.DataFrame) -> None:
    plot_df = shared.head(34).copy()
    if plot_df.empty:
        return
    colors = np.where(plot_df["ortholog_symbol_class"].eq("nonidentical_symbol"), "#c75c2b", "#168b7a")
    fig, ax = plt.subplots(figsize=(10.4, 8.0))
    y = np.arange(len(plot_df))
    ax.barh(y - 0.18, plot_df["dentate_rank_delta_vs_background"], height=0.34, color=colors, label="dentate candidate")
    ax.barh(y + 0.18, plot_df["cerebellar_rank_delta_vs_background"], height=0.34, color="#6d3bbd", label="cerebellar candidate")
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Median within-sample rank delta versus local background")
    ax.set_title("MGI one-to-one full-matrix shared pseudobulk hits")
    ax.legend(frameon=False)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(
    targets: pd.DataFrame,
    expr: pd.DataFrame,
    coverage: pd.DataFrame,
    stats_df: pd.DataFrame,
    shared: pd.DataFrame,
    branch: pd.DataFrame,
    nonidentical_hits: pd.DataFrame,
) -> None:
    contrast_datasets = sorted(expr.loc[expr["eligible_class"], "dataset"].unique())
    same_symbol_hits = int(shared["ortholog_symbol_class"].eq("same_symbol").sum())
    nonidentical_shared = int(shared["ortholog_symbol_class"].eq("nonidentical_symbol").sum())
    lines = [
        "# MGI One-to-One Ortholog Full-Matrix Pseudobulk Screen",
        "",
        "Date built: 2026-06-22",
        "",
        "## Scope",
        "",
        "This screen expands the previous full-matrix same-symbol analysis by resolving mouse matrix rows through the MGI mouse-symbol side of one-to-one human-mouse homology classes. The canonical output gene is the human symbol.",
        "",
        "It remains a rank-based pseudobulk screen, not the final mixed-effect DE model. It is designed to make the current gene-level evidence ortholog-aware before formal modeling.",
        "",
        "## Coverage",
        "",
        f"- MGI one-to-one target genes: {targets['canonical_gene'].nunique():,}.",
        f"- Non-identical human/mouse symbol targets: {int(targets['ortholog_symbol_class'].eq('nonidentical_symbol').sum()):,}.",
        f"- Pseudobulk expression rows: {len(expr):,}.",
        f"- Datasets with expression represented: {coverage.loc[coverage['target_genes_present'].gt(0), 'dataset'].nunique()}/10 primary datasets.",
        f"- Datasets contributing to rank contrasts: {len(contrast_datasets)} ({', '.join(contrast_datasets)}).",
        f"- Genes tested in contrast statistics: {stats_df['canonical_gene'].nunique():,}.",
        f"- Shared-positive ortholog genes: {len(shared):,}.",
        f"- Shared-positive same-symbol genes: {same_symbol_hits:,}.",
        f"- Shared-positive non-identical-symbol genes: {nonidentical_shared:,}.",
        f"- Shared-positive genes passing BH<0.10 in both branches: {int(shared['shared_strict_bh_0_10'].sum()):,}.",
        f"- Branch-specific genes: {len(branch):,}.",
        "",
    ]
    for _, row in coverage.sort_values(["dataset", "sample"]).iterrows():
        note = ""
        if row["target_genes_present"] == 0:
            note = " Ensembl-to-symbol mapping is still needed for this source."
        lines.append(
            f"- `{row['dataset']}` / `{row['sample']}`: "
            f"{int(row['target_genes_present'])}/{int(row['target_genes_total'])} MGI targets, "
            f"{int(row.get('nonidentical_target_genes_present', 0))}/{int(row.get('nonidentical_target_genes_total', 0))} non-identical targets, "
            f"{int(row['n_labeled_observations'])}/{int(row['n_matrix_observations'])} labeled observations "
            f"(`{row['source_species']}`).{note}"
        )
    lines.extend(["", "## Top Shared Hits", ""])
    for _, row in shared.head(30).iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['ortholog_symbol_class']}; mouse `{row['mouse_symbol']}`): "
            f"dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
            f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}, "
            f"BH<0.10 both branches={bool(row['shared_strict_bh_0_10'])}."
        )
    lines.extend(["", "## Non-Identical Symbol Shared Hits", ""])
    if nonidentical_hits.empty:
        lines.append("- No non-identical-symbol shared hits passed the shared-positive rule.")
    else:
        for _, row in nonidentical_hits.head(30).iterrows():
            lines.append(
                f"- `{row['gene']}` / mouse `{row['mouse_symbol']}`: dentate delta "
                f"{row['dentate_rank_delta_vs_background']:.3f}, cerebellar delta "
                f"{row['cerebellar_rank_delta_vs_background']:.3f}."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This pass corrects the main limitation of the earlier same-symbol full-matrix screen: mouse genes with non-identical human ortholog symbols can now contribute.",
            "- `GSE214309` contributes symbol-resolved rows in this pass, although the source file starts with Ensembl-style rows. A dedicated Ensembl-to-symbol bridge would still be useful to rescue any residual Ensembl-only rows.",
            "- The next step is to re-run the dataset-aware meta-model using this MGI ortholog full-matrix expression layer instead of the same-symbol layer.",
            "",
            "## Outputs",
            "",
            f"- Expression table: `{OUT_EXPR.relative_to(ROOT)}`",
            f"- Coverage table: `{OUT_COVERAGE.relative_to(ROOT)}`",
            f"- Statistics table: `{OUT_STATS.relative_to(ROOT)}`",
            f"- Shared hits: `{OUT_SHARED.relative_to(ROOT)}`",
            f"- Branch-specific hits: `{OUT_BRANCH.relative_to(ROOT)}`",
            f"- Non-identical-symbol hits: `{OUT_NONIDENTICAL.relative_to(ROOT)}`",
            f"- Shared-hit plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    targets, human_map, mouse_map = load_targets()
    expr, coverage = collect_expression(human_map, mouse_map)
    expr = same_symbol.add_strict_ranks(expr)
    expr = add_ortholog_columns(expr, targets)
    stats_df, shared, branch = same_symbol.compute_stats(expr, targets)
    stats_df = add_ortholog_columns(stats_df, targets)
    shared = add_ortholog_columns(shared, targets)
    branch = add_ortholog_columns(branch, targets)
    nonidentical_hits = shared.loc[shared["ortholog_symbol_class"].eq("nonidentical_symbol")].copy()
    nonidentical_hits = nonidentical_hits.sort_values(
        ["shared_strict_bh_0_10", "combined_rank_delta", "minimum_branch_detection"],
        ascending=[False, False, False],
    )
    plot_shared(shared)
    expr.to_csv(OUT_EXPR, sep="\t", index=False, compression="gzip")
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    stats_df.to_csv(OUT_STATS, sep="\t", index=False)
    shared.to_csv(OUT_SHARED, sep="\t", index=False)
    branch.to_csv(OUT_BRANCH, sep="\t", index=False)
    nonidentical_hits.to_csv(OUT_NONIDENTICAL, sep="\t", index=False)
    write_report(targets, expr, coverage, stats_df, shared, branch, nonidentical_hits)
    print(f"Wrote {len(expr):,} MGI ortholog full-matrix pseudobulk expression rows")
    print(f"Wrote {len(stats_df):,} MGI ortholog gene statistics")
    print(f"Wrote {len(shared):,} shared-positive hits")
    print(f"Wrote {len(nonidentical_hits):,} non-identical-symbol shared hits")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
