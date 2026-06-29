#!/usr/bin/env python3
"""Compute per-cell marker-panel scores from local granule-cell datasets.

This script intentionally streams only configured marker genes. It provides a
lightweight bridge between the fast dataset audit and a later full Seurat/SCE
workflow with normalization, cell annotation, and trajectory analysis.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import tarfile
from collections import defaultdict
from pathlib import Path

from extract_marker_panel_expression import (
    CANONICAL,
    EXTERNAL,
    PANELS,
    RESULTS,
    clean_gene,
    load_csv_metadata_groups,
    load_gse104323_clusters,
    load_gse214309_groups,
    load_two_column_groups,
    open_text,
    rel,
    resolve_gene_candidates,
    resolve_gene_norm,
    tar_member_text,
)


OUTPUT_PREFIX = os.environ.get("MODULE_SCORE_PREFIX", "per_cell_marker_module")
PER_CELL_OUT = RESULTS / f"{OUTPUT_PREFIX}_scores.tsv.gz"
SUMMARY_OUT = RESULTS / f"{OUTPUT_PREFIX}_score_summary.tsv"

PANEL_BY_GENE: dict[str, list[str]] = defaultdict(list)
for panel_name, genes in PANELS.items():
    for gene_norm in genes:
        PANEL_BY_GENE[gene_norm].append(panel_name)


PER_CELL_FIELDS = [
    "dataset",
    "sample",
    "cell_id",
    "group",
    "species",
    "region",
    "platform",
    "panel",
    "panel_gene_count",
    "genes_found_in_matrix",
    "genes_detected_in_cell",
    "detection_fraction_panel",
    "sum_expression",
    "mean_expression_panel",
    "mean_log1p_expression_panel",
    "mean_log1p_expression_found_genes",
    "source_path",
]

SUMMARY_FIELDS = [
    "dataset",
    "sample",
    "group",
    "species",
    "region",
    "platform",
    "panel",
    "n_cells_or_spots",
    "panel_gene_count",
    "genes_found_in_matrix",
    "mean_genes_detected_per_cell",
    "mean_detection_fraction_panel",
    "mean_sum_expression",
    "mean_expression_panel",
    "mean_log1p_expression_panel",
    "mean_log1p_expression_found_genes",
    "source_path",
]


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def read_lines(path: Path) -> list[str]:
    with open_text(path) as fh:
        return [line.rstrip("\n") for line in fh if line.rstrip("\n")]


def read_tar_lines(tf: tarfile.TarFile, member_name: str) -> list[str]:
    with tar_member_text(tf, member_name) as fh:
        return [line.rstrip("\n") for line in fh if line.rstrip("\n")]


def init_scores(n_cells: int) -> dict[str, dict[str, object]]:
    scores: dict[str, dict[str, object]] = {}
    for panel in PANELS:
        scores[panel] = {
            "sum": [0.0] * n_cells,
            "logsum": [0.0] * n_cells,
            "detected": [0] * n_cells,
            "found_genes": set(),
        }
    return scores


def add_gene_value(scores: dict[str, dict[str, object]], cell_idx: int, gene_norm: str, value: float) -> None:
    for panel in PANEL_BY_GENE.get(gene_norm, []):
        scores[panel]["found_genes"].add(gene_norm)  # type: ignore[index, union-attr]
        if value == 0:
            continue
        scores[panel]["sum"][cell_idx] += value  # type: ignore[index]
        scores[panel]["logsum"][cell_idx] += math.log1p(value)  # type: ignore[index]
        scores[panel]["detected"][cell_idx] += 1  # type: ignore[index]


def mark_gene_found(scores: dict[str, dict[str, object]], gene_norm: str) -> None:
    for panel in PANEL_BY_GENE.get(gene_norm, []):
        scores[panel]["found_genes"].add(gene_norm)  # type: ignore[index, union-attr]


def new_summary_acc() -> dict[str, float]:
    return {
        "n": 0,
        "genes_detected": 0.0,
        "detection_fraction": 0.0,
        "sum_expression": 0.0,
        "mean_expression": 0.0,
        "mean_log1p": 0.0,
        "mean_log1p_found": 0.0,
    }


def write_score_rows(
    writer: csv.DictWriter,
    summary_acc: dict[tuple[str, str, str, str, str, str, str, str, int, int, str], dict[str, float]],
    dataset: str,
    sample: str,
    cells: list[str],
    groups: list[str],
    species: str,
    region: str,
    platform: str,
    source_path: str,
    scores: dict[str, dict[str, object]],
) -> None:
    for cell_idx, cell_id in enumerate(cells):
        group = groups[cell_idx] if cell_idx < len(groups) else "all"
        for panel in sorted(PANELS):
            data = scores[panel]
            panel_gene_count = len(PANELS[panel])
            found_count = len(data["found_genes"])  # type: ignore[arg-type]
            genes_detected = int(data["detected"][cell_idx])  # type: ignore[index]
            sum_expression = float(data["sum"][cell_idx])  # type: ignore[index]
            logsum = float(data["logsum"][cell_idx])  # type: ignore[index]
            detection_fraction = genes_detected / panel_gene_count if panel_gene_count else 0.0
            mean_expression = sum_expression / panel_gene_count if panel_gene_count else 0.0
            mean_log1p = logsum / panel_gene_count if panel_gene_count else 0.0
            mean_log1p_found = logsum / found_count if found_count else 0.0

            writer.writerow(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "cell_id": cell_id,
                    "group": group,
                    "species": species,
                    "region": region,
                    "platform": platform,
                    "panel": panel,
                    "panel_gene_count": panel_gene_count,
                    "genes_found_in_matrix": found_count,
                    "genes_detected_in_cell": genes_detected,
                    "detection_fraction_panel": f"{detection_fraction:.6g}",
                    "sum_expression": f"{sum_expression:.6g}",
                    "mean_expression_panel": f"{mean_expression:.6g}",
                    "mean_log1p_expression_panel": f"{mean_log1p:.6g}",
                    "mean_log1p_expression_found_genes": f"{mean_log1p_found:.6g}",
                    "source_path": source_path,
                }
            )

            key = (
                dataset,
                sample,
                group,
                species,
                region,
                platform,
                panel,
                source_path,
                panel_gene_count,
                found_count,
                source_path,
            )
            acc = summary_acc.setdefault(key, new_summary_acc())
            acc["n"] += 1
            acc["genes_detected"] += genes_detected
            acc["detection_fraction"] += detection_fraction
            acc["sum_expression"] += sum_expression
            acc["mean_expression"] += mean_expression
            acc["mean_log1p"] += mean_log1p
            acc["mean_log1p_found"] += mean_log1p_found


def process_wide_table(
    writer: csv.DictWriter,
    summary_acc: dict,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    species: str,
    region: str,
    platform: str,
    group_map: dict[str, str] | None = None,
    header_has_gene_col: bool = True,
) -> None:
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        cells = [cell.strip().strip('"') for cell in (header[1:] if header_has_gene_col else header)]
        groups = [group_map.get(cell, "unmapped") for cell in cells] if group_map else ["all"] * len(cells)
        scores = init_scores(len(cells))
        for row in reader:
            if not row:
                continue
            gene_norm = resolve_gene_norm(row[0])
            if not gene_norm:
                continue
            mark_gene_found(scores, gene_norm)
            for cell_idx, value in enumerate(row[1:]):
                add_gene_value(scores, cell_idx, gene_norm, parse_float(value))
    write_score_rows(writer, summary_acc, dataset, sample, cells, groups, species, region, platform, rel(path), scores)


def process_wide_table_from_tar(
    writer: csv.DictWriter,
    summary_acc: dict,
    dataset: str,
    sample: str,
    tar_path: Path,
    member_name: str,
    delimiter: str,
    species: str,
    region: str,
    platform: str,
) -> None:
    with tarfile.open(tar_path) as tf:
        with tar_member_text(tf, member_name) as fh:
            reader = csv.reader(fh, delimiter=delimiter)
            header = next(reader)
            cells = [cell.strip().strip('"') for cell in header[1:]]
            scores = init_scores(len(cells))
            for row in reader:
                if not row:
                    continue
                gene_norm = resolve_gene_norm(row[0])
                if not gene_norm:
                    continue
                mark_gene_found(scores, gene_norm)
                for cell_idx, value in enumerate(row[1:]):
                    add_gene_value(scores, cell_idx, gene_norm, parse_float(value))
    source = f"{rel(tar_path)}:{member_name}"
    write_score_rows(writer, summary_acc, dataset, sample, cells, ["all"] * len(cells), species, region, platform, source, scores)


def process_obs_by_gene_table(
    writer: csv.DictWriter,
    summary_acc: dict,
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    species: str,
    region: str,
    platform: str,
    group_map: dict[str, str] | None = None,
) -> None:
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        target_cols: list[tuple[int, str]] = []
        for col_idx, gene_raw in enumerate(header[1:], start=1):
            gene_norm = resolve_gene_norm(gene_raw)
            if gene_norm:
                target_cols.append((col_idx, gene_norm))

        cells: list[str] = []
        groups: list[str] = []
        per_cell_values: list[dict[str, dict[str, float | int | set[str]]]] = []
        for row in reader:
            if not row:
                continue
            cell = row[0].strip().strip('"')
            cells.append(cell)
            groups.append(group_map.get(cell, "unmapped") if group_map else "all")
            scores = init_scores(1)
            for col_idx, gene_norm in target_cols:
                mark_gene_found(scores, gene_norm)
                if col_idx < len(row):
                    add_gene_value(scores, 0, gene_norm, parse_float(row[col_idx]))
            per_cell_values.append(scores)

    scores_all = init_scores(len(cells))
    for cell_idx, cell_scores in enumerate(per_cell_values):
        for panel in PANELS:
            scores_all[panel]["found_genes"].update(cell_scores[panel]["found_genes"])  # type: ignore[index, union-attr]
            scores_all[panel]["sum"][cell_idx] = cell_scores[panel]["sum"][0]  # type: ignore[index]
            scores_all[panel]["logsum"][cell_idx] = cell_scores[panel]["logsum"][0]  # type: ignore[index]
            scores_all[panel]["detected"][cell_idx] = cell_scores[panel]["detected"][0]  # type: ignore[index]
    write_score_rows(writer, summary_acc, dataset, sample, cells, groups, species, region, platform, rel(path), scores_all)


def read_10x_features(path: Path) -> dict[int, str]:
    target_rows: dict[int, str] = {}
    with open_text(path) as fh:
        for row_idx, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            gene_norm = resolve_gene_candidates([part for part in parts[:2] if part])
            if gene_norm:
                target_rows[row_idx] = gene_norm
    return target_rows


def read_10x_features_from_tar(tf: tarfile.TarFile, member_name: str) -> dict[int, str]:
    target_rows: dict[int, str] = {}
    with tar_member_text(tf, member_name) as fh:
        for row_idx, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            gene_norm = resolve_gene_candidates([part for part in parts[:2] if part])
            if gene_norm:
                target_rows[row_idx] = gene_norm
    return target_rows


def process_10x_matrix(
    writer: csv.DictWriter,
    summary_acc: dict,
    dataset: str,
    sample: str,
    matrix: Path,
    features: Path,
    barcodes: Path,
    species: str,
    region: str,
    platform: str,
    groups: list[str] | None = None,
) -> None:
    cells = read_lines(barcodes)
    groups = groups or ["all"] * len(cells)
    target_rows = read_10x_features(features)
    scores = init_scores(len(cells))
    for gene_norm in set(target_rows.values()):
        mark_gene_found(scores, gene_norm)

    shape_seen = False
    with open_text(matrix) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if not shape_seen:
                shape_seen = True
                continue
            feature_idx = int(parts[0])
            gene_norm = target_rows.get(feature_idx)
            if not gene_norm:
                continue
            cell_idx = int(parts[1]) - 1
            add_gene_value(scores, cell_idx, gene_norm, parse_float(parts[2]))
    write_score_rows(writer, summary_acc, dataset, sample, cells, groups, species, region, platform, rel(matrix), scores)


def process_10x_matrix_from_tar(
    writer: csv.DictWriter,
    summary_acc: dict,
    dataset: str,
    sample: str,
    tar_path: Path,
    matrix_member: str,
    features_member: str,
    barcodes_member: str,
    species: str,
    region: str,
    platform: str,
) -> None:
    with tarfile.open(tar_path) as tf:
        cells = read_tar_lines(tf, barcodes_member)
        target_rows = read_10x_features_from_tar(tf, features_member)
        scores = init_scores(len(cells))
        for gene_norm in set(target_rows.values()):
            mark_gene_found(scores, gene_norm)

        shape_seen = False
        with tar_member_text(tf, matrix_member) as fh:
            for line in fh:
                if line.startswith("%"):
                    continue
                parts = line.strip().split()
                if not shape_seen:
                    shape_seen = True
                    continue
                gene_norm = target_rows.get(int(parts[0]))
                if not gene_norm:
                    continue
                add_gene_value(scores, int(parts[1]) - 1, gene_norm, parse_float(parts[2]))
    source = f"{rel(tar_path)}:{matrix_member}"
    write_score_rows(writer, summary_acc, dataset, sample, cells, ["all"] * len(cells), species, region, platform, source, scores)


def load_ordered_metadata_groups(path: Path, n_rows: int, preferred_cols: list[str]) -> list[str]:
    groups: list[str] = []
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = ""
            for col in preferred_cols:
                if row.get(col):
                    label = row[col]
                    break
            groups.append(label or "metadata_group")
    if len(groups) != n_rows:
        return ["metadata_unmatched"] * n_rows
    return groups


def write_summary(summary_acc: dict) -> None:
    rows = []
    for key, acc in sorted(summary_acc.items()):
        dataset, sample, group, species, region, platform, panel, _source_key, panel_gene_count, found_count, source_path = key
        n = int(acc["n"])
        rows.append(
            {
                "dataset": dataset,
                "sample": sample,
                "group": group,
                "species": species,
                "region": region,
                "platform": platform,
                "panel": panel,
                "n_cells_or_spots": n,
                "panel_gene_count": panel_gene_count,
                "genes_found_in_matrix": found_count,
                "mean_genes_detected_per_cell": f"{acc['genes_detected'] / n if n else 0:.6g}",
                "mean_detection_fraction_panel": f"{acc['detection_fraction'] / n if n else 0:.6g}",
                "mean_sum_expression": f"{acc['sum_expression'] / n if n else 0:.6g}",
                "mean_expression_panel": f"{acc['mean_expression'] / n if n else 0:.6g}",
                "mean_log1p_expression_panel": f"{acc['mean_log1p'] / n if n else 0:.6g}",
                "mean_log1p_expression_found_genes": f"{acc['mean_log1p_found'] / n if n else 0:.6g}",
                "source_path": source_path,
            }
        )
    with SUMMARY_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary_acc: dict = {}
    with gzip.open(PER_CELL_OUT, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=PER_CELL_FIELDS)
        writer.writeheader()

        gse104_groups = load_gse104323_clusters(EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz")
        process_wide_table(
            writer,
            summary_acc,
            "GSE104323",
            "10X_all_cells",
            EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
            "\t",
            "mouse",
            "dentate_gyrus",
            "10x",
            gse104_groups,
        )

        gse957_groups = load_two_column_groups(EXTERNAL / "GEO/GSE95752/GSE95752_datasetB_cellid_to_sampleid.txt.gz", "sample")
        process_wide_table(
            writer,
            summary_acc,
            "GSE95752",
            "C1_all_cells",
            EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz",
            "\t",
            "mouse",
            "dentate_gyrus",
            "C1",
            gse957_groups,
        )

        process_wide_table(
            writer,
            summary_acc,
            "GSE214905",
            "patch_RNA_QC_counts",
            EXTERNAL / "GEO/GSE214905/GSE214905_Data-counts.tsv.gz",
            "\t",
            "mouse",
            "dentate_gyrus",
            "patch_seq",
        )

        gse214309_groups = load_gse214309_groups(EXTERNAL / "GEO/GSE214309/GSE214309_series_matrix.txt.gz")
        process_wide_table(
            writer,
            summary_acc,
            "GSE214309",
            "snRNA_counts",
            EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
            ",",
            "mouse",
            "dentate_gyrus",
            "snRNA_seq",
            gse214309_groups,
            header_has_gene_col=False,
        )

        gse292_groups = load_csv_metadata_groups(
            EXTERNAL / "GEO/GSE292261/GSE292261_sample_data_SS2_filtered.csv.gz",
            ["Sample", "louvain", "Leiden", "Run"],
        )
        process_obs_by_gene_table(
            writer,
            summary_acc,
            "GSE292261",
            "SS2_filtered_counts",
            EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
            ",",
            "mouse",
            "dentate_gyrus",
            "Smart_seq2",
            gse292_groups,
        )

        gse122_tar = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
        for member in ["GSM3464549_P0.csv.gz", "GSM3464550_P8a.csv.gz", "GSM3464551_P8b.csv.gz"]:
            process_wide_table_from_tar(
                writer,
                summary_acc,
                "GSE122357",
                member.replace(".csv.gz", ""),
                gse122_tar,
                member,
                ",",
                "mouse",
                "cerebellum",
                "single_cell_RNA_seq",
            )

        gse165 = EXTERNAL / "GEO/GSE165657"
        process_10x_matrix(
            writer,
            summary_acc,
            "GSE165657",
            "Cerebellum_aggr",
            gse165 / "GSE165657_Cerebellum_aggr_matrix.mtx.gz",
            gse165 / "GSE165657_Cerebellum_aggr_genes.tsv.gz",
            gse165 / "GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
            "human",
            "cerebellum",
            "10x",
        )

        gse312 = EXTERNAL / "GEO/GSE312658"
        for sample, prefix in [("Ctrl", "GSM9350909_Ctrl"), ("cKO", "GSM9350910_cKO")]:
            process_10x_matrix(
                writer,
                summary_acc,
                "GSE312658",
                sample,
                gse312 / f"{prefix}_matrix.mtx.gz",
                gse312 / f"{prefix}_features.tsv.gz",
                gse312 / f"{prefix}_barcodes.tsv.gz",
                "mouse",
                "cerebellum",
                "10x",
            )

        gse150 = EXTERNAL / "GEO/GSE150153"
        for sample, prefix in [("NAY6153A1_125", "GSM4524697_NAY6153A1_125"), ("NAY6153A2_678", "GSM4524699_NAY6153A2_678")]:
            barcodes = gse150 / f"{prefix}_barcodes.tsv.gz"
            groups = load_ordered_metadata_groups(gse150 / f"{prefix}_metadata.tsv.gz", len(read_lines(barcodes)), ["hash_ID", "hto_classification"])
            process_10x_matrix(
                writer,
                summary_acc,
                "GSE150153",
                sample,
                gse150 / f"{prefix}_matrix.mtx.gz",
                gse150 / f"{prefix}_genes.tsv.gz",
                barcodes,
                "human",
                "organoid",
                "10x",
                groups,
            )

        gse242_tar = EXTERNAL / "Proteomics/GSE242688/GSE242688_RAW.tar"
        for sample, matrix_member, features_member, barcodes_member in [
            (
                "WT_ZT12_rep1",
                "GSM7767079_WT_ZT12_rep1_matrix.mtx.gz",
                "GSM7767079_WT_ZT12_rep1_features.tsv.gz",
                "GSM7767079_WT_ZT12_rep1_barcodes.tsv.gz",
            ),
            (
                "WT_ZT12_rep2",
                "GSM7767080_WT_ZT12_rep2_matrix.mtx.gz",
                "GSM7767080_WT_ZT12_rep2_features.tsv.gz",
                "GSM7767080_WT_ZT12_rep2_barcodes.tsv.gz",
            ),
        ]:
            process_10x_matrix_from_tar(
                writer,
                summary_acc,
                "GSE242688",
                sample,
                gse242_tar,
                matrix_member,
                features_member,
                barcodes_member,
                "mouse",
                "cerebellum_spatial",
                "visium_spatial",
            )

    write_summary(summary_acc)
    print(f"Wrote {PER_CELL_OUT}")
    print(f"Wrote {SUMMARY_OUT}")


if __name__ == "__main__":
    main()
