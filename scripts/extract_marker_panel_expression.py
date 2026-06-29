#!/usr/bin/env python3
"""Extract focused marker-panel expression summaries from local matrices."""

from __future__ import annotations

import csv
import gzip
import io
import os
import sys
import tarfile
from collections import defaultdict
from pathlib import Path


csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parents[2]
CONFIG = Path(os.environ.get("GRANULE_MARKER_PANEL_CONFIG", ROOT / "Project" / "config" / "granule_marker_panels.tsv"))
RESULTS = ROOT / "Project" / "results"
EXTERNAL = ROOT / "External_Data"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def clean_gene(gene: str) -> str:
    return gene.strip().strip('"').strip("'")


def norm_gene(gene: str) -> str:
    return clean_gene(gene).lower()


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


def load_panels() -> tuple[dict[str, set[str]], dict[str, str]]:
    panels: dict[str, set[str]] = defaultdict(set)
    canonical: dict[str, str] = {}
    with CONFIG.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            panel = row["panel"]
            gene = clean_gene(row["gene"])
            ngene = norm_gene(gene)
            panels[panel].add(ngene)
            canonical.setdefault(ngene, gene)
    return dict(panels), canonical


PANELS, CANONICAL = load_panels()
TARGET_GENES = set(CANONICAL)


def add_feature_aliases(aliases: dict[str, str], parts: list[str]) -> None:
    candidates = [norm_gene(part) for part in parts[:2] if part]
    target = next((candidate for candidate in candidates if candidate in TARGET_GENES), None)
    if not target:
        return
    for candidate in candidates:
        aliases.setdefault(candidate, target)


def build_gene_aliases() -> dict[str, str]:
    aliases = {gene_norm: gene_norm for gene_norm in TARGET_GENES}
    for path in sorted(EXTERNAL.rglob("*")):
        if not path.is_file() or not path.name.endswith(("_features.tsv.gz", "_genes.tsv.gz")):
            continue
        with open_text(path) as fh:
            for line in fh:
                add_feature_aliases(aliases, line.rstrip("\n").split("\t"))

    gse242_tar = EXTERNAL / "Proteomics/GSE242688/GSE242688_RAW.tar"
    if gse242_tar.exists():
        with tarfile.open(gse242_tar) as tf:
            for member in tf.getnames():
                if not member.endswith("_features.tsv.gz"):
                    continue
                with tar_member_text(tf, member) as fh:
                    for line in fh:
                        add_feature_aliases(aliases, line.rstrip("\n").split("\t"))
    return aliases


GENE_ALIASES = build_gene_aliases()


def resolve_gene_norm(gene: str) -> str | None:
    return GENE_ALIASES.get(norm_gene(gene))


def resolve_gene_candidates(candidates: list[str]) -> str | None:
    for candidate in candidates:
        resolved = resolve_gene_norm(candidate)
        if resolved:
            return resolved
    return None


def gene_panels(normed_gene: str) -> list[str]:
    return [panel for panel, genes in PANELS.items() if normed_gene in genes]


def load_gse104323_clusters(path: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cell = row.get("Sample name (24185 single cells)", "")
            cluster = row.get("characteristics: cell cluster", "")
            if cell and cluster:
                groups[cell] = cluster
    return groups


def load_two_column_groups(path: Path, default_label: str) -> dict[str, str]:
    groups: dict[str, str] = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                groups[parts[0]] = parts[1] or default_label
    return groups


def load_csv_metadata_groups(path: Path, preferred_cols: list[str]) -> dict[str, str]:
    groups: dict[str, str] = {}
    with open_text(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            cell = row.get("") or row.get("cell") or row.get("Cell") or row.get("cell_id") or row.get("barcode")
            if not cell:
                continue
            label = ""
            for col in preferred_cols:
                if row.get(col):
                    label = row[col]
                    break
            groups[cell] = label or "metadata_group"
    return groups


def load_gse214309_groups(path: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0] != "!Sample_title":
                continue
            for title in row[1:]:
                parts = [part.strip() for part in clean_gene(title).split(",")]
                if len(parts) < 3:
                    continue
                cell = parts[-1]
                group = "_".join(part.lower().replace(" ", "_") for part in parts[:2])
                groups[cell] = group
            break
    return groups


def empty_gene_stats(dataset: str, sample: str, source: str, group: str, gene_norm: str, observations: int, format_name: str) -> dict[str, str]:
    return {
        "dataset": dataset,
        "sample": sample,
        "group": group,
        "format": format_name,
        "source_path": source,
        "gene": CANONICAL.get(gene_norm, gene_norm),
        "gene_norm": gene_norm,
        "observations": str(observations),
        "nonzero_observations": "0",
        "detection_fraction": "0",
        "total_expression": "0",
        "mean_all_observations": "0",
        "mean_nonzero_observations": "0",
    }


def finalize_row(row: dict[str, str], total: float, nonzero: int, observations: int):
    row["nonzero_observations"] = str(nonzero)
    row["total_expression"] = f"{total:.6g}"
    row["detection_fraction"] = f"{(nonzero / observations) if observations else 0:.6g}"
    row["mean_all_observations"] = f"{(total / observations) if observations else 0:.6g}"
    row["mean_nonzero_observations"] = f"{(total / nonzero) if nonzero else 0:.6g}"


def extract_wide_table(
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    group_map: dict[str, str] | None = None,
    source_override: str | None = None,
    header_has_gene_col: bool = True,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open_text(path) as fh:
        rows.extend(_extract_wide_reader(dataset, sample, fh, delimiter, group_map, source_override or rel(path), header_has_gene_col))
    return rows


def extract_wide_table_from_tar(
    dataset: str,
    sample: str,
    tar_path: Path,
    member_name: str,
    delimiter: str,
    group_map: dict[str, str] | None = None,
    header_has_gene_col: bool = True,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with tarfile.open(tar_path) as tf:
        with tar_member_text(tf, member_name) as fh:
            rows.extend(
                _extract_wide_reader(
                    dataset,
                    sample,
                    fh,
                    delimiter,
                    group_map,
                    f"{rel(tar_path)}:{member_name}",
                    header_has_gene_col,
                )
            )
    return rows


def _extract_wide_reader(
    dataset: str,
    sample: str,
    fh,
    delimiter: str,
    group_map: dict[str, str] | None,
    source: str,
    header_has_gene_col: bool,
) -> list[dict[str, str]]:
    reader = csv.reader(fh, delimiter=delimiter)
    header = next(reader)
    cell_fields = header[1:] if header_has_gene_col else header
    cells = [cell.strip().strip('"') for cell in cell_fields]
    labels = ["all"] * len(cells)
    if group_map:
        labels = [group_map.get(cell, "unmapped") for cell in cells]
    group_counts: dict[str, int] = {"all": len(cells)}
    if group_map:
        for label in labels:
            group_counts[label] = group_counts.get(label, 0) + 1

    result_rows: list[dict[str, str]] = []
    for row in reader:
        if not row:
            continue
        gene_raw = clean_gene(row[0])
        gene_norm = resolve_gene_norm(gene_raw)
        if not gene_norm:
            continue
        totals = {group: 0.0 for group in group_counts}
        nonzeros = {group: 0 for group in group_counts}
        for i, value in enumerate(row[1:]):
            try:
                val = float(value)
            except ValueError:
                val = 0.0
            if val == 0:
                continue
            totals["all"] += val
            nonzeros["all"] += 1
            if group_map and i < len(labels):
                label = labels[i]
                totals[label] += val
                nonzeros[label] += 1
        for group, obs in group_counts.items():
            out = empty_gene_stats(dataset, sample, source, group, gene_norm, obs, "wide_gene_by_observation_table")
            out["source_gene_symbol"] = gene_raw
            finalize_row(out, totals.get(group, 0.0), nonzeros.get(group, 0), obs)
            result_rows.append(out)
    return result_rows


def extract_obs_by_gene_table(
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    group_map: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    source = rel(path)
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        target_cols: list[tuple[int, str, str]] = []
        for idx, gene_raw in enumerate(header[1:], start=1):
            gene_norm = resolve_gene_norm(gene_raw)
            if gene_norm:
                target_cols.append((idx, gene_norm, clean_gene(gene_raw)))

        group_counts: dict[str, int] = {"all": 0}
        totals: dict[tuple[str, str], float] = defaultdict(float)
        nonzeros: dict[tuple[str, str], int] = defaultdict(int)
        source_gene_by_norm: dict[str, str] = {}
        for _, gene_norm, gene_raw in target_cols:
            source_gene_by_norm.setdefault(gene_norm, gene_raw)

        for row in reader:
            if not row:
                continue
            cell = row[0].strip().strip('"')
            label = group_map.get(cell, "unmapped") if group_map else ""
            group_counts["all"] += 1
            if group_map:
                group_counts[label] = group_counts.get(label, 0) + 1
            for col_idx, gene_norm, _ in target_cols:
                if col_idx >= len(row):
                    continue
                try:
                    val = float(row[col_idx])
                except ValueError:
                    val = 0.0
                if val == 0:
                    continue
                totals[("all", gene_norm)] += val
                nonzeros[("all", gene_norm)] += 1
                if group_map:
                    totals[(label, gene_norm)] += val
                    nonzeros[(label, gene_norm)] += 1

    result_rows: list[dict[str, str]] = []
    target_genes = sorted({gene_norm for _, gene_norm, _ in target_cols})
    for gene_norm in target_genes:
        for group, obs in group_counts.items():
            out = empty_gene_stats(dataset, sample, source, group, gene_norm, obs, "wide_observation_by_gene_table")
            out["source_gene_symbol"] = source_gene_by_norm.get(gene_norm, CANONICAL.get(gene_norm, gene_norm))
            finalize_row(out, totals.get((group, gene_norm), 0.0), nonzeros.get((group, gene_norm), 0), obs)
            result_rows.append(out)
    return result_rows


def read_10x_features(path: Path) -> dict[int, tuple[str, str]]:
    features: dict[int, tuple[str, str]] = {}
    with open_text(path) as fh:
        for idx, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            candidates = [p for p in parts[:2] if p]
            ngene = resolve_gene_candidates(candidates)
            if ngene:
                features[idx] = (CANONICAL[ngene], ngene)
    return features


def read_10x_features_from_tar(tf: tarfile.TarFile, member_name: str) -> dict[int, tuple[str, str]]:
    features: dict[int, tuple[str, str]] = {}
    with tar_member_text(tf, member_name) as fh:
        for idx, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            candidates = [p for p in parts[:2] if p]
            ngene = resolve_gene_candidates(candidates)
            if ngene:
                features[idx] = (CANONICAL[ngene], ngene)
    return features


def extract_10x_matrix(dataset: str, sample: str, matrix: Path, features: Path) -> list[dict[str, str]]:
    target_rows = read_10x_features(features)
    result_rows: list[dict[str, str]] = []
    with open_text(matrix) as fh:
        result_rows.extend(_extract_mtx_stream(dataset, sample, fh, target_rows, rel(matrix), "10x_matrix_market"))
    return result_rows


def extract_10x_matrix_from_tar(
    dataset: str,
    sample: str,
    tar_path: Path,
    matrix_member: str,
    features_member: str,
    format_name: str = "10x_matrix_market",
) -> list[dict[str, str]]:
    result_rows: list[dict[str, str]] = []
    with tarfile.open(tar_path) as tf:
        target_rows = read_10x_features_from_tar(tf, features_member)
        with tar_member_text(tf, matrix_member) as fh:
            result_rows.extend(
                _extract_mtx_stream(
                    dataset,
                    sample,
                    fh,
                    target_rows,
                    f"{rel(tar_path)}:{matrix_member}",
                    format_name,
                )
            )
    return result_rows


def _extract_mtx_stream(
    dataset: str,
    sample: str,
    fh,
    target_rows: dict[int, tuple[str, str]],
    source: str,
    format_name: str,
) -> list[dict[str, str]]:
    n_features = n_obs = n_nonzero = None
    totals: dict[str, float] = defaultdict(float)
    nonzeros: dict[str, int] = defaultdict(int)
    shape_seen = False
    for line in fh:
        if line.startswith("%"):
            continue
        parts = line.strip().split()
        if not shape_seen:
            n_features, n_obs, n_nonzero = map(int, parts[:3])
            shape_seen = True
            continue
        feature_idx = int(parts[0])
        if feature_idx not in target_rows:
            continue
        value = float(parts[2])
        if value == 0:
            continue
        _, gene_norm = target_rows[feature_idx]
        totals[gene_norm] += value
        nonzeros[gene_norm] += 1

    obs = int(n_obs or 0)
    result_rows: list[dict[str, str]] = []
    for gene_norm in sorted({gene_norm for _, gene_norm in target_rows.values()}):
        out = empty_gene_stats(dataset, sample, source, "all", gene_norm, obs, format_name)
        out["source_gene_symbol"] = CANONICAL.get(gene_norm, gene_norm)
        finalize_row(out, totals.get(gene_norm, 0.0), nonzeros.get(gene_norm, 0), obs)
        result_rows.append(out)
    return result_rows


def add_panel_rows(gene_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    expanded: list[dict[str, str]] = []
    for row in gene_rows:
        for panel in gene_panels(row["gene_norm"]):
            out = dict(row)
            out["panel"] = panel
            expanded.append(out)
    return expanded


def summarize_panels(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["dataset"], row["sample"], row["group"], row["panel"])].append(row)

    summaries: list[dict[str, str]] = []
    for (dataset, sample, group, panel), panel_rows in sorted(grouped.items()):
        panel_genes = PANELS[panel]
        found = {row["gene_norm"] for row in panel_rows}
        detections = [float(row["detection_fraction"]) for row in panel_rows]
        means = [float(row["mean_all_observations"]) for row in panel_rows]
        summaries.append(
            {
                "dataset": dataset,
                "sample": sample,
                "group": group,
                "panel": panel,
                "panel_gene_count": str(len(panel_genes)),
                "genes_found": str(len(found)),
                "coverage_fraction": f"{len(found) / len(panel_genes):.6g}",
                "mean_detection_fraction_found_genes": f"{(sum(detections) / len(detections)) if detections else 0:.6g}",
                "mean_expression_found_genes": f"{(sum(means) / len(means)) if means else 0:.6g}",
                "found_genes": ",".join(CANONICAL[g] for g in sorted(found)),
                "missing_genes": ",".join(CANONICAL[g] for g in sorted(panel_genes - found)),
            }
        )
    return summaries


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    gene_rows: list[dict[str, str]] = []

    gse104_groups = load_gse104323_clusters(EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz")
    gene_rows += extract_wide_table(
        "GSE104323",
        "10X_all_cells",
        EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz",
        "\t",
        gse104_groups,
    )

    gse957_groups = load_two_column_groups(EXTERNAL / "GEO/GSE95752/GSE95752_datasetB_cellid_to_sampleid.txt.gz", "sample")
    gene_rows += extract_wide_table(
        "GSE95752",
        "C1_all_cells",
        EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz",
        "\t",
        gse957_groups,
    )

    gene_rows += extract_wide_table(
        "GSE214905",
        "patch_RNA_QC_counts",
        EXTERNAL / "GEO/GSE214905/GSE214905_Data-counts.tsv.gz",
        "\t",
    )
    gse214309_groups = load_gse214309_groups(EXTERNAL / "GEO/GSE214309/GSE214309_series_matrix.txt.gz")
    gene_rows += extract_wide_table(
        "GSE214309",
        "snRNA_counts",
        EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
        ",",
        gse214309_groups,
        header_has_gene_col=False,
    )

    gse292_groups = load_csv_metadata_groups(
        EXTERNAL / "GEO/GSE292261/GSE292261_sample_data_SS2_filtered.csv.gz",
        ["Sample", "louvain", "Leiden", "Run"],
    )
    gene_rows += extract_obs_by_gene_table(
        "GSE292261",
        "SS2_filtered_counts",
        EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
        ",",
        gse292_groups,
    )

    gene_rows += extract_10x_matrix(
        "GSE165657",
        "Cerebellum_aggr",
        EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_matrix.mtx.gz",
        EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_genes.tsv.gz",
    )

    for sample, prefix in [("Ctrl", "GSM9350909_Ctrl"), ("cKO", "GSM9350910_cKO")]:
        gene_rows += extract_10x_matrix(
            "GSE312658",
            sample,
            EXTERNAL / f"GEO/GSE312658/{prefix}_matrix.mtx.gz",
            EXTERNAL / f"GEO/GSE312658/{prefix}_features.tsv.gz",
        )

    for sample, prefix in [("NAY6153A1_125", "GSM4524697_NAY6153A1_125"), ("NAY6153A2_678", "GSM4524699_NAY6153A2_678")]:
        gene_rows += extract_10x_matrix(
            "GSE150153",
            sample,
            EXTERNAL / f"GEO/GSE150153/{prefix}_matrix.mtx.gz",
            EXTERNAL / f"GEO/GSE150153/{prefix}_genes.tsv.gz",
        )

    gse122_tar = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
    for member in ["GSM3464549_P0.csv.gz", "GSM3464550_P8a.csv.gz", "GSM3464551_P8b.csv.gz"]:
        gene_rows += extract_wide_table_from_tar(
            "GSE122357",
            member.replace(".csv.gz", ""),
            gse122_tar,
            member,
            ",",
        )

    gse242_tar = EXTERNAL / "Proteomics/GSE242688/GSE242688_RAW.tar"
    for sample, matrix_member, features_member in [
        ("WT_ZT12_rep1", "GSM7767079_WT_ZT12_rep1_matrix.mtx.gz", "GSM7767079_WT_ZT12_rep1_features.tsv.gz"),
        ("WT_ZT12_rep2", "GSM7767080_WT_ZT12_rep2_matrix.mtx.gz", "GSM7767080_WT_ZT12_rep2_features.tsv.gz"),
    ]:
        gene_rows += extract_10x_matrix_from_tar(
            "GSE242688",
            sample,
            gse242_tar,
            matrix_member,
            features_member,
            "10x_matrix_market_spatial",
        )

    panel_gene_rows = add_panel_rows(gene_rows)
    panel_summaries = summarize_panels(panel_gene_rows)

    gene_fields = [
        "dataset",
        "sample",
        "group",
        "panel",
        "gene",
        "source_gene_symbol",
        "format",
        "observations",
        "nonzero_observations",
        "detection_fraction",
        "total_expression",
        "mean_all_observations",
        "mean_nonzero_observations",
        "source_path",
    ]
    summary_fields = [
        "dataset",
        "sample",
        "group",
        "panel",
        "panel_gene_count",
        "genes_found",
        "coverage_fraction",
        "mean_detection_fraction_found_genes",
        "mean_expression_found_genes",
        "found_genes",
        "missing_genes",
    ]
    write_tsv(RESULTS / "marker_gene_expression_summary.tsv", panel_gene_rows, gene_fields)
    write_tsv(RESULTS / "marker_panel_expression_summary.tsv", panel_summaries, summary_fields)
    print(f"Wrote {len(panel_gene_rows)} marker gene rows")
    print(f"Wrote {len(panel_summaries)} panel summary rows")


if __name__ == "__main__":
    main()
