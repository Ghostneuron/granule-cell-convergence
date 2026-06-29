#!/usr/bin/env python3
"""Audit local granule-cell project datasets without loading full matrices."""

from __future__ import annotations

import csv
import gzip
import io
import os
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
EXTERNAL = ROOT / "External_Data"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def byte_size(path: Path) -> int:
    return path.stat().st_size


def human_size(n: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(n)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f}{unit}" if unit != "B" else f"{int(size)}B"
        size /= 1024
    return f"{n}B"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def wide_table_shape(path: Path, delimiter: str, header_has_gene_col: bool = True) -> tuple[int, int]:
    rows = 0
    cols = 0
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        for i, row in enumerate(reader):
            if i == 0:
                cols = max(0, len(row) - 1 if header_has_gene_col else len(row))
            else:
                rows += 1
    return rows, cols


def obs_by_gene_table_shape(path: Path, delimiter: str) -> tuple[int, int]:
    observations = 0
    genes = 0
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        for i, row in enumerate(reader):
            if i == 0:
                genes = max(0, len(row) - 1)
            else:
                observations += 1
    return genes, observations


def mtx_shape(path: Path) -> tuple[int, int, int]:
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith("%"):
                parts = line.strip().split()
                return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"No Matrix Market shape line found in {path}")


def count_lines(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh)


def tar_member_text(tf: tarfile.TarFile, member_name: str):
    raw = tf.extractfile(member_name)
    if raw is None:
        raise FileNotFoundError(member_name)
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")
    return io.TextIOWrapper(raw, newline="")


def tar_wide_table_shape(tar_path: Path, member_name: str, delimiter: str) -> tuple[int, int]:
    rows = 0
    cols = 0
    with tarfile.open(tar_path) as tf:
        with tar_member_text(tf, member_name) as fh:
            reader = csv.reader(fh, delimiter=delimiter)
            for i, row in enumerate(reader):
                if i == 0:
                    cols = max(0, len(row) - 1)
                else:
                    rows += 1
    return rows, cols


def tar_mtx_shape(tar_path: Path, member_name: str) -> tuple[int, int, int]:
    with tarfile.open(tar_path) as tf:
        with tar_member_text(tf, member_name) as fh:
            for line in fh:
                if not line.startswith("%"):
                    parts = line.strip().split()
                    return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"No Matrix Market shape line found in {tar_path}:{member_name}")


def infer_file_kind(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".mtx.gz"):
        return "matrix_market"
    if name.endswith((".csv.gz", ".tsv.gz", ".txt.gz", ".tab.gz")):
        return "compressed_table"
    if name.endswith(".tar"):
        return "tar_archive"
    if name.endswith(".xlsx"):
        return "spreadsheet"
    if name.endswith(".rds.gz"):
        return "compressed_rds"
    if name.endswith(".txt"):
        return "text"
    if name.endswith(".html"):
        return "html"
    return "other"


def discover_files() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(EXTERNAL.rglob("*")):
        if not path.is_file() or path.name == ".DS_Store":
            continue
        accession = next((p for p in path.parts if p.startswith("GSE") or p.startswith("Allen_Institute")), "")
        size = byte_size(path)
        rows.append(
            {
                "accession_or_source": accession,
                "path": rel(path),
                "kind": infer_file_kind(path),
                "bytes": str(size),
                "human_size": human_size(size),
            }
        )
    return rows


def add_10x_matrix(rows: list[dict[str, str]], dataset: str, sample: str, matrix: Path, features: Path | None, barcodes: Path | None):
    genes = ""
    cells = ""
    nnz = ""
    notes = []
    try:
        n_features, n_obs, n_nonzero = mtx_shape(matrix)
        genes = str(n_features)
        cells = str(n_obs)
        nnz = str(n_nonzero)
    except Exception as exc:  # keep audit running across mixed downloads
        notes.append(f"matrix_error={exc}")
    if features and features.exists():
        try:
            notes.append(f"features_lines={count_lines(features)}")
        except Exception as exc:
            notes.append(f"features_error={exc}")
    if barcodes and barcodes.exists():
        try:
            notes.append(f"barcodes_lines={count_lines(barcodes)}")
        except Exception as exc:
            notes.append(f"barcodes_error={exc}")
    rows.append(
        {
            "dataset": dataset,
            "sample": sample,
            "format": "10x_matrix_market",
            "source_path": rel(matrix),
            "features_or_genes": genes,
            "observations": cells,
            "nonzero_entries": nnz,
            "notes": ";".join(notes),
        }
    )


def add_wide_table(
    rows: list[dict[str, str]],
    dataset: str,
    sample: str,
    path: Path,
    delimiter: str,
    note: str = "",
    header_has_gene_col: bool = True,
):
    try:
        genes, cells = wide_table_shape(path, delimiter, header_has_gene_col)
        notes = note
    except Exception as exc:
        genes, cells = "", ""
        notes = f"{note};error={exc}" if note else f"error={exc}"
    rows.append(
        {
            "dataset": dataset,
            "sample": sample,
            "format": "wide_gene_by_observation_table",
            "source_path": rel(path),
            "features_or_genes": str(genes),
            "observations": str(cells),
            "nonzero_entries": "",
            "notes": notes,
        }
    )


def add_obs_by_gene_table(rows: list[dict[str, str]], dataset: str, sample: str, path: Path, delimiter: str, note: str = ""):
    try:
        genes, cells = obs_by_gene_table_shape(path, delimiter)
        notes = note
    except Exception as exc:
        genes, cells = "", ""
        notes = f"{note};error={exc}" if note else f"error={exc}"
    rows.append(
        {
            "dataset": dataset,
            "sample": sample,
            "format": "wide_observation_by_gene_table",
            "source_path": rel(path),
            "features_or_genes": str(genes),
            "observations": str(cells),
            "nonzero_entries": "",
            "notes": notes,
        }
    )


def build_matrix_audit() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    add_wide_table(rows, "GSE104323", "10X_all_cells", EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz", "\t")
    add_wide_table(rows, "GSE95752", "C1_all_cells", EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz", "\t")
    add_wide_table(rows, "GSE214905", "patch_RNA_QC_counts", EXTERNAL / "GEO/GSE214905/GSE214905_Data-counts.tsv.gz", "\t")
    add_wide_table(
        rows,
        "GSE214309",
        "snRNA_counts",
        EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz",
        ",",
        "comma_delimited;ensembl_gene_ids",
        header_has_gene_col=False,
    )
    add_obs_by_gene_table(
        rows,
        "GSE292261",
        "SS2_filtered_counts",
        EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz",
        ",",
        "cell_by_gene_orientation",
    )

    gse165 = EXTERNAL / "GEO/GSE165657"
    add_10x_matrix(
        rows,
        "GSE165657",
        "Cerebellum_aggr",
        gse165 / "GSE165657_Cerebellum_aggr_matrix.mtx.gz",
        gse165 / "GSE165657_Cerebellum_aggr_genes.tsv.gz",
        gse165 / "GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
    )

    gse312 = EXTERNAL / "GEO/GSE312658"
    for sample, prefix in [("Ctrl", "GSM9350909_Ctrl"), ("cKO", "GSM9350910_cKO")]:
        add_10x_matrix(
            rows,
            "GSE312658",
            sample,
            gse312 / f"{prefix}_matrix.mtx.gz",
            gse312 / f"{prefix}_features.tsv.gz",
            gse312 / f"{prefix}_barcodes.tsv.gz",
        )

    gse150 = EXTERNAL / "GEO/GSE150153"
    for sample, prefix in [("NAY6153A1_125", "GSM4524697_NAY6153A1_125"), ("NAY6153A2_678", "GSM4524699_NAY6153A2_678")]:
        add_10x_matrix(
            rows,
            "GSE150153",
            sample,
            gse150 / f"{prefix}_matrix.mtx.gz",
            gse150 / f"{prefix}_genes.tsv.gz",
            gse150 / f"{prefix}_barcodes.tsv.gz",
        )

    gse122_tar = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
    if gse122_tar.exists():
        for member in ["GSM3464549_P0.csv.gz", "GSM3464550_P8a.csv.gz", "GSM3464551_P8b.csv.gz"]:
            try:
                genes, cells = tar_wide_table_shape(gse122_tar, member, ",")
                notes = "inside_tar"
            except Exception as exc:
                genes, cells, notes = "", "", f"inside_tar;error={exc}"
            rows.append(
                {
                    "dataset": "GSE122357",
                    "sample": member.replace(".csv.gz", ""),
                    "format": "wide_gene_by_observation_table",
                    "source_path": f"{rel(gse122_tar)}:{member}",
                    "features_or_genes": str(genes),
                    "observations": str(cells),
                    "nonzero_entries": "",
                    "notes": notes,
                }
            )

    gse242_tar = EXTERNAL / "Proteomics/GSE242688/GSE242688_RAW.tar"
    if gse242_tar.exists():
        for member in [
            "GSM7767079_WT_ZT12_rep1_matrix.mtx.gz",
            "GSM7767080_WT_ZT12_rep2_matrix.mtx.gz",
        ]:
            try:
                genes, spots, nnz = tar_mtx_shape(gse242_tar, member)
                notes = "inside_tar;10x_visium_spatial"
            except Exception as exc:
                genes, spots, nnz, notes = "", "", "", f"inside_tar;error={exc}"
            rows.append(
                {
                    "dataset": "GSE242688",
                    "sample": member.replace("_matrix.mtx.gz", ""),
                    "format": "10x_matrix_market_spatial",
                    "source_path": f"{rel(gse242_tar)}:{member}",
                    "features_or_genes": str(genes),
                    "observations": str(spots),
                    "nonzero_entries": str(nnz),
                    "notes": notes,
                }
            )

    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    file_rows = discover_files()
    matrix_rows = build_matrix_audit()
    write_tsv(
        RESULTS / "dataset_file_inventory.tsv",
        file_rows,
        ["accession_or_source", "path", "kind", "bytes", "human_size"],
    )
    write_tsv(
        RESULTS / "matrix_dimension_audit.tsv",
        matrix_rows,
        ["dataset", "sample", "format", "source_path", "features_or_genes", "observations", "nonzero_entries", "notes"],
    )
    print(f"Wrote {len(file_rows)} file rows")
    print(f"Wrote {len(matrix_rows)} matrix rows")


if __name__ == "__main__":
    main()
