#!/usr/bin/env python3
"""Convert first-pass human DG/hippocampus seed text matrices to sparse objects."""

from __future__ import annotations

import csv
import gzip
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
OUT_ROOT = ROOT / "Project" / "processed" / "human_seed_sparse_objects"
INVENTORY = RESULTS / "human_seed_archive_inventory.tsv"
SUMMARY_TSV = RESULTS / "human_seed_sparse_object_summary.tsv"
SUMMARY_MD = RESULTS / "human_seed_sparse_object_summary.md"


ARCHIVE_BY_DATASET = {
    "GSE185277": ROOT / "External_Data/GEO/GSE185277/GSE185277_RAW.tar",
    "GSE185553": ROOT / "External_Data/GEO/GSE185553/GSE185553_RAW.tar",
}


def clean_id(text: str) -> str:
    text = re.sub(r"\.(dge|deg)?\.?txt\.gz$", "", text)
    text = re.sub(r"[^A-Za-z0-9]+", "_", text)
    return text.strip("_")


def open_member_text(tar: tarfile.TarFile, member_name: str) -> io.TextIOWrapper:
    raw = tar.extractfile(member_name)
    if raw is None:
        raise ValueError(f"Could not read {member_name}")
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace")
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")


def split_header(line: str, delimiter_hint: str) -> list[str]:
    line = line.strip()
    if delimiter_hint == "tab":
        fields = [field.strip('"') for field in line.split("\t")]
        if fields and fields[0].upper() == "GENE":
            return fields[1:]
        return fields
    return [field.strip('"') for field in line.split()]


def split_gene_and_values(line: str, delimiter_hint: str) -> tuple[str, str]:
    line = line.strip()
    if delimiter_hint == "tab":
        gene, values = line.split("\t", 1)
        return gene.strip('"'), values
    gene, values = line.split(None, 1)
    return gene.strip('"'), values


def parse_values(values: str, delimiter_hint: str, expected_n: int) -> np.ndarray:
    sep = "\t" if delimiter_hint == "tab" else " "
    arr = np.fromstring(values, dtype=np.int32, sep=sep)
    if arr.size != expected_n:
        raise ValueError(f"Expected {expected_n} values, parsed {arr.size}")
    return arr


def write_tsv_gz(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def convert_member(dataset: str, member_name: str, delimiter_hint: str, archive: Path) -> dict[str, object]:
    library_id = clean_id(member_name)
    out_dir = OUT_ROOT / dataset / library_id
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = out_dir / "matrix_cells_by_genes.npz"
    barcodes_path = out_dir / "barcodes.tsv.gz"
    genes_path = out_dir / "genes.tsv.gz"
    cell_meta_path = out_dir / "cell_metadata.tsv.gz"
    gene_meta_path = out_dir / "gene_metadata.tsv.gz"

    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    genes: list[str] = []
    gene_counts: list[int] = []
    gene_detected_cells: list[int] = []

    mode = "r:gz" if archive.name.endswith(".tar.gz") else "r:"
    with tarfile.open(archive, mode) as tar:
        with open_member_text(tar, member_name) as fh:
            header_line = fh.readline()
            barcodes = split_header(header_line, delimiter_hint)
            n_cells = len(barcodes)
            cell_counts = np.zeros(n_cells, dtype=np.int64)
            cell_genes = np.zeros(n_cells, dtype=np.int32)

            for gene_idx, line in enumerate(fh):
                if not line.strip():
                    continue
                gene, values = split_gene_and_values(line, delimiter_hint)
                arr = parse_values(values, delimiter_hint, n_cells)
                nz = np.nonzero(arr)[0]
                if nz.size:
                    row_chunks.append(nz.astype(np.int32))
                    col_chunks.append(np.full(nz.size, gene_idx, dtype=np.int32))
                    data_chunks.append(arr[nz].astype(np.int32))
                    cell_counts += arr
                    cell_genes += (arr > 0)
                genes.append(gene)
                gene_counts.append(int(arr.sum()))
                gene_detected_cells.append(int(nz.size))

    if data_chunks:
        row = np.concatenate(row_chunks)
        col = np.concatenate(col_chunks)
        data = np.concatenate(data_chunks)
    else:
        row = np.array([], dtype=np.int32)
        col = np.array([], dtype=np.int32)
        data = np.array([], dtype=np.int32)

    matrix = sparse.coo_matrix((data, (row, col)), shape=(len(barcodes), len(genes))).tocsr()
    sparse.save_npz(matrix_path, matrix, compressed=True)

    cell_ids = [f"{dataset}:{library_id}:{barcode}" for barcode in barcodes]
    write_tsv_gz(barcodes_path, ["cell_id", "barcode"], [[cid, barcode] for cid, barcode in zip(cell_ids, barcodes)])
    write_tsv_gz(genes_path, ["gene"], [[gene] for gene in genes])
    write_tsv_gz(
        cell_meta_path,
        ["cell_id", "dataset", "library_id", "barcode", "source_member", "n_counts", "n_genes"],
        [
            [cid, dataset, library_id, barcode, member_name, int(n_counts), int(n_genes)]
            for cid, barcode, n_counts, n_genes in zip(cell_ids, barcodes, cell_counts, cell_genes)
        ],
    )
    write_tsv_gz(
        gene_meta_path,
        ["gene", "dataset", "library_id", "source_member", "n_counts", "n_cells"],
        [
            [gene, dataset, library_id, member_name, n_counts, n_cells]
            for gene, n_counts, n_cells in zip(genes, gene_counts, gene_detected_cells)
        ],
    )

    return {
        "dataset": dataset,
        "library_id": library_id,
        "source_member": member_name,
        "delimiter_hint": delimiter_hint,
        "n_cells": len(barcodes),
        "n_genes": len(genes),
        "nnz": int(matrix.nnz),
        "total_counts": int(matrix.sum()),
        "matrix_path": str(matrix_path.relative_to(ROOT)),
        "barcodes_path": str(barcodes_path.relative_to(ROOT)),
        "genes_path": str(genes_path.relative_to(ROOT)),
        "cell_metadata_path": str(cell_meta_path.relative_to(ROOT)),
        "gene_metadata_path": str(gene_meta_path.relative_to(ROOT)),
    }


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(INVENTORY, sep="\t")
    expression = inventory.loc[inventory["kind"] == "expression_table"].copy()

    rows: list[dict[str, object]] = []
    for record in expression.to_dict(orient="records"):
        dataset = record["dataset"]
        member = record["member"]
        delimiter_hint = record["delimiter_hint"]
        archive = ARCHIVE_BY_DATASET[dataset]
        print(f"Converting {dataset} {member}", flush=True)
        rows.append(convert_member(dataset, member, delimiter_hint, archive))

    fieldnames = [
        "dataset",
        "library_id",
        "source_member",
        "delimiter_hint",
        "n_cells",
        "n_genes",
        "nnz",
        "total_counts",
        "matrix_path",
        "barcodes_path",
        "genes_path",
        "cell_metadata_path",
        "gene_metadata_path",
    ]
    with SUMMARY_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = pd.DataFrame(rows)
    lines = [
        "# Human seed sparse object build",
        "",
        "Date built: 2026-06-21",
        "",
        "## Summary",
        "",
    ]
    for dataset, group in summary.groupby("dataset"):
        lines.append(
            f"- `{dataset}`: {len(group)} libraries, {int(group['n_cells'].sum())} cells/barcodes, "
            f"{int(group['nnz'].sum())} non-zero counts, {int(group['total_counts'].sum())} total counts."
        )
    lines.extend(
        [
            "",
            "## Output Format",
            "",
            "Each library has a compressed SciPy sparse matrix in cells-by-genes orientation, plus gzipped barcodes, genes, cell metadata, and gene metadata files. This is intentionally neutral: it can be imported into AnnData, Seurat, or SingleCellExperiment later.",
            "",
            "## Immediate Use",
            "",
            "- Use per-library QC first; do not merge libraries until sample labels and donor/age metadata are curated.",
            "- Prefix cell IDs with dataset and library ID to avoid barcode collisions.",
            "- Treat this as the first constructed human dentate/hippocampal branch for marker/QC work, not yet the final harmonized object.",
            "",
            f"Detailed table: `{SUMMARY_TSV.relative_to(ROOT)}`",
            "",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines))
    print(f"Wrote {SUMMARY_TSV}")
    print(f"Wrote {SUMMARY_MD}")


if __name__ == "__main__":
    main()
