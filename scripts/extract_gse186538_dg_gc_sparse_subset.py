#!/usr/bin/env python3
"""Stream-extract GSE186538 human DG granule-cell candidate matrix."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE186538"
RESULTS = ROOT / "Project/results"
OUT_DIR = ROOT / "Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates"

MTX = BASE / "GSE186538_Human_counts.mtx.gz"
GENES = BASE / "GSE186538_Human_genes.txt.gz"
CANDIDATES = RESULTS / "gse186538_human_dg_gc_candidate_cells.tsv"

SUMMARY_TSV = RESULTS / "gse186538_human_dg_gc_sparse_subset_summary.tsv"
SUMMARY_MD = RESULTS / "gse186538_human_dg_gc_sparse_subset_summary.md"


def mtx_dimensions_and_skiprows(path: Path) -> tuple[int, int, int, int]:
    skiprows = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            skiprows += 1
            if line.startswith("%"):
                continue
            parts = line.split()
            return int(parts[0]), int(parts[1]), int(parts[2]), skiprows
    raise ValueError(f"No Matrix Market dimension row found in {path}")


def write_tsv_gz(path: Path, header: list[str], rows) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(CANDIDATES, sep="\t")
    genes = pd.read_csv(GENES, sep="\t", header=None, names=["gene"])
    n_genes_mtx, n_cells_mtx, full_nnz, skiprows = mtx_dimensions_and_skiprows(MTX)

    if n_genes_mtx != len(genes):
        raise ValueError(f"Gene count mismatch: matrix has {n_genes_mtx}, genes file has {len(genes)}")

    selected_cols = candidates["cell_index_1based"].to_numpy(dtype=np.int32)
    if selected_cols.min() < 1 or selected_cols.max() > n_cells_mtx:
        raise ValueError("Candidate cell indices fall outside the Matrix Market column range")

    selected_map = np.full(n_cells_mtx + 1, -1, dtype=np.int32)
    selected_map[selected_cols] = np.arange(len(selected_cols), dtype=np.int32)

    expected_nnz = int(candidates["nFeature_RNA"].sum())
    rows = np.empty(expected_nnz, dtype=np.int32)
    cols = np.empty(expected_nnz, dtype=np.int32)
    data = np.empty(expected_nnz, dtype=np.int32)

    offset = 0
    chunk_size = 5_000_000
    processed_rows = 0

    reader = pd.read_csv(
        MTX,
        sep=" ",
        header=None,
        names=["gene", "cell", "count"],
        skiprows=skiprows,
        compression="gzip",
        chunksize=chunk_size,
        dtype={"gene": np.int32, "cell": np.int32, "count": np.int32},
        engine="c",
    )

    for chunk_idx, chunk in enumerate(reader, start=1):
        processed_rows += len(chunk)
        cell_values = chunk["cell"].to_numpy(dtype=np.int32, copy=False)
        mapped_rows = selected_map[cell_values]
        mask = mapped_rows >= 0
        n_selected = int(mask.sum())

        if n_selected:
            next_offset = offset + n_selected
            if next_offset > len(rows):
                raise RuntimeError(
                    f"Selected entries exceeded expected nFeature_RNA total: {next_offset} > {len(rows)}"
                )
            rows[offset:next_offset] = mapped_rows[mask]
            cols[offset:next_offset] = chunk["gene"].to_numpy(dtype=np.int32, copy=False)[mask] - 1
            data[offset:next_offset] = chunk["count"].to_numpy(dtype=np.int32, copy=False)[mask]
            offset = next_offset

        if chunk_idx == 1 or chunk_idx % 10 == 0:
            print(
                f"chunks={chunk_idx} full_entries={processed_rows} selected_entries={offset}",
                flush=True,
            )

    if processed_rows != full_nnz:
        raise RuntimeError(f"Processed {processed_rows} entries, expected {full_nnz}")
    if offset != expected_nnz:
        raise RuntimeError(f"Selected {offset} entries, expected {expected_nnz} from metadata")

    matrix = sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(len(candidates), len(genes)),
        dtype=np.int32,
    ).tocsr()
    matrix.sum_duplicates()

    matrix_path = OUT_DIR / "matrix_cells_by_genes.npz"
    cell_meta_path = OUT_DIR / "cell_metadata.tsv.gz"
    gene_meta_path = OUT_DIR / "gene_metadata.tsv.gz"
    genes_path = OUT_DIR / "genes.tsv.gz"
    barcodes_path = OUT_DIR / "barcodes.tsv.gz"

    sparse.save_npz(matrix_path, matrix, compressed=True)

    cell_ids = "GSE186538:DG_GC_candidates:" + candidates["cell_name"].astype(str)
    cell_counts = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
    cell_genes = matrix.getnnz(axis=1).astype(np.int32)
    gene_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
    gene_cells = matrix.getnnz(axis=0).astype(np.int32)

    cell_meta = candidates.copy()
    cell_meta.insert(0, "cell_id", cell_ids)
    cell_meta["extracted_n_counts"] = cell_counts
    cell_meta["extracted_n_genes"] = cell_genes
    cell_meta.to_csv(cell_meta_path, sep="\t", index=False, compression="gzip")

    write_tsv_gz(barcodes_path, ["cell_id", "cell_name"], zip(cell_ids, candidates["cell_name"]))
    write_tsv_gz(
        genes_path,
        ["gene_index_1based", "gene"],
        ((idx, gene) for idx, gene in enumerate(genes["gene"], start=1)),
    )

    gene_meta = pd.DataFrame(
        {
            "gene_index_1based": np.arange(1, len(genes) + 1, dtype=np.int32),
            "gene": genes["gene"].astype(str),
            "n_counts": gene_counts,
            "n_cells": gene_cells,
        }
    )
    gene_meta.to_csv(gene_meta_path, sep="\t", index=False, compression="gzip")

    summary_rows = [
        {
            "dataset": "GSE186538",
            "subset": "DG_GC_candidates",
            "source_matrix": str(MTX.relative_to(ROOT)),
            "n_cells": matrix.shape[0],
            "n_genes": matrix.shape[1],
            "nnz": int(matrix.nnz),
            "total_counts": int(matrix.sum()),
            "expected_nnz_from_metadata": expected_nnz,
            "expected_total_counts_from_metadata": int(candidates["nCount_RNA"].sum()),
            "matrix_path": str(matrix_path.relative_to(ROOT)),
            "cell_metadata_path": str(cell_meta_path.relative_to(ROOT)),
            "gene_metadata_path": str(gene_meta_path.relative_to(ROOT)),
        }
    ]
    with SUMMARY_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    cluster_counts = candidates["cluster"].value_counts().to_dict()
    sample_counts = candidates["samplename"].value_counts().to_dict()
    lines = [
        "# GSE186538 human DG GC sparse subset",
        "",
        "Date built: 2026-06-21",
        "",
        "## Summary",
        "",
        f"- Cells: {matrix.shape[0]}",
        f"- Genes: {matrix.shape[1]}",
        f"- Non-zero count entries: {int(matrix.nnz)}",
        f"- Total counts: {int(matrix.sum())}",
        "",
        "## Cluster Counts",
        "",
    ]
    for cluster, n_cells in cluster_counts.items():
        lines.append(f"- `{cluster}`: {n_cells}")
    lines.extend(["", "## Sample Counts", ""])
    for sample, n_cells in sample_counts.items():
        lines.append(f"- `{sample}`: {n_cells}")
    lines.extend(
        [
            "",
            "## Output Format",
            "",
            "The object is saved as a compressed SciPy sparse matrix in cells-by-genes orientation, with matching cell, barcode, gene, and gene-metadata tables.",
            "",
            f"- Matrix: `{matrix_path.relative_to(ROOT)}`",
            f"- Cell metadata: `{cell_meta_path.relative_to(ROOT)}`",
            f"- Gene metadata: `{gene_meta_path.relative_to(ROOT)}`",
            f"- Summary TSV: `{SUMMARY_TSV.relative_to(ROOT)}`",
            "",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines))

    print(f"Wrote {matrix_path}")
    print(f"Wrote {cell_meta_path}")
    print(f"Wrote {gene_meta_path}")
    print(f"Wrote {SUMMARY_TSV}")
    print(f"Wrote {SUMMARY_MD}")


if __name__ == "__main__":
    main()
