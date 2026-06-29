#!/usr/bin/env python3
"""Build a selected human-core gene matrix for GSE268609 by streaming the MEX file."""

from __future__ import annotations

import argparse
import gzip
import re
import time
from array import array
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE268609"
PROCESSED = ROOT / "Project/processed/gse268609_rna_selected"
RESULTS = ROOT / "Project/results"

FEATURES = BASE / "GSE268609_features.tsv.gz"
BARCODES = BASE / "GSE268609_barcodes.tsv.gz"
MATRIX = BASE / "GSE268609_matrix.mtx.gz"
SELECTED_PRESENCE = RESULTS / "gse268609_selected_gene_presence.tsv"
SAMPLE_METADATA = RESULTS / "gse268609_geo_sample_metadata.tsv"

OUT_MATRIX = PROCESSED / "matrix_cells_by_selected_genes.npz"
OUT_CELL_METADATA = PROCESSED / "cell_metadata.tsv.gz"
OUT_VAR = PROCESSED / "var_selected_features.tsv"
OUT_SUMMARY = RESULTS / "gse268609_selected_npz_summary.tsv"
OUT_MD = RESULTS / "gse268609_selected_npz_summary.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-entries",
        type=int,
        default=0,
        help="Stop after this many matrix entries for a smoke test. The default streams the full matrix.",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Parse and report progress without writing NPZ/metadata outputs. Useful with --max-entries.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25_000_000,
        help="Report progress every N matrix entries.",
    )
    return parser.parse_args()


def read_features() -> tuple[bytearray, int, int]:
    feature_types: list[str] = []
    with gzip.open(FEATURES, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            feature_types.append(parts[2] if len(parts) > 2 else "")
    is_gene = bytearray(len(feature_types) + 1)
    n_gene = 0
    for idx, feature_type in enumerate(feature_types, start=1):
        if feature_type == "Gene Expression":
            is_gene[idx] = 1
            n_gene += 1
    return is_gene, len(feature_types), n_gene


def read_barcodes() -> pd.DataFrame:
    rows = []
    with gzip.open(BARCODES, "rt", errors="replace") as fh:
        for idx, line in enumerate(fh, start=1):
            cell_id = line.strip()
            match = re.search(r"-(\d+)$", cell_id)
            sample_id = match.group(1) if match else ""
            rows.append({"matrix_col_1based": idx, "cell_id": cell_id, "sample_id": sample_id})
    cells = pd.DataFrame(rows)
    if SAMPLE_METADATA.exists():
        meta = pd.read_csv(SAMPLE_METADATA, sep="\t", dtype={"sample_id": str})
        rna = meta.loc[meta["library_type"].astype(str).str.upper().eq("RNA")].copy()
        keep = [
            "sample_id",
            "sample_accession",
            "diagnosis",
            "age_at_death_years",
            "pmi_hours",
            "tissue",
            "source_name",
            "instrument_model",
        ]
        keep = [col for col in keep if col in rna.columns]
        cells = cells.merge(rna[keep], on="sample_id", how="left")
    return cells


def read_selected_map() -> tuple[dict[int, int], pd.DataFrame]:
    var = pd.read_csv(SELECTED_PRESENCE, sep="\t")
    present = var.loc[var["present_in_gse268609"].astype(bool)].copy()
    present["source_matrix_row_1based"] = pd.to_numeric(present["source_matrix_row_1based"], errors="coerce").astype(int)
    present["selected_feature_index"] = pd.to_numeric(present["selected_feature_index"], errors="coerce").astype(int)
    row_to_selected = dict(zip(present["source_matrix_row_1based"], present["selected_feature_index"]))
    return row_to_selected, var


def read_matrix_header() -> tuple[str, int, int, int]:
    market_header = ""
    with gzip.open(MATRIX, "rb") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b"%%MatrixMarket"):
                market_header = line.decode("ascii", errors="replace")
                continue
            if line.startswith(b"%"):
                continue
            parts = line.split()
            return market_header, int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"No Matrix Market dimension line found in {MATRIX}")


def stream_matrix(
    is_gene: bytearray,
    row_to_selected: dict[int, int],
    n_cells: int,
    max_entries: int,
    progress_every: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int, bool]:
    n_count_rna = np.zeros(n_cells, dtype=np.int64)
    n_feature_rna = np.zeros(n_cells, dtype=np.int32)
    n_count_selected = np.zeros(n_cells, dtype=np.int64)
    n_feature_selected = np.zeros(n_cells, dtype=np.int32)
    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    indices = array("i")
    data = array("i")

    current_col = 1
    entries_seen = 0
    gene_entries_seen = 0
    selected_entries_seen = 0
    column_sorted = True
    t0 = time.time()
    header_consumed = False

    with gzip.open(MATRIX, "rb") as fh:
        for raw in fh:
            if raw.startswith(b"%"):
                continue
            if not header_consumed:
                header_consumed = True
                continue

            parts = raw.split()
            if len(parts) < 3:
                continue
            row = int(parts[0])
            col = int(parts[1])
            value = int(parts[2])
            entries_seen += 1

            if col < current_col:
                column_sorted = False
            while current_col < col and current_col <= n_cells:
                indptr[current_col] = len(indices)
                current_col += 1

            cell_idx = col - 1
            if row < len(is_gene) and is_gene[row]:
                n_count_rna[cell_idx] += value
                n_feature_rna[cell_idx] += 1
                gene_entries_seen += 1
                selected_idx = row_to_selected.get(row)
                if selected_idx is not None:
                    indices.append(selected_idx)
                    data.append(value)
                    n_count_selected[cell_idx] += value
                    n_feature_selected[cell_idx] += 1
                    selected_entries_seen += 1

            if progress_every and entries_seen % progress_every == 0:
                elapsed = max(time.time() - t0, 1e-9)
                print(
                    "streamed {:,} entries; col {:,}/{:,}; selected nnz {:,}; {:.1f}M entries/min".format(
                        entries_seen,
                        col,
                        n_cells,
                        selected_entries_seen,
                        entries_seen / elapsed / 1_000_000 * 60,
                    ),
                    flush=True,
                )
            if max_entries and entries_seen >= max_entries:
                break

    while current_col <= n_cells:
        indptr[current_col] = len(indices)
        current_col += 1

    np_indices = np.frombuffer(indices, dtype=np.int32)
    np_data = np.frombuffer(data, dtype=np.int32)
    return (
        indptr,
        np_indices,
        np_data,
        n_count_rna,
        n_feature_rna,
        n_count_selected,
        n_feature_selected,
        entries_seen,
        gene_entries_seen,
        column_sorted,
    )


def write_outputs(
    X: sparse.csr_matrix,
    cells: pd.DataFrame,
    var: pd.DataFrame,
    n_count_rna: np.ndarray,
    n_feature_rna: np.ndarray,
    n_count_selected: np.ndarray,
    n_feature_selected: np.ndarray,
    matrix_header: tuple[str, int, int, int],
    entries_seen: int,
    gene_entries_seen: int,
    column_sorted: bool,
) -> None:
    market_header, n_rows, n_cols, n_total_nnz = matrix_header
    sparse.save_npz(OUT_MATRIX, X, compressed=True)

    obs = cells.copy()
    obs["nCount_RNA"] = n_count_rna
    obs["nFeature_RNA"] = n_feature_rna
    obs["nCount_selected"] = n_count_selected
    obs["nFeature_selected"] = n_feature_selected
    obs["selected_count_fraction"] = np.divide(
        n_count_selected,
        n_count_rna,
        out=np.zeros(len(obs), dtype=np.float64),
        where=n_count_rna > 0,
    )
    obs["analysis_include_basic_qc"] = (obs["nCount_RNA"] >= 2000) & (obs["nFeature_RNA"] >= 1000)
    obs.to_csv(OUT_CELL_METADATA, sep="\t", index=False, compression="gzip")

    var_out = var.copy()
    var_out.to_csv(OUT_VAR, sep="\t", index=False)

    row_counts = np.asarray(X.sum(axis=1)).ravel()
    row_genes = X.getnnz(axis=1)
    n_basic_qc = int(obs["analysis_include_basic_qc"].sum())
    summary = pd.DataFrame(
        [
            {
                "dataset": "GSE268609",
                "matrix_path": str(OUT_MATRIX.relative_to(ROOT)),
                "cell_metadata_path": str(OUT_CELL_METADATA.relative_to(ROOT)),
                "var_path": str(OUT_VAR.relative_to(ROOT)),
                "matrix_market_header": market_header,
                "source_matrix_rows": n_rows,
                "source_matrix_cols": n_cols,
                "source_matrix_nnz": n_total_nnz,
                "entries_streamed": entries_seen,
                "gene_expression_entries_streamed": gene_entries_seen,
                "column_sorted": column_sorted,
                "n_cells": X.shape[0],
                "n_selected_genes": X.shape[1],
                "n_present_selected_genes": int(var_out["present_in_gse268609"].astype(bool).sum()),
                "n_missing_selected_genes": int((~var_out["present_in_gse268609"].astype(bool)).sum()),
                "nnz": int(X.nnz),
                "total_selected_counts": int(X.sum()),
                "median_selected_counts_per_cell": float(np.median(row_counts)),
                "median_selected_genes_per_cell": float(np.median(row_genes)),
                "median_total_rna_counts_per_cell": float(np.median(n_count_rna)),
                "median_total_rna_genes_per_cell": float(np.median(n_feature_rna)),
                "n_cells_basic_qc": n_basic_qc,
                "basic_qc_fraction": float(n_basic_qc / X.shape[0]) if X.shape[0] else 0.0,
                "data_dtype": str(X.data.dtype),
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    row = summary.iloc[0]
    lines = [
        "# GSE268609 Selected NPZ Summary",
        "",
        "Date built: 2026-06-21",
        "",
        "## Matrix",
        "",
        f"- Shape: {int(row['n_cells'])} cells/nuclei by {int(row['n_selected_genes'])} selected genes.",
        f"- Present selected genes: {int(row['n_present_selected_genes'])} / {int(row['n_selected_genes'])}.",
        f"- Non-zero selected count entries: {int(row['nnz'])}.",
        f"- Total selected counts: {int(row['total_selected_counts'])}.",
        f"- Cells passing basic RNA QC: {int(row['n_cells_basic_qc'])} ({row['basic_qc_fraction']:.1%}).",
        "",
        "## Source Matrix",
        "",
        f"- Header: `{market_header}`",
        f"- Source shape: {n_rows} rows by {n_cols} columns with {n_total_nnz} non-zero entries.",
        f"- Entries streamed: {entries_seen}; gene-expression entries among streamed rows: {gene_entries_seen}.",
        f"- Matrix entries were column sorted: {column_sorted}.",
        "",
        "## Outputs",
        "",
        f"- Matrix: `{OUT_MATRIX.relative_to(ROOT)}`",
        f"- Cell metadata: `{OUT_CELL_METADATA.relative_to(ROOT)}`",
        f"- Feature metadata: `{OUT_VAR.relative_to(ROOT)}`",
        f"- Summary TSV: `{OUT_SUMMARY.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    matrix_header = read_matrix_header()
    _, n_rows, n_cols, n_total_nnz = matrix_header
    is_gene, n_features, n_gene = read_features()
    cells = read_barcodes()
    row_to_selected, var = read_selected_map()

    if n_rows != n_features:
        raise ValueError(f"Matrix row count {n_rows} does not match feature count {n_features}")
    if n_cols != len(cells):
        raise ValueError(f"Matrix column count {n_cols} does not match barcode count {len(cells)}")

    print(
        "Starting GSE268609 stream: {:,} entries, {:,} cells, {:,} gene rows, {:,} selected present genes".format(
            n_total_nnz, n_cols, n_gene, len(row_to_selected)
        ),
        flush=True,
    )
    (
        indptr,
        indices,
        data,
        n_count_rna,
        n_feature_rna,
        n_count_selected,
        n_feature_selected,
        entries_seen,
        gene_entries_seen,
        column_sorted,
    ) = stream_matrix(
        is_gene=is_gene,
        row_to_selected=row_to_selected,
        n_cells=n_cols,
        max_entries=args.max_entries,
        progress_every=args.progress_every,
    )

    n_selected = int(var["selected_feature_index"].max()) + 1
    print(
        "Finished stream: entries={:,}; gene_entries={:,}; selected_nnz={:,}; column_sorted={}".format(
            entries_seen, gene_entries_seen, len(data), column_sorted
        ),
        flush=True,
    )

    if args.no_save:
        print("No-save mode; not writing matrix outputs.", flush=True)
        return

    X = sparse.csr_matrix((data, indices, indptr), shape=(n_cols, n_selected), dtype=np.int32)
    X.sum_duplicates()
    write_outputs(
        X=X,
        cells=cells,
        var=var,
        n_count_rna=n_count_rna,
        n_feature_rna=n_feature_rna,
        n_count_selected=n_count_selected,
        n_feature_selected=n_feature_selected,
        matrix_header=matrix_header,
        entries_seen=entries_seen,
        gene_entries_seen=gene_entries_seen,
        column_sorted=column_sorted,
    )

    print(f"Wrote {OUT_MATRIX}")
    print(f"Wrote {OUT_CELL_METADATA}")
    print(f"Wrote {OUT_VAR}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_MD}")
    print(f"shape={X.shape}; nnz={X.nnz}; total_selected_counts={int(X.sum())}")


if __name__ == "__main__":
    main()
