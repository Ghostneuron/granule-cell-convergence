#!/usr/bin/env python3
"""Extract selected GSE268609 gene and ATAC peak rows if the full matrix exists.

The selected feature rows are defined by
`gse268609_epigenomic_selective_extraction_manifest.tsv`. When the full
`GSE268609_matrix.mtx.gz` file is not local, this script writes a plan-only
summary and exits cleanly.
"""

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
RESULTS = ROOT / "Project/results"
PROCESSED = ROOT / "Project/processed/gse268609_epigenomic_selected"

MATRIX = BASE / "GSE268609_matrix.mtx.gz"
BARCODES = BASE / "GSE268609_barcodes.tsv.gz"
MANIFEST = RESULTS / "gse268609_epigenomic_selective_extraction_manifest.tsv"

OUT_MATRIX = PROCESSED / "matrix_cells_by_epigenomic_selected_features.npz"
OUT_CELL_METADATA = PROCESSED / "cell_metadata.tsv.gz"
OUT_VAR = PROCESSED / "var_epigenomic_selected_features.tsv"
OUT_SUMMARY = RESULTS / "gse268609_epigenomic_selected_matrix_summary.tsv"
OUT_MD = RESULTS / "gse268609_epigenomic_selected_matrix_extraction.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-entries", type=int, default=0, help="Stop early for a smoke test.")
    parser.add_argument("--progress-every", type=int, default=25_000_000)
    parser.add_argument("--no-save", action="store_true")
    return parser.parse_args()


def read_manifest() -> tuple[dict[int, int], pd.DataFrame]:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    rows = []
    for feature_row, sub in manifest.groupby("feature_row_1based", sort=True):
        first = sub.iloc[0]
        rows.append(
            {
                "selected_feature_index": len(rows),
                "source_matrix_row_1based": int(feature_row),
                "feature_id": first["feature_id"],
                "feature_name": first["feature_name"],
                "feature_type": first["feature_type"],
                "feature_chrom": first["feature_chrom"],
                "feature_start": first["feature_start"],
                "feature_end": first["feature_end"],
                "linked_genes": ";".join(sorted(sub["linked_gene"].astype(str).unique())),
                "best_priorities": ";".join(sorted(sub["best_priority"].astype(str).unique())),
                "model_terms_supported": ";".join(sorted(sub["model_terms_supported"].astype(str).unique())),
                "target_sets": ";".join(sorted(sub["target_sets"].astype(str).unique())),
                "peak_categories": ";".join(
                    sorted(x for x in sub["peak_category"].fillna("").astype(str).unique() if x)
                ),
                "selection_reasons": ";".join(sorted(sub["selection_reason"].astype(str).unique())),
            }
        )
    var = pd.DataFrame(rows)
    row_map = dict(zip(var["source_matrix_row_1based"], var["selected_feature_index"]))
    return row_map, var


def read_barcodes() -> pd.DataFrame:
    rows = []
    with gzip.open(BARCODES, "rt", errors="replace") as handle:
        for idx, line in enumerate(handle, start=1):
            cell_id = line.strip()
            match = re.search(r"-(\d+)$", cell_id)
            rows.append(
                {
                    "matrix_col_1based": idx,
                    "cell_id": cell_id,
                    "sample_id": match.group(1) if match else "",
                }
            )
    return pd.DataFrame(rows)


def read_matrix_header() -> tuple[str, int, int, int]:
    market_header = ""
    with gzip.open(MATRIX, "rb") as handle:
        for raw in handle:
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


def write_plan(var: pd.DataFrame, cells: pd.DataFrame, reason: str) -> None:
    feature_type_counts = var["feature_type"].value_counts().to_dict()
    lines = [
        "# GSE268609 Epigenomic Selected Matrix Extraction",
        "",
        "Date built: 2026-06-26",
        "",
        "## Status",
        "",
        f"- Extraction status: {reason}.",
        f"- Full matrix path: `{MATRIX.relative_to(ROOT)}`.",
        f"- Full matrix present locally: {MATRIX.exists()}.",
        f"- Selected feature rows: {len(var)}.",
        f"- Barcode/cell columns available: {len(cells)}.",
        f"- Selected feature types: {feature_type_counts}.",
        "",
        "## Next Action",
        "",
        "Restore or re-download `GSE268609_matrix.mtx.gz`, then rerun this script to stream only the selected gene and ATAC peak rows.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")
    summary = pd.DataFrame(
        [
            {
                "status": reason,
                "matrix_path": str(MATRIX.relative_to(ROOT)),
                "matrix_present": MATRIX.exists(),
                "n_selected_features": len(var),
                "n_cells": len(cells),
                "n_gene_expression_features": int(var["feature_type"].eq("Gene Expression").sum()),
                "n_peak_features": int(var["feature_type"].eq("Peaks").sum()),
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)


def stream_selected(
    row_map: dict[int, int],
    n_cells: int,
    max_entries: int,
    progress_every: int,
) -> tuple[sparse.csr_matrix, int, int, bool]:
    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    indices = array("i")
    data = array("i")
    current_col = 1
    entries_seen = 0
    selected_entries = 0
    column_sorted = True
    header_consumed = False
    t0 = time.time()

    with gzip.open(MATRIX, "rb") as handle:
        for raw in handle:
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
            selected_idx = row_map.get(row)
            if selected_idx is not None:
                indices.append(selected_idx)
                data.append(value)
                selected_entries += 1
            if progress_every and entries_seen % progress_every == 0:
                elapsed = max(time.time() - t0, 1e-9)
                print(
                    "streamed {:,} entries; selected nnz {:,}; {:.1f}M entries/min".format(
                        entries_seen,
                        selected_entries,
                        entries_seen / elapsed / 1_000_000 * 60,
                    ),
                    flush=True,
                )
            if max_entries and entries_seen >= max_entries:
                break

    while current_col <= n_cells:
        indptr[current_col] = len(indices)
        current_col += 1

    X = sparse.csr_matrix(
        (
            np.frombuffer(data, dtype=np.int32),
            np.frombuffer(indices, dtype=np.int32),
            indptr,
        ),
        shape=(n_cells, len(row_map)),
    )
    return X, entries_seen, selected_entries, column_sorted


def main() -> None:
    args = parse_args()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    row_map, var = read_manifest()
    cells = read_barcodes()

    if not MATRIX.exists():
        write_plan(var, cells, "plan_only_missing_full_matrix")
        print(f"Full matrix is not local; wrote {OUT_MD.relative_to(ROOT)}")
        return

    market_header, n_rows, n_cols, n_entries = read_matrix_header()
    X, entries_seen, selected_entries, column_sorted = stream_selected(
        row_map=row_map,
        n_cells=n_cols,
        max_entries=args.max_entries,
        progress_every=args.progress_every,
    )

    if not args.no_save:
        sparse.save_npz(OUT_MATRIX, X, compressed=True)
        cells.to_csv(OUT_CELL_METADATA, sep="\t", index=False, compression="gzip")
        var.to_csv(OUT_VAR, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "status": "completed" if not args.max_entries else "smoke_test",
                "matrix_header": market_header,
                "source_n_rows": n_rows,
                "source_n_cols": n_cols,
                "source_n_entries": n_entries,
                "entries_streamed": entries_seen,
                "selected_nnz": selected_entries,
                "column_sorted": column_sorted,
                "n_selected_features": len(var),
                "n_cells": len(cells),
                "output_matrix": str(OUT_MATRIX.relative_to(ROOT)) if not args.no_save else "",
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    lines = [
        "# GSE268609 Epigenomic Selected Matrix Extraction",
        "",
        "Date built: 2026-06-26",
        "",
        "## Status",
        "",
        f"- Status: {summary.loc[0, 'status']}.",
        f"- Selected features: {len(var)}.",
        f"- Cells/nuclei: {len(cells)}.",
        f"- Selected nonzero entries: {selected_entries}.",
        f"- Column sorted: {column_sorted}.",
        "",
        "## Outputs",
        "",
        f"- Matrix: `{OUT_MATRIX.relative_to(ROOT)}`",
        f"- Cell metadata: `{OUT_CELL_METADATA.relative_to(ROOT)}`",
        f"- Feature metadata: `{OUT_VAR.relative_to(ROOT)}`",
        f"- Summary: `{OUT_SUMMARY.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
