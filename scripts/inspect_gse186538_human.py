#!/usr/bin/env python3
"""Inspect downloaded GSE186538 human files."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE186538"
RESULTS = ROOT / "Project/results"

META = BASE / "GSE186538_Human_cell_meta.txt.gz"
GENES = BASE / "GSE186538_Human_genes.txt.gz"
MTX = BASE / "GSE186538_Human_counts.mtx.gz"

SUMMARY_TSV = RESULTS / "gse186538_human_file_summary.tsv"
META_PROFILE_TSV = RESULTS / "gse186538_human_metadata_profile.tsv"
SUMMARY_MD = RESULTS / "gse186538_human_file_inspection.md"


def read_mtx_dims(path: Path) -> tuple[str, tuple[int, int, int]]:
    header = ""
    dims = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("%%MatrixMarket"):
                header = line
                continue
            if line.startswith("%"):
                continue
            parts = line.split()
            dims = (int(parts[0]), int(parts[1]), int(parts[2]))
            break
    if dims is None:
        raise ValueError(f"No Matrix Market dimensions found in {path}")
    return header, dims


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(META, sep="\t")
    genes = pd.read_csv(GENES, sep="\t", header=None)
    mtx_header, (mtx_rows, mtx_cols, mtx_nnz) = read_mtx_dims(MTX)

    if mtx_rows == len(genes) and mtx_cols == len(meta):
        orientation = "genes_by_cells"
    elif mtx_cols == len(genes) and mtx_rows == len(meta):
        orientation = "cells_by_genes"
    else:
        orientation = "unmatched_dimensions"

    summary_rows = [
        {
            "dataset": "GSE186538",
            "species": "Homo sapiens",
            "file_role": "cell_metadata",
            "path": str(META.relative_to(ROOT)),
            "rows": len(meta),
            "columns": len(meta.columns),
            "matrix_market_header": "",
            "matrix_rows": "",
            "matrix_columns": "",
            "matrix_nnz": "",
            "orientation_inference": "",
        },
        {
            "dataset": "GSE186538",
            "species": "Homo sapiens",
            "file_role": "genes",
            "path": str(GENES.relative_to(ROOT)),
            "rows": len(genes),
            "columns": len(genes.columns),
            "matrix_market_header": "",
            "matrix_rows": "",
            "matrix_columns": "",
            "matrix_nnz": "",
            "orientation_inference": "",
        },
        {
            "dataset": "GSE186538",
            "species": "Homo sapiens",
            "file_role": "counts_matrix",
            "path": str(MTX.relative_to(ROOT)),
            "rows": "",
            "columns": "",
            "matrix_market_header": mtx_header,
            "matrix_rows": mtx_rows,
            "matrix_columns": mtx_cols,
            "matrix_nnz": mtx_nnz,
            "orientation_inference": orientation,
        },
    ]
    with SUMMARY_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    profile_rows: list[dict[str, object]] = []
    for col in meta.columns:
        series = meta[col]
        values = series.astype("string")
        vc = values.value_counts(dropna=False).head(12)
        profile_rows.append(
            {
                "column": col,
                "n_unique": int(values.nunique(dropna=True)),
                "n_missing": int(values.isna().sum()),
                "top_values": "; ".join(f"{idx}={int(count)}" for idx, count in vc.items()),
            }
        )
    pd.DataFrame(profile_rows).to_csv(META_PROFILE_TSV, sep="\t", index=False)

    likely_columns = [
        col
        for col in meta.columns
        if any(key in col.lower() for key in ["region", "sub", "cluster", "type", "cell", "sample", "donor", "age"])
    ]

    lines = [
        "# GSE186538 human file inspection",
        "",
        "Date inspected: 2026-06-21",
        "",
        "## Files",
        "",
        f"- Cell metadata: {len(meta)} rows x {len(meta.columns)} columns.",
        f"- Genes: {len(genes)} rows x {len(genes.columns)} columns.",
        f"- Count matrix: {mtx_rows} rows x {mtx_cols} columns with {mtx_nnz} non-zero entries.",
        f"- Orientation inference: `{orientation}`.",
        "",
        "## Metadata Columns",
        "",
        ", ".join(meta.columns),
        "",
        "## Likely Annotation Columns",
        "",
        ", ".join(likely_columns) if likely_columns else "No obvious annotation columns detected by name.",
        "",
        "## Outputs",
        "",
        f"- Summary TSV: `{SUMMARY_TSV.relative_to(ROOT)}`",
        f"- Metadata profile TSV: `{META_PROFILE_TSV.relative_to(ROOT)}`",
        "",
    ]
    SUMMARY_MD.write_text("\n".join(lines))
    print(f"Wrote {SUMMARY_TSV}")
    print(f"Wrote {META_PROFILE_TSV}")
    print(f"Wrote {SUMMARY_MD}")
    print(f"orientation={orientation}; meta_rows={len(meta)}; genes={len(genes)}; nnz={mtx_nnz}")


if __name__ == "__main__":
    main()
