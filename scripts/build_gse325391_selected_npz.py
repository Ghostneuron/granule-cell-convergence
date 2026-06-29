#!/usr/bin/env python3
"""Convert the GSE325391 selected-feature CSC bridge into a SciPy NPZ object."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BRIDGE = ROOT / "Project/processed/gse325391_adult_dg_selected"
RESULTS = ROOT / "Project/results"

SHAPE = BRIDGE / "selected_present_genes_csc_shape.tsv"
FEATURES = BRIDGE / "selected_features.tsv"
CELLS = BRIDGE / "cell_metadata.tsv.gz"
I_BIN = BRIDGE / "selected_present_genes_csc_i_int32.bin"
P_BIN = BRIDGE / "selected_present_genes_csc_p_int32.bin"
X_BIN = BRIDGE / "selected_present_genes_csc_x_float64.bin"

OUT_MATRIX = BRIDGE / "matrix_cells_by_selected_genes.npz"
OUT_VAR = BRIDGE / "var_selected_features.tsv"
OUT_SUMMARY = RESULTS / "gse325391_selected_npz_summary.tsv"
OUT_MD = RESULTS / "gse325391_selected_npz_summary.md"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    shape = pd.read_csv(SHAPE, sep="\t").iloc[0]
    features = pd.read_csv(FEATURES, sep="\t")
    cells = pd.read_csv(CELLS, sep="\t", usecols=["cell_id"])

    n_selected = int(shape["n_selected_genes"])
    n_present = int(shape["n_present_selected_genes"])
    n_cells = int(shape["n_cells"])
    nnz = int(shape["nnz_present_selected"])

    i = np.fromfile(I_BIN, dtype="<i4")
    p = np.fromfile(P_BIN, dtype="<i4")
    x = np.fromfile(X_BIN, dtype="<f8")

    if len(i) != nnz or len(x) != nnz:
        raise ValueError(f"NNZ mismatch: shape says {nnz}, i={len(i)}, x={len(x)}")
    if len(p) != n_cells + 1:
        raise ValueError(f"Column pointer length mismatch: expected {n_cells + 1}, got {len(p)}")
    if len(cells) != n_cells:
        raise ValueError(f"Cell metadata mismatch: expected {n_cells}, got {len(cells)}")

    present = features.loc[features["present_in_gse325391"].astype(bool)].copy()
    present = present.sort_values("present_matrix_row_index")
    if len(present) != n_present:
        raise ValueError(f"Present feature mismatch: expected {n_present}, got {len(present)}")

    present_to_selected = present["selected_feature_index"].to_numpy(dtype=np.int32)
    if not np.array_equal(present["present_matrix_row_index"].to_numpy(dtype=np.int32), np.arange(n_present, dtype=np.int32)):
        raise ValueError("Present matrix row indices are not contiguous")

    if len(x):
        max_fractional_error = float(np.max(np.abs(x - np.rint(x))))
    else:
        max_fractional_error = 0.0
    if max_fractional_error > 1e-6:
        data = x.astype(np.float32)
        data_dtype = "float32"
    else:
        data = np.rint(x).astype(np.int32)
        data_dtype = "int32"

    mapped_i = present_to_selected[i]
    genes_by_cells = sparse.csc_matrix((data, mapped_i, p), shape=(n_selected, n_cells))
    cells_by_genes = genes_by_cells.T.tocsr()
    cells_by_genes.sum_duplicates()
    sparse.save_npz(OUT_MATRIX, cells_by_genes, compressed=True)

    var = features.copy()
    keep_cols = [
        "selected_feature_index",
        "feature_index",
        "gene",
        "selection_reason",
        "max_feature_score",
        "GSE185277",
        "GSE185553",
        "GSE186538",
        "is_marker_panel_gene",
        "present_in_gse325391",
        "source_feature",
        "source_row_in_rds_1based",
        "present_matrix_row_index",
    ]
    keep_cols = [col for col in keep_cols if col in var.columns]
    var[keep_cols].to_csv(OUT_VAR, sep="\t", index=False)

    total_counts = int(cells_by_genes.sum())
    row_counts = np.asarray(cells_by_genes.sum(axis=1)).ravel()
    row_genes = cells_by_genes.getnnz(axis=1)
    summary = pd.DataFrame(
        [
            {
                "dataset": "GSE325391",
                "matrix_path": str(OUT_MATRIX.relative_to(ROOT)),
                "cell_metadata_path": str(CELLS.relative_to(ROOT)),
                "var_path": str(OUT_VAR.relative_to(ROOT)),
                "n_cells": cells_by_genes.shape[0],
                "n_selected_genes": cells_by_genes.shape[1],
                "n_present_selected_genes": n_present,
                "n_missing_selected_genes": n_selected - n_present,
                "nnz": int(cells_by_genes.nnz),
                "total_selected_counts": total_counts,
                "median_selected_counts_per_cell": float(np.median(row_counts)),
                "median_selected_genes_per_cell": float(np.median(row_genes)),
                "data_dtype": data_dtype,
                "max_fractional_error_before_cast": max_fractional_error,
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    lines = [
        "# GSE325391 Selected NPZ Summary",
        "",
        "Date built: 2026-06-21",
        "",
        "## Matrix",
        "",
        f"- Shape: {cells_by_genes.shape[0]} cells/nuclei by {cells_by_genes.shape[1]} selected genes.",
        f"- Present selected genes: {n_present} / {n_selected}.",
        f"- Non-zero selected count entries: {cells_by_genes.nnz}.",
        f"- Total selected counts: {total_counts}.",
        f"- Data dtype: `{data_dtype}`.",
        "",
        "## Outputs",
        "",
        f"- Matrix: `{OUT_MATRIX.relative_to(ROOT)}`",
        f"- Cell metadata: `{CELLS.relative_to(ROOT)}`",
        f"- Feature metadata: `{OUT_VAR.relative_to(ROOT)}`",
        f"- Summary TSV: `{OUT_SUMMARY.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))

    print(f"Wrote {OUT_MATRIX}")
    print(f"Wrote {OUT_VAR}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_MD}")
    print(f"shape={cells_by_genes.shape}; nnz={cells_by_genes.nnz}; dtype={data_dtype}")


if __name__ == "__main__":
    main()
