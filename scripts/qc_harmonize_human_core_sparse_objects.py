#!/usr/bin/env python3
"""Create harmonized QC metadata for built human dentate/hippocampal sparse objects."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SEED_SUMMARY = RESULTS / "human_seed_sparse_object_summary.tsv"
GSE186538_SUMMARY = RESULTS / "gse186538_human_dg_gc_sparse_subset_summary.tsv"

OUT_CELLS = RESULTS / "human_core_harmonized_cell_qc.tsv.gz"
OUT_COMPONENT_SUMMARY = RESULTS / "human_core_component_qc_summary.tsv"
OUT_DATASET_SUMMARY = RESULTS / "human_core_dataset_qc_summary.tsv"
OUT_MD = RESULTS / "human_core_qc_harmonization_summary.md"


def read_gene_names(path: Path) -> pd.Series:
    genes = pd.read_csv(path, sep="\t")
    if "gene" not in genes.columns:
        if genes.shape[1] == 1:
            genes.columns = ["gene"]
        else:
            raise ValueError(f"No gene column found in {path}")
    return genes["gene"].astype(str)


def compute_percent_mt(matrix_path: Path, genes_path: Path, n_counts: np.ndarray) -> tuple[np.ndarray, int]:
    matrix = sparse.load_npz(matrix_path).tocsr()
    genes = read_gene_names(genes_path)
    if matrix.shape[1] != len(genes):
        raise ValueError(f"Gene length mismatch for {matrix_path}: matrix {matrix.shape[1]}, genes {len(genes)}")

    mito_mask = genes.str.upper().str.startswith("MT-").to_numpy()
    n_mito_genes = int(mito_mask.sum())
    if n_mito_genes:
        mito_counts = np.asarray(matrix[:, mito_mask].sum(axis=1)).ravel().astype(float)
    else:
        mito_counts = np.zeros(matrix.shape[0], dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        percent_mt = np.where(n_counts > 0, mito_counts / n_counts * 100.0, 0.0)
    return percent_mt, n_mito_genes


def sample_hint_from_library(library_id: str) -> str:
    match = re.search(r"(sample\d+|Sample\d+)", library_id)
    if match:
        return match.group(1)
    match = re.search(r"(GSM\d+)", library_id)
    if match:
        return match.group(1)
    return library_id


def age_hint_from_library(library_id: str) -> str:
    match = re.search(r"(\d+yrs|\d+yr|\d+mos|\d+mo|\d+days|\d+d)", library_id, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def preliminary_pass(n_counts: pd.Series, n_genes: pd.Series, percent_mt: pd.Series) -> pd.Series:
    return (n_counts >= 500) & (n_genes >= 300) & (percent_mt <= 20)


def harmonize_seed_row(row: pd.Series) -> pd.DataFrame:
    matrix_path = ROOT / row["matrix_path"]
    genes_path = ROOT / row["genes_path"]
    cell_meta_path = ROOT / row["cell_metadata_path"]
    cell_meta = pd.read_csv(cell_meta_path, sep="\t")
    n_counts = cell_meta["n_counts"].to_numpy(dtype=np.int64)
    percent_mt, n_mito_genes = compute_percent_mt(matrix_path, genes_path, n_counts)

    out = pd.DataFrame(
        {
            "cell_id": cell_meta["cell_id"],
            "dataset": row["dataset"],
            "component_id": row["library_id"],
            "component_type": "library",
            "source_member": row["source_member"],
            "source_cell_name": "",
            "barcode": cell_meta["barcode"],
            "sample_hint": sample_hint_from_library(row["library_id"]),
            "age_hint": age_hint_from_library(row["library_id"]),
            "region": "hippocampus_dentate_reference",
            "cluster": "",
            "n_counts": cell_meta["n_counts"].astype(np.int64),
            "n_genes": cell_meta["n_genes"].astype(np.int64),
            "percent_mt": percent_mt,
            "source_percent_mt": np.nan,
            "n_mito_genes_detected_in_feature_space": n_mito_genes,
            "matrix_path": row["matrix_path"],
            "cell_metadata_path": row["cell_metadata_path"],
        }
    )
    out["preliminary_qc_pass"] = preliminary_pass(out["n_counts"], out["n_genes"], out["percent_mt"])
    return out


def harmonize_gse186538_row(row: pd.Series) -> pd.DataFrame:
    matrix_path = ROOT / row["matrix_path"]
    cell_meta_path = ROOT / row["cell_metadata_path"]
    gene_meta_path = ROOT / row["gene_metadata_path"]
    cell_meta = pd.read_csv(cell_meta_path, sep="\t")

    n_counts = cell_meta["extracted_n_counts"].to_numpy(dtype=np.int64)
    percent_mt, n_mito_genes = compute_percent_mt(matrix_path, gene_meta_path, n_counts)

    out = pd.DataFrame(
        {
            "cell_id": cell_meta["cell_id"],
            "dataset": row["dataset"],
            "component_id": row["subset"],
            "component_type": "curated_subset",
            "source_member": row["source_matrix"],
            "source_cell_name": cell_meta["cell_name"],
            "barcode": cell_meta["cell_name"].astype(str).str.rsplit("_", n=1).str[-1],
            "sample_hint": cell_meta["samplename"],
            "age_hint": "",
            "region": cell_meta["region"],
            "cluster": cell_meta["cluster"],
            "n_counts": cell_meta["extracted_n_counts"].astype(np.int64),
            "n_genes": cell_meta["extracted_n_genes"].astype(np.int64),
            "percent_mt": percent_mt,
            "source_percent_mt": cell_meta["percent.mt"].astype(float),
            "n_mito_genes_detected_in_feature_space": n_mito_genes,
            "matrix_path": row["matrix_path"],
            "cell_metadata_path": row["cell_metadata_path"],
        }
    )
    out["preliminary_qc_pass"] = preliminary_pass(out["n_counts"], out["n_genes"], out["percent_mt"])
    return out


def summarize_group(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    grouped = df.groupby(group_cols, dropna=False)
    rows = []
    for keys, group in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_cells": len(group),
                "preliminary_qc_pass_cells": int(group["preliminary_qc_pass"].sum()),
                "median_counts": float(group["n_counts"].median()),
                "median_genes": float(group["n_genes"].median()),
                "median_percent_mt": float(group["percent_mt"].median()),
                "p05_counts": float(group["n_counts"].quantile(0.05)),
                "p95_counts": float(group["n_counts"].quantile(0.95)),
                "p05_genes": float(group["n_genes"].quantile(0.05)),
                "p95_genes": float(group["n_genes"].quantile(0.95)),
                "p95_percent_mt": float(group["percent_mt"].quantile(0.95)),
                "high_mito_gt20pct_cells": int((group["percent_mt"] > 20).sum()),
                "low_gene_lt300_cells": int((group["n_genes"] < 300).sum()),
                "low_count_lt500_cells": int((group["n_counts"] < 500).sum()),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    seed_summary = pd.read_csv(SEED_SUMMARY, sep="\t")
    gse186538_summary = pd.read_csv(GSE186538_SUMMARY, sep="\t")

    harmonized_parts: list[pd.DataFrame] = []
    for _, row in seed_summary.iterrows():
        print(f"QC seed {row['dataset']} {row['library_id']}", flush=True)
        harmonized_parts.append(harmonize_seed_row(row))

    for _, row in gse186538_summary.iterrows():
        print(f"QC {row['dataset']} {row['subset']}", flush=True)
        harmonized_parts.append(harmonize_gse186538_row(row))

    cells = pd.concat(harmonized_parts, ignore_index=True)
    cells.to_csv(OUT_CELLS, sep="\t", index=False, compression="gzip")

    component_summary = summarize_group(cells, ["dataset", "component_id", "component_type"])
    dataset_summary = summarize_group(cells, ["dataset"])
    component_summary.to_csv(OUT_COMPONENT_SUMMARY, sep="\t", index=False)
    dataset_summary.to_csv(OUT_DATASET_SUMMARY, sep="\t", index=False)

    lines = [
        "# Human core QC harmonization",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This QC pass covers the currently built human dentate/hippocampal core objects: `GSE185277`, `GSE185553`, and the `GSE186538` DG GC sparse subset.",
        "",
        "## Dataset Summary",
        "",
    ]
    for _, row in dataset_summary.sort_values("dataset").iterrows():
        pass_rate = row["preliminary_qc_pass_cells"] / row["n_cells"] * 100 if row["n_cells"] else 0
        lines.append(
            f"- `{row['dataset']}`: {int(row['n_cells'])} cells, "
            f"median {row['median_counts']:.0f} counts, median {row['median_genes']:.0f} genes, "
            f"median {row['median_percent_mt']:.2f}% MT, preliminary pass {pass_rate:.1f}%."
        )
    lines.extend(
        [
            "",
            "## Preliminary QC Rule",
            "",
            "A cell is marked as preliminary pass if `n_counts >= 500`, `n_genes >= 300`, and `percent_mt <= 20`. This is a diagnostic flag, not a final filtering rule.",
            "",
            "## Outputs",
            "",
            f"- Harmonized cell QC table: `{OUT_CELLS.relative_to(ROOT)}`",
            f"- Component summary: `{OUT_COMPONENT_SUMMARY.relative_to(ROOT)}`",
            f"- Dataset summary: `{OUT_DATASET_SUMMARY.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote {OUT_CELLS}")
    print(f"Wrote {OUT_COMPONENT_SUMMARY}")
    print(f"Wrote {OUT_DATASET_SUMMARY}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
