#!/usr/bin/env python3
"""Inspect GSE268609 RNA companion files before large matrix extraction."""

from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE268609"
PROCESSED = ROOT / "Project/processed/gse268609_rna_selected"
RESULTS = ROOT / "Project/results"

FEATURES = BASE / "GSE268609_features.tsv.gz"
BARCODES = BASE / "GSE268609_barcodes.tsv.gz"
MATRIX = BASE / "GSE268609_matrix.mtx.gz"
SELECTED_VAR = ROOT / "Project/processed/human_core_normalized_reduced_object/var.tsv"
SAMPLE_METADATA = RESULTS / "gse268609_geo_sample_metadata.tsv"

OUT_FEATURE_TYPES = RESULTS / "gse268609_feature_type_summary.tsv"
OUT_BARCODE_SUFFIX = RESULTS / "gse268609_barcode_sample_suffix_summary.tsv"
OUT_SELECTED_PRESENCE = RESULTS / "gse268609_selected_gene_presence.tsv"
OUT_MATRIX_HEADER = RESULTS / "gse268609_matrix_header.tsv"
OUT_MD = RESULTS / "gse268609_rna_file_inspection.md"


def norm_gene(value: object) -> str:
    return str(value).strip().upper()


def read_features() -> pd.DataFrame:
    rows = []
    with gzip.open(FEATURES, "rt", errors="replace") as fh:
        for idx, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            while len(parts) < 6:
                parts.append("")
            rows.append(
                {
                    "matrix_row_1based": idx + 1,
                    "feature_id": parts[0],
                    "feature_name": parts[1],
                    "feature_type": parts[2],
                    "chromosome": parts[3],
                    "start": parts[4],
                    "end": parts[5],
                }
            )
    return pd.DataFrame(rows)


def read_barcodes() -> pd.DataFrame:
    rows = []
    suffix_counter: Counter[str] = Counter()
    with gzip.open(BARCODES, "rt", errors="replace") as fh:
        for idx, line in enumerate(fh):
            barcode = line.strip()
            match = re.search(r"-(\d+)$", barcode)
            suffix = match.group(1) if match else ""
            suffix_counter[suffix] += 1
            rows.append({"matrix_col_1based": idx + 1, "cell_id": barcode, "sample_id": suffix})
    out = pd.DataFrame(rows)
    suffix = (
        pd.DataFrame(
            [{"sample_id": key, "n_barcodes": value} for key, value in suffix_counter.items()]
        )
        .assign(sample_id_int=lambda d: pd.to_numeric(d["sample_id"], errors="coerce"))
        .sort_values("sample_id_int")
        .drop(columns=["sample_id_int"])
    )
    return out, suffix


def read_matrix_header() -> dict[str, object]:
    header = {
        "matrix_path": str(MATRIX.relative_to(ROOT)),
        "file_exists": MATRIX.exists(),
        "actual_bytes": MATRIX.stat().st_size if MATRIX.exists() else 0,
        "matrix_market_header": "",
        "n_rows": None,
        "n_cols": None,
        "n_nonzero_entries": None,
        "header_read_error": "",
    }
    if not MATRIX.exists():
        return header
    try:
        with gzip.open(MATRIX, "rt", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("%%MatrixMarket"):
                    header["matrix_market_header"] = line
                    continue
                if line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    header["n_rows"] = int(parts[0])
                    header["n_cols"] = int(parts[1])
                    header["n_nonzero_entries"] = int(parts[2])
                break
    except Exception as exc:  # pragma: no cover - records partial gzip/header failures for project provenance.
        header["header_read_error"] = repr(exc)
    return header


def selected_gene_presence(features: pd.DataFrame) -> pd.DataFrame:
    selected = pd.read_csv(SELECTED_VAR, sep="\t")
    gene_features = features.loc[features["feature_type"].eq("Gene Expression")].copy()
    gene_features["gene_norm"] = gene_features["feature_name"].map(norm_gene)

    # If duplicate gene symbols exist, keep the first matrix row and record multiplicity.
    multiplicity = gene_features.groupby("gene_norm").size().rename("n_matching_gene_rows").reset_index()
    first_gene = (
        gene_features.sort_values("matrix_row_1based")
        .drop_duplicates("gene_norm")
        [["gene_norm", "feature_id", "feature_name", "matrix_row_1based"]]
    )
    out = selected.copy()
    out["selected_feature_index"] = out["feature_index"]
    out["gene_norm"] = out["gene"].map(norm_gene)
    out = out.merge(first_gene, on="gene_norm", how="left")
    out = out.merge(multiplicity, on="gene_norm", how="left")
    out["present_in_gse268609"] = out["matrix_row_1based"].notna()
    out["n_matching_gene_rows"] = out["n_matching_gene_rows"].fillna(0).astype(int)
    out = out.rename(
        columns={
            "feature_id": "gse268609_feature_id",
            "feature_name": "gse268609_feature_name",
            "matrix_row_1based": "source_matrix_row_1based",
        }
    )
    preferred = [
        "selected_feature_index",
        "feature_index",
        "gene",
        "selection_reason",
        "max_feature_score",
        "GSE185277",
        "GSE185553",
        "GSE186538",
        "is_marker_panel_gene",
        "present_in_gse268609",
        "gse268609_feature_id",
        "gse268609_feature_name",
        "source_matrix_row_1based",
        "n_matching_gene_rows",
    ]
    return out[[col for col in preferred if col in out.columns]]


def write_md(
    features: pd.DataFrame,
    suffix: pd.DataFrame,
    selected_presence: pd.DataFrame,
    matrix_header: dict[str, object],
) -> None:
    feature_type_summary = features["feature_type"].value_counts().rename_axis("feature_type").reset_index(name="n_features")
    n_gene = int(feature_type_summary.loc[feature_type_summary["feature_type"].eq("Gene Expression"), "n_features"].sum())
    n_peak = int(feature_type_summary.loc[feature_type_summary["feature_type"].eq("Peaks"), "n_features"].sum())
    n_present = int(selected_presence["present_in_gse268609"].sum())
    n_selected = len(selected_presence)
    rows_match = matrix_header.get("n_rows") == len(features)
    cols_match = matrix_header.get("n_cols") == int(suffix["n_barcodes"].sum()) if len(suffix) else False

    lines = [
        "# GSE268609 RNA File Inspection",
        "",
        "Date inspected: 2026-06-21",
        "",
        "## Companion Files",
        "",
        f"- Features: {len(features)} total ({n_gene} gene-expression rows, {n_peak} peak rows).",
        f"- Barcodes: {int(suffix['n_barcodes'].sum())} total across {len(suffix)} sample suffixes.",
        f"- Selected human-core genes present as gene-expression rows: {n_present} / {n_selected}.",
        "",
        "## Matrix Header",
        "",
        f"- Header: `{matrix_header.get('matrix_market_header', '')}`",
        f"- Dimensions: {matrix_header.get('n_rows')} rows by {matrix_header.get('n_cols')} columns.",
        f"- Non-zero entries: {matrix_header.get('n_nonzero_entries')}.",
        f"- Matrix rows match feature count: {rows_match}.",
        f"- Matrix columns match barcode count: {cols_match}.",
        f"- Header read error: `{matrix_header.get('header_read_error', '')}`",
        "",
        "## Extraction Decision",
        "",
        "- Treat this as a primary candidate only through its RNA gene-expression rows for transcriptomic comparison.",
        "- ATAC peak rows are valuable for later regulatory analysis, but they should not enter the first cross-dataset morphology-gene module projection.",
        "- Full selected-gene extraction should be attempted only after confirming the non-zero count and available memory/disk, because the combined multiome matrix is much larger than the previous adult dentate anchor.",
        "",
        "## Outputs",
        "",
        f"- Feature type summary: `{OUT_FEATURE_TYPES.relative_to(ROOT)}`",
        f"- Barcode suffix summary: `{OUT_BARCODE_SUFFIX.relative_to(ROOT)}`",
        f"- Selected gene presence: `{OUT_SELECTED_PRESENCE.relative_to(ROOT)}`",
        f"- Matrix header: `{OUT_MATRIX_HEADER.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    features = read_features()
    barcodes, suffix = read_barcodes()
    matrix_header = read_matrix_header()
    selected_presence = selected_gene_presence(features)

    feature_summary = features["feature_type"].value_counts().rename_axis("feature_type").reset_index(name="n_features")
    feature_summary.to_csv(OUT_FEATURE_TYPES, sep="\t", index=False)

    if SAMPLE_METADATA.exists() and len(suffix):
        sample_meta = pd.read_csv(SAMPLE_METADATA, sep="\t", dtype={"sample_id": str})
        rna_meta = sample_meta.loc[sample_meta["library_type"].astype(str).str.upper().eq("RNA")].copy()
        suffix = suffix.merge(
            rna_meta[["sample_id", "sample_accession", "diagnosis", "age_at_death_years", "pmi_hours", "tissue"]],
            on="sample_id",
            how="left",
        )
    suffix.to_csv(OUT_BARCODE_SUFFIX, sep="\t", index=False)
    selected_presence.to_csv(OUT_SELECTED_PRESENCE, sep="\t", index=False)

    matrix_header = {
        **matrix_header,
        "n_features_file": len(features),
        "n_barcodes_file": len(barcodes),
        "rows_match_features": matrix_header.get("n_rows") == len(features),
        "cols_match_barcodes": matrix_header.get("n_cols") == len(barcodes),
    }
    pd.DataFrame([matrix_header]).to_csv(OUT_MATRIX_HEADER, sep="\t", index=False)
    write_md(features, suffix, selected_presence, matrix_header)

    print(f"Wrote {OUT_FEATURE_TYPES}")
    print(f"Wrote {OUT_BARCODE_SUFFIX}")
    print(f"Wrote {OUT_SELECTED_PRESENCE}")
    print(f"Wrote {OUT_MATRIX_HEADER}")
    print(f"Wrote {OUT_MD}")
    print(
        "features={}; barcodes={}; selected_present={}/{}; matrix_nnz={}".format(
            len(features),
            len(barcodes),
            int(selected_presence["present_in_gse268609"].sum()),
            len(selected_presence),
            matrix_header.get("n_nonzero_entries"),
        )
    )


if __name__ == "__main__":
    main()
