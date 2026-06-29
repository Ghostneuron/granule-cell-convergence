#!/usr/bin/env python3
"""Extract selected GSE322785 human H5 features and provisional marker calls.

This script keeps the GSE322785 extension conservative: it does not claim
source-author cell labels. Instead, it builds a reusable selected matrix and a
first-pass marker-score annotation scaffold for the downloaded human cerebellar
multiome H5 files.
"""

from __future__ import annotations

import argparse
import gzip
import time
from array import array
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE322785"
RESULTS = ROOT / "Project/results"
PROCESSED = ROOT / "Project/processed/gse322785_human_h5_selected"

PLAN = RESULTS / "gse322785_cerebellar_multiome_download_plan.tsv"
MANIFEST = RESULTS / "gse322785_human_h5_epigenomic_selective_manifest.tsv"

OUT_SUMMARY = RESULTS / "gse322785_human_h5_selected_matrix_summary.tsv"
OUT_PANEL_COVERAGE = RESULTS / "gse322785_human_h5_marker_panel_coverage.tsv"
OUT_CALL_SUMMARY = RESULTS / "gse322785_human_h5_marker_celltype_summary.tsv"
OUT_HIGH_CONF = RESULTS / "gse322785_human_h5_marker_high_confidence_barcodes.tsv.gz"
OUT_MD = RESULTS / "gse322785_human_h5_marker_annotation.md"


MARKER_PANELS: dict[str, list[str]] = {
    "cerebellar_granule": [
        "ATOH1",
        "BARHL1",
        "BARHL2",
        "CBLN1",
        "CBLN3",
        "GABRA6",
        "GRIN2C",
        "PAX6",
        "ROR1",
        "SLC17A7",
        "ZIC1",
        "ZIC2",
    ],
    "purkinje": ["CALB1", "CA8", "PCP2", "FOXP2", "GRID2", "ITPR1", "PVALB", "RORA"],
    "inhibitory_interneuron": ["GAD1", "GAD2", "SLC6A1", "DLX1", "DLX2", "PAX2"],
    "astrocyte_bergmann": ["AQP4", "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "FABP7", "SOX9"],
    "oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "CLDN11"],
    "opc": ["PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2"],
    "microglia": ["CX3CR1", "P2RY12", "AIF1", "C1QA", "TYROBP"],
    "vascular": ["CLDN5", "PECAM1", "VWF", "KDR", "FLT1"],
    "excitatory_neuron": ["SLC17A7", "SLC17A6", "CAMK2A", "RBFOX3", "SNAP25", "SYT1"],
    "neuronal_synaptic": ["RBFOX3", "SNAP25", "SYT1", "SYN1", "MAP2", "TUBB3"],
    "dentate_like_check": ["PROX1", "MEX3A", "C1QL3", "ITPKA", "GLIS3", "EGR3", "CALB2"],
    "morphogenesis": ["GAP43", "STMN2", "STMN3", "DPYSL2", "DPYSL3", "NCAM1", "L1CAM", "CNTN2"],
}

BACKGROUND_PANELS = [
    "astrocyte_bergmann",
    "oligodendrocyte",
    "opc",
    "microglia",
    "vascular",
]

CALL_PANEL_ORDER = [
    "cerebellar_granule",
    "purkinje",
    "inhibitory_interneuron",
    "astrocyte_bergmann",
    "oligodendrocyte",
    "opc",
    "microglia",
    "vascular",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-barcodes", type=int, default=20_000)
    parser.add_argument("--max-barcodes", type=int, default=0, help="Smoke-test limit per sample.")
    parser.add_argument("--no-save-matrix", action="store_true")
    return parser.parse_args()


def norm_gene(gene: object) -> str:
    return str(gene).strip().strip('"').strip("'").upper()


def read_human_samples() -> pd.DataFrame:
    plan = pd.read_csv(PLAN, sep="\t")
    samples = plan.loc[plan["download_tier"].eq("tier1_human_cerebellar_h5")].copy()
    samples["local_path"] = samples["file_name"].map(lambda x: str(BASE / str(x)))
    samples["local_exists"] = samples["local_path"].map(lambda x: Path(x).exists())
    return samples.loc[samples["local_exists"]].copy()


def read_feature_table(h5: h5py.File) -> pd.DataFrame:
    features = h5["matrix/features"]
    return pd.DataFrame(
        {
            "feature_index_0based": np.arange(len(features["name"]), dtype=np.int64),
            "feature_id": [x.decode() for x in features["id"][:]],
            "feature_name": [x.decode() for x in features["name"][:]],
            "feature_type": [x.decode() for x in features["feature_type"][:]],
            "genome": [x.decode() for x in features["genome"][:]],
            "interval": [x.decode() for x in features["interval"][:]],
        }
    )


def build_var(sample_accession: str, features: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, list[int]], pd.DataFrame]:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    sample_manifest = manifest.loc[manifest["sample_accession"].eq(sample_accession)].copy()

    var_rows: list[dict[str, object]] = []
    for feature_idx, sub in sample_manifest.groupby("feature_index_0based", sort=True):
        first = sub.iloc[0]
        var_rows.append(
            {
                "source_feature_index_0based": int(feature_idx),
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
                "peak_categories": ";".join(sorted(x for x in sub["peak_category"].fillna("").astype(str).unique() if x)),
                "selection_reasons": ";".join(sorted(sub["selection_reason"].astype(str).unique())),
                "selected_for_epigenomic_target": True,
                "selected_for_marker_panel": False,
                "marker_panels": "",
            }
        )

    gene_features = features.loc[features["feature_type"].eq("Gene Expression")].copy()
    gene_features["gene_upper"] = gene_features["feature_name"].map(norm_gene)
    gene_to_rows: dict[str, list[int]] = defaultdict(list)
    for row in gene_features.itertuples(index=False):
        gene_to_rows[str(row.gene_upper)].append(int(row.feature_index_0based))

    row_to_panels: dict[int, list[int]] = defaultdict(list)
    existing = {int(row["source_feature_index_0based"]): row for row in var_rows}
    for panel_idx, (panel, genes) in enumerate(MARKER_PANELS.items()):
        for gene in genes:
            for feature_idx in gene_to_rows.get(norm_gene(gene), []):
                row_to_panels[feature_idx].append(panel_idx)
                if feature_idx in existing:
                    existing[feature_idx]["selected_for_marker_panel"] = True
                    panels = set(str(existing[feature_idx]["marker_panels"]).split(";")) if existing[feature_idx]["marker_panels"] else set()
                    panels.add(panel)
                    existing[feature_idx]["marker_panels"] = ";".join(sorted(panels))
                    continue
                f = features.loc[features["feature_index_0based"].eq(feature_idx)].iloc[0]
                rec = {
                    "source_feature_index_0based": feature_idx,
                    "feature_id": f["feature_id"],
                    "feature_name": f["feature_name"],
                    "feature_type": f["feature_type"],
                    "feature_chrom": "",
                    "feature_start": "",
                    "feature_end": "",
                    "linked_genes": norm_gene(f["feature_name"]),
                    "best_priorities": "",
                    "model_terms_supported": "",
                    "target_sets": "",
                    "peak_categories": "",
                    "selection_reasons": "marker_panel_gene",
                    "selected_for_epigenomic_target": False,
                    "selected_for_marker_panel": True,
                    "marker_panels": panel,
                }
                existing[feature_idx] = rec

    var = pd.DataFrame(sorted(existing.values(), key=lambda x: int(x["source_feature_index_0based"]))).reset_index(drop=True)
    var.insert(0, "selected_feature_index", np.arange(len(var), dtype=np.int32))

    coverage_rows = []
    present_genes = set(gene_features["gene_upper"])
    for panel, genes in MARKER_PANELS.items():
        requested = [norm_gene(gene) for gene in genes]
        found = [gene for gene in requested if gene in present_genes]
        coverage_rows.append(
            {
                "sample_accession": sample_accession,
                "panel": panel,
                "panel_gene_count": len(requested),
                "genes_found_in_h5": len(found),
                "genes_found": ",".join(found),
                "genes_missing": ",".join([gene for gene in requested if gene not in present_genes]),
            }
        )
    return var, dict(row_to_panels), pd.DataFrame(coverage_rows)


def rank01(values: np.ndarray, include: np.ndarray) -> np.ndarray:
    ranks = np.zeros(len(values), dtype=np.float32)
    idx = np.flatnonzero(include)
    if len(idx) == 0:
        return ranks
    vals = values[idx]
    order = np.argsort(vals, kind="mergesort")
    sorted_vals = vals[order]
    local = np.empty(len(vals), dtype=np.float32)
    start = 0
    while start < len(vals):
        end = start + 1
        while end < len(vals) and sorted_vals[end] == sorted_vals[start]:
            end += 1
        pct = (start + end) / 2 / max(len(vals), 1)
        local[order[start:end]] = pct
        start = end
    ranks[idx] = local
    return ranks


def marker_calls(obs: pd.DataFrame) -> pd.DataFrame:
    include = obs["analysis_include_basic_qc"].to_numpy(dtype=bool)
    for panel in MARKER_PANELS:
        obs[f"rank_{panel}"] = rank01(obs[f"score_{panel}"].to_numpy(dtype=np.float32), include)

    call_scores = np.vstack([obs[f"score_{panel}"].to_numpy(dtype=np.float32) for panel in CALL_PANEL_ORDER]).T
    call_detected = np.vstack([obs[f"detected_{panel}"].to_numpy(dtype=np.int16) for panel in CALL_PANEL_ORDER]).T
    best_idx = call_scores.argmax(axis=1)
    second = np.partition(call_scores, -2, axis=1)[:, -2]
    best_scores = call_scores[np.arange(len(obs)), best_idx]
    best_detected = call_detected[np.arange(len(obs)), best_idx]
    best_panels = np.array(CALL_PANEL_ORDER, dtype=object)[best_idx]
    best_ranks = np.array([obs[f"rank_{panel}"].to_numpy(dtype=np.float32) for panel in CALL_PANEL_ORDER]).T[
        np.arange(len(obs)), best_idx
    ]

    neuronal = np.maximum.reduce(
        [
            obs["score_cerebellar_granule"].to_numpy(dtype=np.float32),
            obs["score_purkinje"].to_numpy(dtype=np.float32),
            obs["score_inhibitory_interneuron"].to_numpy(dtype=np.float32),
            obs["score_excitatory_neuron"].to_numpy(dtype=np.float32),
            obs["score_neuronal_synaptic"].to_numpy(dtype=np.float32),
        ]
    )
    background = np.maximum.reduce([obs[f"score_{panel}"].to_numpy(dtype=np.float32) for panel in BACKGROUND_PANELS])
    granule = obs["score_cerebellar_granule"].to_numpy(dtype=np.float32)
    purkinje = obs["score_purkinje"].to_numpy(dtype=np.float32)

    calls = np.full(len(obs), "low_information_or_low_qc", dtype=object)
    confidence = np.full(len(obs), "low", dtype=object)
    reason = np.full(len(obs), "below basic RNA/ATAC QC or no marker signal", dtype=object)
    margin = best_scores - second

    high = include & (best_detected >= 2) & (best_ranks >= 0.80) & (margin > 0)
    for panel in CALL_PANEL_ORDER:
        panel_mask = high & (best_panels == panel)
        calls[panel_mask] = f"{panel}_candidate"
        confidence[panel_mask] = np.where(margin[panel_mask] >= 0.05, "high", "medium")
        reason[panel_mask] = "highest marker panel among provisional cerebellar cortex panels"

    ambiguous_neuronal = include & calls.astype(str).__eq__("low_information_or_low_qc") & (neuronal > background) & (neuronal > 0)
    calls[ambiguous_neuronal] = "ambiguous_neuronal"
    confidence[ambiguous_neuronal] = "low"
    reason[ambiguous_neuronal] = "neuronal signal present but no high-confidence panel winner"

    ambiguous_background = include & calls.astype(str).__eq__("low_information_or_low_qc") & (background > 0)
    calls[ambiguous_background] = "ambiguous_non_neuronal_or_niche"
    confidence[ambiguous_background] = "low"
    reason[ambiguous_background] = "non-neuronal marker signal present but no high-confidence panel winner"

    obs["marker_call"] = calls
    obs["marker_confidence"] = confidence
    obs["marker_call_reason"] = reason
    obs["marker_best_panel"] = best_panels
    obs["marker_best_score"] = best_scores
    obs["marker_second_score"] = second
    obs["marker_score_margin"] = margin
    obs["neuronal_minus_background_score"] = neuronal - background
    obs["granule_minus_purkinje_score"] = granule - purkinje
    obs["cerebellar_granule_rank"] = obs["rank_cerebellar_granule"]
    obs["purkinje_rank"] = obs["rank_purkinje"]
    return obs


def process_sample(sample: pd.Series, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_accession = str(sample["sample_accession"])
    donor_id = str(sample["donor_id"])
    path = Path(sample["local_path"])
    out_dir = PROCESSED / donor_id
    out_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][:])
        n_features, n_barcodes_full = shape
        n_barcodes = min(n_barcodes_full, args.max_barcodes) if args.max_barcodes else n_barcodes_full
        features = read_feature_table(h5)
        var, row_to_panels, coverage = build_var(sample_accession, features)

        feature_type = features["feature_type"].to_numpy()
        is_gene = feature_type == "Gene Expression"
        is_peak = feature_type == "Peaks"
        selected_lookup = np.full(n_features, -1, dtype=np.int32)
        selected_lookup[var["source_feature_index_0based"].astype(int).to_numpy()] = var["selected_feature_index"].astype(int).to_numpy()
        marker_row_mask = np.zeros(n_features, dtype=bool)
        for feature_idx in row_to_panels:
            marker_row_mask[feature_idx] = True

        panel_names = list(MARKER_PANELS)
        panel_logsum = np.zeros((len(panel_names), n_barcodes), dtype=np.float32)
        panel_detected = np.zeros((len(panel_names), n_barcodes), dtype=np.uint8)
        n_count_rna = np.zeros(n_barcodes, dtype=np.int32)
        n_feature_rna = np.zeros(n_barcodes, dtype=np.int32)
        n_count_atac = np.zeros(n_barcodes, dtype=np.int32)
        n_feature_atac = np.zeros(n_barcodes, dtype=np.int32)

        out_indptr = np.zeros(n_barcodes + 1, dtype=np.int64)
        out_indices = array("i")
        out_data = array("i")

        h5_indptr = matrix["indptr"]
        h5_indices = matrix["indices"]
        h5_data = matrix["data"]
        t0 = time.time()

        for col_start in range(0, n_barcodes, args.chunk_barcodes):
            col_end = min(n_barcodes, col_start + args.chunk_barcodes)
            ptr = h5_indptr[col_start : col_end + 1]
            start = int(ptr[0])
            end = int(ptr[-1])
            rows = h5_indices[start:end].astype(np.int64)
            vals = h5_data[start:end].astype(np.int32)
            counts = np.diff(ptr).astype(np.int64)
            cols = np.repeat(np.arange(col_start, col_end, dtype=np.int64), counts)

            gene_mask = is_gene[rows]
            if gene_mask.any():
                np.add.at(n_count_rna, cols[gene_mask], vals[gene_mask])
                np.add.at(n_feature_rna, cols[gene_mask], 1)

            peak_mask = is_peak[rows]
            if peak_mask.any():
                np.add.at(n_count_atac, cols[peak_mask], vals[peak_mask])
                np.add.at(n_feature_atac, cols[peak_mask], 1)

            selected = selected_lookup[rows]
            selected_mask = selected >= 0
            if selected_mask.any():
                selected_cols = cols[selected_mask]
                selected_counts = np.bincount(selected_cols - col_start, minlength=col_end - col_start)
                out_indices.extend(selected[selected_mask].astype(np.int32).tolist())
                out_data.extend(vals[selected_mask].astype(np.int32).tolist())
                out_indptr[col_start + 1 : col_end + 1] = len(out_indices) - int(selected_counts.sum()) + np.cumsum(selected_counts)
            else:
                out_indptr[col_start + 1 : col_end + 1] = len(out_indices)

            marker_mask = marker_row_mask[rows]
            if marker_mask.any():
                for feature_idx, col_idx, value in zip(rows[marker_mask], cols[marker_mask], vals[marker_mask]):
                    log_value = np.log1p(float(value))
                    for panel_idx in row_to_panels[int(feature_idx)]:
                        panel_logsum[panel_idx, int(col_idx)] += log_value
                        panel_detected[panel_idx, int(col_idx)] += 1

        barcodes = [x.decode() for x in matrix["barcodes"][:n_barcodes]]

    obs = pd.DataFrame(
        {
            "sample_accession": sample_accession,
            "donor_id": donor_id,
            "barcode_index_0based": np.arange(n_barcodes, dtype=np.int32),
            "barcode": barcodes,
            "nCount_RNA": n_count_rna,
            "nFeature_RNA": n_feature_rna,
            "nCount_ATAC": n_count_atac,
            "nFeature_ATAC": n_feature_atac,
        }
    )
    obs["analysis_include_basic_qc"] = (
        (obs["nFeature_RNA"] >= 100)
        & (obs["nFeature_ATAC"] >= 100)
        & ((obs["nCount_RNA"] + obs["nCount_ATAC"]) >= 500)
    )
    for panel_idx, panel in enumerate(panel_names):
        denom = max(len(MARKER_PANELS[panel]), 1)
        obs[f"score_{panel}"] = panel_logsum[panel_idx] / denom
        obs[f"detected_{panel}"] = panel_detected[panel_idx]
    obs = marker_calls(obs)

    if not args.no_save_matrix:
        X = sparse.csr_matrix(
            (
                np.frombuffer(out_data, dtype=np.int32),
                np.frombuffer(out_indices, dtype=np.int32),
                out_indptr,
            ),
            shape=(n_barcodes, len(var)),
        )
        sparse.save_npz(out_dir / "matrix_barcodes_by_selected_features.npz", X, compressed=True)

    obs.to_csv(out_dir / "cell_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    var.to_csv(out_dir / "var_selected_features.tsv", sep="\t", index=False)

    include = obs["analysis_include_basic_qc"]
    call_summary = (
        obs.groupby(["sample_accession", "donor_id", "marker_call", "marker_confidence"], dropna=False)
        .size()
        .reset_index(name="n_barcodes")
    )
    qc_call_summary = (
        obs.loc[include]
        .groupby(["sample_accession", "donor_id", "marker_call", "marker_confidence"], dropna=False)
        .size()
        .reset_index(name="n_basic_qc_barcodes")
    )
    call_summary = call_summary.merge(
        qc_call_summary,
        on=["sample_accession", "donor_id", "marker_call", "marker_confidence"],
        how="left",
    ).fillna({"n_basic_qc_barcodes": 0})

    high_conf = obs.loc[
        include
        & obs["marker_confidence"].isin(["high", "medium"])
        & ~obs["marker_call"].eq("low_information_or_low_qc")
    ].copy()
    high_conf = high_conf[
        [
            "sample_accession",
            "donor_id",
            "barcode_index_0based",
            "barcode",
            "marker_call",
            "marker_confidence",
            "marker_best_panel",
            "marker_best_score",
            "marker_score_margin",
            "nCount_RNA",
            "nFeature_RNA",
            "nCount_ATAC",
            "nFeature_ATAC",
            "score_cerebellar_granule",
            "score_purkinje",
            "score_astrocyte_bergmann",
            "score_oligodendrocyte",
            "score_microglia",
            "score_vascular",
            "granule_minus_purkinje_score",
            "neuronal_minus_background_score",
        ]
    ]

    elapsed = time.time() - t0
    matrix_path = out_dir / "matrix_barcodes_by_selected_features.npz"
    summary = pd.DataFrame(
        [
            {
                "sample_accession": sample_accession,
                "donor_id": donor_id,
                "source_h5": str(path.relative_to(ROOT)),
                "n_h5_features": n_features,
                "n_h5_barcodes": n_barcodes_full,
                "n_processed_barcodes": n_barcodes,
                "n_selected_features": len(var),
                "n_selected_gene_features": int(var["feature_type"].eq("Gene Expression").sum()),
                "n_selected_peak_features": int(var["feature_type"].eq("Peaks").sum()),
                "selected_matrix_nnz": len(out_indices),
                "n_basic_qc_barcodes": int(include.sum()),
                "n_high_or_medium_confidence_calls": len(high_conf),
                "matrix_path": str(matrix_path.relative_to(ROOT)) if not args.no_save_matrix else "",
                "cell_metadata_path": str((out_dir / "cell_metadata.tsv.gz").relative_to(ROOT)),
                "var_path": str((out_dir / "var_selected_features.tsv").relative_to(ROOT)),
                "elapsed_seconds": round(elapsed, 2),
            }
        ]
    )
    return summary, coverage, call_summary, high_conf


def write_markdown(summary: pd.DataFrame, call_summary: pd.DataFrame) -> None:
    n_samples = summary.shape[0]
    n_barcodes = int(summary["n_processed_barcodes"].sum())
    n_qc = int(summary["n_basic_qc_barcodes"].sum())
    n_conf = int(summary["n_high_or_medium_confidence_calls"].sum())
    n_features = int(summary["n_selected_features"].max()) if not summary.empty else 0
    n_nnz = int(summary["selected_matrix_nnz"].sum()) if not summary.empty else 0

    top_calls = (
        call_summary.groupby("marker_call", as_index=False)["n_basic_qc_barcodes"]
        .sum()
        .sort_values("n_basic_qc_barcodes", ascending=False)
        .head(8)
    )
    top_lines = [f"- `{row.marker_call}`: {int(row.n_basic_qc_barcodes)} QC barcodes." for row in top_calls.itertuples()]

    lines = [
        "# GSE322785 Human H5 Selected Matrix and Marker Annotation",
        "",
        "Date built: 2026-06-26",
        "",
        "## Scope",
        "",
        f"- Human H5 files processed: {n_samples}.",
        f"- Barcodes processed: {n_barcodes}.",
        f"- Basic-QC barcodes: {n_qc}.",
        f"- High/medium-confidence provisional marker calls: {n_conf}.",
        f"- Selected feature columns per sample: up to {n_features}.",
        f"- Selected matrix nonzero entries across samples: {n_nnz}.",
        "",
        "## Top Provisional Calls Among Basic-QC Barcodes",
        "",
        *top_lines,
        "",
        "## Interpretation",
        "",
        "These are marker-score calls, not source-author cell labels. They are sufficient to prioritize label-transfer and selected count extraction, but any manuscript-level chromatin claim still needs verified clustering or label transfer.",
        "",
        "## Outputs",
        "",
        f"- Summary: `{OUT_SUMMARY.relative_to(ROOT)}`",
        f"- Marker panel coverage: `{OUT_PANEL_COVERAGE.relative_to(ROOT)}`",
        f"- Marker cell-type summary: `{OUT_CALL_SUMMARY.relative_to(ROOT)}`",
        f"- High-confidence barcode table: `{OUT_HIGH_CONF.relative_to(ROOT)}`",
        f"- Per-donor selected matrices and metadata: `{PROCESSED.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    summaries = []
    coverages = []
    call_summaries = []
    high_conf_tables = []
    for sample in read_human_samples().to_dict("records"):
        summary, coverage, call_summary, high_conf = process_sample(pd.Series(sample), args)
        summaries.append(summary)
        coverages.append(coverage)
        call_summaries.append(call_summary)
        high_conf_tables.append(high_conf)

    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    coverage = pd.concat(coverages, ignore_index=True) if coverages else pd.DataFrame()
    call_summary = pd.concat(call_summaries, ignore_index=True) if call_summaries else pd.DataFrame()
    high_conf = pd.concat(high_conf_tables, ignore_index=True) if high_conf_tables else pd.DataFrame()

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    coverage.to_csv(OUT_PANEL_COVERAGE, sep="\t", index=False)
    call_summary.to_csv(OUT_CALL_SUMMARY, sep="\t", index=False)
    high_conf.to_csv(OUT_HIGH_CONF, sep="\t", index=False, compression="gzip")
    write_markdown(summary, call_summary)


if __name__ == "__main__":
    main()
