#!/usr/bin/env python3
"""Score GSE322785 selected RNA/ATAC features by provisional marker call."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SUMMARY = RESULTS / "gse322785_human_h5_selected_matrix_summary.tsv"

OUT_FEATURE = RESULTS / "gse322785_human_h5_epigenomic_marker_group_feature_scores.tsv.gz"
OUT_MODULE = RESULTS / "gse322785_human_h5_epigenomic_marker_group_module_scores.tsv"
OUT_CONTRASTS = RESULTS / "gse322785_human_h5_epigenomic_marker_group_contrasts.tsv"
OUT_MD = RESULTS / "gse322785_human_h5_epigenomic_marker_group_scoring.md"

MIN_GROUP_BARCODES = 20
FOCAL_GROUPS = [
    "cerebellar_granule_candidate",
    "purkinje_candidate",
    "astrocyte_bergmann_candidate",
    "oligodendrocyte_candidate",
    "opc_candidate",
    "microglia_candidate",
    "vascular_candidate",
    "ambiguous_neuronal",
]


def split_semicolon(value: object) -> list[str]:
    text = str(value)
    if not text or text == "nan":
        return []
    return [part for part in text.split(";") if part]


def feature_scores_for_group(
    sample_accession: str,
    donor_id: str,
    marker_call: str,
    marker_confidence: str,
    X_group: sparse.csr_matrix,
    var: pd.DataFrame,
) -> pd.DataFrame:
    n = X_group.shape[0]
    sums = np.asarray(X_group.sum(axis=0)).ravel()
    detected = X_group.getnnz(axis=0)
    mean_count = sums / max(n, 1)
    detection_fraction = detected / max(n, 1)
    score = np.log1p(mean_count) + detection_fraction
    out = var[
        [
            "selected_feature_index",
            "source_feature_index_0based",
            "feature_id",
            "feature_name",
            "feature_type",
            "linked_genes",
            "best_priorities",
            "model_terms_supported",
            "target_sets",
            "peak_categories",
            "selection_reasons",
            "selected_for_epigenomic_target",
            "selected_for_marker_panel",
        ]
    ].copy()
    out.insert(0, "sample_accession", sample_accession)
    out.insert(1, "donor_id", donor_id)
    out.insert(2, "marker_call", marker_call)
    out.insert(3, "marker_confidence", marker_confidence)
    out.insert(4, "n_barcodes", n)
    out["mean_count"] = mean_count
    out["detected_barcodes"] = detected
    out["detection_fraction"] = detection_fraction
    out["feature_score_log1p_mean_plus_detection"] = score
    return out


def explode_module_scores(feature_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rec in feature_scores.to_dict("records"):
        model_terms = split_semicolon(rec["model_terms_supported"]) or ["marker_or_unassigned"]
        target_sets = split_semicolon(rec["target_sets"]) or ["marker_or_unassigned"]
        peak_categories = split_semicolon(rec["peak_categories"]) or ["gene_expression_or_marker"]
        linked_genes = split_semicolon(rec["linked_genes"]) or [str(rec["feature_name"])]
        for term in model_terms:
            for target_set in target_sets:
                for peak_category in peak_categories:
                    rows.append(
                        {
                            "sample_accession": rec["sample_accession"],
                            "donor_id": rec["donor_id"],
                            "marker_call": rec["marker_call"],
                            "marker_confidence": rec["marker_confidence"],
                            "feature_type": rec["feature_type"],
                            "model_term_supported": term,
                            "target_set": target_set,
                            "peak_category": peak_category,
                            "n_barcodes": rec["n_barcodes"],
                            "feature_name": rec["feature_name"],
                            "linked_genes": ";".join(linked_genes),
                            "mean_count": rec["mean_count"],
                            "detection_fraction": rec["detection_fraction"],
                            "feature_score": rec["feature_score_log1p_mean_plus_detection"],
                        }
                    )
    exploded = pd.DataFrame(rows)
    if exploded.empty:
        return exploded
    module = (
        exploded.groupby(
            [
                "sample_accession",
                "donor_id",
                "marker_call",
                "marker_confidence",
                "feature_type",
                "model_term_supported",
                "target_set",
                "peak_category",
            ],
            dropna=False,
        )
        .agg(
            n_barcodes=("n_barcodes", "first"),
            n_features=("feature_name", "nunique"),
            n_linked_gene_strings=("linked_genes", "nunique"),
            mean_count=("mean_count", "mean"),
            mean_detection_fraction=("detection_fraction", "mean"),
            mean_feature_score=("feature_score", "mean"),
            median_feature_score=("feature_score", "median"),
        )
        .reset_index()
    )
    return module


def build_contrasts(module: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comparators = [
        "purkinje_candidate",
        "astrocyte_bergmann_candidate",
        "oligodendrocyte_candidate",
        "ambiguous_neuronal",
    ]
    key_cols = ["feature_type", "model_term_supported", "target_set", "peak_category"]
    for donor_id, sub in module.groupby("donor_id"):
        granule = sub.loc[sub["marker_call"].eq("cerebellar_granule_candidate")]
        if granule.empty:
            continue
        for comparator in comparators:
            comp = sub.loc[sub["marker_call"].eq(comparator)]
            if comp.empty:
                continue
            merged = granule.merge(comp, on=key_cols, suffixes=("_granule", "_comparator"))
            for rec in merged.to_dict("records"):
                rows.append(
                    {
                        "donor_id": donor_id,
                        "comparator": comparator,
                        "feature_type": rec["feature_type"],
                        "model_term_supported": rec["model_term_supported"],
                        "target_set": rec["target_set"],
                        "peak_category": rec["peak_category"],
                        "n_barcodes_granule": rec["n_barcodes_granule"],
                        "n_barcodes_comparator": rec["n_barcodes_comparator"],
                        "mean_feature_score_granule": rec["mean_feature_score_granule"],
                        "mean_feature_score_comparator": rec["mean_feature_score_comparator"],
                        "delta_granule_minus_comparator": rec["mean_feature_score_granule"]
                        - rec["mean_feature_score_comparator"],
                        "mean_detection_granule": rec["mean_detection_fraction_granule"],
                        "mean_detection_comparator": rec["mean_detection_fraction_comparator"],
                        "delta_detection_granule_minus_comparator": rec["mean_detection_fraction_granule"]
                        - rec["mean_detection_fraction_comparator"],
                    }
                )
    contrasts = pd.DataFrame(rows)
    if contrasts.empty:
        return contrasts
    pooled = (
        contrasts.groupby(["comparator", "feature_type", "model_term_supported", "target_set", "peak_category"], dropna=False)
        .agg(
            n_donors=("donor_id", "nunique"),
            mean_delta_score=("delta_granule_minus_comparator", "mean"),
            median_delta_score=("delta_granule_minus_comparator", "median"),
            mean_delta_detection=("delta_detection_granule_minus_comparator", "mean"),
            median_delta_detection=("delta_detection_granule_minus_comparator", "median"),
        )
        .reset_index()
    )
    pooled.insert(0, "donor_id", "pooled")
    pooled["n_barcodes_granule"] = ""
    pooled["n_barcodes_comparator"] = ""
    pooled["mean_feature_score_granule"] = ""
    pooled["mean_feature_score_comparator"] = ""
    pooled["delta_granule_minus_comparator"] = pooled["mean_delta_score"]
    pooled["mean_detection_granule"] = ""
    pooled["mean_detection_comparator"] = ""
    pooled["delta_detection_granule_minus_comparator"] = pooled["mean_delta_detection"]
    pooled = pooled[
        [
            "donor_id",
            "comparator",
            "feature_type",
            "model_term_supported",
            "target_set",
            "peak_category",
            "n_barcodes_granule",
            "n_barcodes_comparator",
            "mean_feature_score_granule",
            "mean_feature_score_comparator",
            "delta_granule_minus_comparator",
            "mean_detection_granule",
            "mean_detection_comparator",
            "delta_detection_granule_minus_comparator",
            "n_donors",
            "mean_delta_score",
            "median_delta_score",
            "mean_delta_detection",
            "median_delta_detection",
        ]
    ]
    contrasts["n_donors"] = ""
    contrasts["mean_delta_score"] = ""
    contrasts["median_delta_score"] = ""
    contrasts["mean_delta_detection"] = ""
    contrasts["median_delta_detection"] = ""
    return pd.concat([contrasts, pooled], ignore_index=True, sort=False)


def process_one(row: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    X = sparse.load_npz(ROOT / row["matrix_path"]).tocsr()
    obs = pd.read_csv(ROOT / row["cell_metadata_path"], sep="\t")
    var = pd.read_csv(ROOT / row["var_path"], sep="\t")

    keep = obs["analysis_include_basic_qc"].astype(bool) & obs["marker_call"].isin(FOCAL_GROUPS)
    keep &= obs["marker_confidence"].isin(["high", "medium", "low"])
    obs_keep = obs.loc[keep].copy()

    feature_tables = []
    for (marker_call, marker_confidence), sub in obs_keep.groupby(["marker_call", "marker_confidence"], dropna=False):
        if len(sub) < MIN_GROUP_BARCODES:
            continue
        idx = sub["barcode_index_0based"].to_numpy(dtype=np.int64)
        feature_tables.append(
            feature_scores_for_group(
                str(row["sample_accession"]),
                str(row["donor_id"]),
                str(marker_call),
                str(marker_confidence),
                X[idx, :],
                var,
            )
        )
    feature_scores = pd.concat(feature_tables, ignore_index=True) if feature_tables else pd.DataFrame()
    module_scores = explode_module_scores(feature_scores) if not feature_scores.empty else pd.DataFrame()
    return feature_scores, module_scores


def write_markdown(feature_scores: pd.DataFrame, module: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    n_groups = module[["donor_id", "marker_call", "marker_confidence"]].drop_duplicates().shape[0] if not module.empty else 0
    n_feature_rows = len(feature_scores)
    n_module_rows = len(module)
    n_contrast_rows = len(contrasts)
    pooled = contrasts.loc[contrasts["donor_id"].eq("pooled")].copy() if not contrasts.empty else pd.DataFrame()
    top = pooled.sort_values("mean_delta_score", ascending=False).head(8) if not pooled.empty else pd.DataFrame()
    top_lines = [
        "- `{}` vs `{}` / `{}` / `{}`: delta {:.4g}.".format(
            row.comparator,
            row.feature_type,
            row.model_term_supported,
            row.target_set,
            float(row.mean_delta_score),
        )
        for row in top.itertuples()
    ]

    lines = [
        "# GSE322785 Human H5 Epigenomic Marker-Group Scoring",
        "",
        "Date built: 2026-06-26",
        "",
        "## Scope",
        "",
        f"- Marker groups scored: {n_groups}.",
        f"- Feature-score rows: {n_feature_rows}.",
        f"- Module-score rows: {n_module_rows}.",
        f"- Granule-versus-comparator contrast rows: {n_contrast_rows}.",
        "",
        "## Top Pooled Granule-Positive Contrasts",
        "",
        *(top_lines or ["- No pooled contrasts available."]),
        "",
        "## Interpretation",
        "",
        "This is a provisional marker-group scoring layer. It tests whether selected target genes and nearby ATAC peak features can be summarized by candidate cell groups, but it should not be interpreted as source-author cell-type annotation or causal chromatin evidence.",
        "",
        "## Outputs",
        "",
        f"- Feature scores: `{OUT_FEATURE.relative_to(ROOT)}`",
        f"- Module scores: `{OUT_MODULE.relative_to(ROOT)}`",
        f"- Granule/comparator contrasts: `{OUT_CONTRASTS.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    summary = pd.read_csv(SUMMARY, sep="\t")
    feature_tables = []
    module_tables = []
    for row in summary.to_dict("records"):
        feature_scores, module_scores = process_one(pd.Series(row))
        if not feature_scores.empty:
            feature_tables.append(feature_scores)
        if not module_scores.empty:
            module_tables.append(module_scores)

    feature_scores = pd.concat(feature_tables, ignore_index=True) if feature_tables else pd.DataFrame()
    module = pd.concat(module_tables, ignore_index=True) if module_tables else pd.DataFrame()
    contrasts = build_contrasts(module) if not module.empty else pd.DataFrame()

    feature_scores.to_csv(OUT_FEATURE, sep="\t", index=False, compression="gzip")
    module.to_csv(OUT_MODULE, sep="\t", index=False)
    contrasts.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    write_markdown(feature_scores, module, contrasts)


if __name__ == "__main__":
    main()
