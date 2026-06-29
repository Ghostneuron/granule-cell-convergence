#!/usr/bin/env python3
"""Validate GSE322785 provisional marker calls with selected-gene clustering."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SELECTED_SUMMARY = RESULTS / "gse322785_human_h5_selected_matrix_summary.tsv"

OUT_ASSIGNMENTS = RESULTS / "gse322785_human_h5_cluster_validation_barcode_assignments.tsv.gz"
OUT_CLUSTER_SUMMARY = RESULTS / "gse322785_human_h5_cluster_validation_summary.tsv"
OUT_CALL_ENRICHMENT = RESULTS / "gse322785_human_h5_cluster_validation_marker_call_enrichment.tsv"
OUT_SUPPORT = RESULTS / "gse322785_human_h5_cluster_validation_marker_support.tsv"
OUT_METRICS = RESULTS / "gse322785_human_h5_cluster_validation_metrics.tsv"
OUT_MD = RESULTS / "gse322785_human_h5_cluster_validation.md"

N_CLUSTERS = 12
N_SVD = 20
RANDOM_STATE = 17

MARKER_SCORE_COLUMNS = [
    "score_cerebellar_granule",
    "score_purkinje",
    "score_inhibitory_interneuron",
    "score_astrocyte_bergmann",
    "score_oligodendrocyte",
    "score_opc",
    "score_microglia",
    "score_vascular",
    "score_neuronal_synaptic",
    "score_dentate_like_check",
    "score_morphogenesis",
]


def normalize_selected_gene_matrix(X: sparse.csr_matrix) -> sparse.csr_matrix:
    X = X.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    scale = np.divide(10_000.0, totals, out=np.zeros_like(totals), where=totals > 0)
    X = X.multiply(scale[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X


def marker_label_for_metrics(obs: pd.DataFrame) -> pd.Series:
    label = obs["marker_call"].astype(str).copy()
    label.loc[obs["marker_confidence"].eq("low")] = "low_confidence_or_ambiguous"
    return label


def process_donor(row: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    donor_id = str(row["donor_id"])
    sample_accession = str(row["sample_accession"])
    X = sparse.load_npz(ROOT / row["matrix_path"]).tocsr()
    obs = pd.read_csv(ROOT / row["cell_metadata_path"], sep="\t")
    var = pd.read_csv(ROOT / row["var_path"], sep="\t")

    qc_mask = obs["analysis_include_basic_qc"].astype(bool).to_numpy()
    qc_idx = np.flatnonzero(qc_mask)
    gene_cols = np.flatnonzero(var["feature_type"].eq("Gene Expression").to_numpy())
    X_gene = normalize_selected_gene_matrix(X[qc_idx, :][:, gene_cols])

    n_components = min(N_SVD, X_gene.shape[1] - 1, X_gene.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=RANDOM_STATE)
    emb = svd.fit_transform(X_gene)
    emb_scaled = StandardScaler().fit_transform(emb)

    n_clusters = min(N_CLUSTERS, max(2, X_gene.shape[0] // 200))
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=RANDOM_STATE,
        n_init=5,
        batch_size=4096,
        max_iter=200,
        reassignment_ratio=0.01,
    )
    clusters = kmeans.fit_predict(emb_scaled)

    obs_qc = obs.loc[qc_mask].copy().reset_index(drop=True)
    obs_qc["cluster_id"] = [f"{donor_id}_C{int(x):02d}" for x in clusters]
    obs_qc["cluster_numeric"] = clusters
    obs_qc["svd1"] = emb[:, 0]
    obs_qc["svd2"] = emb[:, 1] if emb.shape[1] > 1 else 0.0

    assign_cols = [
        "sample_accession",
        "donor_id",
        "barcode_index_0based",
        "barcode",
        "cluster_id",
        "cluster_numeric",
        "marker_call",
        "marker_confidence",
        "marker_best_panel",
        "marker_best_score",
        "marker_score_margin",
        "nCount_RNA",
        "nFeature_RNA",
        "nCount_ATAC",
        "nFeature_ATAC",
        "svd1",
        "svd2",
    ]
    assignments = obs_qc[assign_cols].copy()

    cluster_rows = []
    enrichment_rows = []
    for cluster_id, sub in obs_qc.groupby("cluster_id", sort=True):
        calls = sub["marker_call"].value_counts()
        top_call = calls.index[0]
        top_count = int(calls.iloc[0])
        n_cluster = len(sub)
        high_medium = int(sub["marker_confidence"].isin(["high", "medium"]).sum())
        row_out = {
            "sample_accession": sample_accession,
            "donor_id": donor_id,
            "cluster_id": cluster_id,
            "n_basic_qc_barcodes": n_cluster,
            "n_high_or_medium_marker_calls": high_medium,
            "top_marker_call": top_call,
            "top_marker_call_count": top_count,
            "top_marker_call_fraction": top_count / max(n_cluster, 1),
            "mean_nCount_RNA": sub["nCount_RNA"].mean(),
            "mean_nFeature_RNA": sub["nFeature_RNA"].mean(),
            "mean_nCount_ATAC": sub["nCount_ATAC"].mean(),
            "mean_nFeature_ATAC": sub["nFeature_ATAC"].mean(),
        }
        for col in MARKER_SCORE_COLUMNS:
            row_out[f"median_{col}"] = sub[col].median()
        cluster_rows.append(row_out)

        for marker_call, n_call in calls.items():
            enrichment_rows.append(
                {
                    "sample_accession": sample_accession,
                    "donor_id": donor_id,
                    "cluster_id": cluster_id,
                    "marker_call": marker_call,
                    "n_marker_call_barcodes_in_cluster": int(n_call),
                    "n_cluster_barcodes": n_cluster,
                    "fraction_of_cluster": int(n_call) / max(n_cluster, 1),
                }
            )

    cluster_summary = pd.DataFrame(cluster_rows)
    enrichment = pd.DataFrame(enrichment_rows)

    total_by_call = obs_qc["marker_call"].value_counts().to_dict()
    support_rows = []
    for marker_call, sub in enrichment.groupby("marker_call", sort=True):
        best = sub.sort_values("n_marker_call_barcodes_in_cluster", ascending=False).iloc[0]
        support_rows.append(
            {
                "sample_accession": sample_accession,
                "donor_id": donor_id,
                "marker_call": marker_call,
                "n_marker_call_barcodes": int(total_by_call.get(marker_call, 0)),
                "dominant_cluster_id": best["cluster_id"],
                "n_in_dominant_cluster": int(best["n_marker_call_barcodes_in_cluster"]),
                "fraction_in_dominant_cluster": int(best["n_marker_call_barcodes_in_cluster"])
                / max(int(total_by_call.get(marker_call, 0)), 1),
                "dominant_cluster_fraction_marker_call": float(best["fraction_of_cluster"]),
            }
        )
    support = pd.DataFrame(support_rows)

    metric_labels = marker_label_for_metrics(obs_qc)
    metrics = {
        "sample_accession": sample_accession,
        "donor_id": donor_id,
        "n_basic_qc_barcodes": len(obs_qc),
        "n_gene_features_used": len(gene_cols),
        "n_svd_components": n_components,
        "n_clusters": n_clusters,
        "svd_explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
        "adjusted_rand_marker_call_vs_cluster": adjusted_rand_score(metric_labels, clusters),
        "normalized_mutual_info_marker_call_vs_cluster": normalized_mutual_info_score(metric_labels, clusters),
    }
    return assignments, cluster_summary, enrichment, support, metrics


def write_markdown(cluster_summary: pd.DataFrame, support: pd.DataFrame, metrics: pd.DataFrame) -> None:
    n_clusters = cluster_summary.shape[0]
    n_barcodes = int(metrics["n_basic_qc_barcodes"].sum()) if not metrics.empty else 0
    mean_ari = metrics["adjusted_rand_marker_call_vs_cluster"].mean() if not metrics.empty else 0
    mean_nmi = metrics["normalized_mutual_info_marker_call_vs_cluster"].mean() if not metrics.empty else 0
    focal = support.loc[
        support["marker_call"].isin(
            [
                "cerebellar_granule_candidate",
                "purkinje_candidate",
                "astrocyte_bergmann_candidate",
                "oligodendrocyte_candidate",
            ]
        )
    ].copy()
    focal["line"] = focal.apply(
        lambda r: "- `{}` in `{}`: {} barcodes; dominant cluster `{}` captures {:.1%}.".format(
            r["marker_call"],
            r["donor_id"],
            int(r["n_marker_call_barcodes"]),
            r["dominant_cluster_id"],
            float(r["fraction_in_dominant_cluster"]),
        ),
        axis=1,
    )
    lines = [
        "# GSE322785 Human H5 Cluster Validation",
        "",
        "Date built: 2026-06-26",
        "",
        "## Scope",
        "",
        f"- Basic-QC barcodes clustered: {n_barcodes}.",
        f"- Donor-specific clusters: {n_clusters}.",
        f"- Mean adjusted Rand index between provisional marker calls and clusters: {mean_ari:.3f}.",
        f"- Mean normalized mutual information between provisional marker calls and clusters: {mean_nmi:.3f}.",
        "",
        "## Focal Marker-Call Concentration",
        "",
        *(focal["line"].tolist() or ["- No focal marker calls available."]),
        "",
        "## Interpretation",
        "",
        "Selected-gene SVD/k-means clustering provides an internal validation layer for the provisional marker calls. This analysis can support prioritization of marker-group epigenomic contrasts, but it remains weaker than source-author taxonomy or full multimodal clustering.",
        "",
        "## Outputs",
        "",
        f"- Barcode assignments: `{OUT_ASSIGNMENTS.relative_to(ROOT)}`",
        f"- Cluster summary: `{OUT_CLUSTER_SUMMARY.relative_to(ROOT)}`",
        f"- Marker-call enrichment: `{OUT_CALL_ENRICHMENT.relative_to(ROOT)}`",
        f"- Marker support: `{OUT_SUPPORT.relative_to(ROOT)}`",
        f"- Metrics: `{OUT_METRICS.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    summary = pd.read_csv(SELECTED_SUMMARY, sep="\t")
    assignment_tables = []
    cluster_tables = []
    enrichment_tables = []
    support_tables = []
    metric_rows = []
    for row in summary.to_dict("records"):
        assignments, cluster_summary, enrichment, support, metrics = process_donor(pd.Series(row))
        assignment_tables.append(assignments)
        cluster_tables.append(cluster_summary)
        enrichment_tables.append(enrichment)
        support_tables.append(support)
        metric_rows.append(metrics)

    assignments = pd.concat(assignment_tables, ignore_index=True)
    cluster_summary = pd.concat(cluster_tables, ignore_index=True)
    enrichment = pd.concat(enrichment_tables, ignore_index=True)
    support = pd.concat(support_tables, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)

    assignments.to_csv(OUT_ASSIGNMENTS, sep="\t", index=False, compression="gzip")
    cluster_summary.to_csv(OUT_CLUSTER_SUMMARY, sep="\t", index=False)
    enrichment.to_csv(OUT_CALL_ENRICHMENT, sep="\t", index=False)
    support.to_csv(OUT_SUPPORT, sep="\t", index=False)
    metrics.to_csv(OUT_METRICS, sep="\t", index=False)
    write_markdown(cluster_summary, support, metrics)


if __name__ == "__main__":
    main()
