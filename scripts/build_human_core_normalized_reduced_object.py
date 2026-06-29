#!/usr/bin/env python3
"""Build a normalized reduced sparse object for the current human core.

The output is intentionally portable because scanpy/anndata are not available
in the local Python environment: a SciPy sparse matrix plus obs/var TSV files.
The feature set combines high-information genes selected from component gene
metadata with all refined marker-panel genes.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD

from validate_human_core_marker_programs import load_panels


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
PROCESSED = ROOT / "Project/processed"
OUT_DIR = PROCESSED / "human_core_normalized_reduced_object"

SEED_SUMMARY = RESULTS / "human_seed_sparse_object_summary.tsv"
GSE186538_SUMMARY = RESULTS / "gse186538_human_dg_gc_sparse_subset_summary.tsv"
MARKER_SCORES = RESULTS / "human_core_marker_validation_scores.tsv.gz"
ENRICHED_METADATA = RESULTS / "human_core_enriched_cell_metadata.tsv.gz"

OUT_X = OUT_DIR / "X_log1p_cp10k_selected_genes.npz"
OUT_OBS = OUT_DIR / "obs.tsv.gz"
OUT_VAR = OUT_DIR / "var.tsv"
OUT_README = OUT_DIR / "README.md"

OUT_FEATURE_SELECTION = RESULTS / "human_core_normalized_feature_selection.tsv"
OUT_DATASET_SUMMARY = RESULTS / "human_core_normalized_object_dataset_summary.tsv"
OUT_SVD = RESULTS / "human_core_normalized_svd_embedding.tsv.gz"
OUT_SVD_VARIANCE = RESULTS / "human_core_normalized_svd_variance.tsv"
OUT_SVD_CENTROIDS = RESULTS / "human_core_normalized_svd_group_centroids.tsv"
OUT_SVD_PLOT = RESULTS / "human_core_normalized_svd_pc12.png"
OUT_MD = RESULTS / "human_core_normalized_object_summary.md"

TOP_GENES_PER_DATASET = 1500
MAX_NON_MARKER_FEATURES = 4500
N_SVD_COMPONENTS = 30
RANDOM_STATE = 13


def norm_gene(value: object) -> str:
    return str(value).strip().strip('"').strip("'").upper()


def is_excluded_feature(gene: str) -> bool:
    gene = norm_gene(gene)
    return (
        gene.startswith("MT-")
        or gene.startswith("RPL")
        or gene.startswith("RPS")
        or gene in {"MALAT1", "NEAT1", "XIST", "Y_RNA"}
    )


def read_gene_names(path: Path) -> list[str]:
    genes = pd.read_csv(path, sep="\t")
    gene_col = "gene" if "gene" in genes.columns else genes.columns[0]
    return [norm_gene(item) for item in genes[gene_col].astype(str)]


def read_gene_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    out = pd.DataFrame(
        {
            "gene": df[gene_col].map(norm_gene),
            "n_counts": pd.to_numeric(df.get("n_counts", 0), errors="coerce").fillna(0.0),
            "n_cells_detected": pd.to_numeric(df.get("n_cells", 0), errors="coerce").fillna(0.0),
        }
    )
    return out.groupby("gene", as_index=False).sum()


def component_table() -> pd.DataFrame:
    seed = pd.read_csv(SEED_SUMMARY, sep="\t").rename(columns={"library_id": "component_id"})
    seed = seed[
        [
            "dataset",
            "component_id",
            "n_cells",
            "matrix_path",
            "gene_metadata_path",
        ]
    ].copy()

    gse = pd.read_csv(GSE186538_SUMMARY, sep="\t").rename(columns={"subset": "component_id"})
    gse = gse[
        [
            "dataset",
            "component_id",
            "n_cells",
            "matrix_path",
            "gene_metadata_path",
        ]
    ].copy()
    return pd.concat([seed, gse], ignore_index=True)


def load_obs() -> pd.DataFrame:
    scores = pd.read_csv(MARKER_SCORES, sep="\t", low_memory=False)
    enriched = pd.read_csv(ENRICHED_METADATA, sep="\t", low_memory=False)
    extra_cols = [
        "cell_id",
        "gsm",
        "geo_sample_title",
        "specimen_id",
        "geo_tissue",
        "geo_development_stage",
        "geo_age",
        "geo_age_years",
        "geo_sex",
        "sample_match_method",
    ]
    extra_cols = [col for col in extra_cols if col in enriched.columns]
    obs = scores.merge(enriched[extra_cols], on="cell_id", how="left", validate="one_to_one")
    obs["preliminary_qc_pass"] = obs["preliminary_qc_pass"].astype(str).str.lower().isin({"true", "1"})
    obs["analysis_include"] = obs["preliminary_qc_pass"]
    obs["analysis_label"] = obs["marker_call"].map(
        {
            "curated_human_dg_gc_reference": "curated_human_dg_gc_anchor",
            "marker_supported_human_dg_like": "marker_supported_human_dg_like",
            "immature_neuron_or_neurogenic_candidate": "immature_neurogenic_candidate",
            "neuronal_non_dg_or_ambiguous": "hippocampal_neuronal_ambiguous",
            "likely_non_neuronal_background": "likely_non_neuronal_background",
            "cerebellar_marker_high_warning": "cerebellar_marker_high_warning",
            "low_information_or_low_qc": "low_information_or_low_qc",
        }
    ).fillna("other_or_unmapped")
    return obs


def marker_genes() -> set[str]:
    panels = load_panels()
    return {norm_gene(gene) for genes in panels.values() for gene in genes}


def feature_stats(components: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for _, comp in components.iterrows():
        genes = read_gene_metadata(ROOT / comp["gene_metadata_path"])
        genes["dataset"] = comp["dataset"]
        genes["component_id"] = comp["component_id"]
        genes["component_n_cells"] = int(comp["n_cells"])
        rows.append(genes)
    per_component = pd.concat(rows, ignore_index=True)
    dataset_cells = components.groupby("dataset")["n_cells"].sum().to_dict()
    stats = (
        per_component.groupby(["dataset", "gene"], as_index=False)
        .agg(
            total_counts=("n_counts", "sum"),
            detected_cells=("n_cells_detected", "sum"),
            components_with_gene=("component_id", "nunique"),
        )
    )
    stats["dataset_n_cells"] = stats["dataset"].map(dataset_cells).astype(float)
    stats["detected_fraction"] = stats["detected_cells"] / stats["dataset_n_cells"].replace(0, np.nan)
    stats["mean_counts"] = stats["total_counts"] / stats["dataset_n_cells"].replace(0, np.nan)
    stats["feature_score"] = np.log1p(stats["total_counts"]) * np.sqrt(stats["detected_fraction"].clip(lower=0))
    stats["excluded_feature"] = stats["gene"].map(is_excluded_feature)
    return stats


def select_features(stats: pd.DataFrame, marker_gene_set: set[str]) -> pd.DataFrame:
    selected: dict[str, dict[str, object]] = {}

    for gene in marker_gene_set:
        selected[gene] = {
            "gene": gene,
            "selection_reason": "marker_panel",
            "max_feature_score": float(stats.loc[stats["gene"] == gene, "feature_score"].max() or 0.0),
        }

    candidate_rows = []
    for dataset, sub in stats.groupby("dataset"):
        detection_threshold = max(20, 0.0025 * float(sub["dataset_n_cells"].iloc[0]))
        candidates = sub[
            (~sub["excluded_feature"])
            & (sub["detected_cells"] >= detection_threshold)
            & (sub["total_counts"] >= 100)
        ].sort_values("feature_score", ascending=False)
        for _, row in candidates.head(TOP_GENES_PER_DATASET).iterrows():
            candidate_rows.append(row)

    candidate_df = pd.DataFrame(candidate_rows)
    if len(candidate_df):
        max_scores = candidate_df.groupby("gene", as_index=False)["feature_score"].max().sort_values("feature_score", ascending=False)
        non_marker = max_scores[~max_scores["gene"].isin(marker_gene_set)].head(MAX_NON_MARKER_FEATURES)
        for _, row in non_marker.iterrows():
            selected[row["gene"]] = {
                "gene": row["gene"],
                "selection_reason": "high_information_gene",
                "max_feature_score": float(row["feature_score"]),
            }

    selected_df = pd.DataFrame(selected.values())
    if selected_df.empty:
        raise RuntimeError("No selected features")

    coverage = (
        stats[stats["gene"].isin(selected_df["gene"])]
        .pivot_table(index="gene", columns="dataset", values="detected_fraction", aggfunc="max", fill_value=0.0)
        .reset_index()
    )
    selected_df = selected_df.merge(coverage, on="gene", how="left")
    score_max = stats.groupby("gene", as_index=False)["feature_score"].max().rename(columns={"feature_score": "max_feature_score_from_stats"})
    selected_df = selected_df.merge(score_max, on="gene", how="left")
    selected_df["max_feature_score"] = selected_df["max_feature_score_from_stats"].fillna(selected_df["max_feature_score"])
    selected_df = selected_df.drop(columns=["max_feature_score_from_stats"])
    selected_df["is_marker_panel_gene"] = selected_df["gene"].isin(marker_gene_set)
    return selected_df.sort_values(["is_marker_panel_gene", "max_feature_score", "gene"], ascending=[False, False, True]).reset_index(drop=True)


def build_matrix(obs: pd.DataFrame, selected_features: pd.DataFrame) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    feature_genes = selected_features["gene"].tolist()
    feature_to_idx = {gene: idx for idx, gene in enumerate(feature_genes)}
    matrix_parts: list[sparse.csr_matrix] = []
    obs_parts: list[pd.DataFrame] = []
    presence_rows: list[dict[str, object]] = []

    for matrix_path, group in obs.groupby("matrix_path", sort=False):
        include_mask = group["analysis_include"].to_numpy(dtype=bool)
        if not include_mask.any():
            continue

        matrix_abs = ROOT / matrix_path
        gene_meta_path = ROOT / group["cell_metadata_path"].iloc[0].replace("cell_metadata.tsv.gz", "gene_metadata.tsv.gz")
        if not gene_meta_path.exists():
            comp = component_table()
            match = comp.loc[comp["matrix_path"] == matrix_path, "gene_metadata_path"]
            if not len(match):
                raise FileNotFoundError(f"No gene metadata for {matrix_path}")
            gene_meta_path = ROOT / match.iloc[0]

        genes = read_gene_names(gene_meta_path)
        gene_to_col: dict[str, int] = {}
        duplicate_count = 0
        for col_idx, gene in enumerate(genes):
            if gene in gene_to_col:
                duplicate_count += 1
                continue
            gene_to_col[gene] = col_idx

        present = [gene for gene in feature_genes if gene in gene_to_col]
        source_cols = [gene_to_col[gene] for gene in present]
        global_cols = np.array([feature_to_idx[gene] for gene in present], dtype=np.int32)
        row_idx = np.flatnonzero(include_mask)

        print(f"Building normalized block {matrix_path}: {len(row_idx)} cells x {len(present)} features", flush=True)
        matrix = sparse.load_npz(matrix_abs).tocsr()
        sub = matrix[row_idx, :][:, source_cols].astype(np.float32).tocsr()

        n_counts = pd.to_numeric(group.loc[group["analysis_include"], "n_counts"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
        scale = np.divide(10000.0, n_counts, out=np.zeros_like(n_counts), where=n_counts > 0)
        sub = sub.multiply(scale[:, None]).tocsr()
        sub.data = np.log1p(sub.data).astype(np.float32)
        remapped_indices = global_cols[sub.indices]
        block = sparse.csr_matrix((sub.data, remapped_indices, sub.indptr.copy()), shape=(sub.shape[0], len(feature_genes)))
        block.sort_indices()

        matrix_parts.append(block)
        obs_parts.append(group.loc[group["analysis_include"]].copy())
        presence_rows.append(
            {
                "dataset": group["dataset"].iloc[0],
                "component_id": group["component_id"].iloc[0],
                "matrix_path": matrix_path,
                "analysis_cells": int(len(row_idx)),
                "selected_features_present": int(len(present)),
                "selected_features_absent": int(len(feature_genes) - len(present)),
                "duplicate_feature_names_in_source": int(duplicate_count),
            }
        )

    if not matrix_parts:
        raise RuntimeError("No matrix blocks were built")

    X = sparse.vstack(matrix_parts, format="csr")
    obs_out = pd.concat(obs_parts, ignore_index=True)
    presence = pd.DataFrame(presence_rows)
    return X, obs_out, presence


def summarize_dataset(obs: pd.DataFrame) -> pd.DataFrame:
    summary = (
        obs.groupby(["dataset", "analysis_label", "marker_call", "marker_state"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            median_counts=("n_counts", "median"),
            median_genes=("n_genes", "median"),
            median_percent_mt=("percent_mt", "median"),
            median_dentate_identity=("score_dentate_identity", "median"),
            median_cerebellar_identity=("score_cerebellar_identity", "median"),
            median_structural_program=("structural_program_mean", "median"),
        )
        .reset_index()
    )
    totals = summary.groupby("dataset")["n_cells"].transform("sum")
    summary["fraction_of_dataset_analysis_cells"] = summary["n_cells"] / totals
    return summary.sort_values(["dataset", "n_cells"], ascending=[True, False])


def run_svd(X: sparse.csr_matrix, obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    n_components = min(N_SVD_COMPONENTS, X.shape[1] - 1, X.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=RANDOM_STATE, n_iter=7)
    embedding = svd.fit_transform(X)
    pc_cols = [f"svd_{i + 1}" for i in range(n_components)]
    emb = pd.DataFrame(embedding, columns=pc_cols)
    emb.insert(0, "cell_id", obs["cell_id"].to_numpy())
    emb.insert(1, "dataset", obs["dataset"].to_numpy())
    emb.insert(2, "analysis_label", obs["analysis_label"].to_numpy())
    emb.insert(3, "marker_call", obs["marker_call"].to_numpy())
    emb.insert(4, "component_id", obs["component_id"].to_numpy())

    variance = pd.DataFrame(
        {
            "component": pc_cols,
            "explained_variance": svd.explained_variance_,
            "explained_variance_ratio": svd.explained_variance_ratio_,
            "singular_value": svd.singular_values_,
        }
    )

    centroid_cols = ["dataset", "analysis_label"]
    centroids = emb.groupby(centroid_cols, dropna=False)[pc_cols].median().reset_index()
    counts = emb.groupby(centroid_cols, dropna=False).size().reset_index(name="n_cells")
    centroids = centroids.merge(counts, on=centroid_cols, how="left")
    return emb, variance, centroids


def plot_svd(embedding: pd.DataFrame) -> None:
    if len(embedding) > 20000:
        sample_parts = []
        for _, group in embedding.groupby(["dataset", "analysis_label"], sort=False):
            n_group = min(len(group), max(200, int(20000 * len(group) / len(embedding))))
            sample_parts.append(group.sample(n_group, random_state=RANDOM_STATE))
        sample = pd.concat(sample_parts, ignore_index=True)
        if len(sample) > 22000:
            sample = sample.sample(22000, random_state=RANDOM_STATE)
    else:
        sample = embedding

    label_order = [
        "curated_human_dg_gc_anchor",
        "marker_supported_human_dg_like",
        "immature_neurogenic_candidate",
        "hippocampal_neuronal_ambiguous",
        "likely_non_neuronal_background",
        "cerebellar_marker_high_warning",
        "low_information_or_low_qc",
    ]
    colors = {
        "curated_human_dg_gc_anchor": "#227c70",
        "marker_supported_human_dg_like": "#2a9d8f",
        "immature_neurogenic_candidate": "#d95f02",
        "hippocampal_neuronal_ambiguous": "#457b9d",
        "likely_non_neuronal_background": "#6c757d",
        "cerebellar_marker_high_warning": "#b56576",
        "low_information_or_low_qc": "#adb5bd",
    }

    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    for label in label_order:
        sub = sample[sample["analysis_label"] == label]
        if len(sub) == 0:
            continue
        ax.scatter(sub["svd_1"], sub["svd_2"], s=4, alpha=0.45, linewidths=0, color=colors.get(label, "#333333"), label=label)
    ax.set_xlabel("SVD 1")
    ax.set_ylabel("SVD 2")
    ax.set_title("Human core normalized reduced object")
    ax.legend(frameon=False, fontsize=8, markerscale=3, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_SVD_PLOT, dpi=180)
    plt.close(fig)


def write_readme(X: sparse.csr_matrix, obs: pd.DataFrame, var: pd.DataFrame) -> None:
    lines = [
        "# Human Core Normalized Reduced Object",
        "",
        "Date built: 2026-06-21",
        "",
        "This directory stores a portable reduced single-cell/nucleus object for the constructed human dentate/hippocampal core.",
        "",
        "## Files",
        "",
        "- `X_log1p_cp10k_selected_genes.npz`: cells-by-genes CSR matrix with log1p(CP10K) values.",
        "- `obs.tsv.gz`: cell metadata in the same row order as the matrix.",
        "- `var.tsv`: selected feature metadata in the same column order as the matrix.",
        "",
        "## Shape",
        "",
        f"- Cells/nuclei: {X.shape[0]}",
        f"- Selected genes: {X.shape[1]}",
        f"- Non-zero matrix entries: {X.nnz}",
        "",
        "The object uses preliminary QC-pass cells only. The full all-cell metadata remains in `Project/results/human_core_enriched_cell_metadata.tsv.gz`.",
        "",
    ]
    OUT_README.write_text("\n".join(lines))


def write_summary(
    X: sparse.csr_matrix,
    obs_all: pd.DataFrame,
    obs_analysis: pd.DataFrame,
    var: pd.DataFrame,
    dataset_summary: pd.DataFrame,
    variance: pd.DataFrame,
    presence: pd.DataFrame,
) -> None:
    all_counts = obs_all.groupby("dataset").size()
    analysis_counts = obs_analysis.groupby("dataset").size()
    lines = [
        "# Human Core Normalized Reduced Object",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This checkpoint builds a portable normalized object from the marker-validated human core. It includes preliminary QC-pass cells from `GSE185277`, `GSE185553`, and the `GSE186538` DG GC subset.",
        "",
        "## Object Shape",
        "",
        f"- Analysis cells/nuclei: {X.shape[0]} of {len(obs_all)} total human-core cells.",
        f"- Selected genes: {X.shape[1]} ({int(var['is_marker_panel_gene'].sum())} marker-panel genes plus high-information genes).",
        f"- Non-zero log-normalized entries: {X.nnz}.",
        "",
        "## Dataset Inclusion",
        "",
    ]
    for dataset in sorted(all_counts.index):
        total = int(all_counts.loc[dataset])
        included = int(analysis_counts.get(dataset, 0))
        lines.append(f"- `{dataset}`: {included} / {total} cells included ({included / total * 100:.1f}%).")

    lines.extend(
        [
            "",
            "## First Embedding Check",
            "",
            f"- SVD components computed: {len(variance)}.",
            f"- First component explained variance ratio: {variance['explained_variance_ratio'].iloc[0]:.4f}.",
            f"- First five components cumulative explained variance ratio: {variance['explained_variance_ratio'].head(5).sum():.4f}.",
            "",
            "## Interpretation",
            "",
            "- This is an analysis-ready reduced object, not yet a final integrated atlas. It is suitable for first-pass PCA/SVD, label QC, and dataset-aware module checks.",
            "- `GSE186538` remains the curated human DG GC anchor; `GSE185277` and `GSE185553` provide candidate-rich human hippocampal contexts.",
            "- This object is the substrate for tuned-label and dataset-aware module testing; `GSE325391` and the RNA branch of `GSE268609` have since been added as selected-feature human dentate/hippocampal bridge objects.",
            "",
            "## Outputs",
            "",
            f"- Normalized matrix: `{OUT_X.relative_to(ROOT)}`",
            f"- Observation metadata: `{OUT_OBS.relative_to(ROOT)}`",
            f"- Feature metadata: `{OUT_VAR.relative_to(ROOT)}`",
            f"- Dataset summary: `{OUT_DATASET_SUMMARY.relative_to(ROOT)}`",
            f"- SVD embedding: `{OUT_SVD.relative_to(ROOT)}`",
            f"- SVD plot: `{OUT_SVD_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    components = component_table()
    obs_all = load_obs()

    print("Selecting features", flush=True)
    marker_gene_set = marker_genes()
    stats = feature_stats(components)
    var = select_features(stats, marker_gene_set)
    var.insert(0, "feature_index", range(len(var)))
    var.to_csv(OUT_VAR, sep="\t", index=False)
    var.to_csv(OUT_FEATURE_SELECTION, sep="\t", index=False)

    print("Building normalized matrix", flush=True)
    X, obs_analysis, presence = build_matrix(obs_all, var)
    sparse.save_npz(OUT_X, X, compressed=True)
    obs_analysis.to_csv(OUT_OBS, sep="\t", index=False, compression="gzip")

    dataset_summary = summarize_dataset(obs_analysis)
    dataset_summary.to_csv(OUT_DATASET_SUMMARY, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)

    print("Running SVD", flush=True)
    embedding, variance, centroids = run_svd(X, obs_analysis)
    embedding.to_csv(OUT_SVD, sep="\t", index=False, compression="gzip", float_format="%.6g")
    variance.to_csv(OUT_SVD_VARIANCE, sep="\t", index=False, float_format="%.6g")
    centroids.to_csv(OUT_SVD_CENTROIDS, sep="\t", index=False, float_format="%.6g")
    presence.to_csv(RESULTS / "human_core_normalized_component_feature_presence.tsv", sep="\t", index=False)
    plot_svd(embedding)
    write_readme(X, obs_analysis, var)
    write_summary(X, obs_all, obs_analysis, var, dataset_summary, variance, presence)

    print(f"Wrote {OUT_X}")
    print(f"Wrote {OUT_OBS}")
    print(f"Wrote {OUT_VAR}")
    print(f"Wrote {OUT_SVD}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
