#!/usr/bin/env python3
"""Full-transcriptome/available-feature diffusion trajectories for primary core.

This is the manuscript-level upgrade after the endpoint-gene graph audit. It
uses all locally available genes for true full matrices, highly variable genes
for geometry, and a common PCA -> kNN -> diffusion/shortest-path pseudotime
workflow. Very large human bridge datasets that only have selected-feature
objects in the local workspace are retained but explicitly marked as selected
feature support rather than strict full-transcriptome trajectories.
"""

from __future__ import annotations

import csv
import gc
import gzip
import io
import math
import os
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.linalg import eigsh
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from sklearn.utils.sparsefuncs import mean_variance_axis


csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_2005_paper_endpoint_pseudotime_audit import PAPER_ENDPOINT_MODULES  # noqa: E402


CALLS = RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"
CORE = RESULTS / "integrated_primary_core_datasets.tsv"
NICHE_GENE_SETS = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
GSE325391_LABELS = RESULTS / "gse325391_human_core_label_projection.tsv.gz"
GSE268609_LABELS = RESULTS / "gse268609_human_core_label_projection.tsv.gz"

OUT_CELL = RESULTS / "primary_core_full_transcriptome_diffusion_cell_scores.tsv.gz"
OUT_DATASET_SUMMARY = RESULTS / "primary_core_full_transcriptome_diffusion_dataset_summary.tsv"
OUT_GROUP_SUMMARY = RESULTS / "primary_core_full_transcriptome_diffusion_group_summary.tsv"
OUT_MODULE_CORR = RESULTS / "primary_core_full_transcriptome_diffusion_module_correlations.tsv"
OUT_FIG_IMPACT = RESULTS / "primary_core_full_transcriptome_diffusion_fig1_5_impact.tsv"
OUT_PLOT = RESULTS / "primary_core_full_transcriptome_diffusion_overview.png"
OUT_MD = RESULTS / "primary_core_full_transcriptome_diffusion.md"

MAX_CELLS = 8000
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 18

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


GSE104323_ORDER = {
    "RGL_young": 0.0,
    "RGL": 1.0,
    "nIPC": 2.0,
    "nIPC-perin": 3.0,
    "Neuroblast": 4.0,
    "Immature-GC": 5.0,
    "GC-juv": 6.0,
    "GC-adult": 7.0,
}
GSE292261_ORDER = {f"DG_P{age}": float(age) for age in [5, 7, 10, 15, 28]}
GSE122357_STAGE = {
    "GSM3464549_P0": ("GSM3464549_P0.csv.gz", "P0", 0.0),
    "GSM3464550_P8a": ("GSM3464550_P8a.csv.gz", "P8a", 8.1),
    "GSM3464551_P8b": ("GSM3464551_P8b.csv.gz", "P8b", 8.2),
}


PROGENITOR_ROOT_GENES = {
    "SOX2",
    "NES",
    "HES1",
    "HES5",
    "ASCL1",
    "ATOH1",
    "EOMES",
    "DCX",
    "NEUROD1",
    "MKI67",
    "TOP2A",
    "PCNA",
}
MATURE_GRANULE_GENES = {
    "PROX1",
    "RBFOX3",
    "GABRA6",
    "CALB2",
    "SYN1",
    "SYN2",
    "SYP",
    "SNAP25",
    "MAP2",
    "CAMK2A",
    "GRIN2B",
}
PLOT_MODULES = [
    "proliferation_brdU_proxy",
    "immature_progenitor_state",
    "neuronal_differentiation_maturation",
    "tgf_smad_pai1_response",
    "bdnf_erk_response",
    "secreted_stop_candidate_axis",
    "downstream_neurite_morphology",
    "downstream_synaptic_excitability",
]


def rel(path: Path | str) -> str:
    path = Path(path)
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').strip("'").upper()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def read_sparse_npz(path: Path) -> sparse.csr_matrix:
    obj = np.load(path, allow_pickle=False)
    if {"data", "indices", "indptr", "shape"}.issubset(obj.files):
        return sparse.csr_matrix((obj["data"], obj["indices"], obj["indptr"]), shape=tuple(obj["shape"]))
    return sparse.load_npz(path).tocsr()


def read_gene_vector(path: Path, column: str = "gene") -> list[str]:
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in {column, "gene", "gene_name"}, low_memory=False)
    if column in df:
        return df[column].astype(str).tolist()
    if "gene" in df:
        return df["gene"].astype(str).tolist()
    return df.iloc[:, 0].astype(str).tolist()


def read_10x_features(path: Path) -> list[str]:
    genes: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                genes.append(parts[1])
            elif parts:
                genes.append(parts[0])
    return genes


def read_barcodes(path: Path) -> list[str]:
    with open_text(path) as fh:
        return [line.rstrip("\n").split("\t")[0] for line in fh if line.strip()]


def load_calls() -> pd.DataFrame:
    cols = ["dataset", "sample", "cell_id", "group", "candidate_call", "call_confidence"]
    calls = pd.read_csv(CALLS, sep="\t", usecols=cols, low_memory=False, dtype=str)
    calls["cell_id"] = calls["cell_id"].astype(str)
    return calls


def broad_class(call: str) -> str:
    if call == "candidate_dentate_granule":
        return "dentate_candidate"
    if call == "candidate_cerebellar_granule":
        return "cerebellar_candidate"
    if call == "dentate_like_low_support":
        return "dentate_low_support"
    if call == "cerebellum_dentate_panel_high_warning":
        return "cerebellar_lineage_warning"
    if call == "known_non_dentate_reference":
        return "non_dentate_background"
    return "other_or_ambiguous"


def make_obs_from_calls(
    calls: pd.DataFrame,
    dataset: str,
    keep: Callable[[pd.DataFrame], pd.Series],
    axis_type: str,
    order_map: dict[str, float] | None = None,
) -> pd.DataFrame:
    sub = calls.loc[calls["dataset"].eq(dataset)].copy()
    sub = sub.loc[keep(sub)].copy()
    sub["axis_type"] = axis_type
    sub["axis_label"] = sub["group"]
    sub["axis_order"] = sub["group"].map(order_map or {})
    sub["broad_class"] = sub["candidate_call"].map(broad_class)
    sub["source_cell_id"] = sub["cell_id"].astype(str)
    sub["root_cell"] = False
    return sub[
        [
            "dataset",
            "sample",
            "cell_id",
            "source_cell_id",
            "axis_type",
            "axis_label",
            "axis_order",
            "candidate_call",
            "call_confidence",
            "broad_class",
            "root_cell",
        ]
    ].reset_index(drop=True)


def stratified_obs_sample(obs: pd.DataFrame, max_cells: int = MAX_CELLS) -> pd.DataFrame:
    if len(obs) <= max_cells:
        return obs.reset_index(drop=True)
    rng = np.random.default_rng(7)
    strata = obs[["axis_label", "broad_class", "sample"]].fillna("missing").astype(str).agg("|".join, axis=1)
    counts = strata.value_counts()
    allocation = {}
    base = max(1, max_cells // max(1, len(counts)))
    remaining = max_cells
    for key, count in counts.items():
        take = min(int(count), base)
        allocation[key] = take
        remaining -= take
    if remaining > 0:
        for key, count in counts.sort_values(ascending=False).items():
            room = int(count) - allocation[key]
            if room <= 0:
                continue
            add = min(room, remaining)
            allocation[key] += add
            remaining -= add
            if remaining <= 0:
                break
    chosen: list[int] = []
    for key, n_take in allocation.items():
        idx = np.flatnonzero(strata.to_numpy() == key)
        if len(idx) <= n_take:
            chosen.extend(idx.tolist())
        else:
            chosen.extend(rng.choice(idx, size=n_take, replace=False).tolist())
    sampled = obs.iloc[sorted(chosen)].copy()
    sampled["trajectory_sampling"] = f"stratified_{len(sampled)}_of_{len(obs)}"
    return sampled.reset_index(drop=True)


def finalize_obs(obs: pd.DataFrame, dataset: str, source_scope: str, evidence_grade: str) -> pd.DataFrame:
    obs = obs.copy()
    obs["dataset"] = dataset
    obs["source_scope"] = source_scope
    obs["trajectory_evidence_grade"] = evidence_grade
    if "trajectory_sampling" not in obs:
        obs["trajectory_sampling"] = f"all_{len(obs)}"
    return obs.reset_index(drop=True)


def load_wide_gene_by_cell(
    path: Path,
    obs: pd.DataFrame,
    delimiter: str,
    *,
    tar_member: str | None = None,
    header_has_gene_col: bool = True,
) -> tuple[pd.DataFrame, sparse.csr_matrix, list[str]]:
    selected = set(obs["source_cell_id"].astype(str))
    if tar_member:
        tf = tarfile.open(path)
        raw = tf.extractfile(tar_member)
        if raw is None:
            tf.close()
            raise FileNotFoundError(tar_member)
        fh = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")
        close_tf = tf
    else:
        fh = open_text(path)
        close_tf = None

    genes: list[str] = []
    rows: list[np.ndarray] = []
    try:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader)
        cell_fields = header[1:] if header_has_gene_col else header
        cells = [str(cell).strip().strip('"') for cell in cell_fields]
        positions = [(idx, cell) for idx, cell in enumerate(cells) if cell in selected]
        obs_order = obs.set_index("source_cell_id").loc[[cell for _, cell in positions]].reset_index(drop=True)
        col_indices = [idx + (1 if header_has_gene_col else 0) for idx, _ in positions]
        for row in reader:
            if not row:
                continue
            gene = str(row[0]).strip().strip('"') if header_has_gene_col else f"feature_{len(genes)+1}"
            vals = np.zeros(len(col_indices), dtype=np.float32)
            any_nonzero = False
            for out_idx, col_idx in enumerate(col_indices):
                if col_idx >= len(row):
                    continue
                try:
                    value = float(row[col_idx])
                except ValueError:
                    value = 0.0
                if value != 0:
                    any_nonzero = True
                    vals[out_idx] = value
            if any_nonzero:
                genes.append(gene)
                rows.append(vals)
    finally:
        fh.close()
        if close_tf is not None:
            close_tf.close()
    if not rows:
        raise RuntimeError(f"No expression rows retained for {path}")
    X = sparse.csr_matrix(np.vstack(rows).T)
    return obs_order.reset_index(drop=True), X, genes


def load_obs_by_gene_table(path: Path, obs: pd.DataFrame, delimiter: str) -> tuple[pd.DataFrame, sparse.csr_matrix, list[str]]:
    selected = set(obs["source_cell_id"].astype(str))
    df = pd.read_csv(path, sep=delimiter, low_memory=False)
    cell_col = df.columns[0]
    df[cell_col] = df[cell_col].astype(str)
    df = df.loc[df[cell_col].isin(selected)].copy()
    obs_order = obs.set_index("source_cell_id").loc[df[cell_col].tolist()].reset_index(drop=True)
    genes = [str(col) for col in df.columns[1:]]
    X = sparse.csr_matrix(df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32))
    keep = np.asarray(X.sum(axis=0)).ravel() > 0
    return obs_order.reset_index(drop=True), X[:, keep].tocsr(), [g for g, k in zip(genes, keep) if k]


def load_10x_matrix(
    matrix: Path,
    features: Path,
    barcodes: Path,
    obs: pd.DataFrame,
) -> tuple[pd.DataFrame, sparse.csr_matrix, list[str]]:
    genes = read_10x_features(features)
    barcodes_list = read_barcodes(barcodes)
    selected = set(obs["source_cell_id"].astype(str))
    keep_idx = [idx for idx, barcode in enumerate(barcodes_list) if barcode in selected]
    selected_barcodes = [barcodes_list[idx] for idx in keep_idx]
    obs_order = obs.set_index("source_cell_id").loc[selected_barcodes].reset_index(drop=True)
    X = mmread(matrix).tocsr().T
    X = X[keep_idx, :].tocsr()
    keep_genes = np.asarray(X.sum(axis=0)).ravel() > 0
    return obs_order, X[:, keep_genes].tocsr(), [g for g, keep in zip(genes, keep_genes) if keep]


def subset_sparse_by_obs(
    matrix_path: Path,
    genes_path: Path,
    metadata_path: Path,
    obs: pd.DataFrame,
    *,
    gene_col: str = "gene",
) -> tuple[pd.DataFrame, sparse.csr_matrix, list[str]]:
    meta = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    if "cell_id" not in meta:
        raise RuntimeError(f"{metadata_path} lacks cell_id")
    meta["cell_id"] = meta["cell_id"].astype(str)
    selected = set(obs["source_cell_id"].astype(str))
    keep_idx = [idx for idx, cid in enumerate(meta["cell_id"].tolist()) if cid in selected]
    selected_ids = [meta["cell_id"].iloc[idx] for idx in keep_idx]
    obs_order = obs.set_index("source_cell_id").loc[selected_ids].reset_index(drop=True)
    X = read_sparse_npz(matrix_path).tocsr()[keep_idx, :]
    if genes_path.suffix == ".gz":
        gene_df = pd.read_csv(genes_path, sep="\t", low_memory=False)
    else:
        gene_df = pd.read_csv(genes_path, sep="\t", low_memory=False)
    if gene_col not in gene_df:
        gene_col = "gene" if "gene" in gene_df else gene_df.columns[0]
    genes = gene_df[gene_col].astype(str).tolist()
    keep_genes = np.asarray(X.sum(axis=0)).ravel() > 0
    return obs_order, X[:, keep_genes].tocsr(), [g for g, keep in zip(genes, keep_genes) if keep]


def module_gene_sets() -> dict[str, set[str]]:
    modules: dict[str, set[str]] = {}
    for module in PAPER_ENDPOINT_MODULES:
        modules[module["module_id"]] = {canon(gene) for gene in module["genes"]}
    if NICHE_GENE_SETS.exists():
        niche = pd.read_csv(NICHE_GENE_SETS, sep="\t")
        for module_id, sub in niche.groupby("module_id", sort=False):
            modules[module_id] = {canon(gene) for gene in sub["canonical_gene"]}
    modules["root_progenitor_immature"] = set(PROGENITOR_ROOT_GENES)
    modules["root_mature_granule"] = set(MATURE_GRANULE_GENES)
    return modules


def normalize_log1p(X: sparse.csr_matrix) -> sparse.csr_matrix:
    X = X.tocsr().astype(np.float32)
    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
    lib[lib <= 0] = 1.0
    scale = 10000.0 / lib
    Xn = X.multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def select_hvg(Xlog: sparse.csr_matrix, genes: list[str], n_hvg: int = N_HVG) -> tuple[sparse.csr_matrix, list[str]]:
    if Xlog.shape[1] <= n_hvg:
        return Xlog, genes
    detected = np.asarray((Xlog > 0).sum(axis=0)).ravel()
    mean, var = mean_variance_axis(Xlog, axis=0)
    valid = (detected >= 3) & np.isfinite(var) & (var > 0)
    if valid.sum() < min(50, Xlog.shape[1]):
        valid = np.isfinite(var) & (var > 0)
    idx = np.flatnonzero(valid)
    if len(idx) > n_hvg:
        ranked = idx[np.argsort(var[idx])[-n_hvg:]]
    else:
        ranked = idx
    ranked = np.sort(ranked)
    return Xlog[:, ranked].tocsr(), [genes[i] for i in ranked]


def zscore_dense(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float32)
    sd = a.std(axis=0)
    sd[sd == 0] = 1.0
    return (a - a.mean(axis=0)) / sd


def module_scores(Xlog: sparse.csr_matrix, genes: list[str], modules: dict[str, set[str]]) -> pd.DataFrame:
    gene_to_cols: dict[str, list[int]] = defaultdict(list)
    for idx, gene in enumerate(genes):
        gene_to_cols[canon(gene)].append(idx)
    rows: dict[str, np.ndarray | int] = {}
    for module_id, wanted in modules.items():
        cols = sorted({col for gene in wanted for col in gene_to_cols.get(gene, [])})
        rows[f"{module_id}_n_genes"] = len(cols)
        if not cols:
            rows[module_id] = np.full(Xlog.shape[0], np.nan, dtype=np.float32)
            continue
        vals = Xlog[:, cols].toarray().astype(np.float32)
        vals = zscore_dense(vals)
        raw = vals.mean(axis=1)
        rows[module_id] = pd.Series(raw).rank(method="average", pct=True).to_numpy(dtype=np.float32)
    return pd.DataFrame(rows)


def infer_root_from_modules(scores: pd.DataFrame, min_roots: int = 20) -> tuple[np.ndarray, str]:
    progenitor = pd.to_numeric(scores.get("root_progenitor_immature"), errors="coerce")
    mature = pd.to_numeric(scores.get("root_mature_granule"), errors="coerce")
    root_score = progenitor.fillna(0.5) - mature.fillna(0.5)
    if root_score.nunique(dropna=True) < 3:
        n_roots = max(3, min(min_roots, len(root_score) // 20 if len(root_score) >= 100 else len(root_score)))
        mask = np.zeros(len(root_score), dtype=bool)
        mask[:n_roots] = True
        return mask, "fallback_first_cells_no_marker_root"
    q = max(0.90, 1.0 - max(min_roots, int(len(root_score) * 0.05)) / max(1, len(root_score)))
    threshold = root_score.quantile(q)
    mask = root_score.ge(threshold).to_numpy()
    if mask.sum() < 3:
        idx = np.argsort(root_score.to_numpy())[-max(3, min_roots):]
        mask = np.zeros(len(root_score), dtype=bool)
        mask[idx] = True
    return mask, "inferred_high_progenitor_low_mature_root"


def trajectory_geometry(
    Xhvg: sparse.csr_matrix,
    root_mask: np.ndarray,
    source_pseudotime: pd.Series | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    n_cells = Xhvg.shape[0]
    n_components = max(2, min(N_PCS, n_cells - 1, Xhvg.shape[1] - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=11)
    pcs = svd.fit_transform(Xhvg)
    pcs = zscore_dense(pcs)

    k = max(2, min(N_NEIGHBORS, n_cells - 1))
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nn.fit(pcs)
    distances, indices = nn.kneighbors(pcs)
    rows = np.repeat(np.arange(n_cells), k)
    cols = indices[:, 1:].reshape(-1)
    d = distances[:, 1:].reshape(-1)
    sigma = np.median(d[d > 0]) if np.any(d > 0) else 1.0
    weights = np.exp(-((d / sigma) ** 2)).astype(np.float32)
    affinity = sparse.coo_matrix((weights, (rows, cols)), shape=(n_cells, n_cells)).tocsr()
    affinity = affinity.maximum(affinity.T)

    roots = np.where(root_mask)[0]
    if roots.size == 0 and source_pseudotime is not None:
        s = pd.to_numeric(source_pseudotime, errors="coerce")
        roots = np.flatnonzero(s.le(s.quantile(0.05)).fillna(False).to_numpy())
    if roots.size == 0:
        roots = np.array([0], dtype=int)
    n = affinity.shape[0]
    root_node = n
    # Convert affinity to sparse edge costs without materializing non-edges.
    cost = affinity.tocoo()
    edge_cost = sparse.coo_matrix(
        (np.maximum(1e-6, 1.0 - cost.data), (cost.row, cost.col)),
        shape=cost.shape,
    ).tocsr()
    graph = sparse.vstack(
        [
            sparse.hstack([edge_cost, sparse.csr_matrix((n, 1))]),
            sparse.csr_matrix((1, n + 1)),
        ]
    ).tolil()
    for root in roots:
        graph[root_node, root] = 1e-6
        graph[root, root_node] = 1e-6
    dist = dijkstra(graph.tocsr(), directed=False, indices=root_node)[:n]
    finite = np.isfinite(dist)
    pseudotime = np.zeros(n, dtype=np.float32)
    if finite.any() and dist[finite].max() > 0:
        pseudotime[finite] = (dist[finite] / dist[finite].max()).astype(np.float32)
        pseudotime[~finite] = 1.0

    diff = np.zeros((n_cells, 2), dtype=np.float32)
    try:
        transition = normalize(affinity, norm="l1", axis=1)
        vals, vecs = eigsh(transition.T @ transition, k=min(4, n_cells - 1), which="LM")
        order = np.argsort(vals)[::-1]
        vecs = vecs[:, order]
        if vecs.shape[1] >= 3:
            diff[:, 0] = vecs[:, 1]
            diff[:, 1] = vecs[:, 2]
        elif vecs.shape[1] >= 2:
            diff[:, 0] = vecs[:, 1]
    except Exception:
        diff[:, 0] = pcs[:, 0]
        diff[:, 1] = pcs[:, 1] if pcs.shape[1] > 1 else 0.0
    return pseudotime, pcs[:, :2], diff, float(svd.explained_variance_ratio_.sum())


def run_dataset(
    dataset: str,
    obs: pd.DataFrame,
    loader: Callable[[pd.DataFrame], tuple[pd.DataFrame, sparse.csr_matrix, list[str]]],
    *,
    source_scope: str,
    evidence_grade: str,
    root_strategy: str,
    source_pseudotime_col: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    obs = stratified_obs_sample(obs)
    obs, X, genes = loader(obs)
    obs = finalize_obs(obs, dataset, source_scope, evidence_grade)
    obs["n_total_genes_available"] = len(genes)
    obs["n_cells_before_sampling"] = obs["trajectory_sampling"].str.extract(r"_of_(\d+)$")[0].fillna(len(obs)).astype(int)

    Xlog = normalize_log1p(X)
    scores = module_scores(Xlog, genes, module_gene_sets())
    if obs["root_cell"].astype(bool).sum() >= 3:
        root_mask = obs["root_cell"].astype(bool).to_numpy()
        actual_root_strategy = root_strategy
    elif source_pseudotime_col and source_pseudotime_col in obs:
        s = pd.to_numeric(obs[source_pseudotime_col], errors="coerce")
        root_mask = s.le(s.quantile(0.05)).fillna(False).to_numpy()
        actual_root_strategy = f"source_{source_pseudotime_col}_lowest_5pct"
    else:
        root_mask, inferred = infer_root_from_modules(scores)
        actual_root_strategy = inferred
    Xhvg, hvg_genes = select_hvg(Xlog, genes)
    pseudotime, pcs, diff, explained = trajectory_geometry(
        Xhvg,
        root_mask,
        obs[source_pseudotime_col] if source_pseudotime_col and source_pseudotime_col in obs else None,
    )

    out = pd.concat([obs.reset_index(drop=True), scores.reset_index(drop=True)], axis=1)
    out["full_transcriptome_pseudotime"] = pseudotime
    out["pca1"] = pcs[:, 0]
    out["pca2"] = pcs[:, 1]
    out["diffusion1"] = diff[:, 0]
    out["diffusion2"] = diff[:, 1]
    out["root_used"] = root_mask
    out["root_strategy"] = actual_root_strategy
    out["n_hvg_used"] = len(hvg_genes)
    out["svd_variance_explained"] = explained

    summary = {
        "dataset": dataset,
        "source_scope": source_scope,
        "trajectory_evidence_grade": evidence_grade,
        "root_strategy": actual_root_strategy,
        "n_cells": len(out),
        "n_roots": int(root_mask.sum()),
        "n_total_genes_available": len(genes),
        "n_hvg_used": len(hvg_genes),
        "sampling": out["trajectory_sampling"].iloc[0],
        "svd_variance_explained": explained,
    }
    if out["axis_order"].notna().sum() >= 10 and out["axis_order"].nunique(dropna=True) >= 3:
        rho, p = stats.spearmanr(out["full_transcriptome_pseudotime"], out["axis_order"], nan_policy="omit")
        summary["rho_vs_known_order"] = float(rho)
        summary["p_vs_known_order"] = float(p)
    else:
        summary["rho_vs_known_order"] = np.nan
        summary["p_vs_known_order"] = np.nan
    if source_pseudotime_col and source_pseudotime_col in out and out[source_pseudotime_col].notna().sum() >= 10:
        xy = out[["full_transcriptome_pseudotime", source_pseudotime_col]].dropna()
        if len(xy) >= 10 and xy[source_pseudotime_col].nunique() >= 3:
            rho, p = stats.spearmanr(xy["full_transcriptome_pseudotime"], xy[source_pseudotime_col])
            summary["rho_vs_source_pseudotime"] = float(rho)
            summary["p_vs_source_pseudotime"] = float(p)
        else:
            summary["rho_vs_source_pseudotime"] = np.nan
            summary["p_vs_source_pseudotime"] = np.nan
    else:
        summary["rho_vs_source_pseudotime"] = np.nan
        summary["p_vs_source_pseudotime"] = np.nan
    return out, summary


def dataset_jobs(calls: pd.DataFrame) -> list[tuple[str, pd.DataFrame, Callable, dict[str, str]]]:
    jobs = []

    obs = make_obs_from_calls(
        calls,
        "GSE104323",
        lambda d: d["group"].isin(GSE104323_ORDER),
        "adult_dentate_lineage_state",
        GSE104323_ORDER,
    )
    obs["root_cell"] = obs["axis_label"].isin(["RGL_young", "RGL"])
    jobs.append(
        (
            "GSE104323",
            obs,
            lambda o: load_wide_gene_by_cell(
                EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz", o, "\t"
            ),
            {
                "source_scope": "full_transcriptome_raw_matrix",
                "evidence_grade": "primary_lineage_pseudotime",
                "root_strategy": "curated_RGL_RGL_young_root",
            },
        )
    )

    obs = make_obs_from_calls(
        calls,
        "GSE95752",
        lambda d: d["candidate_call"].isin(["candidate_dentate_granule", "dentate_like_low_support"]),
        "dentate_maturation_validation_cluster",
        None,
    )
    jobs.append(
        (
            "GSE95752",
            obs,
            lambda o: load_wide_gene_by_cell(EXTERNAL / "GEO/GSE95752/GSE95752_C1_expression_data.tab.gz", o, "\t"),
            {
                "source_scope": "full_transcriptome_raw_matrix",
                "evidence_grade": "supporting_intrinsic_diffusion",
                "root_strategy": "marker_inferred_immature_root",
            },
        )
    )

    obs = make_obs_from_calls(
        calls,
        "GSE292261",
        lambda d: d["candidate_call"].isin(["candidate_dentate_granule", "dentate_like_low_support"]),
        "postnatal_dentate_age",
        GSE292261_ORDER,
    )
    obs["root_cell"] = obs["axis_label"].eq("DG_P5")
    jobs.append(
        (
            "GSE292261",
            obs,
            lambda o: load_obs_by_gene_table(
                EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz", o, ","
            ),
            {
                "source_scope": "full_transcriptome_raw_matrix",
                "evidence_grade": "primary_stage_pseudotime",
                "root_strategy": "earliest_postnatal_P5_root",
            },
        )
    )

    obs = make_obs_from_calls(
        calls,
        "GSE214309",
        lambda d: d["candidate_call"].isin(["candidate_dentate_granule", "dentate_like_low_support"]),
        "adult_dentate_immature_mature_activity_state",
        None,
    )
    state_order = {
        "immature_1hr": 0.0,
        "immature_4hr": 0.1,
        "immatureactive_1hr": 0.4,
        "immatureactive_4hr": 0.5,
        "mature_1hr": 1.0,
        "mature_4hr": 1.1,
        "matureactive_1hr": 1.4,
        "matureactive_4hr": 1.5,
    }
    obs["axis_order"] = obs["axis_label"].map(state_order)
    obs["root_cell"] = obs["axis_label"].isin(["immature_1hr", "immature_4hr"])
    jobs.append(
        (
            "GSE214309",
            obs,
            lambda o: load_wide_gene_by_cell(
                EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz", o, ",", header_has_gene_col=True
            ),
            {
                "source_scope": "full_transcriptome_ensembl_raw_matrix",
                "evidence_grade": "state_axis_pseudotime",
                "root_strategy": "immature_DGC_state_root",
            },
        )
    )

    for sample, (member, stage, order) in GSE122357_STAGE.items():
        pass
    obs = make_obs_from_calls(
        calls,
        "GSE122357",
        lambda d: d["candidate_call"].isin(["candidate_cerebellar_granule", "cerebellum_dentate_panel_high_warning"]),
        "postnatal_cerebellar_age",
        None,
    )
    obs["axis_label"] = obs["sample"].map({k: v[1] for k, v in GSE122357_STAGE.items()})
    obs["axis_order"] = obs["sample"].map({k: v[2] for k, v in GSE122357_STAGE.items()})
    obs["root_cell"] = obs["axis_label"].eq("P0")

    def load_gse122357(o: pd.DataFrame):
        obs_parts, x_parts, gene_union = [], [], []
        all_genes: list[str] = []
        for sample, (member, _, _) in GSE122357_STAGE.items():
            sub = o.loc[o["sample"].eq(sample)].copy()
            if sub.empty:
                continue
            sub_obs, sub_x, genes = load_wide_gene_by_cell(
                EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar", sub, ",", tar_member=member
            )
            sub_obs["source_cell_id"] = sub_obs["cell_id"]
            sub_obs["cell_id"] = sample + ":" + sub_obs["cell_id"].astype(str)
            obs_parts.append(sub_obs)
            x_parts.append(sub_x)
            gene_union.extend(genes)
            all_genes.append(genes)
        union = sorted(set(gene_union))
        colmap = {gene: idx for idx, gene in enumerate(union)}
        aligned = []
        for x, genes in zip(x_parts, all_genes):
            cols = [colmap[g] for g in genes]
            coo = x.tocoo()
            aligned.append(sparse.coo_matrix((coo.data, (coo.row, np.array(cols)[coo.col])), shape=(x.shape[0], len(union))).tocsr())
        return pd.concat(obs_parts, ignore_index=True), sparse.vstack(aligned).tocsr(), union

    jobs.append(
        (
            "GSE122357",
            obs,
            load_gse122357,
            {
                "source_scope": "full_transcriptome_raw_matrix",
                "evidence_grade": "primary_stage_pseudotime",
                "root_strategy": "earliest_postnatal_P0_root",
            },
        )
    )

    obs = make_obs_from_calls(
        calls,
        "GSE165657",
        lambda d: d["candidate_call"].isin(["candidate_cerebellar_granule", "cerebellum_dentate_panel_high_warning"]),
        "human_cerebellum_intrinsic_state",
        None,
    )
    jobs.append(
        (
            "GSE165657",
            obs,
            lambda o: load_10x_matrix(
                EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_matrix.mtx.gz",
                EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_genes.tsv.gz",
                EXTERNAL / "GEO/GSE165657/GSE165657_Cerebellum_aggr_barcodes.tsv.gz",
                o,
            ),
            {
                "source_scope": "full_transcriptome_10x_matrix",
                "evidence_grade": "supporting_intrinsic_diffusion",
                "root_strategy": "marker_inferred_immature_root",
            },
        )
    )

    obs = make_obs_from_calls(
        calls,
        "GSE312658",
        lambda d: d["candidate_call"].isin(["candidate_cerebellar_granule", "cerebellum_dentate_panel_high_warning"]),
        "mouse_cerebellum_ctrl_cKO_intrinsic_state",
        None,
    )

    def load_gse312658(o: pd.DataFrame):
        parts_obs, parts_x, genes_list, union_genes = [], [], [], []
        for sample, prefix in [("Ctrl", "GSM9350909_Ctrl"), ("cKO", "GSM9350910_cKO")]:
            sub = o.loc[o["sample"].eq(sample)].copy()
            if sub.empty:
                continue
            sub_obs, sub_x, genes = load_10x_matrix(
                EXTERNAL / f"GEO/GSE312658/{prefix}_matrix.mtx.gz",
                EXTERNAL / f"GEO/GSE312658/{prefix}_features.tsv.gz",
                EXTERNAL / f"GEO/GSE312658/{prefix}_barcodes.tsv.gz",
                sub,
            )
            sub_obs["cell_id"] = sample + ":" + sub_obs["cell_id"].astype(str)
            parts_obs.append(sub_obs)
            parts_x.append(sub_x)
            genes_list.append(genes)
            union_genes.extend(genes)
        union = sorted(set(union_genes))
        colmap = {gene: idx for idx, gene in enumerate(union)}
        aligned = []
        for x, genes in zip(parts_x, genes_list):
            cols = np.array([colmap[g] for g in genes])
            coo = x.tocoo()
            aligned.append(sparse.coo_matrix((coo.data, (coo.row, cols[coo.col])), shape=(x.shape[0], len(union))).tocsr())
        return pd.concat(parts_obs, ignore_index=True), sparse.vstack(aligned).tocsr(), union

    jobs.append(
        (
            "GSE312658",
            obs,
            load_gse312658,
            {
                "source_scope": "full_transcriptome_10x_matrix",
                "evidence_grade": "perturbation_context_diffusion",
                "root_strategy": "marker_inferred_immature_root",
            },
        )
    )

    meta = pd.read_csv(
        PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/cell_metadata.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    obs = pd.DataFrame(
        {
            "dataset": "GSE186538",
            "sample": meta["samplename"].astype(str),
            "cell_id": meta["cell_id"].astype(str),
            "source_cell_id": meta["cell_id"].astype(str),
            "axis_type": "human_DG_GC_taxonomy_intrinsic_state",
            "axis_label": meta["cluster"].astype(str),
            "axis_order": np.nan,
            "candidate_call": "candidate_dentate_granule",
            "call_confidence": "taxonomy_anchor",
            "broad_class": "dentate_candidate",
            "root_cell": False,
        }
    )
    jobs.append(
        (
            "GSE186538",
            obs,
            lambda o: subset_sparse_by_obs(
                PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/matrix_cells_by_genes.npz",
                PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/gene_metadata.tsv.gz",
                PROCESSED / "human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/cell_metadata.tsv.gz",
                o,
            ),
            {
                "source_scope": "full_transcriptome_DG_GC_sparse_subset",
                "evidence_grade": "human_DG_intrinsic_diffusion",
                "root_strategy": "marker_inferred_immature_root",
            },
        )
    )

    labels325 = pd.read_csv(GSE325391_LABELS, sep="\t", low_memory=False)
    labels325 = labels325.loc[labels325["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])].copy()
    labels325["source_cell_id"] = labels325["cell_id"].astype(str)
    labels325["root_cell"] = labels325["source_anchor_label"].eq("adult_human_dg_differentiating_anchor")
    obs = pd.DataFrame(
        {
            "dataset": "GSE325391",
            "sample": labels325["sample"].astype(str),
            "cell_id": labels325["cell_id"].astype(str),
            "source_cell_id": labels325["source_cell_id"].astype(str),
            "axis_type": "adult_human_dentate_source_sling_pseudotime",
            "axis_label": labels325["sub_cell_type"].astype(str),
            "axis_order": pd.to_numeric(labels325.get("slingAvgPseudotime"), errors="coerce"),
            "candidate_call": labels325["source_anchor_label"].astype(str),
            "call_confidence": labels325["source_label_confidence"].astype(str),
            "broad_class": "dentate_candidate",
            "root_cell": labels325["root_cell"].fillna(False),
            "source_sling_pseudotime": pd.to_numeric(labels325.get("slingAvgPseudotime"), errors="coerce"),
        }
    )
    jobs.append(
        (
            "GSE325391",
            obs,
            lambda o: subset_sparse_by_obs(
                PROCESSED / "gse325391_adult_dg_selected/matrix_cells_by_selected_genes.npz",
                PROCESSED / "gse325391_adult_dg_selected/var_selected_features.tsv",
                PROCESSED / "gse325391_adult_dg_selected/cell_metadata.tsv.gz",
                o,
            ),
            {
                "source_scope": "selected_feature_bridge_2169_genes",
                "evidence_grade": "selected_feature_source_pseudotime_validation",
                "root_strategy": "source_differentiating_or_low_sling_root",
                "source_pseudotime_col": "source_sling_pseudotime",
            },
        )
    )

    labels268 = pd.read_csv(GSE268609_LABELS, sep="\t", low_memory=False, dtype={"sample_id": str})
    keep_labels = {
        "human_dg_like_high_confidence",
        "human_dg_like_candidate",
        "immature_neurogenic_candidate",
        "hippocampal_neuronal_ambiguous",
    }
    labels268 = labels268.loc[
        labels268["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
        & labels268["projected_label"].isin(keep_labels)
    ].copy()
    obs = pd.DataFrame(
        {
            "dataset": "GSE268609",
            "sample": labels268["sample_id"].astype(str),
            "cell_id": labels268["cell_id"].astype(str),
            "source_cell_id": labels268["cell_id"].astype(str),
            "axis_type": "human_hippocampus_DG_projected_intrinsic_state",
            "axis_label": labels268["projected_label"].astype(str),
            "axis_order": np.nan,
            "candidate_call": labels268["projected_label"].astype(str),
            "call_confidence": labels268["projection_confidence"].astype(str),
            "broad_class": np.where(
                labels268["projected_label"].eq("immature_neurogenic_candidate"),
                "dentate_immature_candidate",
                "dentate_candidate_or_hippocampal_neuronal",
            ),
            "root_cell": labels268["projected_label"].eq("immature_neurogenic_candidate"),
        }
    )
    jobs.append(
        (
            "GSE268609",
            obs,
            lambda o: subset_sparse_by_obs(
                PROCESSED / "gse268609_rna_selected/matrix_cells_by_selected_genes.npz",
                PROCESSED / "gse268609_rna_selected/var_selected_features.tsv",
                PROCESSED / "gse268609_rna_selected/cell_metadata.tsv.gz",
                o,
            ),
            {
                "source_scope": "selected_feature_bridge_2169_genes",
                "evidence_grade": "selected_feature_human_context_diffusion",
                "root_strategy": "projected_immature_neurogenic_root",
            },
        )
    )
    return jobs


def build_group_summary(cells: pd.DataFrame) -> pd.DataFrame:
    score_cols = [col for col in cells.columns if col in module_gene_sets()]
    group_cols = ["dataset", "source_scope", "axis_type", "axis_label", "candidate_call", "broad_class"]
    summary = (
        cells.groupby(group_cols, dropna=False, sort=False)
        .agg(
            n_cells=("cell_id", "count"),
            median_axis_order=("axis_order", "median"),
            median_pseudotime=("full_transcriptome_pseudotime", "median"),
            mean_pseudotime=("full_transcriptome_pseudotime", "mean"),
            root_fraction=("root_used", "mean"),
            **{f"median_{col}": (col, "median") for col in score_cols},
        )
        .reset_index()
    )
    return summary.sort_values(["dataset", "median_axis_order", "axis_label", "candidate_call"], na_position="last")


def build_module_correlations(cells: pd.DataFrame) -> pd.DataFrame:
    modules = module_gene_sets()
    rows = []
    for dataset, sub in cells.groupby("dataset", sort=False):
        for module_id in modules:
            if module_id not in sub:
                continue
            xy = sub[["full_transcriptome_pseudotime", module_id]].dropna()
            n_genes_col = f"{module_id}_n_genes"
            n_genes = int(pd.to_numeric(sub[n_genes_col], errors="coerce").max()) if n_genes_col in sub else 0
            if len(xy) >= 20 and xy[module_id].nunique() >= 3 and xy["full_transcriptome_pseudotime"].nunique() >= 3:
                rho, p = stats.spearmanr(xy["full_transcriptome_pseudotime"], xy[module_id])
            else:
                rho, p = np.nan, np.nan
            rows.append(
                {
                    "dataset": dataset,
                    "module_id": module_id,
                    "n_module_genes_detected": n_genes,
                    "spearman_rho_vs_full_transcriptome_pseudotime": rho,
                    "spearman_p": p,
                    "n_cells": len(xy),
                }
            )
    return pd.DataFrame(rows)


def build_fig_impact(dataset_summary: pd.DataFrame, module_corr: pd.DataFrame) -> pd.DataFrame:
    def rho(ds: str, module: str) -> float:
        sub = module_corr.loc[module_corr["dataset"].eq(ds) & module_corr["module_id"].eq(module)]
        if sub.empty:
            return np.nan
        return float(sub["spearman_rho_vs_full_transcriptome_pseudotime"].iloc[0])

    mature_support = np.nanmedian(
        [
            rho("GSE104323", "neuronal_differentiation_maturation"),
            rho("GSE292261", "neuronal_differentiation_maturation"),
            rho("GSE122357", "neuronal_differentiation_maturation"),
        ]
    )
    tgf_support = np.nanmedian(
        [
            rho("GSE104323", "tgf_smad_pai1_response"),
            rho("GSE122357", "tgf_smad_pai1_response"),
            rho("GSE214309", "tgf_smad_pai1_response"),
        ]
    )
    construction_support = np.nanmedian(
        [
            rho("GSE104323", "downstream_neurite_morphology"),
            rho("GSE122357", "downstream_neurite_morphology"),
            rho("GSE165657", "downstream_neurite_morphology"),
            rho("GSE186538", "downstream_neurite_morphology"),
        ]
    )
    rows = [
        {
            "figure": "Figure 1",
            "current_result": "primary-core design and convergence hypothesis",
            "trajectory_effect": "no structural change",
            "recommended_action": "Add note that trajectory support is now available as an Aim 2 refinement.",
            "basis": f"{dataset_summary['dataset'].nunique()} primary datasets processed in trajectory layer.",
        },
        {
            "figure": "Figure 2",
            "current_result": "ortholog candidate tiers",
            "trajectory_effect": "candidate tiers unchanged for now",
            "recommended_action": "Do not rerank genes until trajectory-aware gene tests are built; use module overlays as supporting evidence.",
            "basis": "The diffusion layer tests module timing, not gene-level tier statistics.",
        },
        {
            "figure": "Figure 3",
            "current_result": "specificity/niche/circuit model",
            "trajectory_effect": "refines interpretation",
            "recommended_action": "Mention that downstream modules are stage-windowed and not uniquely granule-specific.",
            "basis": f"median neurite-morphology rho across key trajectory datasets = {construction_support:.3f}.",
        },
        {
            "figure": "Figure 4",
            "current_result": "transcriptomic configuration model",
            "trajectory_effect": "strengthens but makes it stage-aware",
            "recommended_action": "Add a trajectory panel or caption sentence that configuration is identity-coupled and maturation-window dependent.",
            "basis": f"median neuronal-maturation rho across developmental anchors = {mature_support:.3f}.",
        },
        {
            "figure": "Figure 5",
            "current_result": "final Aim2/Aim3 model",
            "trajectory_effect": "should be revised/expanded",
            "recommended_action": "Add pseudotime as an explicit layer between regional fate/niche and terminal morphology.",
            "basis": f"median TGF/SMAD trajectory rho across selected anchors = {tgf_support:.3f}.",
        },
    ]
    return pd.DataFrame(rows)


def plot_overview(dataset_summary: pd.DataFrame, module_corr: pd.DataFrame) -> None:
    plot = module_corr.loc[module_corr["module_id"].isin(PLOT_MODULES)].copy()
    if plot.empty:
        return
    pivot = plot.pivot_table(
        index="dataset",
        columns="module_id",
        values="spearman_rho_vs_full_transcriptome_pseudotime",
        aggfunc="first",
    )
    pivot = pivot.reindex(dataset_summary["dataset"].tolist())
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), gridspec_kw={"width_ratios": [1.2, 1.0]})
    im = axes[0].imshow(pivot.to_numpy(dtype=float), vmin=-0.75, vmax=0.75, cmap="RdBu_r", aspect="auto")
    axes[0].set_yticks(np.arange(len(pivot.index)))
    axes[0].set_yticklabels(pivot.index)
    axes[0].set_xticks(np.arange(len(pivot.columns)))
    axes[0].set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    axes[0].set_title("Module correlation with full-transcriptome pseudotime")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="Spearman rho")

    axes[1].barh(dataset_summary["dataset"], dataset_summary["n_cells"], color="#4f6f8f")
    axes[1].set_xlabel("Cells/nuclei in trajectory")
    axes[1].set_title("Trajectory cell counts")
    axes[1].grid(axis="x", color="#dddddd", linewidth=0.6)
    axes[1].invert_yaxis()
    fig.suptitle("Primary-core full-transcriptome/available-feature diffusion trajectories", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def fmt(value: object, digits: int = 3) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(val):
        return "NA"
    return f"{val:.{digits}f}"


def write_report(
    cells: pd.DataFrame,
    dataset_summary: pd.DataFrame,
    module_corr: pd.DataFrame,
    fig_impact: pd.DataFrame,
) -> None:
    strict_full = dataset_summary["source_scope"].str.contains("full_transcriptome|full_transcriptome_DG", regex=True).sum()
    selected = dataset_summary["source_scope"].str.contains("selected_feature").sum()
    lines = [
        "# Primary-Core Full-Transcriptome Diffusion Trajectories",
        "",
        "## Scope",
        "",
        "This layer upgrades the endpoint-gene graph audit to a highly-variable-gene diffusion/pseudotime workflow across the strict 10-dataset primary core.",
        "",
        "For datasets with local full matrices, geometry uses all available genes before HVG selection. `GSE325391` and `GSE268609` are retained as selected-feature bridge trajectories because the local analysis objects contain 2,169 selected genes rather than full cell-by-gene matrices ready for in-memory trajectory analysis.",
        "",
        "The method is dependency-light: log1p(CP10K) normalization, HVG selection, TruncatedSVD, kNN graph, diffusion components, and shortest-path pseudotime from curated or inferred roots.",
        "",
        "## Coverage",
        "",
        f"- Primary datasets processed: {dataset_summary['dataset'].nunique()}/10.",
        f"- Cells/nuclei/spots in trajectory tables after deterministic sampling: {len(cells):,}.",
        f"- Strict local full-transcriptome or full-DG-subset datasets: {strict_full}.",
        f"- Selected-feature bridge datasets: {selected}.",
        "",
        "## Dataset Notes",
        "",
    ]
    for _, row in dataset_summary.iterrows():
        order_text = ""
        if pd.notna(row.get("rho_vs_known_order")):
            order_text = f"; rho vs known/stage order {fmt(row['rho_vs_known_order'])}"
        source_text = ""
        if pd.notna(row.get("rho_vs_source_pseudotime")):
            source_text = f"; rho vs source pseudotime {fmt(row['rho_vs_source_pseudotime'])}"
        lines.append(
            f"- `{row['dataset']}`: {int(row['n_cells']):,} cells, {int(row['n_total_genes_available']):,} genes available, "
            f"{int(row['n_hvg_used']):,} HVGs, root `{row['root_strategy']}` ({row['trajectory_evidence_grade']}){order_text}{source_text}."
        )
    lines.extend(["", "## Main Interpretation", ""])
    lines.append(
        "The diffusion layer does not overturn the existing Fig1-5 structure. It makes the mechanism more precise: the shared construction/configuration signal is stage-windowed, and age labels alone should not be treated as pseudotime."
    )
    lines.append(
        "The clearest manuscript-level change is Figure 5: pseudotime should become an explicit layer between regional fate/niche input and final granule-cell morphology. Figure 4 can keep the configuration model, but should describe it as maturation-window dependent."
    )
    lines.extend(["", "## Figure Impact", ""])
    for _, row in fig_impact.iterrows():
        lines.append(f"- `{row['figure']}`: {row['trajectory_effect']}. {row['recommended_action']}")
    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- This is not Scanpy/Monocle, because those packages are not installed locally; it is an equivalent lightweight diffusion/kNN trajectory layer.",
            "- Pseudotime direction is root-dependent. Curated developmental roots are stronger than marker-inferred roots.",
            "- `GSE325391` and `GSE268609` should not be called strict full-transcriptome trajectories until their complete source matrices are converted into trajectory-ready sparse objects.",
            "- RNA trajectories still do not directly measure BrdU, p21/p27 protein, pERK, pSMAD, or secreted protein bioactivity.",
            "",
            "## Outputs",
            "",
            f"- Cell scores: `{rel(OUT_CELL)}` ({len(cells):,} rows).",
            f"- Dataset summary: `{rel(OUT_DATASET_SUMMARY)}` ({len(dataset_summary):,} rows).",
            f"- Group summary: `{rel(OUT_GROUP_SUMMARY)}`.",
            f"- Module correlations: `{rel(OUT_MODULE_CORR)}`.",
            f"- Fig1-5 impact table: `{rel(OUT_FIG_IMPACT)}`.",
            f"- Overview plot: `{rel(OUT_PLOT)}`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    calls = load_calls()
    all_cells: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []
    for dataset, obs, loader, kwargs in dataset_jobs(calls):
        print(f"Processing {dataset}: {len(obs):,} starting cells", flush=True)
        cells, summary = run_dataset(dataset, obs, loader, **kwargs)
        all_cells.append(cells)
        summaries.append(summary)
        print(
            f"  -> {len(cells):,} cells; {summary['n_total_genes_available']:,} genes; "
            f"{summary['n_hvg_used']:,} HVGs; root={summary['root_strategy']}",
            flush=True,
        )
        del cells
        gc.collect()

    cells = pd.concat(all_cells, ignore_index=True, sort=False)
    dataset_summary = pd.DataFrame(summaries)
    group_summary = build_group_summary(cells)
    module_corr = build_module_correlations(cells)
    fig_impact = build_fig_impact(dataset_summary, module_corr)

    cells.to_csv(OUT_CELL, sep="\t", index=False, compression="gzip")
    dataset_summary.to_csv(OUT_DATASET_SUMMARY, sep="\t", index=False)
    group_summary.to_csv(OUT_GROUP_SUMMARY, sep="\t", index=False)
    module_corr.to_csv(OUT_MODULE_CORR, sep="\t", index=False)
    fig_impact.to_csv(OUT_FIG_IMPACT, sep="\t", index=False)
    plot_overview(dataset_summary, module_corr)
    write_report(cells, dataset_summary, module_corr, fig_impact)
    print(f"Wrote {rel(OUT_MD)}")
    print(dataset_summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
