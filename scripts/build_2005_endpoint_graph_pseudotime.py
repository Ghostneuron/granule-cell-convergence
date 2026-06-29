#!/usr/bin/env python3
"""Lightweight cell-level graph pseudotime for 2005 paper endpoint genes.

This script builds a dependency-light pseudotime layer without Scanpy:
selected paper-endpoint genes -> log1p -> z-score -> PCA -> kNN graph ->
shortest-path distance from biologically defined root cells.

It is intentionally a first-pass trajectory audit, not a final full-transcriptome
trajectory model.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.csgraph import dijkstra
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_2005_paper_endpoint_pseudotime_audit import (  # noqa: E402
    MODULE_COLORS,
    PAPER_ENDPOINT_MODULES,
    PLOT_MODULES,
    module_gene_table,
)


GSE104323_EXPR = EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz"
GSE104323_META = EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz"
GSE292261_EXPR = EXTERNAL / "GEO/GSE292261/GSE292261_counts_SS2_filtered_raw.csv.gz"
GSE292261_META = EXTERNAL / "GEO/GSE292261/GSE292261_sample_data_SS2_filtered.csv.gz"
GSE122357_TAR = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
REFINED_CALLS = RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"

OUT_CELL_SCORES = RESULTS / "primary_core_2005_endpoint_graph_pseudotime_cell_scores.tsv.gz"
OUT_GROUP_SUMMARY = RESULTS / "primary_core_2005_endpoint_graph_pseudotime_group_summary.tsv"
OUT_BIN_SUMMARY = RESULTS / "primary_core_2005_endpoint_graph_pseudotime_bin_summary.tsv"
OUT_CORRELATIONS = RESULTS / "primary_core_2005_endpoint_graph_pseudotime_module_correlations.tsv"
OUT_PLOT = RESULTS / "primary_core_2005_endpoint_graph_pseudotime.png"
OUT_MD = RESULTS / "primary_core_2005_endpoint_graph_pseudotime.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


GSE104323_ORDER = {
    "RGL_young": 0,
    "RGL": 1,
    "nIPC": 2,
    "nIPC-perin": 3,
    "Neuroblast": 4,
    "Immature-GC": 5,
    "GC-juv": 6,
    "GC-adult": 7,
}

GSE292261_ORDER = {f"P{age}": float(age) for age in [5, 7, 10, 15, 28]}
GSE122357_ORDER = {"P0": 0.0, "P8a": 8.1, "P8b": 8.2}
GSE122357_MEMBERS = {
    "GSM3464549_P0": ("GSM3464549_P0.csv.gz", "P0"),
    "GSM3464550_P8a": ("GSM3464550_P8a.csv.gz", "P8a"),
    "GSM3464551_P8b": ("GSM3464551_P8b.csv.gz", "P8b"),
}


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').strip("'").upper()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def target_genes() -> set[str]:
    genes = set(module_gene_table()["canonical_gene"])
    # Add a few granule-lineage anchors that stabilize graph geometry.
    genes.update({"PROX1", "ZBTB20", "CALB2", "GABRA6", "PAX6", "NEUROD6"})
    return genes


def load_candidate_calls(dataset: str, candidate_call: str) -> pd.DataFrame:
    cols = ["dataset", "sample", "cell_id", "group", "candidate_call"]
    calls = pd.read_csv(REFINED_CALLS, sep="\t", usecols=cols, dtype=str, low_memory=False)
    return calls.loc[calls["dataset"].eq(dataset) & calls["candidate_call"].eq(candidate_call)].copy()


def matrix_from_gene_rows(
    *,
    path: Path | None,
    target: set[str],
    selected_cells: dict[str, dict[str, object]],
    delimiter: str,
    tar_member: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if tar_member is None:
        fh_cm = open_text(path)  # type: ignore[arg-type]
        close_tf = None
    else:
        close_tf = tarfile.open(path)  # type: ignore[arg-type]
        raw = close_tf.extractfile(tar_member)
        if raw is None:
            raise FileNotFoundError(tar_member)
        fh_cm = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")

    gene_values: dict[str, np.ndarray] = {}
    obs: pd.DataFrame
    try:
        reader = csv.reader(fh_cm, delimiter=delimiter)
        header = next(reader)
        cell_ids = [str(cell).strip().strip('"') for cell in header[1:]]
        selected_positions = [(idx + 1, cell) for idx, cell in enumerate(cell_ids) if cell in selected_cells]
        obs_rows = []
        for _, cell in selected_positions:
            obs_rows.append({"cell_id": cell, **selected_cells[cell]})
        obs = pd.DataFrame(obs_rows)

        for row in reader:
            if not row:
                continue
            gene = canon(row[0])
            if gene not in target:
                continue
            vals = np.zeros(len(selected_positions), dtype=float)
            for out_idx, (col_idx, _) in enumerate(selected_positions):
                try:
                    vals[out_idx] = float(row[col_idx])
                except (IndexError, ValueError):
                    vals[out_idx] = 0.0
            gene_values[gene] = vals
    finally:
        fh_cm.close()
        if close_tf is not None:
            close_tf.close()

    genes = sorted(gene_values)
    x = pd.DataFrame({gene: gene_values[gene] for gene in genes})
    return obs.reset_index(drop=True), x.reset_index(drop=True)


def load_gse104323_cells(target: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(GSE104323_META, sep="\t", dtype=str)
    cell_col = "Sample name (24185 single cells)"
    cluster_col = "characteristics: cell cluster"
    meta = meta.loc[meta[cluster_col].isin(GSE104323_ORDER)].copy()
    selected: dict[str, dict[str, object]] = {}
    for _, row in meta.iterrows():
        group = str(row[cluster_col])
        cell = str(row[cell_col])
        selected[cell] = {
            "dataset": "GSE104323",
            "axis_type": "adult_dentate_lineage_state",
            "axis_label": group,
            "axis_order": GSE104323_ORDER[group],
            "comparison_group": "curated_dentate_granule_lineage",
            "root_cell": group in {"RGL_young", "RGL"},
        }
    return matrix_from_gene_rows(path=GSE104323_EXPR, target=target, selected_cells=selected, delimiter="\t")


def load_gse292261_cells(target: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    calls = load_candidate_calls("GSE292261", "candidate_dentate_granule")
    candidate_cells = set(calls["cell_id"].astype(str))
    meta = pd.read_csv(GSE292261_META, dtype=str)
    first_col = meta.columns[0]
    meta = meta.rename(columns={first_col: "cell_id"})
    meta["axis_label"] = meta["Sample"].astype(str).str.replace("DG_", "", regex=False)
    meta = meta.loc[meta["cell_id"].isin(candidate_cells) & meta["axis_label"].isin(GSE292261_ORDER)].copy()

    header = pd.read_csv(GSE292261_EXPR, nrows=0).columns.tolist()
    gene_cols = [col for col in header[1:] if canon(col) in target]
    df = pd.read_csv(GSE292261_EXPR, usecols=[header[0], *gene_cols], low_memory=False).rename(columns={header[0]: "cell_id"})
    df["cell_id"] = df["cell_id"].astype(str)
    df = df.merge(meta[["cell_id", "axis_label"]], on="cell_id", how="inner")

    obs = pd.DataFrame(
        {
            "cell_id": df["cell_id"].astype(str),
            "dataset": "GSE292261",
            "axis_type": "postnatal_dentate_age",
            "axis_label": df["axis_label"],
            "axis_order": df["axis_label"].map(GSE292261_ORDER),
            "comparison_group": "candidate_dentate_granule_only",
            "root_cell": df["axis_label"].eq("P5"),
        }
    )
    x = df[gene_cols].copy()
    x.columns = [canon(col) for col in x.columns]
    x = x.T.groupby(level=0).sum().T
    return obs.reset_index(drop=True), x.reset_index(drop=True)


def load_gse122357_cells(target: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    calls = load_candidate_calls("GSE122357", "candidate_cerebellar_granule")
    selected_by_sample: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for _, row in calls.iterrows():
        sample = str(row["sample"])
        if sample not in GSE122357_MEMBERS:
            continue
        stage = GSE122357_MEMBERS[sample][1]
        cell = str(row["cell_id"])
        selected_by_sample[sample][cell] = {
            "dataset": "GSE122357",
            "axis_type": "postnatal_cerebellar_age",
            "axis_label": stage,
            "axis_order": GSE122357_ORDER[stage],
            "comparison_group": "candidate_cerebellar_granule_only",
            "root_cell": stage == "P0",
        }

    obs_pieces: list[pd.DataFrame] = []
    x_pieces: list[pd.DataFrame] = []
    for sample, (member, _) in GSE122357_MEMBERS.items():
        obs, x = matrix_from_gene_rows(
            path=GSE122357_TAR,
            target=target,
            selected_cells=selected_by_sample[sample],
            delimiter=",",
            tar_member=member,
        )
        obs["cell_id"] = sample + ":" + obs["cell_id"].astype(str)
        obs_pieces.append(obs)
        x_pieces.append(x)
    all_obs = pd.concat(obs_pieces, ignore_index=True)
    all_x = pd.concat(x_pieces, ignore_index=True).fillna(0.0)
    return all_obs.reset_index(drop=True), all_x.reset_index(drop=True)


def align_matrix(x: pd.DataFrame) -> pd.DataFrame:
    x = x.copy()
    x.columns = [canon(col) for col in x.columns]
    x = x.T.groupby(level=0).sum().T
    keep = [col for col in x.columns if pd.to_numeric(x[col], errors="coerce").fillna(0).sum() > 0]
    return x[sorted(keep)].apply(pd.to_numeric, errors="coerce").fillna(0.0)


def graph_pseudotime(x: pd.DataFrame, root_mask: pd.Series, n_neighbors: int = 15) -> tuple[np.ndarray, np.ndarray]:
    x = align_matrix(x)
    if x.shape[0] < 5 or x.shape[1] < 3:
        raise RuntimeError(f"Too few cells/genes for graph pseudotime: {x.shape}")
    logx = np.log1p(x.to_numpy(dtype=float))
    std = logx.std(axis=0)
    keep = std > 1e-8
    logx = logx[:, keep]
    z = (logx - logx.mean(axis=0)) / logx.std(axis=0)
    n_components = max(2, min(12, z.shape[0] - 1, z.shape[1]))
    emb = PCA(n_components=n_components, random_state=1).fit_transform(z)
    k = max(2, min(n_neighbors, emb.shape[0] - 1))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    distances, indices = nn.kneighbors(emb)
    rows = np.repeat(np.arange(emb.shape[0]), k)
    cols = indices[:, 1:].reshape(-1)
    vals = distances[:, 1:].reshape(-1)
    adj = sparse.coo_matrix((vals, (rows, cols)), shape=(emb.shape[0], emb.shape[0])).tocsr()
    adj = (adj + adj.T) * 0.5

    roots = np.where(root_mask.to_numpy(dtype=bool))[0]
    if roots.size == 0:
        raise RuntimeError("No root cells found for pseudotime.")
    n = adj.shape[0]
    root_node = n
    extra_rows = np.concatenate([np.full(roots.size, root_node), roots])
    extra_cols = np.concatenate([roots, np.full(roots.size, root_node)])
    extra_vals = np.full(roots.size * 2, 1e-6, dtype=float)
    graph = sparse.vstack(
        [
            sparse.hstack([adj, sparse.csr_matrix((n, 1))]),
            sparse.csr_matrix((1, n + 1)),
        ]
    ).tolil()
    for r, c, v in zip(extra_rows, extra_cols, extra_vals):
        graph[r, c] = v
    dist = dijkstra(graph.tocsr(), directed=False, indices=root_node)[:n]
    finite = np.isfinite(dist)
    if finite.sum() == 0:
        pseudotime = np.zeros(n)
    else:
        maxv = dist[finite].max()
        pseudotime = np.zeros(n)
        pseudotime[finite] = dist[finite] / maxv if maxv > 0 else 0.0
        pseudotime[~finite] = 1.0
    return pseudotime, emb[:, :2]


def add_module_scores(obs: pd.DataFrame, x: pd.DataFrame) -> pd.DataFrame:
    x = align_matrix(x)
    logx = np.log1p(x)
    z = (logx - logx.mean(axis=0)) / logx.std(axis=0).replace(0, np.nan)
    z = z.fillna(0.0)
    out = obs.copy()
    for module in PAPER_ENDPOINT_MODULES:
        genes = [canon(gene) for gene in module["genes"] if canon(gene) in z.columns]
        if genes:
            raw = z[genes].mean(axis=1)
            out[module["module_id"]] = raw.rank(method="average", pct=True)
            out[f"{module['module_id']}_n_genes"] = len(genes)
        else:
            out[module["module_id"]] = np.nan
            out[f"{module['module_id']}_n_genes"] = 0
    return out


def process_dataset(name: str, loader) -> pd.DataFrame:
    obs, x = loader(target_genes())
    x = align_matrix(x)
    pseudotime, emb = graph_pseudotime(x, obs["root_cell"])
    scores = add_module_scores(obs, x)
    scores["graph_pseudotime"] = pseudotime
    scores["pca1"] = emb[:, 0]
    scores["pca2"] = emb[:, 1]
    scores["n_genes_used"] = x.shape[1]
    return scores


def finite_median(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def build_group_summary(cells: pd.DataFrame) -> pd.DataFrame:
    module_ids = [module["module_id"] for module in PAPER_ENDPOINT_MODULES]
    agg = {
        "cell_id": "count",
        "graph_pseudotime": ["median", "mean"],
        "axis_order": "median",
        **{module_id: "median" for module_id in module_ids},
    }
    summary = cells.groupby(["dataset", "comparison_group", "axis_type", "axis_label"], sort=False).agg(agg)
    summary.columns = ["_".join(col).rstrip("_") if isinstance(col, tuple) else col for col in summary.columns]
    summary = summary.reset_index().rename(columns={"cell_id_count": "n_cells"})
    return summary.sort_values(["dataset", "comparison_group", "axis_order_median"])


def build_bin_summary(cells: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    module_ids = [module["module_id"] for module in PAPER_ENDPOINT_MODULES]
    for (dataset, group), sub in cells.groupby(["dataset", "comparison_group"], sort=False):
        sub = sub.copy()
        ranks = sub["graph_pseudotime"].rank(method="first")
        try:
            sub["pseudotime_bin"] = pd.qcut(ranks, q=min(n_bins, len(sub)), labels=False, duplicates="drop")
        except ValueError:
            continue
        for bin_id, bsub in sub.groupby("pseudotime_bin", sort=True):
            row = {
                "dataset": dataset,
                "comparison_group": group,
                "pseudotime_bin": int(bin_id),
                "n_cells": int(len(bsub)),
                "median_graph_pseudotime": finite_median(bsub["graph_pseudotime"]),
                "dominant_axis_label": bsub["axis_label"].value_counts().idxmax(),
            }
            for module_id in module_ids:
                row[module_id] = finite_median(bsub[module_id])
            rows.append(row)
    return pd.DataFrame(rows)


def build_correlations(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (dataset, group), sub in cells.groupby(["dataset", "comparison_group"], sort=False):
        for module in PAPER_ENDPOINT_MODULES:
            module_id = module["module_id"]
            xy = sub[["graph_pseudotime", module_id]].dropna()
            if len(xy) < 10 or xy["graph_pseudotime"].nunique() < 3 or xy[module_id].nunique() < 3:
                rho = p = np.nan
            else:
                rho, p = stats.spearmanr(xy["graph_pseudotime"], xy[module_id])
            rows.append(
                {
                    "dataset": dataset,
                    "comparison_group": group,
                    "module_id": module_id,
                    "module_label": module["module_label"],
                    "paper_readout": module["paper_readout"],
                    "trajectory_role": module["trajectory_role"],
                    "spearman_rho_vs_graph_pseudotime": float(rho) if np.isfinite(rho) else np.nan,
                    "spearman_p": float(p) if np.isfinite(p) else np.nan,
                    "n_cells": int(len(xy)),
                }
            )
    return pd.DataFrame(rows).sort_values(["dataset", "comparison_group", "module_id"])


def plot_bins(bin_summary: pd.DataFrame) -> None:
    panels = [
        ("GSE104323", "curated_dentate_granule_lineage", "GSE104323 dentate lineage"),
        ("GSE292261", "candidate_dentate_granule_only", "GSE292261 postnatal dentate candidates"),
        ("GSE122357", "candidate_cerebellar_granule_only", "GSE122357 cerebellar candidates"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
    for ax, (dataset, group, title) in zip(axes, panels):
        sub = bin_summary.loc[bin_summary["dataset"].eq(dataset) & bin_summary["comparison_group"].eq(group)].copy()
        if sub.empty:
            ax.axis("off")
            continue
        for module_id in PLOT_MODULES:
            if module_id not in sub:
                continue
            ax.plot(
                sub["median_graph_pseudotime"],
                sub[module_id],
                marker="o",
                linewidth=1.8,
                markersize=3.5,
                color=MODULE_COLORS.get(module_id, "#333333"),
                label=module_id,
            )
        ax.set_title(title)
        ax.set_xlabel("Graph pseudotime")
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    axes[0].set_ylabel("Median module score")
    handles, labels = axes[0].get_legend_handles_labels()
    label_lookup = {module["module_id"]: module["module_label"] for module in PAPER_ENDPOINT_MODULES}
    labels = [label_lookup.get(label, label) for label in labels]
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Cell-level graph pseudotime for 2005 endpoint modules", fontsize=14)
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


def corr_text(corr: pd.DataFrame, dataset: str, module_id: str) -> str:
    sub = corr.loc[corr["dataset"].eq(dataset) & corr["module_id"].eq(module_id)]
    if sub.empty:
        return "NA"
    row = sub.iloc[0]
    return f"rho {fmt(row['spearman_rho_vs_graph_pseudotime'])}, p {fmt(row['spearman_p'])}"


def write_markdown(cells: pd.DataFrame, group_summary: pd.DataFrame, correlations: pd.DataFrame) -> None:
    lines = [
        "# 2005 Endpoint Graph Pseudotime",
        "",
        "## Purpose",
        "",
        "This is the first cell-level trajectory layer for the 2005 paper endpoints. It uses paper-endpoint genes, PCA, kNN graph construction, and shortest-path distance from biologically defined root cells. It is not yet a full-transcriptome Scanpy/Monocle trajectory.",
        "",
        "## Roots",
        "",
        "- `GSE104323`: RGL_young/RGL dentate lineage cells.",
        "- `GSE292261`: P5 candidate dentate granule cells.",
        "- `GSE122357`: P0 candidate cerebellar granule cells.",
        "",
        "## Main Findings",
        "",
        f"- `GSE104323` is the cleanest trajectory: proliferation versus graph pseudotime is {corr_text(correlations, 'GSE104323', 'proliferation_brdU_proxy')}; neuronal maturation is {corr_text(correlations, 'GSE104323', 'neuronal_differentiation_maturation')}.",
        f"- `GSE292261` shows postnatal dentate timing but not a simple age-only line: proliferation is {corr_text(correlations, 'GSE292261', 'proliferation_brdU_proxy')}; TGF/SMAD is {corr_text(correlations, 'GSE292261', 'tgf_smad_pai1_response')}.",
        f"- `GSE122357` shows a cerebellar candidate trajectory from P0 roots, but P8a/P8b differ: TGF/SMAD is {corr_text(correlations, 'GSE122357', 'tgf_smad_pai1_response')}; secreted stop-factor axis is {corr_text(correlations, 'GSE122357', 'secreted_stop_candidate_axis')}.",
        "",
        "## Interpretation",
        "",
        "The graph layer supports the user's concern: the 2005 readouts are stage/trajectory dependent. The cleanest evidence is the adult dentate lineage, where neuronal maturation rises along graph pseudotime and proliferation behaves as a trajectory-windowed state rather than a simple late-versus-early marker. The postnatal dentate and cerebellar datasets show trajectory windows rather than simple monotonic age effects, so future manuscript claims should avoid treating age as pseudotime.",
        "",
        "The next stronger analysis should use full-transcriptome highly variable genes and a dedicated trajectory package or an equivalent diffusion/PAGA workflow, then overlay the 2005 endpoint modules, TGF-beta/BDNF, and conditioned-medium secretome candidates.",
        "",
        "## Outputs",
        "",
        f"- Cell scores: `{OUT_CELL_SCORES.relative_to(ROOT)}` ({len(cells):,} cells).",
        f"- Group summary: `{OUT_GROUP_SUMMARY.relative_to(ROOT)}` ({len(group_summary):,} rows).",
        f"- Module correlations: `{OUT_CORRELATIONS.relative_to(ROOT)}` ({len(correlations):,} rows).",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    pieces = [
        process_dataset("GSE104323", load_gse104323_cells),
        process_dataset("GSE292261", load_gse292261_cells),
        process_dataset("GSE122357", load_gse122357_cells),
    ]
    cells = pd.concat(pieces, ignore_index=True)
    group_summary = build_group_summary(cells)
    bin_summary = build_bin_summary(cells)
    correlations = build_correlations(cells)

    cells.to_csv(OUT_CELL_SCORES, sep="\t", index=False, compression="gzip")
    group_summary.to_csv(OUT_GROUP_SUMMARY, sep="\t", index=False)
    bin_summary.to_csv(OUT_BIN_SUMMARY, sep="\t", index=False)
    correlations.to_csv(OUT_CORRELATIONS, sep="\t", index=False)
    plot_bins(bin_summary)
    write_markdown(cells, group_summary, correlations)

    print(f"Wrote {OUT_MD}")
    print(f"Cells: {len(cells):,}")
    print(f"Group rows: {len(group_summary):,}")
    print(f"Correlation rows: {len(correlations):,}")
    print(
        correlations[
            [
                "dataset",
                "module_id",
                "spearman_rho_vs_graph_pseudotime",
                "spearman_p",
                "n_cells",
            ]
        ]
        .round(3)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
