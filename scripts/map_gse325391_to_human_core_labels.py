#!/usr/bin/env python3
"""Map GSE325391 adult dentate nuclei into the tuned human-core label convention."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from tune_human_core_labels_and_test_modules import compute_normalized_panel_scores


ROOT = Path(__file__).resolve().parents[2]
BRIDGE = ROOT / "Project/processed/gse325391_adult_dg_selected"
RESULTS = ROOT / "Project/results"

X_IN = BRIDGE / "matrix_cells_by_selected_genes.npz"
OBS_IN = BRIDGE / "cell_metadata.tsv.gz"
VAR_IN = BRIDGE / "var_selected_features.tsv"
UMAP_IN = BRIDGE / "harmony_umap.tsv.gz"

OUT_LABELS = RESULTS / "gse325391_human_core_label_projection.tsv.gz"
OUT_LABEL_SUMMARY = RESULTS / "gse325391_human_core_label_projection_summary.tsv"
OUT_REPLICATE_SUMMARY = RESULTS / "gse325391_human_core_replicate_module_summary.tsv"
OUT_TESTS = RESULTS / "gse325391_human_core_module_tests.tsv"
OUT_HEATMAP = RESULTS / "gse325391_human_core_module_heatmap.png"
OUT_DELTA_PLOT = RESULTS / "gse325391_human_core_diffn_vs_matn_deltas.png"
OUT_UMAP = RESULTS / "gse325391_human_core_label_umap.png"
OUT_MD = RESULTS / "gse325391_human_core_label_projection.md"

MODULES = [
    "norm_dentate_identity",
    "norm_human_dg_immature",
    "norm_human_dg_mature",
    "norm_cerebellar_identity",
    "norm_shared_granule_neuronal",
    "norm_morphogenesis_cytoskeleton",
    "norm_axon_guidance_synapse",
    "norm_background_max",
]


def normalize_counts(X: sparse.csr_matrix, obs: pd.DataFrame) -> sparse.csr_matrix:
    n_counts = pd.to_numeric(obs["nCount_RNA"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
    scale = np.divide(10000.0, n_counts, out=np.zeros_like(n_counts), where=n_counts > 0)
    Xn = X.astype(np.float32).multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data).astype(np.float32)
    return Xn


def add_ranks(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    score_cols = [col for col in out.columns if col.startswith("norm_") and not col.startswith("norm_detected_")]
    score_cols = [col for col in score_cols if col != "norm_background_dominant_panel"]
    for col in score_cols:
        out[f"{col}_rank"] = out[col].rank(method="average", pct=True).astype(np.float32)
    return out


def source_anchor_label(row: pd.Series) -> tuple[str, str, str]:
    if str(row.get("scDblFinder_class", "")).lower() == "doublet":
        return (
            "adult_dg_doublet_flag",
            "low",
            "scDblFinder marks this nucleus as a doublet",
        )

    background_dominates = float(row["norm_background_max"]) > max(
        float(row["norm_shared_granule_neuronal"]),
        float(row["norm_human_dg_mature"]),
        float(row["norm_human_dg_immature"]),
    )
    if background_dominates and float(row["norm_background_max_rank"]) >= 0.98:
        return (
            "adult_dg_background_warning",
            "low",
            f"background-dominant panel: {row['norm_background_dominant_panel']}",
        )

    cell_type = str(row.get("cell_type", ""))
    sub_cell_type = str(row.get("sub_cell_type", ""))
    combined = f"{cell_type};{sub_cell_type}".lower()
    if "diffn" in combined:
        return (
            "adult_human_dg_differentiating_anchor",
            "high",
            "source adult DG annotation marks this nucleus as differentiating granule-like",
        )
    if "matn" in combined:
        return (
            "adult_human_dg_mature_anchor",
            "high",
            "source adult DG annotation marks this nucleus as mature granule-like",
        )
    return (
        "adult_human_dg_anchor_unresolved_state",
        "medium",
        "source broad cell type is ExN_GC but maturation state is not parsed",
    )


def convention_projection(row: pd.Series) -> str:
    label = row["source_anchor_label"]
    if label == "adult_dg_doublet_flag":
        return "low_information_or_low_qc"
    if label == "adult_dg_background_warning":
        return "non_neuronal_background"
    if label == "adult_human_dg_differentiating_anchor":
        return "immature_neurogenic_candidate"
    if label == "adult_human_dg_mature_anchor":
        return "curated_human_dg_gc_anchor"
    return "human_dg_like_candidate"


def summarize_labels(obs: pd.DataFrame) -> pd.DataFrame:
    module_aggs = {f"median_{module.replace('norm_', '')}": (module, "median") for module in MODULES}
    return (
        obs.groupby(["group", "cell_type", "sub_cell_type", "source_anchor_label", "human_core_convention_label"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            n_samples=("sample", "nunique"),
            median_counts=("nCount_RNA", "median"),
            median_genes=("nFeature_RNA", "median"),
            median_percent_mito=("percent_mito", "median"),
            median_structural=("norm_structural_program_mean", "median"),
            **module_aggs,
        )
        .reset_index()
        .sort_values(["group", "source_anchor_label", "n_cells"], ascending=[True, True, False])
    )


def replicate_summary(obs: pd.DataFrame) -> pd.DataFrame:
    return (
        obs.groupby(["sample", "group", "cell_type", "source_anchor_label"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            median_counts=("nCount_RNA", "median"),
            **{f"median_{module.replace('norm_', '')}": (module, "median") for module in MODULES},
        )
        .reset_index()
    )


def paired_tests(rep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    a_label = "adult_human_dg_differentiating_anchor"
    b_label = "adult_human_dg_mature_anchor"
    for module in MODULES:
        col = f"median_{module.replace('norm_', '')}"
        pivot = rep.pivot_table(index="sample", columns="source_anchor_label", values=col, aggfunc="median")
        if a_label not in pivot or b_label not in pivot:
            continue
        paired = pivot[[a_label, b_label]].dropna()
        if len(paired) < 3:
            continue
        diff = paired[a_label] - paired[b_label]
        try:
            p_value = stats.wilcoxon(diff).pvalue
        except ValueError:
            p_value = np.nan
        rows.append(
            {
                "comparison": "differentiating_vs_mature_within_sample",
                "module": module,
                "n_sample_units": len(paired),
                "median_delta": float(np.median(diff)),
                "mean_delta": float(np.mean(diff)),
                "p_value": p_value,
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        valid = out["p_value"].notna()
        ranked = out.loc[valid, "p_value"].rank(method="first").to_numpy()
        m = int(valid.sum())
        adjusted = np.minimum.accumulate((out.loc[valid, "p_value"].to_numpy()[np.argsort(ranked)] * m / np.arange(1, m + 1))[::-1])[::-1]
        # Simpler stable BH assignment by sorted p values.
        valid_idx = out.index[valid].to_numpy()
        order = out.loc[valid_idx, "p_value"].to_numpy().argsort()
        p_sorted = out.loc[valid_idx[order], "p_value"].to_numpy()
        bh = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
        p_adj = pd.Series(np.nan, index=out.index, dtype=float)
        p_adj.loc[valid_idx[order]] = np.minimum(bh, 1.0)
        out["p_adj_bh"] = p_adj
    else:
        out["p_adj_bh"] = []
    return out


def plot_heatmap(summary: pd.DataFrame) -> None:
    heat = (
        summary.groupby(["source_anchor_label", "cell_type"], dropna=False)
        .agg(**{module: (module.replace("norm_", "median_"), "median") for module in MODULES if module.replace("norm_", "median_") in summary.columns})
        .reset_index()
    )
    if heat.empty:
        return
    heat["row_label"] = heat["source_anchor_label"] + " | " + heat["cell_type"].astype(str)
    matrix = heat[[module for module in MODULES if module in heat.columns]].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(10.5, max(3.5, 0.35 * len(heat))))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels(heat["row_label"], fontsize=8)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels([col.replace("norm_", "") for col in MODULES if col in heat.columns], rotation=35, ha="right", fontsize=8)
    ax.set_title("GSE325391 module medians by projected label")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("median normalized module score")
    fig.tight_layout()
    fig.savefig(OUT_HEATMAP, dpi=180)
    plt.close(fig)


def plot_deltas(tests: pd.DataFrame) -> None:
    if tests.empty:
        return
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    labels = tests["module"].str.replace("norm_", "", regex=False)
    colors = ["#1b7f6b" if value >= 0 else "#a74d2a" for value in tests["median_delta"]]
    ax.bar(np.arange(len(tests)), tests["median_delta"], color=colors)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_xticks(np.arange(len(tests)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("median paired delta")
    ax.set_title("Differentiating minus mature DG anchors by sample")
    fig.tight_layout()
    fig.savefig(OUT_DELTA_PLOT, dpi=180)
    plt.close(fig)


def plot_umap(obs: pd.DataFrame) -> None:
    if not UMAP_IN.exists():
        return
    umap = pd.read_csv(UMAP_IN, sep="\t")
    merged = umap.merge(obs[["cell_id", "source_anchor_label", "cell_type", "group"]], on="cell_id", how="inner")
    coord_cols = [col for col in merged.columns if col not in {"cell_id", "cell_name", "source_anchor_label", "cell_type", "group"}]
    if len(coord_cols) < 2:
        return
    rng = np.random.default_rng(13)
    if len(merged) > 40000:
        merged = merged.iloc[rng.choice(len(merged), size=40000, replace=False)].copy()
    labels = sorted(merged["source_anchor_label"].unique())
    palette = {
        "adult_human_dg_mature_anchor": "#2f6fbb",
        "adult_human_dg_differentiating_anchor": "#c26b2e",
        "adult_dg_doublet_flag": "#6f6f6f",
        "adult_dg_background_warning": "#9b2f3d",
        "adult_human_dg_anchor_unresolved_state": "#6d4c9f",
    }
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for label in labels:
        sub = merged[merged["source_anchor_label"] == label]
        ax.scatter(sub[coord_cols[0]], sub[coord_cols[1]], s=2, alpha=0.45, c=palette.get(label, "#777777"), label=label)
    ax.set_xlabel(coord_cols[0])
    ax.set_ylabel(coord_cols[1])
    ax.set_title("GSE325391 harmony UMAP label projection")
    ax.legend(markerscale=5, fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(OUT_UMAP, dpi=180)
    plt.close(fig)


def write_md(obs: pd.DataFrame, tests: pd.DataFrame) -> None:
    counts = obs["source_anchor_label"].value_counts()
    convention_counts = obs["human_core_convention_label"].value_counts()
    lines = [
        "# GSE325391 Human-Core Label Projection",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This pass maps the adult human dentate `GSE325391` Seurat object into the current human-core selected feature space and tuned label convention.",
        "",
        "## Source-Aware Labels",
        "",
    ]
    for label, count in counts.items():
        lines.append(f"- `{label}`: {int(count)} cells.")
    lines.extend(["", "## Human-Core Convention Projection", ""])
    for label, count in convention_counts.items():
        lines.append(f"- `{label}`: {int(count)} cells.")
    lines.extend(["", "## Differentiating Versus Mature DG Module Tests", ""])
    if tests.empty:
        lines.append("No paired tests were run.")
    else:
        for _, row in tests.sort_values("p_adj_bh").iterrows():
            lines.append(
                f"- `{row['module']}`: median differentiating-minus-mature delta {row['median_delta']:.4f} "
                f"across {int(row['n_sample_units'])} samples (BH-adjusted p={row['p_adj_bh']:.3g})."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `GSE325391` is now a direct adult human dentate granule-cell anchor rather than only a candidate resource.",
            "- Source mature and differentiating labels are preserved; the projected convention allows comparison with `GSE186538`, `GSE185277`, and `GSE185553` without erasing the stronger source annotation.",
            "",
            "## Outputs",
            "",
            f"- Cell labels and scores: `{OUT_LABELS.relative_to(ROOT)}`",
            f"- Label summary: `{OUT_LABEL_SUMMARY.relative_to(ROOT)}`",
            f"- Replicate summary: `{OUT_REPLICATE_SUMMARY.relative_to(ROOT)}`",
            f"- Module tests: `{OUT_TESTS.relative_to(ROOT)}`",
            f"- Heatmap: `{OUT_HEATMAP.relative_to(ROOT)}`",
            f"- Delta plot: `{OUT_DELTA_PLOT.relative_to(ROOT)}`",
            f"- UMAP plot: `{OUT_UMAP.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    print("Loading selected GSE325391 object", flush=True)
    X = sparse.load_npz(X_IN).tocsr()
    obs = pd.read_csv(OBS_IN, sep="\t", low_memory=False)
    var = pd.read_csv(VAR_IN, sep="\t")
    if X.shape != (len(obs), len(var)):
        raise ValueError(f"Shape mismatch: X {X.shape}, obs {len(obs)}, var {len(var)}")

    print("Normalizing and scoring panels", flush=True)
    Xn = normalize_counts(X, obs)
    scores = compute_normalized_panel_scores(Xn, var)
    scores = add_ranks(scores)
    obs = pd.concat([obs.reset_index(drop=True), scores.reset_index(drop=True)], axis=1)

    labels = obs.apply(source_anchor_label, axis=1, result_type="expand")
    labels.columns = ["source_anchor_label", "source_label_confidence", "source_label_reason"]
    obs = pd.concat([obs, labels], axis=1)
    obs["human_core_convention_label"] = obs.apply(convention_projection, axis=1)
    obs["analysis_include"] = (
        obs["scDblFinder_class"].astype(str).str.lower().eq("singlet")
        & (pd.to_numeric(obs["nFeature_RNA"], errors="coerce") >= 300)
        & (pd.to_numeric(obs["nCount_RNA"], errors="coerce") >= 500)
        & (pd.to_numeric(obs["percent_mito"], errors="coerce") <= 20)
    )

    print("Summarizing and plotting", flush=True)
    summary = summarize_labels(obs)
    rep = replicate_summary(obs[obs["analysis_include"]].copy())
    tests = paired_tests(rep)

    keep_cols = [
        "cell_id",
        "cell_name",
        "sample",
        "group",
        "Run",
        "scDblFinder_class",
        "broad_cell_type",
        "cell_type",
        "sub_cell_type",
        "nCount_RNA",
        "nFeature_RNA",
        "percent_mito",
        "slingAvgPseudotime",
        "source_anchor_label",
        "source_label_confidence",
        "source_label_reason",
        "human_core_convention_label",
        "analysis_include",
    ]
    score_cols = [col for col in obs.columns if col.startswith("norm_")]
    obs[keep_cols + score_cols].to_csv(OUT_LABELS, sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT_LABEL_SUMMARY, sep="\t", index=False)
    rep.to_csv(OUT_REPLICATE_SUMMARY, sep="\t", index=False)
    tests.to_csv(OUT_TESTS, sep="\t", index=False)
    plot_heatmap(summary)
    plot_deltas(tests)
    plot_umap(obs)
    write_md(obs, tests)

    print(f"Wrote {OUT_LABELS}")
    print(f"Wrote {OUT_LABEL_SUMMARY}")
    print(f"Wrote {OUT_REPLICATE_SUMMARY}")
    print(f"Wrote {OUT_TESTS}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
