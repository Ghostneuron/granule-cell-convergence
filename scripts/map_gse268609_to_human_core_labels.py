#!/usr/bin/env python3
"""Project GSE268609 RNA nuclei into the human-core marker/module convention."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from validate_human_core_marker_programs import BACKGROUND_PANELS, load_panels


ROOT = Path(__file__).resolve().parents[2]
BRIDGE = ROOT / "Project/processed/gse268609_rna_selected"
RESULTS = ROOT / "Project/results"

X_IN = BRIDGE / "matrix_cells_by_selected_genes.npz"
OBS_IN = BRIDGE / "cell_metadata.tsv.gz"
VAR_IN = BRIDGE / "var_selected_features.tsv"

OUT_LABELS = RESULTS / "gse268609_human_core_label_projection.tsv.gz"
OUT_LABEL_SUMMARY = RESULTS / "gse268609_human_core_label_projection_summary.tsv"
OUT_SAMPLE_SUMMARY = RESULTS / "gse268609_human_core_sample_module_summary.tsv"
OUT_TESTS = RESULTS / "gse268609_human_core_diagnosis_module_tests.tsv"
OUT_HEATMAP = RESULTS / "gse268609_human_core_label_module_heatmap.png"
OUT_DELTA_PLOT = RESULTS / "gse268609_human_core_diagnosis_module_deltas.png"
OUT_MD = RESULTS / "gse268609_human_core_label_projection.md"

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

RANK_COLS = [
    "norm_dentate_identity",
    "norm_cerebellar_identity",
    "norm_human_dg_immature",
    "norm_human_dg_mature",
    "norm_shared_granule_neuronal",
    "norm_morphogenesis_cytoskeleton",
    "norm_axon_guidance_synapse",
    "norm_structural_program_mean",
    "norm_background_max",
]


def norm_gene(value: object) -> str:
    return str(value).strip().strip('"').strip("'").upper()


def percentile_rank(series: pd.Series) -> pd.Series:
    if len(series) <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def compute_panel_scores_from_counts(X: sparse.csr_matrix, obs: pd.DataFrame, var: pd.DataFrame) -> pd.DataFrame:
    panels = load_panels()
    gene_to_idx = {norm_gene(gene): idx for idx, gene in enumerate(var["gene"])}
    n_counts = pd.to_numeric(obs["nCount_RNA"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
    scale = np.divide(10000.0, n_counts, out=np.zeros_like(n_counts), where=n_counts > 0)

    score_cols: dict[str, np.ndarray] = {}
    detected_cols: dict[str, np.ndarray] = {}
    for panel, genes in panels.items():
        idxs = [gene_to_idx[norm_gene(gene)] for gene in genes if norm_gene(gene) in gene_to_idx]
        if idxs:
            sub = X[:, idxs].astype(np.float32, copy=True)
            sub = sub.multiply(scale[:, None]).tocsr()
            sub.data = np.log1p(sub.data).astype(np.float32)
            values = np.asarray(sub.mean(axis=1)).ravel().astype(np.float32)
            detected = np.asarray((X[:, idxs] > 0).sum(axis=1)).ravel().astype(np.int16)
        else:
            values = np.zeros(X.shape[0], dtype=np.float32)
            detected = np.zeros(X.shape[0], dtype=np.int16)
        score_cols[f"norm_{panel}"] = values
        detected_cols[f"norm_detected_{panel}"] = detected

    scores = pd.DataFrame({**score_cols, **detected_cols})
    background_cols = [f"norm_{panel}" for panel in BACKGROUND_PANELS if f"norm_{panel}" in scores.columns]
    scores["norm_background_max"] = scores[background_cols].max(axis=1)
    scores["norm_background_dominant_panel"] = scores[background_cols].idxmax(axis=1).str.replace("norm_", "", regex=False)
    scores["norm_structural_program_mean"] = scores[
        ["norm_shared_granule_neuronal", "norm_morphogenesis_cytoskeleton", "norm_axon_guidance_synapse"]
    ].mean(axis=1)
    scores["norm_identity_contrast_dentate_minus_cerebellar"] = scores["norm_dentate_identity"] - scores["norm_cerebellar_identity"]
    return scores


def add_ranks(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in RANK_COLS:
        if col in out:
            out[f"{col}_global_rank"] = percentile_rank(out[col]).astype(np.float32)
            out[f"{col}_sample_rank"] = out.groupby("sample_id", dropna=False)[col].transform(percentile_rank).astype(np.float32)
    return out


def projection_label(row: pd.Series) -> tuple[str, str, str]:
    if not bool(row.get("analysis_include_basic_qc", False)):
        return ("low_information_or_low_qc", "low", "below the source RNA-count/gene-count QC thresholds")

    background_dominates = float(row["norm_background_max"]) > max(
        float(row["norm_shared_granule_neuronal"]),
        float(row["norm_dentate_identity"]),
        float(row["norm_human_dg_immature"]),
        float(row["norm_human_dg_mature"]),
    )
    if background_dominates and float(row["norm_background_max_sample_rank"]) >= 0.85:
        return ("non_neuronal_background", "low", f"background-dominant panel: {row['norm_background_dominant_panel']}")

    dentate = float(row["norm_dentate_identity_sample_rank"])
    shared = float(row["norm_shared_granule_neuronal_sample_rank"])
    immature = float(row["norm_human_dg_immature_sample_rank"])
    mature = float(row["norm_human_dg_mature_sample_rank"])
    structural = float(row["norm_structural_program_mean_sample_rank"])
    background = float(row["norm_background_max_sample_rank"])
    contrast = float(row["norm_identity_contrast_dentate_minus_cerebellar"])

    if dentate >= 0.80 and shared >= 0.65 and immature >= 0.75 and background < 0.90:
        return (
            "immature_neurogenic_candidate",
            "medium",
            "dentate and shared neuronal signal with high immature/neurogenic module support",
        )

    if dentate >= 0.88 and shared >= 0.70 and contrast >= 0 and background < 0.90:
        state = "immature" if immature >= mature else "mature"
        return (
            "human_dg_like_high_confidence",
            "medium",
            f"{state}-shifted dentate identity with shared neuronal support; source cell type not yet available",
        )

    if dentate >= 0.72 and shared >= 0.55 and background < 0.95:
        return (
            "human_dg_like_candidate",
            "low",
            "marker/module projection supports DG-like identity, but without source taxonomy annotation",
        )

    if shared >= 0.70 and structural >= 0.75 and dentate < 0.70:
        return (
            "broad_neuronal_structural_warning",
            "low",
            "strong neuronal/structural program without strict dentate identity",
        )

    if shared >= 0.55 and background < 0.95:
        return (
            "hippocampal_neuronal_ambiguous",
            "low",
            "neuronal signal present but strict DG support is incomplete",
        )

    return ("low_information_or_low_qc", "low", "QC-pass nucleus with weak selected-panel support")


def summarize_labels(obs: pd.DataFrame) -> pd.DataFrame:
    module_aggs = {f"median_{module.replace('norm_', '')}": (module, "median") for module in MODULES}
    summary = (
        obs.groupby(["diagnosis", "projected_label", "projection_confidence"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            n_samples=("sample_id", "nunique"),
            median_counts=("nCount_RNA", "median"),
            median_genes=("nFeature_RNA", "median"),
            median_selected_counts=("nCount_selected", "median"),
            median_selected_genes=("nFeature_selected", "median"),
            **module_aggs,
        )
        .reset_index()
    )
    totals = summary.groupby("diagnosis")["n_cells"].transform("sum")
    summary["fraction_of_diagnosis"] = summary["n_cells"] / totals
    return summary.sort_values(["diagnosis", "n_cells"], ascending=[True, False])


def sample_summary(obs: pd.DataFrame) -> pd.DataFrame:
    module_aggs = {f"median_{module.replace('norm_', '')}": (module, "median") for module in MODULES}
    return (
        obs.groupby(["sample_id", "diagnosis", "projected_label"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            age_at_death_years=("age_at_death_years", "first"),
            pmi_hours=("pmi_hours", "first"),
            median_counts=("nCount_RNA", "median"),
            median_genes=("nFeature_RNA", "median"),
            **module_aggs,
        )
        .reset_index()
    )


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    idx = p_values.index[valid].to_numpy()
    if len(idx) == 0:
        return out
    order = p_values.loc[idx].to_numpy().argsort()
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    m = len(p_sorted)
    adjusted = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def diagnosis_tests(sample: pd.DataFrame) -> pd.DataFrame:
    included_labels = {
        "human_dg_like_high_confidence",
        "human_dg_like_candidate",
        "immature_neurogenic_candidate",
        "hippocampal_neuronal_ambiguous",
    }
    sub = sample.loc[sample["projected_label"].isin(included_labels)].copy()
    if sub.empty:
        return pd.DataFrame()
    by_sample = (
        sub.groupby(["sample_id", "diagnosis"], dropna=False)
        .agg(**{f"median_{module.replace('norm_', '')}": (f"median_{module.replace('norm_', '')}", "median") for module in MODULES})
        .reset_index()
    )
    comparisons = [("HA", "YA"), ("MCI", "HA"), ("AD", "HA"), ("SA", "HA"), ("AD", "YA")]
    rows = []
    for a, b in comparisons:
        for module in MODULES:
            col = f"median_{module.replace('norm_', '')}"
            av = by_sample.loc[by_sample["diagnosis"].eq(a), col].dropna().to_numpy(dtype=float)
            bv = by_sample.loc[by_sample["diagnosis"].eq(b), col].dropna().to_numpy(dtype=float)
            if len(av) < 2 or len(bv) < 2:
                continue
            p_value = stats.mannwhitneyu(av, bv, alternative="two-sided").pvalue
            rows.append(
                {
                    "comparison": f"{a}_vs_{b}",
                    "module": module,
                    "n_a": len(av),
                    "n_b": len(bv),
                    "median_a": float(np.median(av)),
                    "median_b": float(np.median(bv)),
                    "median_delta_a_minus_b": float(np.median(av) - np.median(bv)),
                    "p_value": float(p_value),
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj_bh"] = bh_adjust(out["p_value"])
    return out


def plot_heatmap(summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    module_cols = [f"median_{module.replace('norm_', '')}" for module in MODULES]
    heat = (
        summary.groupby(["projected_label"], dropna=False)
        .agg(**{col: (col, "median") for col in module_cols if col in summary.columns})
        .reset_index()
        .sort_values("projected_label")
    )
    if heat.empty:
        return
    matrix = heat[[col for col in module_cols if col in heat.columns]].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(10.5, max(3.5, 0.38 * len(heat))))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels(heat["projected_label"], fontsize=8)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels([col.replace("median_", "") for col in module_cols if col in heat.columns], rotation=35, ha="right", fontsize=8)
    ax.set_title("GSE268609 module medians by projected label")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("median normalized module score")
    fig.tight_layout()
    fig.savefig(OUT_HEATMAP, dpi=180)
    plt.close(fig)


def plot_deltas(tests: pd.DataFrame) -> None:
    if tests.empty:
        return
    plot_df = tests.loc[tests["comparison"].isin(["AD_vs_HA", "MCI_vs_HA", "HA_vs_YA"])].copy()
    plot_df = plot_df.loc[plot_df["module"].isin(MODULES[:-1])]
    if plot_df.empty:
        return
    plot_df["label"] = plot_df["comparison"] + "\n" + plot_df["module"].str.replace("norm_", "", regex=False)
    fig, ax = plt.subplots(figsize=(12, 4.8))
    colors = ["#1b7f6b" if value >= 0 else "#a74d2a" for value in plot_df["median_delta_a_minus_b"]]
    ax.bar(np.arange(len(plot_df)), plot_df["median_delta_a_minus_b"], color=colors)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_xticks(np.arange(len(plot_df)))
    ax.set_xticklabels(plot_df["label"], rotation=55, ha="right", fontsize=7)
    ax.set_ylabel("median sample-level delta")
    ax.set_title("GSE268609 exploratory diagnosis module deltas")
    fig.tight_layout()
    fig.savefig(OUT_DELTA_PLOT, dpi=180)
    plt.close(fig)


def write_md(obs: pd.DataFrame, summary: pd.DataFrame, tests: pd.DataFrame) -> None:
    label_counts = obs["projected_label"].value_counts()
    diagnosis_counts = obs.groupby("diagnosis")["sample_id"].nunique().sort_index()
    lines = [
        "# GSE268609 Human-Core Label Projection",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This pass maps the `GSE268609` human dentate/hippocampal RNA branch into the current human-core selected feature space.",
        "",
        "## Samples",
        "",
    ]
    for diagnosis, count in diagnosis_counts.items():
        lines.append(f"- `{diagnosis}`: {int(count)} RNA samples.")
    lines.extend(["", "## Projected Labels", ""])
    for label, count in label_counts.items():
        lines.append(f"- `{label}`: {int(count)} cells.")
    lines.extend(["", "## Exploratory Diagnosis Tests", ""])
    if tests.empty:
        lines.append("No diagnosis-level tests were run.")
    else:
        for _, row in tests.sort_values(["p_adj_bh", "comparison"]).head(12).iterrows():
            lines.append(
                f"- `{row['comparison']}` / `{row['module']}`: median delta {row['median_delta_a_minus_b']:.4f} "
                f"(n={int(row['n_a'])} vs {int(row['n_b'])}, BH-adjusted p={row['p_adj_bh']:.3g})."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `GSE268609` is suitable as a primary human dentate/hippocampal RNA candidate after RNA-row extraction, but labels here are marker/module projections rather than source taxonomy labels.",
            "- The source Seurat object or metadata with final cell-type annotations would further strengthen this dataset by distinguishing granule cells, immature neurons, interneurons, glia, and doublets.",
            "",
            "## Outputs",
            "",
            f"- Cell labels and scores: `{OUT_LABELS.relative_to(ROOT)}`",
            f"- Label summary: `{OUT_LABEL_SUMMARY.relative_to(ROOT)}`",
            f"- Sample/module summary: `{OUT_SAMPLE_SUMMARY.relative_to(ROOT)}`",
            f"- Diagnosis tests: `{OUT_TESTS.relative_to(ROOT)}`",
            f"- Heatmap: `{OUT_HEATMAP.relative_to(ROOT)}`",
            f"- Delta plot: `{OUT_DELTA_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    print("Loading selected GSE268609 RNA object", flush=True)
    X = sparse.load_npz(X_IN).tocsr()
    obs = pd.read_csv(OBS_IN, sep="\t", low_memory=False, dtype={"sample_id": str})
    var = pd.read_csv(VAR_IN, sep="\t")
    if X.shape != (len(obs), len(var)):
        raise ValueError(f"Shape mismatch: X {X.shape}, obs {len(obs)}, var {len(var)}")
    obs["analysis_include_basic_qc"] = obs["analysis_include_basic_qc"].astype(str).str.lower().isin(["true", "1", "yes"])

    print("Scoring human-core panels from counts", flush=True)
    scores = compute_panel_scores_from_counts(X, obs, var)
    obs = pd.concat([obs.reset_index(drop=True), scores.reset_index(drop=True)], axis=1)
    obs = add_ranks(obs)

    print("Assigning projected labels", flush=True)
    labels = obs.apply(projection_label, axis=1, result_type="expand")
    labels.columns = ["projected_label", "projection_confidence", "projection_reason"]
    obs = pd.concat([obs, labels], axis=1)
    obs["human_core_convention_label"] = obs["projected_label"]
    obs["analysis_include"] = obs["analysis_include_basic_qc"].astype(bool) & ~obs["projected_label"].isin(
        ["low_information_or_low_qc", "non_neuronal_background"]
    )

    print("Summarizing and plotting", flush=True)
    summary = summarize_labels(obs)
    sample = sample_summary(obs.loc[obs["analysis_include"]].copy())
    tests = diagnosis_tests(sample)

    keep_cols = [
        "cell_id",
        "sample_id",
        "sample_accession",
        "diagnosis",
        "age_at_death_years",
        "pmi_hours",
        "tissue",
        "nCount_RNA",
        "nFeature_RNA",
        "nCount_selected",
        "nFeature_selected",
        "selected_count_fraction",
        "analysis_include_basic_qc",
        "projected_label",
        "projection_confidence",
        "projection_reason",
        "human_core_convention_label",
        "analysis_include",
    ]
    score_cols = [col for col in obs.columns if col.startswith("norm_")]
    obs[keep_cols + score_cols].to_csv(OUT_LABELS, sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT_LABEL_SUMMARY, sep="\t", index=False)
    sample.to_csv(OUT_SAMPLE_SUMMARY, sep="\t", index=False)
    tests.to_csv(OUT_TESTS, sep="\t", index=False)
    plot_heatmap(summary)
    plot_deltas(tests)
    write_md(obs, summary, tests)

    print(f"Wrote {OUT_LABELS}")
    print(f"Wrote {OUT_LABEL_SUMMARY}")
    print(f"Wrote {OUT_SAMPLE_SUMMARY}")
    print(f"Wrote {OUT_TESTS}")
    print(f"Wrote {OUT_MD}")
    print(f"labels={len(obs)}; analysis_include={int(obs['analysis_include'].sum())}")


if __name__ == "__main__":
    main()
