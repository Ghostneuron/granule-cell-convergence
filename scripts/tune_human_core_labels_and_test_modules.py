#!/usr/bin/env python3
"""Tune human-core labels and run dataset-aware module tests."""

from __future__ import annotations

import csv
import math
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from validate_human_core_marker_programs import BACKGROUND_PANELS, load_panels


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "Project/processed"
RESULTS = ROOT / "Project/results"
OBJECT_DIR = PROCESSED / "human_core_normalized_reduced_object"

X_IN = OBJECT_DIR / "X_log1p_cp10k_selected_genes.npz"
OBS_IN = OBJECT_DIR / "obs.tsv.gz"
VAR_IN = OBJECT_DIR / "var.tsv"
SVD_IN = RESULTS / "human_core_normalized_svd_embedding.tsv.gz"

OUT_LABELS = RESULTS / "human_core_tuned_labels.tsv.gz"
OUT_LABEL_SUMMARY = RESULTS / "human_core_tuned_label_summary.tsv"
OUT_REPLICATE_SUMMARY = RESULTS / "human_core_tuned_label_replicate_summary.tsv"
OUT_MODULE_TESTS = RESULTS / "human_core_dataset_aware_module_tests.tsv"
OUT_HEATMAP = RESULTS / "human_core_tuned_label_module_heatmap.png"
OUT_DELTA_PLOT = RESULTS / "human_core_dataset_aware_module_deltas.png"
OUT_MD = RESULTS / "human_core_label_tuning_and_module_tests.md"

PRIMARY_MODULES = [
    "norm_dentate_identity",
    "norm_human_dg_immature",
    "norm_human_dg_mature",
    "norm_cerebellar_identity",
    "norm_shared_granule_neuronal",
    "norm_morphogenesis_cytoskeleton",
    "norm_axon_guidance_synapse",
    "norm_background_max",
]

COMPARISONS = [
    ("human_dg_like_high_confidence", "non_neuronal_background", "dg_high_vs_background"),
    ("human_dg_like_high_confidence", "hippocampal_neuronal_ambiguous", "dg_high_vs_neuronal_ambiguous"),
    ("immature_neurogenic_candidate", "non_neuronal_background", "immature_vs_background"),
    ("broad_neuronal_structural_warning", "hippocampal_neuronal_ambiguous", "warning_vs_neuronal_ambiguous"),
]

MIN_CELLS_PER_LABEL_UNIT = 20


def norm_gene(value: object) -> str:
    return str(value).strip().strip('"').strip("'").upper()


def percentile_rank(series: pd.Series) -> pd.Series:
    if len(series) <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def load_object() -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    X = sparse.load_npz(X_IN).tocsr()
    obs = pd.read_csv(OBS_IN, sep="\t", low_memory=False)
    var = pd.read_csv(VAR_IN, sep="\t")
    if X.shape != (len(obs), len(var)):
        raise ValueError(f"Shape mismatch: X {X.shape}, obs {len(obs)}, var {len(var)}")
    return X, obs, var


def compute_normalized_panel_scores(X: sparse.csr_matrix, var: pd.DataFrame) -> pd.DataFrame:
    panels = load_panels()
    gene_to_idx = {norm_gene(gene): idx for idx, gene in enumerate(var["gene"])}
    score_cols: dict[str, np.ndarray] = {}
    detected_cols: dict[str, np.ndarray] = {}

    for panel, genes in panels.items():
        panel_genes = [norm_gene(gene) for gene in genes]
        idxs = [gene_to_idx[gene] for gene in panel_genes if gene in gene_to_idx]
        if idxs:
            sub = X[:, idxs]
            values = np.asarray(sub.mean(axis=1)).ravel()
            detected = np.asarray((sub > 0).sum(axis=1)).ravel()
        else:
            values = np.zeros(X.shape[0], dtype=float)
            detected = np.zeros(X.shape[0], dtype=int)
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


def replicate_unit(row: pd.Series) -> str:
    if row["dataset"] == "GSE186538":
        return str(row.get("sample_hint", "") or row.get("specimen_id", "") or row["component_id"])
    specimen = str(row.get("specimen_id", ""))
    if specimen and specimen.lower() != "nan":
        return f"specimen_{specimen}"
    return str(row["component_id"])


def tune_one(row: pd.Series) -> tuple[str, str, str]:
    marker_call = row["marker_call"]
    dataset = row["dataset"]

    if marker_call == "curated_human_dg_gc_reference" and dataset == "GSE186538":
        return (
            "curated_human_dg_gc_anchor",
            "anchor",
            "source taxonomy labels this QC-pass subset as DG GC PROX1",
        )

    neuronal_rank = max(float(row["norm_excitatory_neuron_rank"]), float(row["norm_shared_granule_neuronal_rank"]))
    background_dominates = float(row["norm_background_max"]) > max(
        float(row["norm_excitatory_neuron"]), float(row["norm_shared_granule_neuronal"])
    )

    if marker_call == "likely_non_neuronal_background" or (background_dominates and float(row["norm_background_max_rank"]) >= 0.70):
        return (
            "non_neuronal_background",
            "medium",
            f"background-dominant panel: {row['norm_background_dominant_panel']}",
        )

    if marker_call == "marker_supported_human_dg_like":
        high = (
            float(row["norm_dentate_identity_rank"]) >= 0.50
            and float(row["norm_shared_granule_neuronal_rank"]) >= 0.35
            and float(row["norm_identity_contrast_dentate_minus_cerebellar"]) >= 0
            and float(row["norm_background_max_rank"]) < 0.90
        )
        if high:
            state = "immature" if row["norm_human_dg_immature"] >= row["norm_human_dg_mature"] else "mature"
            return (
                "human_dg_like_high_confidence",
                "medium",
                f"{state}-shifted dentate identity with shared neuronal support",
            )
        return (
            "human_dg_like_candidate",
            "low",
            "original marker-supported DG-like call but normalized ranks are below strict high-confidence thresholds",
        )

    if marker_call == "immature_neuron_or_neurogenic_candidate":
        if float(row["norm_human_dg_immature_rank"]) >= 0.55 and neuronal_rank >= 0.35:
            return (
                "immature_neurogenic_candidate",
                "medium",
                "immature/neurogenic module high with neuronal support",
            )
        return (
            "immature_neurogenic_candidate_low_support",
            "low",
            "immature marker call retained but normalized support is modest",
        )

    if marker_call == "cerebellar_marker_high_warning":
        return (
            "broad_neuronal_structural_warning",
            "low",
            "hippocampal cell with cerebellar-panel signal; keep as broad structural/neurogenic warning, not cerebellar identity",
        )

    if marker_call == "neuronal_non_dg_or_ambiguous":
        return (
            "hippocampal_neuronal_ambiguous",
            "low",
            "neuronal signal present but strict DG support is incomplete",
        )

    if marker_call == "low_information_or_low_qc":
        return (
            "low_information_or_low_qc",
            "low",
            "passed preliminary QC but remains low-information in marker/object-level pass",
        )

    return ("other_or_unmapped", "low", "no tuned rule matched")


def add_tuned_labels(obs: pd.DataFrame, norm_scores: pd.DataFrame, embedding: pd.DataFrame) -> pd.DataFrame:
    out = pd.concat([obs.reset_index(drop=True), norm_scores.reset_index(drop=True)], axis=1)
    out["replicate_unit"] = out.apply(replicate_unit, axis=1)

    rank_cols = [
        "norm_dentate_identity",
        "norm_cerebellar_identity",
        "norm_human_dg_immature",
        "norm_human_dg_mature",
        "norm_shared_granule_neuronal",
        "norm_excitatory_neuron",
        "norm_background_max",
        "norm_structural_program_mean",
    ]
    for col in rank_cols:
        out[f"{col}_rank"] = out.groupby(["dataset", "replicate_unit"], dropna=False)[col].transform(percentile_rank)

    tuned = out.apply(tune_one, axis=1)
    out["tuned_label"] = [item[0] for item in tuned]
    out["tuned_confidence"] = [item[1] for item in tuned]
    out["tuned_reason"] = [item[2] for item in tuned]

    emb_cols = ["cell_id"] + [col for col in embedding.columns if col.startswith("svd_")][:10]
    out = out.merge(embedding[emb_cols], on="cell_id", how="left", validate="one_to_one")
    svd_cols = [col for col in out.columns if col.startswith("svd_")]
    anchor = out.loc[out["tuned_label"] == "curated_human_dg_gc_anchor", svd_cols].median().to_numpy(dtype=float)
    svd_values = out[svd_cols].to_numpy(dtype=float)
    out["svd_distance_to_dg_anchor"] = np.sqrt(((svd_values - anchor) ** 2).sum(axis=1))
    return out


def label_summary(labels: pd.DataFrame) -> pd.DataFrame:
    summary = (
        labels.groupby(["dataset", "tuned_label", "tuned_confidence"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            n_replicate_units=("replicate_unit", "nunique"),
            median_counts=("n_counts", "median"),
            median_genes=("n_genes", "median"),
            median_norm_dentate_identity=("norm_dentate_identity", "median"),
            median_norm_cerebellar_identity=("norm_cerebellar_identity", "median"),
            median_norm_human_dg_immature=("norm_human_dg_immature", "median"),
            median_norm_human_dg_mature=("norm_human_dg_mature", "median"),
            median_norm_shared_granule_neuronal=("norm_shared_granule_neuronal", "median"),
            median_norm_structural_program=("norm_structural_program_mean", "median"),
            median_norm_background_max=("norm_background_max", "median"),
            median_svd_distance_to_dg_anchor=("svd_distance_to_dg_anchor", "median"),
        )
        .reset_index()
    )
    totals = summary.groupby("dataset")["n_cells"].transform("sum")
    summary["fraction_of_dataset"] = summary["n_cells"] / totals
    return summary.sort_values(["dataset", "n_cells"], ascending=[True, False])


def replicate_summary(labels: pd.DataFrame) -> pd.DataFrame:
    agg = (
        labels.groupby(["dataset", "replicate_unit", "tuned_label"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            **{f"mean_{module}": (module, "mean") for module in PRIMARY_MODULES},
        )
        .reset_index()
    )
    return agg.sort_values(["dataset", "replicate_unit", "tuned_label"])


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    valid = p.notna()
    out = pd.Series(np.nan, index=p.index)
    if valid.sum() == 0:
        return out
    ranks = p[valid].rank(method="first").astype(int)
    m = valid.sum()
    adjusted = p[valid] * m / ranks
    order = p[valid].sort_values(ascending=False).index
    running = 1.0
    for idx in order:
        running = min(running, adjusted.loc[idx])
        out.loc[idx] = min(running, 1.0)
    return out


def module_tests(rep: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset in sorted(rep["dataset"].unique()):
        dataset_rep = rep[rep["dataset"] == dataset]
        for target, reference, comparison in COMPARISONS:
            wide = dataset_rep[dataset_rep["tuned_label"].isin([target, reference])]
            if wide.empty:
                continue
            for module in PRIMARY_MODULES:
                table = wide.pivot_table(
                    index="replicate_unit",
                    columns="tuned_label",
                    values=f"mean_{module}",
                    aggfunc="first",
                )
                n_target = wide.loc[wide["tuned_label"] == target].set_index("replicate_unit")["n_cells"]
                n_reference = wide.loc[wide["tuned_label"] == reference].set_index("replicate_unit")["n_cells"]
                table = table.dropna(subset=[target, reference], how="any")
                if table.empty:
                    continue
                valid_units = [
                    unit
                    for unit in table.index
                    if n_target.get(unit, 0) >= MIN_CELLS_PER_LABEL_UNIT and n_reference.get(unit, 0) >= MIN_CELLS_PER_LABEL_UNIT
                ]
                table = table.loc[valid_units]
                if table.empty:
                    continue
                delta = table[target] - table[reference]
                p_value = np.nan
                statistic = np.nan
                if len(delta) >= 3 and (delta != 0).any():
                    try:
                        statistic, p_value = stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided")
                    except ValueError:
                        statistic, p_value = np.nan, np.nan
                rows.append(
                    {
                        "dataset": dataset,
                        "comparison": comparison,
                        "target_label": target,
                        "reference_label": reference,
                        "module": module,
                        "n_replicate_units": len(delta),
                        "median_delta_target_minus_reference": float(delta.median()),
                        "mean_delta_target_minus_reference": float(delta.mean()),
                        "fraction_units_delta_gt0": float((delta > 0).mean()),
                        "wilcoxon_statistic": statistic,
                        "p_value": p_value,
                        "replicate_units": ";".join(map(str, delta.index)),
                    }
                )
    tests = pd.DataFrame(rows)
    if not tests.empty:
        tests["p_adj_bh"] = benjamini_hochberg(tests["p_value"])
    return tests.sort_values(["dataset", "comparison", "module"]) if not tests.empty else tests


def plot_heatmap(summary: pd.DataFrame) -> None:
    modules = [
        "median_norm_dentate_identity",
        "median_norm_human_dg_immature",
        "median_norm_human_dg_mature",
        "median_norm_shared_granule_neuronal",
        "median_norm_structural_program",
        "median_norm_background_max",
    ]
    heat = summary.copy()
    heat["row_label"] = heat["dataset"] + " | " + heat["tuned_label"]
    heat = heat.sort_values(["dataset", "tuned_label"])
    matrix = heat[modules].to_numpy(dtype=float)
    fig_h = max(4.5, 0.32 * len(heat))
    fig, ax = plt.subplots(figsize=(9.5, fig_h))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels(heat["row_label"], fontsize=7)
    ax.set_xticks(np.arange(len(modules)))
    ax.set_xticklabels([m.replace("median_norm_", "") for m in modules], rotation=35, ha="right", fontsize=8)
    ax.set_title("Human core tuned-label module medians")
    cbar = fig.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label("median normalized module score")
    fig.tight_layout()
    fig.savefig(OUT_HEATMAP, dpi=180)
    plt.close(fig)


def plot_deltas(tests: pd.DataFrame) -> None:
    if tests.empty:
        return
    plot_modules = [
        "norm_dentate_identity",
        "norm_human_dg_immature",
        "norm_human_dg_mature",
        "norm_shared_granule_neuronal",
        "norm_morphogenesis_cytoskeleton",
        "norm_axon_guidance_synapse",
        "norm_background_max",
    ]
    plot = tests[
        tests["comparison"].isin(["dg_high_vs_background", "dg_high_vs_neuronal_ambiguous"])
        & tests["module"].isin(plot_modules)
    ].copy()
    if plot.empty:
        return
    plot["label"] = plot["dataset"] + "\n" + plot["comparison"].str.replace("_", " ")
    labels = plot["label"].drop_duplicates().tolist()
    x = np.arange(len(labels))
    width = 0.1
    colors = ["#227c70", "#d95f02", "#7570b3", "#1f78b4", "#b15928", "#6a994e", "#6c757d"]
    fig, ax = plt.subplots(figsize=(13.5, 5.4))
    for idx, module in enumerate(plot_modules):
        sub = plot[plot["module"] == module].set_index("label").reindex(labels)
        ax.bar(
            x + (idx - len(plot_modules) / 2) * width + width / 2,
            sub["median_delta_target_minus_reference"],
            width,
            label=module.replace("norm_", ""),
            color=colors[idx],
        )
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=12, ha="right")
    ax.set_ylabel("median replicate delta")
    ax.set_title("Dataset-aware module deltas for DG-like labels")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DELTA_PLOT, dpi=180)
    plt.close(fig)


def write_summary_md(labels: pd.DataFrame, summary: pd.DataFrame, tests: pd.DataFrame) -> None:
    counts = summary.groupby("tuned_label")["n_cells"].sum().sort_values(ascending=False)
    lines = [
        "# Human Core Label Tuning And Module Tests",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This pass tunes labels inside the normalized reduced human-core object and tests marker/structural modules using replicate-level summaries. `GSE186538` remains the strict human DG granule-cell anchor; `GSE185277` and `GSE185553` provide candidate-rich hippocampal contexts.",
        "",
        "## Tuned Label Counts",
        "",
    ]
    for label, n_cells in counts.items():
        lines.append(f"- `{label}`: {int(n_cells)} cells.")

    lines.extend(["", "## Dataset-Aware Test Highlights", ""])
    if tests.empty:
        lines.append("- No paired replicate tests were available.")
    else:
        highlight = tests[
            tests["comparison"].isin(["dg_high_vs_background", "dg_high_vs_neuronal_ambiguous"])
            & tests["module"].isin(["norm_dentate_identity", "norm_shared_granule_neuronal", "norm_morphogenesis_cytoskeleton", "norm_background_max"])
        ].copy()
        highlight = highlight.sort_values(["dataset", "comparison", "module"])
        for _, row in highlight.iterrows():
            lines.append(
                f"- `{row['dataset']}` `{row['comparison']}` `{row['module']}`: "
                f"median delta {row['median_delta_target_minus_reference']:.4f} across {int(row['n_replicate_units'])} replicate units "
                f"(BH-adjusted p={row['p_adj_bh']:.3g})."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The tuned labels intentionally keep `GSE186538` as the strict anchor and avoid calling hippocampal warning cells cerebellar granule cells.",
            "- The useful signal is whether DG-like hippocampal candidates show higher dentate/shared-structural modules than matched background or ambiguous cells within the same dataset/replicate structure.",
            "- These tuned labels have now been used for the `GSE325391` adult DG projection, but they should still be treated as computational labels rather than final manual annotation.",
            "",
            "## Outputs",
            "",
            f"- Tuned labels: `{OUT_LABELS.relative_to(ROOT)}`",
            f"- Label summary: `{OUT_LABEL_SUMMARY.relative_to(ROOT)}`",
            f"- Replicate summary: `{OUT_REPLICATE_SUMMARY.relative_to(ROOT)}`",
            f"- Dataset-aware tests: `{OUT_MODULE_TESTS.relative_to(ROOT)}`",
            f"- Heatmap: `{OUT_HEATMAP.relative_to(ROOT)}`",
            f"- Delta plot: `{OUT_DELTA_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    print("Loading normalized object", flush=True)
    X, obs, var = load_object()
    embedding = pd.read_csv(SVD_IN, sep="\t", low_memory=False)

    print("Computing normalized panel scores", flush=True)
    norm_scores = compute_normalized_panel_scores(X, var)

    print("Tuning labels", flush=True)
    labels = add_tuned_labels(obs, norm_scores, embedding)
    summary = label_summary(labels)
    rep = replicate_summary(labels)
    tests = module_tests(rep)

    keep_cols = [
        "cell_id",
        "dataset",
        "component_id",
        "replicate_unit",
        "analysis_label",
        "marker_call",
        "marker_state",
        "tuned_label",
        "tuned_confidence",
        "tuned_reason",
        "n_counts",
        "n_genes",
        "percent_mt",
        "geo_development_stage",
        "geo_age",
        "geo_age_years",
        "geo_sex",
        "norm_dentate_identity",
        "norm_cerebellar_identity",
        "norm_human_dg_immature",
        "norm_human_dg_mature",
        "norm_shared_granule_neuronal",
        "norm_morphogenesis_cytoskeleton",
        "norm_axon_guidance_synapse",
        "norm_structural_program_mean",
        "norm_background_max",
        "norm_background_dominant_panel",
        "norm_identity_contrast_dentate_minus_cerebellar",
        "svd_distance_to_dg_anchor",
    ]
    rank_cols = [col for col in labels.columns if col.endswith("_rank") and col.startswith("norm_")]
    keep_cols.extend(rank_cols)
    labels.to_csv(OUT_LABELS, sep="\t", index=False, columns=[col for col in keep_cols if col in labels.columns], compression="gzip", float_format="%.6g")
    summary.to_csv(OUT_LABEL_SUMMARY, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)
    rep.to_csv(OUT_REPLICATE_SUMMARY, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)
    tests.to_csv(OUT_MODULE_TESTS, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)

    print("Plotting", flush=True)
    plot_heatmap(summary)
    plot_deltas(tests)
    write_summary_md(labels, summary, tests)

    print(f"Wrote {OUT_LABELS}")
    print(f"Wrote {OUT_LABEL_SUMMARY}")
    print(f"Wrote {OUT_REPLICATE_SUMMARY}")
    print(f"Wrote {OUT_MODULE_TESTS}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
