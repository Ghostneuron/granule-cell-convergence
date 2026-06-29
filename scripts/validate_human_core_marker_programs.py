#!/usr/bin/env python3
"""Validate marker programs in the built human dentate/hippocampal core.

This is a marker-QC layer, not a final cell-type annotation. It scores the
currently built human sparse objects with the refined granule panels plus
human-focused maturation and background panels.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "Project/config/granule_marker_panels_refined.tsv"
RESULTS = ROOT / "Project/results"

QC_CELLS = RESULTS / "human_core_harmonized_cell_qc.tsv.gz"
SEED_SUMMARY = RESULTS / "human_seed_sparse_object_summary.tsv"
GSE186538_SUMMARY = RESULTS / "gse186538_human_dg_gc_sparse_subset_summary.tsv"

OUT_SCORES = RESULTS / "human_core_marker_validation_scores.tsv.gz"
OUT_PANEL_SUMMARY = RESULTS / "human_core_marker_validation_panel_summary.tsv"
OUT_CALL_SUMMARY = RESULTS / "human_core_marker_validation_call_summary.tsv"
OUT_GENE_COVERAGE = RESULTS / "human_core_marker_validation_gene_coverage.tsv"
OUT_PLOT = RESULTS / "human_core_marker_validation_panel_medians.png"
OUT_MD = RESULTS / "human_core_marker_validation_interpretation.md"


EXTRA_PANELS = {
    "human_dg_immature": [
        "PROX1",
        "DCX",
        "NEUROD1",
        "NEUROD2",
        "SOX11",
        "SOX4",
        "MEX3A",
        "CALB2",
        "STMN2",
        "TUBB3",
        "BCL11B",
    ],
    "human_dg_mature": [
        "PROX1",
        "CALB1",
        "SLC17A7",
        "GRIN1",
        "GRIN2A",
        "CAMK2A",
        "RIMS1",
        "SYT1",
        "SNAP25",
        "RBFOX3",
    ],
    "excitatory_neuron": [
        "SLC17A7",
        "SLC17A6",
        "CAMK2A",
        "SATB2",
        "BCL11B",
        "RBFOX3",
        "SNAP25",
        "SYT1",
    ],
    "inhibitory_neuron": ["GAD1", "GAD2", "SLC6A1", "DLX1", "DLX2"],
    "astrocyte_background": ["AQP4", "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3"],
    "oligodendrocyte_background": ["MBP", "PLP1", "MOG", "MOBP", "CLDN11"],
    "opc_background": ["PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2"],
    "microglia_background": ["CX3CR1", "P2RY12", "AIF1", "C1QA", "TYROBP"],
    "vascular_background": ["CLDN5", "PECAM1", "VWF", "KDR", "FLT1"],
    "ependymal_choroid_background": ["FOXJ1", "TTR", "AQP1", "PIFO", "DNAH5"],
}

BACKGROUND_PANELS = [
    "astrocyte_background",
    "oligodendrocyte_background",
    "opc_background",
    "microglia_background",
    "vascular_background",
    "ependymal_choroid_background",
]

FOCAL_PLOT_PANELS = [
    "dentate_identity",
    "human_dg_immature",
    "human_dg_mature",
    "cerebellar_identity",
    "shared_granule_neuronal",
    "morphogenesis_cytoskeleton",
]


def norm_gene(gene: str) -> str:
    return str(gene).strip().strip('"').strip("'").upper()


def add_gene(panel_map: dict[str, list[str]], panel: str, gene: str) -> None:
    gene_norm = norm_gene(gene)
    if gene_norm and gene_norm not in panel_map[panel]:
        panel_map[panel].append(gene_norm)


def load_panels() -> dict[str, list[str]]:
    panels: dict[str, list[str]] = defaultdict(list)
    config = pd.read_csv(CONFIG, sep="\t")
    for _, row in config.iterrows():
        add_gene(panels, row["panel"], row["gene"])
    for panel, genes in EXTRA_PANELS.items():
        for gene in genes:
            add_gene(panels, panel, gene)
    return dict(sorted(panels.items()))


def read_gene_names(path: Path) -> list[str]:
    genes = pd.read_csv(path, sep="\t")
    gene_col = "gene" if "gene" in genes.columns else genes.columns[0]
    return [norm_gene(gene) for gene in genes[gene_col].astype(str)]


def build_matrix_to_gene_path() -> dict[str, str]:
    mapping: dict[str, str] = {}
    seed = pd.read_csv(SEED_SUMMARY, sep="\t")
    for _, row in seed.iterrows():
        mapping[row["matrix_path"]] = row["gene_metadata_path"]

    gse186538 = pd.read_csv(GSE186538_SUMMARY, sep="\t")
    for _, row in gse186538.iterrows():
        mapping[row["matrix_path"]] = row["gene_metadata_path"]
    return mapping


def score_matrix(matrix: sparse.csr_matrix, genes: list[str], panels: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_to_cols: dict[str, list[int]] = defaultdict(list)
    for idx, gene in enumerate(genes):
        gene_to_cols[gene].append(idx)

    all_target_cols = sorted({col for panel_genes in panels.values() for gene in panel_genes for col in gene_to_cols.get(gene, [])})
    selected_pos = {col: pos for pos, col in enumerate(all_target_cols)}
    selected = matrix[:, all_target_cols].toarray() if all_target_cols else np.zeros((matrix.shape[0], 0), dtype=float)

    score_cols: dict[str, np.ndarray] = {}
    detected_cols: dict[str, np.ndarray] = {}
    coverage_rows: list[dict[str, object]] = []

    for panel, panel_genes in panels.items():
        found_genes = [gene for gene in panel_genes if gene in gene_to_cols]
        sum_expr = np.zeros(matrix.shape[0], dtype=float)
        logsum = np.zeros(matrix.shape[0], dtype=float)
        detected = np.zeros(matrix.shape[0], dtype=np.int16)

        for gene in found_genes:
            col_positions = [selected_pos[col] for col in gene_to_cols[gene]]
            gene_expr = selected[:, col_positions].sum(axis=1)
            sum_expr += gene_expr
            logsum += np.log1p(gene_expr)
            detected += (gene_expr > 0).astype(np.int16)

        panel_gene_count = len(panel_genes)
        found_count = len(found_genes)
        score = logsum / panel_gene_count if panel_gene_count else logsum
        score_cols[f"score_{panel}"] = score
        detected_cols[f"detected_{panel}"] = detected
        coverage_rows.append(
            {
                "panel": panel,
                "panel_gene_count": panel_gene_count,
                "genes_found_in_matrix": found_count,
                "genes_found": ",".join(found_genes),
                "genes_missing": ",".join([gene for gene in panel_genes if gene not in gene_to_cols]),
            }
        )

    return pd.DataFrame({**score_cols, **detected_cols}), pd.DataFrame(coverage_rows)


def percentile_rank(series: pd.Series) -> pd.Series:
    if len(series) <= 1:
        return pd.Series(np.ones(len(series)), index=series.index)
    return series.rank(method="average", pct=True)


def call_marker_status(row: pd.Series) -> tuple[str, str, str]:
    if str(row["dataset"]) == "GSE186538" and str(row["cluster"]).startswith("DG GC"):
        state = "immature_shifted" if row["score_human_dg_immature"] > row["score_human_dg_mature"] else "mature_shifted"
        return ("curated_human_dg_gc_reference", state, "source taxonomy labels this subset as DG GC PROX1")

    if not bool(row["preliminary_qc_pass"]):
        return ("low_information_or_low_qc", "not_assigned", "below preliminary QC diagnostic threshold")

    neuronal = max(row["score_shared_granule_neuronal"], row["score_excitatory_neuron"])
    background = row["background_max_score"]
    if background > neuronal and row["background_max_rank"] >= 0.60:
        return ("likely_non_neuronal_background", row["background_dominant_panel"], "background panel is stronger than neuronal panels")

    if row["score_cerebellar_identity"] > row["score_dentate_identity"] and row["cerebellar_identity_rank"] >= 0.80:
        return ("cerebellar_marker_high_warning", "not_assigned", "cerebellar identity panel exceeds dentate panel in hippocampal source")

    if (
        row["score_dentate_identity"] > row["score_cerebellar_identity"]
        and row["dentate_identity_rank"] >= 0.55
        and row["shared_granule_neuronal_rank"] >= 0.35
    ):
        state = "immature_shifted" if row["score_human_dg_immature"] >= row["score_human_dg_mature"] else "mature_shifted"
        return ("marker_supported_human_dg_like", state, "dentate identity exceeds cerebellar identity with shared neuronal support")

    if row["human_dg_immature_rank"] >= 0.70 and row["shared_granule_neuronal_rank"] >= 0.35:
        return ("immature_neuron_or_neurogenic_candidate", "immature_shifted", "immature DG/neurogenic panel is high with neuronal support")

    if max(row["shared_granule_neuronal_rank"], row["excitatory_neuron_rank"]) >= 0.50:
        return ("neuronal_non_dg_or_ambiguous", "not_assigned", "neuronal signal is present but dentate-granule support is incomplete")

    return ("low_information_or_low_qc", "not_assigned", "weak marker support in this sparse marker pass")


def add_marker_calls(df: pd.DataFrame) -> pd.DataFrame:
    score_cols = [col for col in df.columns if col.startswith("score_")]
    background_score_cols = [f"score_{panel}" for panel in BACKGROUND_PANELS]
    df["background_max_score"] = df[background_score_cols].max(axis=1)
    df["background_dominant_panel"] = df[background_score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    df["identity_contrast_dentate_minus_cerebellar"] = df["score_dentate_identity"] - df["score_cerebellar_identity"]
    df["structural_program_mean"] = df[
        ["score_shared_granule_neuronal", "score_morphogenesis_cytoskeleton", "score_axon_guidance_synapse"]
    ].mean(axis=1)

    rank_inputs = [
        "score_dentate_identity",
        "score_cerebellar_identity",
        "score_human_dg_immature",
        "score_shared_granule_neuronal",
        "score_excitatory_neuron",
        "background_max_score",
    ]
    rank_names = {
        "score_dentate_identity": "dentate_identity_rank",
        "score_cerebellar_identity": "cerebellar_identity_rank",
        "score_human_dg_immature": "human_dg_immature_rank",
        "score_shared_granule_neuronal": "shared_granule_neuronal_rank",
        "score_excitatory_neuron": "excitatory_neuron_rank",
        "background_max_score": "background_max_rank",
    }
    for col in rank_inputs:
        df[rank_names[col]] = df.groupby(["dataset", "component_id"], dropna=False)[col].transform(percentile_rank)

    call_tuples = df.apply(call_marker_status, axis=1)
    df["marker_call"] = [item[0] for item in call_tuples]
    df["marker_state"] = [item[1] for item in call_tuples]
    df["marker_call_reason"] = [item[2] for item in call_tuples]

    keep_cols = [
        "cell_id",
        "dataset",
        "component_id",
        "component_type",
        "sample_hint",
        "age_hint",
        "region",
        "cluster",
        "n_counts",
        "n_genes",
        "percent_mt",
        "preliminary_qc_pass",
        "marker_call",
        "marker_state",
        "marker_call_reason",
        "identity_contrast_dentate_minus_cerebellar",
        "structural_program_mean",
        "background_max_score",
        "background_dominant_panel",
    ]
    keep_cols.extend(sorted(score_cols))
    keep_cols.extend(sorted([col for col in df.columns if col.startswith("detected_")]))
    keep_cols.extend(
        [
            "dentate_identity_rank",
            "cerebellar_identity_rank",
            "human_dg_immature_rank",
            "shared_granule_neuronal_rank",
            "excitatory_neuron_rank",
            "background_max_rank",
            "matrix_path",
            "cell_metadata_path",
        ]
    )
    return df[keep_cols]


def summarize_panels(scores: pd.DataFrame, panels: dict[str, list[str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["dataset", "component_id", "region", "cluster"]
    for keys, group in scores.groupby(group_cols, dropna=False):
        base = dict(zip(group_cols, keys))
        base["n_cells"] = len(group)
        base["preliminary_qc_pass_cells"] = int(group["preliminary_qc_pass"].sum())
        for panel, genes in panels.items():
            score = group[f"score_{panel}"]
            detected = group[f"detected_{panel}"]
            row = base.copy()
            row.update(
                {
                    "panel": panel,
                    "panel_gene_count": len(genes),
                    "median_score": float(score.median()),
                    "mean_score": float(score.mean()),
                    "p75_score": float(score.quantile(0.75)),
                    "fraction_detecting_any_panel_gene": float((detected > 0).mean()),
                    "median_detected_genes": float(detected.median()),
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_calls(scores: pd.DataFrame) -> pd.DataFrame:
    summary = (
        scores.groupby(["dataset", "component_id", "region", "cluster", "marker_call", "marker_state"], dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            median_dentate_identity=("score_dentate_identity", "median"),
            median_cerebellar_identity=("score_cerebellar_identity", "median"),
            median_identity_contrast=("identity_contrast_dentate_minus_cerebellar", "median"),
            median_shared_granule_neuronal=("score_shared_granule_neuronal", "median"),
            median_structural_program=("structural_program_mean", "median"),
            median_background_max=("background_max_score", "median"),
        )
        .reset_index()
    )
    total = summary.groupby(["dataset", "component_id", "region", "cluster"], dropna=False)["n_cells"].transform("sum")
    summary["fraction_of_component_or_cluster"] = summary["n_cells"] / total
    return summary.sort_values(["dataset", "component_id", "cluster", "marker_call", "marker_state"])


def summarize_coverage(coverage_parts: list[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(coverage_parts, ignore_index=True).sort_values(["dataset", "component_id", "panel"])


def plot_panel_medians(panel_summary: pd.DataFrame) -> None:
    focal = panel_summary[panel_summary["panel"].isin(FOCAL_PLOT_PANELS)].copy()
    dataset_panel = (
        focal.groupby(["dataset", "panel"], dropna=False)
        .apply(lambda g: np.average(g["median_score"], weights=g["n_cells"]), include_groups=False)
        .reset_index(name="weighted_component_median_score")
    )
    pivot = dataset_panel.pivot(index="dataset", columns="panel", values="weighted_component_median_score").reindex(
        columns=FOCAL_PLOT_PANELS
    )

    fig, ax = plt.subplots(figsize=(10, 3.8))
    x = np.arange(len(pivot.index))
    width = 0.12
    colors = ["#227c70", "#d95f02", "#7570b3", "#6c757d", "#1f78b4", "#b15928"]
    for idx, panel in enumerate(FOCAL_PLOT_PANELS):
        ax.bar(x + (idx - len(FOCAL_PLOT_PANELS) / 2) * width + width / 2, pivot[panel], width, label=panel, color=colors[idx])
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=0)
    ax.set_ylabel("Weighted median log1p marker score")
    ax.set_title("Human core marker-validation panel medians")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left", bbox_to_anchor=(0, 1.22))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180)
    plt.close(fig)


def write_interpretation(scores: pd.DataFrame, call_summary: pd.DataFrame, panel_summary: pd.DataFrame) -> None:
    dataset_counts = scores.groupby("dataset").size().to_dict()
    dataset_call = (
        scores.groupby(["dataset", "marker_call"]).size().reset_index(name="n_cells").sort_values(["dataset", "n_cells"], ascending=[True, False])
    )
    call_lines = []
    for dataset, sub in dataset_call.groupby("dataset"):
        total = int(dataset_counts[dataset])
        top = sub.head(5)
        parts = [f"{row.marker_call}: {int(row.n_cells)} ({row.n_cells / total * 100:.1f}%)" for row in top.itertuples()]
        call_lines.append(f"- `{dataset}`: " + "; ".join(parts) + ".")

    focal_panel = panel_summary[panel_summary["panel"].isin(["dentate_identity", "human_dg_immature", "human_dg_mature", "cerebellar_identity"])]
    dataset_panel = (
        focal_panel.groupby(["dataset", "panel"], dropna=False)
        .apply(lambda g: np.average(g["median_score"], weights=g["n_cells"]), include_groups=False)
        .reset_index(name="weighted_component_median_score")
    )
    panel_lines = []
    for dataset, sub in dataset_panel.groupby("dataset"):
        values = {row.panel: row.weighted_component_median_score for row in sub.itertuples()}
        panel_lines.append(
            f"- `{dataset}`: dentate {values.get('dentate_identity', 0):.3f}, "
            f"immature {values.get('human_dg_immature', 0):.3f}, mature {values.get('human_dg_mature', 0):.3f}, "
            f"cerebellar {values.get('cerebellar_identity', 0):.3f}."
        )

    text = [
        "# Human Core Marker Validation",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This pass scores the built human core objects with refined granule-cell identity panels, human DG maturation panels, and non-neuronal background panels. Calls are triage labels for object construction, not final biological annotation.",
        "",
        "## Dataset-Level Call Mix",
        "",
        *call_lines,
        "",
        "## Focal Panel Medians",
        "",
        *panel_lines,
        "",
        "## Interpretation",
        "",
        "- `GSE186538` now functions as the explicit human DG granule-cell anchor because the extracted cells carry source taxonomy labels beginning with `DG GC`.",
        "- `GSE185277` and `GSE185553` are retained as human hippocampal/dentate construction references, but their marker-supported subsets should be treated as candidates until richer sample and cell-state annotations are added.",
        "- The useful analysis contrast remains strict identity versus shared structural/neurogenic programs. The marker-validation layer is meant to prevent human hippocampal background cells from being mistaken for DG granule cells in the next integration step.",
        "",
        "## Outputs",
        "",
        f"- Per-cell scores and triage calls: `{OUT_SCORES.relative_to(ROOT)}`",
        f"- Panel summary: `{OUT_PANEL_SUMMARY.relative_to(ROOT)}`",
        f"- Call summary: `{OUT_CALL_SUMMARY.relative_to(ROOT)}`",
        f"- Gene coverage: `{OUT_GENE_COVERAGE.relative_to(ROOT)}`",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(text))


def main() -> None:
    panels = load_panels()
    matrix_to_gene_path = build_matrix_to_gene_path()
    qc = pd.read_csv(QC_CELLS, sep="\t", low_memory=False)

    score_parts: list[pd.DataFrame] = []
    coverage_parts: list[pd.DataFrame] = []

    for matrix_path, meta in qc.groupby("matrix_path", sort=False):
        matrix_abs = ROOT / matrix_path
        gene_abs = ROOT / matrix_to_gene_path[matrix_path]
        print(f"Scoring {matrix_path}", flush=True)
        matrix = sparse.load_npz(matrix_abs).tocsr()
        if matrix.shape[0] != len(meta):
            raise ValueError(f"Cell count mismatch for {matrix_path}: matrix {matrix.shape[0]}, metadata {len(meta)}")
        genes = read_gene_names(gene_abs)
        if matrix.shape[1] != len(genes):
            raise ValueError(f"Gene count mismatch for {matrix_path}: matrix {matrix.shape[1]}, genes {len(genes)}")

        marker_scores, coverage = score_matrix(matrix, genes, panels)
        marker_scores.index = meta.index
        component = pd.concat([meta.reset_index(drop=True), marker_scores.reset_index(drop=True)], axis=1)
        score_parts.append(component)

        coverage["dataset"] = meta["dataset"].iloc[0]
        coverage["component_id"] = meta["component_id"].iloc[0]
        coverage["matrix_path"] = matrix_path
        coverage_parts.append(coverage)

    scores = pd.concat(score_parts, ignore_index=True)
    scores = add_marker_calls(scores)
    panel_summary = summarize_panels(scores, panels)
    call_summary = summarize_calls(scores)
    gene_coverage = summarize_coverage(coverage_parts)

    scores.to_csv(OUT_SCORES, sep="\t", index=False, compression="gzip", float_format="%.6g")
    panel_summary.to_csv(OUT_PANEL_SUMMARY, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)
    call_summary.to_csv(OUT_CALL_SUMMARY, sep="\t", index=False, float_format="%.6g", quoting=csv.QUOTE_MINIMAL)
    gene_coverage.to_csv(OUT_GENE_COVERAGE, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    plot_panel_medians(panel_summary)
    write_interpretation(scores, call_summary, panel_summary)

    print(f"Wrote {OUT_SCORES}")
    print(f"Wrote {OUT_PANEL_SUMMARY}")
    print(f"Wrote {OUT_CALL_SUMMARY}")
    print(f"Wrote {OUT_GENE_COVERAGE}")
    print(f"Wrote {OUT_PLOT}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
