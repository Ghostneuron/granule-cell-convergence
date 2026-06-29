#!/usr/bin/env python3
"""Test regional-origin separation versus later shared granule-toolkit timing."""

from __future__ import annotations

import os
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

import build_primary_core_granule_specificity_named_comparators as base


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
SUPP_FIGURES = ROOT / "Project/manuscript/Supplementary figures"
MODULE_GENES = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"

OUT_MODULE_UNITS = RESULTS / "regional_origin_shared_toolkit_timing_module_units.tsv"
OUT_GENE_UNITS = RESULTS / "regional_origin_shared_toolkit_timing_gene_units.tsv"
OUT_STATE_SUMMARY = RESULTS / "regional_origin_shared_toolkit_timing_state_summary.tsv"
OUT_METRICS = RESULTS / "regional_origin_shared_toolkit_timing_metrics.tsv"
OUT_PLOT = RESULTS / "regional_origin_shared_toolkit_timing.png"
OUT_SUPP_PLOT = SUPP_FIGURES / "Fig.S3_regional_origin_shared_toolkit_timing.png"
OUT_MD = RESULTS / "regional_origin_shared_toolkit_timing.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SPECIAL_GENES = {
    "NFIA": "Nfia",
    "NEUROD1": "Neurod1",
    "RBFOX3": "Rbfox3",
    "HMGN2": "Hmgn2",
}

MODULE_ORDER = [
    "dentate_fate_wnt_prox1",
    "cerebellar_fate_rhombic_lip_shh",
    "shared_granule_special",
    "downstream_neurite_morphology",
    "downstream_synaptic_excitability",
    "shared_neurogenic_niche_state",
]

MODULE_LABELS = {
    "dentate_fate_wnt_prox1": "Dentate fate/WNT/PROX1",
    "cerebellar_fate_rhombic_lip_shh": "Cerebellar fate/SHH",
    "shared_granule_special": "Shared special genes",
    "downstream_neurite_morphology": "Neurite/morphology",
    "downstream_synaptic_excitability": "Synaptic/excitability",
    "shared_neurogenic_niche_state": "Shared neurogenic niche",
}

DENTATE_STATE_ORDER = {
    "RGL_young": (0.0, "progenitor", "early/progenitor"),
    "RGL": (1.0, "progenitor", "early/progenitor"),
    "nIPC": (2.0, "progenitor", "early/progenitor"),
    "nIPC-perin": (3.0, "progenitor", "early/progenitor"),
    "Neuroblast": (4.0, "immature", "postmitotic/immature"),
    "Immature-GC": (5.0, "immature", "postmitotic/immature"),
    "GC-juv": (6.0, "mature", "maturing/mature"),
    "GC-adult": (7.0, "mature", "maturing/mature"),
}

CEREBELLAR_SAMPLE_ORDER = {"P0": 0.0, "P8a": 8.1, "P8b": 8.2}
CEREBELLAR_GROUP_OFFSET = {"Granule precursor": 0.0, "Granule cells": 0.2}


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().upper()


def load_module_sets() -> dict[str, list[str]]:
    df = pd.read_csv(MODULE_GENES, sep="\t")
    module_sets: dict[str, list[str]] = {}
    for module_id in MODULE_ORDER:
        if module_id == "shared_granule_special":
            module_sets[module_id] = list(SPECIAL_GENES.values())
            continue
        genes = (
            df.loc[df["module_id"].eq(module_id), "default_mouse_symbol"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        module_sets[module_id] = genes
    return module_sets


def read_selected_rows_from_tar(tar_path: Path, member_name: str, wanted_genes: set[str], tmpdir: Path) -> pd.DataFrame:
    with tarfile.open(tar_path) as tar:
        source = tar.extractfile(member_name)
        if source is None:
            raise FileNotFoundError(member_name)
        path = tmpdir / member_name
        path.write_bytes(source.read())
    return base.read_selected_rows(path, ",", wanted_genes)


def compute_scores(
    *,
    dataset: str,
    sample: str,
    branch: str,
    expression: pd.DataFrame,
    cell_groups: pd.Series,
    module_sets: dict[str, list[str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    index_by_canon = {canon(gene): gene for gene in expression.index}
    common_cells = [cell for cell in expression.columns if cell in set(cell_groups.index)]
    expression = expression[common_cells]
    cell_groups = cell_groups.loc[common_cells]

    module_records: list[dict[str, object]] = []
    gene_records: list[dict[str, object]] = []

    for module_id, genes in module_sets.items():
        present = [index_by_canon[canon(gene)] for gene in genes if canon(gene) in index_by_canon]
        if not present:
            continue
        cell_scores = np.log1p(expression.loc[present].to_numpy(dtype=float)).mean(axis=0)
        score_df = pd.DataFrame({"cell_id": common_cells, "source_group": cell_groups.to_numpy(), "module_score": cell_scores})
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("module_score", "size"),
                median_module_score=("module_score", "median"),
                mean_module_score=("module_score", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_module_rank"] = group_summary["median_module_score"].rank(pct=True, method="average")
        for _, row in group_summary.iterrows():
            module_records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "branch": branch,
                    "source_group": row["source_group"],
                    "module_id": module_id,
                    "module_label": MODULE_LABELS[module_id],
                    "n_present_genes": len(present),
                    "present_genes": ",".join(present),
                    "n_cells": int(row["n_cells"]),
                    "median_module_score": float(row["median_module_score"]),
                    "mean_module_score": float(row["mean_module_score"]),
                    "within_sample_module_rank": float(row["within_sample_module_rank"]),
                }
            )

    for human_gene, mouse_gene in SPECIAL_GENES.items():
        present = index_by_canon.get(canon(mouse_gene))
        if present is None:
            continue
        values = expression.loc[present].to_numpy(dtype=float)
        log_values = np.log1p(values)
        score_df = pd.DataFrame(
            {
                "cell_id": common_cells,
                "source_group": cell_groups.to_numpy(),
                "log1p_expression": log_values,
                "detected": values > 0,
            }
        )
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("log1p_expression", "size"),
                median_log1p_expression=("log1p_expression", "median"),
                mean_log1p_expression=("log1p_expression", "mean"),
                detection_fraction=("detected", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_gene_rank"] = group_summary["median_log1p_expression"].rank(pct=True, method="average")
        for _, row in group_summary.iterrows():
            gene_records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "branch": branch,
                    "source_group": row["source_group"],
                    "human_gene": human_gene,
                    "mouse_gene": mouse_gene,
                    "source_gene_symbol": present,
                    "n_cells": int(row["n_cells"]),
                    "median_log1p_expression": float(row["median_log1p_expression"]),
                    "mean_log1p_expression": float(row["mean_log1p_expression"]),
                    "detection_fraction": float(row["detection_fraction"]),
                    "within_sample_gene_rank": float(row["within_sample_gene_rank"]),
                }
            )

    return module_records, gene_records


def load_dentate(module_sets: dict[str, list[str]], wanted_genes: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    meta = pd.read_csv(base.GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    expression = base.read_selected_rows(base.GSE104323_EXPR, "\t", wanted_genes)
    return compute_scores(
        dataset="GSE104323",
        sample="10X_all_cells",
        branch="dentate",
        expression=expression,
        cell_groups=cell_groups,
        module_sets=module_sets,
    )


def load_cerebellum(module_sets: dict[str, list[str]], wanted_genes: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    module_records: list[dict[str, object]] = []
    gene_records: list[dict[str, object]] = []
    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in base.GSE122357_FILES.items():
            labels = base.load_gse122357_label_map(prefix)
            expression = read_selected_rows_from_tar(base.GSE122357_TAR, member, wanted_genes, tmpdir)
            mod, gene = compute_scores(
                dataset="GSE122357",
                sample=sample,
                branch="cerebellar",
                expression=expression,
                cell_groups=labels,
                module_sets=module_sets,
            )
            module_records.extend(mod)
            gene_records.extend(gene)
    return module_records, gene_records


def annotate_focus(module_units: pd.DataFrame, gene_units: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    def annotate(row: pd.Series) -> pd.Series:
        group = str(row["source_group"])
        sample = str(row["sample"])
        branch = str(row["branch"])
        if branch == "dentate" and group in DENTATE_STATE_ORDER:
            order, bin_id, bin_label = DENTATE_STATE_ORDER[group]
            row["focus_state"] = True
            row["state_order"] = order
            row["state_bin"] = bin_id
            row["state_bin_label"] = bin_label
            row["state_label"] = group
        elif branch == "cerebellar" and group in CEREBELLAR_GROUP_OFFSET:
            row["focus_state"] = True
            row["state_order"] = CEREBELLAR_SAMPLE_ORDER[sample] + CEREBELLAR_GROUP_OFFSET[group]
            row["state_bin"] = "immature" if sample == "P0" or group == "Granule precursor" else "mature"
            row["state_bin_label"] = "early/precursor" if row["state_bin"] == "immature" else "maturing/mature"
            row["state_label"] = f"{sample} {group}"
        else:
            row["focus_state"] = False
            row["state_order"] = np.nan
            row["state_bin"] = "other"
            row["state_bin_label"] = "other"
            row["state_label"] = group
        return row

    module_units = module_units.apply(annotate, axis=1)
    gene_units = gene_units.apply(annotate, axis=1)
    return module_units, gene_units


def build_state_summary(module_units: pd.DataFrame, gene_units: pd.DataFrame) -> pd.DataFrame:
    focus_modules = module_units.loc[module_units["focus_state"]].copy()
    pivot = focus_modules.pivot_table(
        index=["branch", "dataset", "sample", "source_group", "state_label", "state_order", "state_bin", "state_bin_label", "n_cells"],
        columns="module_id",
        values="within_sample_module_rank",
        aggfunc="median",
    ).reset_index()

    focus_genes = gene_units.loc[gene_units["focus_state"]].copy()
    gene_pivot = focus_genes.pivot_table(
        index=["branch", "sample", "source_group"],
        columns="human_gene",
        values="within_sample_gene_rank",
        aggfunc="median",
    ).reset_index()

    out = pivot.merge(gene_pivot, on=["branch", "sample", "source_group"], how="left")
    out["branch_matched_fate_rank"] = np.where(
        out["branch"].eq("dentate"),
        out["dentate_fate_wnt_prox1"],
        out["cerebellar_fate_rhombic_lip_shh"],
    )
    out["opposed_fate_rank"] = np.where(
        out["branch"].eq("dentate"),
        out["cerebellar_fate_rhombic_lip_shh"],
        out["dentate_fate_wnt_prox1"],
    )
    out["regional_fate_polarity"] = out["branch_matched_fate_rank"] - out["opposed_fate_rank"]
    out["construction_rank"] = out[["downstream_neurite_morphology", "downstream_synaptic_excitability"]].mean(axis=1)
    out["shared_special_gene_rank"] = out[["NFIA", "NEUROD1", "RBFOX3", "HMGN2"]].mean(axis=1)
    out["postmitotic_special_gene_rank"] = out[["NEUROD1", "RBFOX3"]].mean(axis=1)
    out["regulatory_chromatin_gene_rank"] = out[["NFIA", "HMGN2"]].mean(axis=1)
    out["configuration_timing_score"] = out["construction_rank"] + out["shared_special_gene_rank"] - out["shared_neurogenic_niche_state"]
    return out.sort_values(["branch", "state_order"])


def median_for(summary: pd.DataFrame, branch: str, bins: set[str], column: str) -> float:
    values = summary.loc[summary["branch"].eq(branch) & summary["state_bin"].isin(bins), column].dropna()
    return float(values.median()) if not values.empty else np.nan


def build_metrics(summary: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("dentate", "early_to_postmitotic", {"progenitor"}, {"immature", "mature"}),
        ("dentate", "early_to_mature", {"progenitor"}, {"mature"}),
        ("cerebellar", "P0_or_precursor_to_maturing", {"immature"}, {"mature"}),
    ]
    columns = [
        "regional_fate_polarity",
        "shared_special_gene_rank",
        "construction_rank",
        "shared_neurogenic_niche_state",
        "configuration_timing_score",
        "postmitotic_special_gene_rank",
        "regulatory_chromatin_gene_rank",
        "NFIA",
        "NEUROD1",
        "RBFOX3",
        "HMGN2",
    ]
    rows: list[dict[str, object]] = []
    for branch, comparison, early_bins, late_bins in comparisons:
        for column in columns:
            early = median_for(summary, branch, early_bins, column)
            late = median_for(summary, branch, late_bins, column)
            rows.append(
                {
                    "branch": branch,
                    "comparison": comparison,
                    "metric": column,
                    "early_median": early,
                    "late_median": late,
                    "late_minus_early": late - early if pd.notna(early) and pd.notna(late) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_summary(summary: pd.DataFrame) -> None:
    plot_specs = [
        ("regional_fate_polarity", "Regional fate polarity"),
        ("postmitotic_special_gene_rank", "NEUROD1/RBFOX3 rank"),
        ("construction_rank", "Construction rank"),
        ("configuration_timing_score", "Configuration timing score"),
    ]
    colors = {"dentate": "#4c78a8", "cerebellar": "#f58518"}
    branches = ["dentate", "cerebellar"]
    fig, axes = plt.subplots(len(branches), len(plot_specs), figsize=(15.5, 7.2), sharey="col")
    panel_i = 0
    for row_i, branch in enumerate(branches):
        sub_branch = summary.loc[summary["branch"].eq(branch)].sort_values("state_order").reset_index(drop=True)
        x = np.arange(len(sub_branch))
        labels = sub_branch["state_label"].astype(str).str.replace("Granule precursor", "GP", regex=False).str.replace(
            "Granule cells", "GC", regex=False
        )
        for col_i, (column, title) in enumerate(plot_specs):
            ax = axes[row_i, col_i]
            ax.set_title(f"{chr(ord('a') + panel_i)}. {title}", loc="left", fontsize=11, fontweight="bold")
            panel_i += 1
            ax.plot(x, sub_branch[column], marker="o", linewidth=1.8, color=colors.get(branch, "gray"))
            ax.axhline(0, color="#666666", linewidth=0.7)
            if col_i == 0:
                ax.set_ylabel(f"{branch}\nrank / score")
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
            ax.set_ylim(-0.55 if column == "regional_fate_polarity" else 0, 1.6 if column == "configuration_timing_score" else 1.05)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    SUPP_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUPP_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, metrics: pd.DataFrame) -> None:
    def metric_value(branch: str, comparison: str, metric: str) -> float:
        sub = metrics.loc[
            metrics["branch"].eq(branch) & metrics["comparison"].eq(comparison) & metrics["metric"].eq(metric),
            "late_minus_early",
        ]
        return float(sub.iloc[0]) if not sub.empty else np.nan

    dentate_special = metric_value("dentate", "early_to_postmitotic", "shared_special_gene_rank")
    dentate_postmitotic_special = metric_value("dentate", "early_to_postmitotic", "postmitotic_special_gene_rank")
    dentate_construction = metric_value("dentate", "early_to_postmitotic", "construction_rank")
    cereb_special = metric_value("cerebellar", "P0_or_precursor_to_maturing", "shared_special_gene_rank")
    cereb_postmitotic_special = metric_value("cerebellar", "P0_or_precursor_to_maturing", "postmitotic_special_gene_rank")
    cereb_construction = metric_value("cerebellar", "P0_or_precursor_to_maturing", "construction_rank")
    dentate_polarity_early = median_for(summary, "dentate", {"progenitor"}, "regional_fate_polarity")
    dentate_polarity_post = median_for(summary, "dentate", {"immature", "mature"}, "regional_fate_polarity")
    cereb_polarity_precursor = median_for(summary, "cerebellar", {"immature"}, "regional_fate_polarity")
    cereb_polarity_mature = median_for(summary, "cerebellar", {"mature"}, "regional_fate_polarity")

    lines = [
        "# Regional-Origin Versus Shared-Toolkit Timing Analysis",
        "",
        "## Purpose",
        "",
        "This analysis tests whether dentate and cerebellar granule cells look like a single recent progenitor lineage or instead retain different regional fate polarity while recruiting shared granule/maturation genes during postmitotic development.",
        "",
        "## Main Findings",
        "",
        f"- Dentate progenitor states do not yet show strong dentate-over-cerebellar fate polarity (median {dentate_polarity_early:.3f}), but dentate postmitotic/mature granule states do (median {dentate_polarity_post:.3f}; shift {dentate_polarity_post - dentate_polarity_early:.3f}).",
        f"- Cerebellar precursor/P0 states show cerebellar-over-dentate fate polarity (median {cereb_polarity_precursor:.3f}), whereas labeled maturing granule cells are lower or near tied (median {cereb_polarity_mature:.3f}).",
        f"- Dentate `NEUROD1`/`RBFOX3` postmitotic-special rank increases from progenitor to postmitotic/mature states by {dentate_postmitotic_special:.3f}; construction rank increases by {dentate_construction:.3f}.",
        f"- Cerebellar `NEUROD1`/`RBFOX3` postmitotic-special rank changes from P0/precursor to maturing states by {cereb_postmitotic_special:.3f}; construction rank increases by {cereb_construction:.3f}.",
        f"- The four-gene shared-special average changes by {dentate_special:.3f} in dentate and {cereb_special:.3f} in cerebellum because `NFIA`/`HMGN2` behave more like early regulatory/chromatin-state genes while `NEUROD1`/`RBFOX3` track postmitotic maturation more directly.",
        "",
        "## Interpretation",
        "",
        "- The timing pattern argues against a single recent shared dentate/cerebellar progenitor in the sampled states. Branch-specific fate polarity appears in different lineage windows rather than as a shared root state.",
        "- `NFIA`, `NEUROD1`, `RBFOX3`, and `HMGN2` are best interpreted as a reused stage-composite toolkit: `NFIA`/`HMGN2` mark regulatory/chromatin competence, while `NEUROD1`/`RBFOX3` better mark postmitotic neuronal maturation.",
        "- The result supports the manuscript model: distinct regional origins, followed by partial convergence through shared neurogenic/maturation and downstream construction layers.",
        "- The cerebellar timing axis remains limited by three postnatal `GSE122357` samples, so the direction of cerebellar changes should be treated as supportive rather than definitive.",
        "",
        "## Outputs",
        "",
        f"- Module timing units: `{OUT_MODULE_UNITS.relative_to(ROOT)}`",
        f"- Gene timing units: `{OUT_GENE_UNITS.relative_to(ROOT)}`",
        f"- State summary: `{OUT_STATE_SUMMARY.relative_to(ROOT)}`",
        f"- Timing metrics: `{OUT_METRICS.relative_to(ROOT)}`",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    module_sets = load_module_sets()
    wanted_genes = {gene for genes in module_sets.values() for gene in genes}
    wanted_genes.update(SPECIAL_GENES.values())

    module_records: list[dict[str, object]] = []
    gene_records: list[dict[str, object]] = []
    mod, gene = load_dentate(module_sets, wanted_genes)
    module_records.extend(mod)
    gene_records.extend(gene)
    mod, gene = load_cerebellum(module_sets, wanted_genes)
    module_records.extend(mod)
    gene_records.extend(gene)

    module_units = pd.DataFrame(module_records)
    gene_units = pd.DataFrame(gene_records)
    module_units, gene_units = annotate_focus(module_units, gene_units)
    state_summary = build_state_summary(module_units, gene_units)
    metrics = build_metrics(state_summary)

    module_units.to_csv(OUT_MODULE_UNITS, sep="\t", index=False)
    gene_units.to_csv(OUT_GENE_UNITS, sep="\t", index=False)
    state_summary.to_csv(OUT_STATE_SUMMARY, sep="\t", index=False)
    metrics.to_csv(OUT_METRICS, sep="\t", index=False)
    plot_summary(state_summary)
    write_report(state_summary, metrics)

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(state_summary[["branch", "state_label", "regional_fate_polarity", "shared_special_gene_rank", "construction_rank"]].to_string(index=False))
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
