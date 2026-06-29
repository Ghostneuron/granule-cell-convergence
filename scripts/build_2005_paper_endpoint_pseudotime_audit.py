#!/usr/bin/env python3
"""Stage and pseudotime audit of the 2005 Development paper endpoints.

The 2005 paper tested conditioned-medium effects on proliferation, p21/p27,
neuronal differentiation markers, SMAD/TGF-beta signaling, BDNF/ERK signaling,
and apoptosis/survival. This script asks whether the same endpoint classes are
stage-dependent in the current primary-core sequencing resources.

This is an RNA trajectory audit. It cannot directly measure BrdU incorporation,
protein abundance, ERK phosphorylation, SMAD nuclear translocation, or TUNEL.
Those paper readouts are represented by transcript modules.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_primary_core_aim2b_stage_resolved_tgf_bdnf import (  # noqa: E402
    add_within_axis_gene_ranks,
    load_gse104323_lineage_units,
    load_gse122357_cerebellar_units,
    load_gse214309_state_units,
    load_gse292261_stage_units,
)


OUT_GENE_UNITS = RESULTS / "primary_core_2005_endpoint_pseudotime_gene_units.tsv"
OUT_MODULE_UNITS = RESULTS / "primary_core_2005_endpoint_pseudotime_module_units.tsv"
OUT_TRAJECTORY = RESULTS / "primary_core_2005_endpoint_pseudotime_trajectory_scores.tsv"
OUT_CORRELATIONS = RESULTS / "primary_core_2005_endpoint_pseudotime_correlations.tsv"
OUT_TRANSITIONS = RESULTS / "primary_core_2005_endpoint_pseudotime_transitions.tsv"
OUT_PLOT = RESULTS / "primary_core_2005_endpoint_pseudotime_audit.png"
OUT_MD = RESULTS / "primary_core_2005_endpoint_pseudotime_audit.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


PAPER_ENDPOINT_MODULES = [
    {
        "module_id": "proliferation_brdU_proxy",
        "module_label": "BrdU/proliferation proxy",
        "paper_readout": "BrdU incorporation and cell proliferation",
        "trajectory_role": "early_or_proliferative",
        "orientation": -1,
        "genes": [
            "MKI67",
            "PCNA",
            "TOP2A",
            "MCM2",
            "MCM3",
            "MCM4",
            "MCM5",
            "MCM6",
            "CCND1",
            "CCND2",
            "CCNE1",
            "CCNB1",
            "AURKB",
            "E2F1",
            "MYBL2",
        ],
    },
    {
        "module_id": "cell_cycle_exit_p21_p27",
        "module_label": "p21/p27 cell-cycle exit",
        "paper_readout": "p21 and p27 induction",
        "trajectory_role": "stop_or_maturation",
        "orientation": 1,
        "genes": ["CDKN1A", "CDKN1B", "CDKN1C", "BTG1", "BTG2", "TOB1", "TOB2", "GADD45A", "GADD45B"],
    },
    {
        "module_id": "immature_progenitor_state",
        "module_label": "Nestin/Mash1/Math1 immature state",
        "paper_readout": "Nestin, MASH1/ASCL1, MATH1/ATOH1, immature neuron context",
        "trajectory_role": "early_or_proliferative",
        "orientation": -1,
        "genes": ["NES", "SOX2", "HES1", "HES5", "ASCL1", "ATOH1", "DCX", "EOMES", "NEUROD1", "NEUROD2", "SOX11"],
    },
    {
        "module_id": "neuronal_differentiation_maturation",
        "module_label": "TuJ1/MAP2/synapsin maturation",
        "paper_readout": "TuJ1, MAP2a/b, synapsin, neuronal maturity/differentiation",
        "trajectory_role": "stop_or_maturation",
        "orientation": 1,
        "genes": ["TUBB3", "MAP2", "RBFOX3", "NEFL", "NEFM", "SYN1", "SYN2", "SYP", "SNAP25", "DLG4", "CAMK2A"],
    },
    {
        "module_id": "tgf_smad_pai1_response",
        "module_label": "TGF-beta/SMAD/PAI1 response",
        "paper_readout": "TGF-beta receptor, SMAD2/3/4, p3TP/PAI1 response",
        "trajectory_role": "stop_or_maturation",
        "orientation": 1,
        "genes": [
            "TGFB2",
            "TGFB3",
            "TGFBR1",
            "TGFBR2",
            "SMAD2",
            "SMAD3",
            "SMAD4",
            "SERPINE1",
            "TGFBI",
            "ID1",
            "ID2",
            "ID3",
            "CDKN1A",
            "CDKN1B",
        ],
    },
    {
        "module_id": "bdnf_erk_response",
        "module_label": "BDNF/TrkB/ERK response",
        "paper_readout": "BDNF, TRKB, ERK1/2, immediate early MAPK response",
        "trajectory_role": "context_dependent_maturation",
        "orientation": 1,
        "genes": [
            "BDNF",
            "NTRK2",
            "SHC1",
            "GRB2",
            "SOS1",
            "MAP2K1",
            "MAP2K2",
            "MAPK1",
            "MAPK3",
            "CREB1",
            "EGR1",
            "FOS",
            "JUN",
            "DUSP1",
            "DUSP6",
            "ELK1",
        ],
    },
    {
        "module_id": "apoptosis_execution",
        "module_label": "TUNEL/apoptosis proxy",
        "paper_readout": "TUNEL and Hoechst apoptotic morphology",
        "trajectory_role": "apoptosis_or_stress",
        "orientation": 0,
        "genes": ["BAX", "BAK1", "CASP3", "CASP7", "CASP8", "CASP9", "APAF1", "BID", "BAD", "BBC3", "PMAIP1"],
    },
    {
        "module_id": "survival_neuroprotection",
        "module_label": "Survival/neuroprotection proxy",
        "paper_readout": "neuroprotective conditioned-medium fractions",
        "trajectory_role": "survival_context",
        "orientation": 0,
        "genes": ["BCL2", "BCL2L1", "BCL2L2", "MCL1", "BDNF", "NTRK2", "AKT1", "AKT2", "CREB1"],
    },
    {
        "module_id": "secreted_stop_candidate_axis",
        "module_label": "Secreted stop-factor candidate axis",
        "paper_readout": "proteinaceous conditioned-medium stop/differentiation factors",
        "trajectory_role": "stop_or_maturation",
        "orientation": 1,
        "genes": [
            "TGFB2",
            "BDNF",
            "SFRP1",
            "SFRP2",
            "BMP6",
            "RELN",
            "GDF11",
            "INHBB",
            "TGFB3",
            "SEMA3A",
            "SLIT2",
            "NTF3",
            "CXCL12",
        ],
    },
]


PLOT_MODULES = [
    "proliferation_brdU_proxy",
    "cell_cycle_exit_p21_p27",
    "neuronal_differentiation_maturation",
    "tgf_smad_pai1_response",
    "bdnf_erk_response",
    "secreted_stop_candidate_axis",
]

MODULE_COLORS = {
    "proliferation_brdU_proxy": "#b2182b",
    "cell_cycle_exit_p21_p27": "#ef8a62",
    "neuronal_differentiation_maturation": "#2166ac",
    "tgf_smad_pai1_response": "#762a83",
    "bdnf_erk_response": "#1b7837",
    "secreted_stop_candidate_axis": "#8c510a",
}


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def finite_mean(values: pd.Series | list[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.mean(arr))


def finite_median(values: pd.Series | list[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def module_gene_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module in PAPER_ENDPOINT_MODULES:
        for order, gene in enumerate(module["genes"], start=1):
            rows.append(
                {
                    "module_id": module["module_id"],
                    "module_label": module["module_label"],
                    "paper_readout": module["paper_readout"],
                    "trajectory_role": module["trajectory_role"],
                    "orientation": module["orientation"],
                    "canonical_gene": canon(gene),
                    "gene_order": order,
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["module_id", "canonical_gene"])


def target_genes() -> set[str]:
    return set(module_gene_table()["canonical_gene"])


def build_gene_units() -> pd.DataFrame:
    targets = target_genes()
    pieces = [
        load_gse104323_lineage_units(targets),
        load_gse292261_stage_units(targets),
        load_gse122357_cerebellar_units(targets),
        load_gse214309_state_units(targets),
    ]
    gene_units = pd.concat([piece for piece in pieces if not piece.empty], ignore_index=True)
    gene_units = add_within_axis_gene_ranks(gene_units)
    gene_units = gene_units.sort_values(["dataset", "comparison_group", "axis_order", "canonical_gene"])
    return gene_units


def build_module_units(gene_units: pd.DataFrame, module_genes: pd.DataFrame) -> pd.DataFrame:
    expr = gene_units.merge(module_genes, on="canonical_gene", how="inner")
    group_cols = [
        "dataset",
        "species",
        "region",
        "axis_type",
        "axis_label",
        "axis_order",
        "comparison_group",
        "module_id",
        "module_label",
        "paper_readout",
        "trajectory_role",
        "orientation",
    ]
    units = (
        expr.groupby(group_cols, sort=False)
        .agg(
            n_cells=("n_cells", "max"),
            n_genes_present=("canonical_gene", "nunique"),
            genes_present=("canonical_gene", lambda values: ",".join(sorted(set(values)))),
            module_score=("gene_rank_within_axis", "median"),
            mean_module_rank=("gene_rank_within_axis", "mean"),
            median_detection_fraction=("detection_fraction", "median"),
            median_mean_log1p_expression=("mean_log1p_expression", "median"),
        )
        .reset_index()
    )
    defined = module_genes.groupby("module_id")["canonical_gene"].nunique().to_dict()
    units["n_genes_defined"] = units["module_id"].map(defined).astype(int)
    units["module_gene_coverage"] = units["n_genes_present"] / units["n_genes_defined"]
    return units.sort_values(["dataset", "comparison_group", "axis_order", "module_id"])


def add_assigned_pseudotime(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["assigned_pseudotime"] = np.nan
    for _, idx in out.groupby(["dataset", "comparison_group"], sort=False).groups.items():
        vals = pd.to_numeric(out.loc[idx, "axis_order"], errors="coerce")
        minv = vals.min()
        maxv = vals.max()
        if np.isfinite(minv) and np.isfinite(maxv) and maxv > minv:
            out.loc[idx, "assigned_pseudotime"] = (vals - minv) / (maxv - minv)
        else:
            out.loc[idx, "assigned_pseudotime"] = 0.0
    return out


def build_trajectory_scores(module_units: pd.DataFrame) -> pd.DataFrame:
    idx_cols = ["dataset", "species", "region", "axis_type", "axis_label", "axis_order", "comparison_group"]
    pivot = module_units.pivot_table(
        index=idx_cols,
        columns="module_id",
        values="module_score",
        aggfunc="median",
    ).reset_index()
    pivot.columns.name = None
    pivot = add_assigned_pseudotime(pivot)

    positive = [
        "cell_cycle_exit_p21_p27",
        "neuronal_differentiation_maturation",
        "tgf_smad_pai1_response",
        "bdnf_erk_response",
        "secreted_stop_candidate_axis",
    ]
    negative = ["proliferation_brdU_proxy", "immature_progenitor_state"]
    pivot["stop_maturation_raw"] = pivot[[c for c in positive if c in pivot]].mean(axis=1, skipna=True) - pivot[
        [c for c in negative if c in pivot]
    ].mean(axis=1, skipna=True)
    pivot["anti_proliferation_balance"] = pivot[
        [c for c in ["cell_cycle_exit_p21_p27", "tgf_smad_pai1_response", "secreted_stop_candidate_axis"] if c in pivot]
    ].mean(axis=1, skipna=True) - pivot[["proliferation_brdU_proxy"]].mean(axis=1, skipna=True)
    pivot["differentiation_over_progenitor"] = pivot[
        [c for c in ["neuronal_differentiation_maturation", "bdnf_erk_response"] if c in pivot]
    ].mean(axis=1, skipna=True) - pivot[["immature_progenitor_state"]].mean(axis=1, skipna=True)

    pivot["data_driven_maturation_pseudotime"] = np.nan
    for _, idx in pivot.groupby(["dataset", "comparison_group"], sort=False).groups.items():
        vals = pd.to_numeric(pivot.loc[idx, "stop_maturation_raw"], errors="coerce")
        minv = vals.min()
        maxv = vals.max()
        if np.isfinite(minv) and np.isfinite(maxv) and maxv > minv:
            pivot.loc[idx, "data_driven_maturation_pseudotime"] = (vals - minv) / (maxv - minv)
        else:
            pivot.loc[idx, "data_driven_maturation_pseudotime"] = 0.5
    return pivot.sort_values(["dataset", "comparison_group", "axis_order"])


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    xy = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(xy) < 3 or xy["x"].nunique() < 2 or xy["y"].nunique() < 2:
        return np.nan, np.nan, int(len(xy))
    rho, p = stats.spearmanr(xy["x"], xy["y"])
    return float(rho), float(p), int(len(xy))


def build_correlations(module_units: pd.DataFrame, trajectory: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    module_units = add_assigned_pseudotime(module_units)
    for (dataset, group, module_id), sub in module_units.groupby(["dataset", "comparison_group", "module_id"], sort=False):
        rho, p, n = spearman(sub["assigned_pseudotime"], sub["module_score"])
        meta = sub.iloc[0]
        rows.append(
            {
                "dataset": dataset,
                "region": meta["region"],
                "axis_type": meta["axis_type"],
                "comparison_group": group,
                "metric_type": "module",
                "metric_id": module_id,
                "metric_label": meta["module_label"],
                "paper_readout": meta["paper_readout"],
                "trajectory_role": meta["trajectory_role"],
                "orientation": meta["orientation"],
                "spearman_rho_vs_assigned_pseudotime": rho,
                "spearman_p": p,
                "n_axis_units": n,
            }
        )
    for (dataset, group), sub in trajectory.groupby(["dataset", "comparison_group"], sort=False):
        meta = sub.iloc[0]
        for metric_id, label in [
            ("data_driven_maturation_pseudotime", "data-driven maturation pseudotime"),
            ("anti_proliferation_balance", "anti-proliferation balance"),
            ("differentiation_over_progenitor", "differentiation over progenitor"),
        ]:
            rho, p, n = spearman(sub["assigned_pseudotime"], sub[metric_id])
            rows.append(
                {
                    "dataset": dataset,
                    "region": meta["region"],
                    "axis_type": meta["axis_type"],
                    "comparison_group": group,
                    "metric_type": "trajectory",
                    "metric_id": metric_id,
                    "metric_label": label,
                    "paper_readout": "composite RNA trajectory score",
                    "trajectory_role": "composite",
                    "orientation": 1,
                    "spearman_rho_vs_assigned_pseudotime": rho,
                    "spearman_p": p,
                    "n_axis_units": n,
                }
            )
    return pd.DataFrame(rows).sort_values(["dataset", "comparison_group", "metric_type", "metric_id"])


def build_transitions(module_units: pd.DataFrame, trajectory: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (dataset, group, module_id), sub in module_units.groupby(["dataset", "comparison_group", "module_id"], sort=False):
        sub = sub.sort_values("axis_order")
        if len(sub) < 2:
            continue
        first = sub.iloc[0]
        last = sub.iloc[-1]
        rows.append(
            {
                "dataset": dataset,
                "region": first["region"],
                "axis_type": first["axis_type"],
                "comparison_group": group,
                "metric_type": "module",
                "metric_id": module_id,
                "metric_label": first["module_label"],
                "paper_readout": first["paper_readout"],
                "start_label": first["axis_label"],
                "end_label": last["axis_label"],
                "start_score": float(first["module_score"]),
                "end_score": float(last["module_score"]),
                "delta_end_minus_start": float(last["module_score"] - first["module_score"]),
            }
        )
    for (dataset, group), sub in trajectory.groupby(["dataset", "comparison_group"], sort=False):
        sub = sub.sort_values("axis_order")
        if len(sub) < 2:
            continue
        first = sub.iloc[0]
        last = sub.iloc[-1]
        for metric_id, label in [
            ("data_driven_maturation_pseudotime", "data-driven maturation pseudotime"),
            ("anti_proliferation_balance", "anti-proliferation balance"),
            ("differentiation_over_progenitor", "differentiation over progenitor"),
        ]:
            rows.append(
                {
                    "dataset": dataset,
                    "region": first["region"],
                    "axis_type": first["axis_type"],
                    "comparison_group": group,
                    "metric_type": "trajectory",
                    "metric_id": metric_id,
                    "metric_label": label,
                    "paper_readout": "composite RNA trajectory score",
                    "start_label": first["axis_label"],
                    "end_label": last["axis_label"],
                    "start_score": float(first[metric_id]),
                    "end_score": float(last[metric_id]),
                    "delta_end_minus_start": float(last[metric_id] - first[metric_id]),
                }
            )
    return pd.DataFrame(rows).sort_values(["dataset", "comparison_group", "metric_type", "metric_id"])


def plot_trajectories(module_units: pd.DataFrame, trajectory: pd.DataFrame) -> None:
    panels = [
        ("GSE104323", "curated_dentate_granule_lineage", "Adult dentate lineage"),
        ("GSE292261", "candidate_dentate_granule_only", "Postnatal dentate candidates"),
        ("GSE122357", "candidate_cerebellar_granule_only", "Postnatal cerebellar candidates"),
        ("GSE214309", "adult_DGC_maturation_activity_state", "Adult dentate maturation/activity"),
    ]
    label_map = {
        "RGL_young": "RGL\nyoung",
        "Immature-GC": "Immature\nGC",
        "GC-juv": "GC\njuv",
        "GC-adult": "GC\nadult",
        "immature_1hr": "immature\n1 hr",
        "immatureactive_1hr": "imm active\n1 hr",
        "mature_1hr": "mature\n1 hr",
        "matureactive_1hr": "mat active\n1 hr",
        "immature_4hr": "immature\n4 hr",
        "immatureactive_4hr": "imm active\n4 hr",
        "mature_4hr": "mature\n4 hr",
        "matureactive_4hr": "mat active\n4 hr",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for ax, (dataset, group, title) in zip(axes, panels):
        sub = module_units.loc[module_units["dataset"].eq(dataset) & module_units["comparison_group"].eq(group)].copy()
        if sub.empty:
            ax.axis("off")
            continue
        sub = sub.sort_values("axis_order")
        for module_id in PLOT_MODULES:
            ms = sub.loc[sub["module_id"].eq(module_id)].copy()
            if ms.empty:
                continue
            x_labels = [label_map.get(str(label), str(label)) for label in ms["axis_label"]]
            ax.plot(
                x_labels,
                ms["module_score"],
                marker="o",
                linewidth=1.8,
                markersize=3.5,
                color=MODULE_COLORS[module_id],
                label=ms["module_label"].iloc[0],
            )
        ts = trajectory.loc[trajectory["dataset"].eq(dataset) & trajectory["comparison_group"].eq(group)].sort_values(
            "axis_order"
        )
        if not ts.empty:
            x_labels = [label_map.get(str(label), str(label)) for label in ts["axis_label"]]
            ax.plot(
                x_labels,
                ts["data_driven_maturation_pseudotime"],
                color="#000000",
                linewidth=2.5,
                marker="s",
                markersize=3.5,
                label="RNA maturation pseudotime",
            )
        ax.axhline(0.5, color="#777777", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.set_title(title)
        ax.set_ylim(-0.05, 1.05)
        ax.set_ylabel("Rank score / pseudotime")
        ax.tick_params(axis="x", rotation=0, labelsize=9)
        ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle("2005 paper endpoint modules across stage and RNA pseudotime", fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.17, hspace=0.42, wspace=0.22)
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


def corr_line(correlations: pd.DataFrame, dataset: str, group: str, metric_id: str) -> str:
    sub = correlations.loc[
        correlations["dataset"].eq(dataset)
        & correlations["comparison_group"].eq(group)
        & correlations["metric_id"].eq(metric_id)
    ]
    if sub.empty:
        return "NA"
    row = sub.iloc[0]
    return f"rho {fmt(row['spearman_rho_vs_assigned_pseudotime'])}, p {fmt(row['spearman_p'])}"


def transition_delta(transitions: pd.DataFrame, dataset: str, group: str, metric_id: str) -> str:
    sub = transitions.loc[
        transitions["dataset"].eq(dataset)
        & transitions["comparison_group"].eq(group)
        & transitions["metric_id"].eq(metric_id)
    ]
    if sub.empty:
        return "NA"
    row = sub.iloc[0]
    return f"{row['start_label']} to {row['end_label']}: {fmt(row['delta_end_minus_start'])}"


def write_markdown(
    gene_units: pd.DataFrame,
    module_units: pd.DataFrame,
    trajectory: pd.DataFrame,
    correlations: pd.DataFrame,
    transitions: pd.DataFrame,
) -> None:
    lines = [
        "# 2005 Paper Endpoint Stage/Pseudotime Audit",
        "",
        "## Question",
        "",
        "Are the other experimental readouts in the 2005 Development paper stage-dependent, and do they support a trajectory/pseudotime interpretation?",
        "",
        "## Answer",
        "",
        "Yes. The paper's readouts are better interpreted as trajectory-dependent than as one static factor effect. RNA modules corresponding to BrdU/cell-cycle, p21/p27 cell-cycle exit, neuronal differentiation/maturation, TGF-beta/SMAD/PAI1 response, BDNF/ERK response, apoptosis/survival, and secreted stop-factor candidates vary across ordered dentate and cerebellar axes.",
        "",
        "Important limitation: this analysis uses RNA proxies. It cannot directly measure BrdU incorporation, p21/p27 protein, MAP2/synapsin protein, ERK phosphorylation, SMAD nuclear translocation, PAI1 reporter activity, or apoptosis by TUNEL/Hoechst.",
        "",
        "## Key Stage-Dependence Results",
        "",
        f"- Adult dentate lineage (`GSE104323`): RNA maturation pseudotime versus lineage order is {corr_line(correlations, 'GSE104323', 'curated_dentate_granule_lineage', 'data_driven_maturation_pseudotime')}. Proliferation delta is {transition_delta(transitions, 'GSE104323', 'curated_dentate_granule_lineage', 'proliferation_brdU_proxy')}; differentiation delta is {transition_delta(transitions, 'GSE104323', 'curated_dentate_granule_lineage', 'neuronal_differentiation_maturation')}.",
        f"- Postnatal dentate candidates (`GSE292261`): RNA maturation pseudotime versus age is {corr_line(correlations, 'GSE292261', 'candidate_dentate_granule_only', 'data_driven_maturation_pseudotime')}. Proliferation delta is {transition_delta(transitions, 'GSE292261', 'candidate_dentate_granule_only', 'proliferation_brdU_proxy')}; TGF/SMAD delta is {transition_delta(transitions, 'GSE292261', 'candidate_dentate_granule_only', 'tgf_smad_pai1_response')}.",
        f"- Postnatal cerebellar candidates (`GSE122357`): RNA maturation pseudotime versus age is {corr_line(correlations, 'GSE122357', 'candidate_cerebellar_granule_only', 'data_driven_maturation_pseudotime')}. Proliferation delta is {transition_delta(transitions, 'GSE122357', 'candidate_cerebellar_granule_only', 'proliferation_brdU_proxy')}; TGF/SMAD delta is {transition_delta(transitions, 'GSE122357', 'candidate_cerebellar_granule_only', 'tgf_smad_pai1_response')}.",
        f"- Adult dentate activity/maturation (`GSE214309`): RNA maturation pseudotime versus ordered state is {corr_line(correlations, 'GSE214309', 'adult_DGC_maturation_activity_state', 'data_driven_maturation_pseudotime')}. This axis is activity/time-state ordered, not a clean developmental lineage.",
        "",
        "## Interpretation For The 2005 Paper",
        "",
        "- The anti-proliferation result maps to a trajectory shift: proliferative modules and immature/progenitor modules tend to separate from cell-cycle-exit, maturation, and stop-factor modules.",
        "- The p21/p27 result should be treated as a stage-sensitive cell-cycle-exit module, not only as an acute response marker.",
        "- The MAP2/TuJ1/synapsin differentiation result is strongly compatible with a maturation trajectory.",
        "- The TGF-beta/SMAD and BDNF/ERK results are partly stage-dependent, but RNA cannot substitute for pSMAD nuclear translocation or pERK signaling assays.",
        "- The apoptosis/neuroprotection fraction in the paper is the least safely inferred from RNA alone and should remain an experimental validation target.",
        "",
        "## Why True Pseudotime Is Needed Next",
        "",
        "This audit uses ordered stages and module-defined RNA pseudotime. A stronger next step is cell-level graph pseudotime, ideally using `GSE104323` dentate lineage cells and `GSE292261` postnatal dentate cells first, then `GSE122357` cerebellar P0/P8 cells. That would test whether TGF/BDNF/secretome modules rise at a specific trajectory window after proliferation and before mature neuronal markers.",
        "",
        "## Outputs",
        "",
        f"- Gene units: `{OUT_GENE_UNITS.relative_to(ROOT)}` ({len(gene_units):,} rows).",
        f"- Module units: `{OUT_MODULE_UNITS.relative_to(ROOT)}` ({len(module_units):,} rows).",
        f"- Trajectory scores: `{OUT_TRAJECTORY.relative_to(ROOT)}` ({len(trajectory):,} rows).",
        f"- Correlations: `{OUT_CORRELATIONS.relative_to(ROOT)}` ({len(correlations):,} rows).",
        f"- Transitions: `{OUT_TRANSITIONS.relative_to(ROOT)}` ({len(transitions):,} rows).",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    module_genes = module_gene_table()
    gene_units = build_gene_units()
    module_units = build_module_units(gene_units, module_genes)
    trajectory = build_trajectory_scores(module_units)
    correlations = build_correlations(module_units, trajectory)
    transitions = build_transitions(module_units, trajectory)

    gene_units.to_csv(OUT_GENE_UNITS, sep="\t", index=False)
    module_units.to_csv(OUT_MODULE_UNITS, sep="\t", index=False)
    trajectory.to_csv(OUT_TRAJECTORY, sep="\t", index=False)
    correlations.to_csv(OUT_CORRELATIONS, sep="\t", index=False)
    transitions.to_csv(OUT_TRANSITIONS, sep="\t", index=False)
    plot_trajectories(module_units, trajectory)
    write_markdown(gene_units, module_units, trajectory, correlations, transitions)

    print(f"Wrote {OUT_MD}")
    print(f"Gene units: {len(gene_units):,}")
    print(f"Module units: {len(module_units):,}")
    print(f"Trajectory rows: {len(trajectory):,}")
    print("Selected correlations:")
    keep = correlations.loc[
        correlations["metric_id"].isin(
            [
                "data_driven_maturation_pseudotime",
                "proliferation_brdU_proxy",
                "cell_cycle_exit_p21_p27",
                "neuronal_differentiation_maturation",
                "tgf_smad_pai1_response",
            ]
        )
    ]
    print(
        keep[
            [
                "dataset",
                "comparison_group",
                "metric_id",
                "spearman_rho_vs_assigned_pseudotime",
                "spearman_p",
                "n_axis_units",
            ]
        ]
        .round(3)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
