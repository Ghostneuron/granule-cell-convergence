#!/usr/bin/env python3
"""Build a first candidate-gene packet from the constructed human bridge.

This is deliberately conservative: it summarizes refined panel genes in the
human dentate/hippocampal bridge objects and links that evidence to the current
rank-level backbone integration. It is not a full differential-expression or
cross-species ortholog model.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from integrate_human_bridge_with_backbone import (
    broad_class,
    classify_gse268609,
    classify_gse325391,
    classify_human_core,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "Project/config"
RESULTS = ROOT / "Project/results"
PROCESSED = ROOT / "Project/processed"

PANELS = CONFIG / "granule_marker_panels_refined.tsv"

HUMAN_CORE_X = PROCESSED / "human_core_normalized_reduced_object/X_log1p_cp10k_selected_genes.npz"
HUMAN_CORE_OBS = PROCESSED / "human_core_normalized_reduced_object/obs.tsv.gz"
HUMAN_CORE_VAR = PROCESSED / "human_core_normalized_reduced_object/var.tsv"
HUMAN_CORE_LABELS = RESULTS / "human_core_tuned_labels.tsv.gz"

GSE325391_X = PROCESSED / "gse325391_adult_dg_selected/matrix_cells_by_selected_genes.npz"
GSE325391_OBS = PROCESSED / "gse325391_adult_dg_selected/cell_metadata.tsv.gz"
GSE325391_VAR = PROCESSED / "gse325391_adult_dg_selected/var_selected_features.tsv"
GSE325391_LABELS = RESULTS / "gse325391_human_core_label_projection.tsv.gz"

GSE268609_X = PROCESSED / "gse268609_rna_selected/matrix_cells_by_selected_genes.npz"
GSE268609_OBS = PROCESSED / "gse268609_rna_selected/cell_metadata.tsv.gz"
GSE268609_VAR = PROCESSED / "gse268609_rna_selected/var_selected_features.tsv"
GSE268609_LABELS = RESULTS / "gse268609_human_core_label_projection.tsv.gz"

INTEGRATED_STATS = RESULTS / "human_bridge_backbone_rank_statistics.tsv"
INTEGRATED_SUMMARY = RESULTS / "human_bridge_backbone_rank_source_summary.tsv"

OUT_GENE_SUMMARY = RESULTS / "human_bridge_marker_gene_summary.tsv"
OUT_PACKET = RESULTS / "human_bridge_candidate_gene_packet.tsv"
OUT_STRUCTURAL_PLOT = RESULTS / "human_bridge_structural_executor_candidates.png"
OUT_MD = RESULTS / "human_bridge_candidate_gene_packet.md"

MIN_DETECTION_SUPPORT = 0.05
HIGH_DETECTION_SUPPORT = 0.20


def gene_key(gene: object) -> str:
    return str(gene).strip().upper()


def load_panels() -> pd.DataFrame:
    panels = pd.read_csv(PANELS, sep="\t")
    panels["gene_key"] = panels["gene"].map(gene_key)
    panels["candidate_role"] = panels["panel"].map(
        {
            "dentate_identity": "regional_dentate_identity",
            "cerebellar_identity": "regional_cerebellar_identity",
            "shared_granule_neuronal": "shared_granule_neuronal_state",
            "morphogenesis_cytoskeleton": "shared_structural_executor",
            "axon_guidance_synapse": "shared_structural_executor",
            "metabolic_proteomic_validation": "supporting_metabolic_validation",
        }
    ).fillna("supporting_or_unclassified")
    return panels


def gene_indices(var: pd.DataFrame, panels: pd.DataFrame) -> tuple[list[int], pd.DataFrame]:
    var = var.copy()
    gene_col = "gene" if "gene" in var.columns else "feature_name"
    var["gene_key"] = var[gene_col].map(gene_key)
    first_index = {}
    for idx, key in enumerate(var["gene_key"]):
        if key and key not in first_index:
            first_index[key] = idx
    panel_rows = []
    indices = []
    for _, row in panels.iterrows():
        idx = first_index.get(row["gene_key"])
        out = row.to_dict()
        out["matrix_col_index"] = idx
        out["present_in_matrix"] = idx is not None
        panel_rows.append(out)
        if idx is not None:
            indices.append(idx)
    return indices, pd.DataFrame(panel_rows)


def matrix_group_summary(
    *,
    source_layer: str,
    expression_scale: str,
    X_path: Path,
    var_path: Path,
    meta: pd.DataFrame,
    panels: pd.DataFrame,
    already_log_scaled: bool,
) -> pd.DataFrame:
    X = sparse.load_npz(X_path)
    var = pd.read_csv(var_path, sep="\t")
    indices, panel_presence = gene_indices(var, panels)
    if X.shape[0] != len(meta):
        raise ValueError(f"{source_layer}: matrix rows {X.shape[0]} != metadata rows {len(meta)}")

    present = panel_presence.loc[panel_presence["present_in_matrix"]].copy()
    present = present.reset_index(drop=True)
    if not len(present):
        return pd.DataFrame()

    X_panel = X[:, present["matrix_col_index"].astype(int).to_numpy()].tocsr()
    rows = []
    keep = meta["analysis_include"].astype(bool).to_numpy()
    for broad in sorted(meta.loc[keep, "broad_class"].dropna().unique()):
        mask = keep & meta["broad_class"].eq(broad).to_numpy()
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        sub = X_panel[mask, :].tocsc()
        nonzero = np.diff(sub.indptr)
        value_sum = np.asarray(sub.sum(axis=0)).ravel()
        if already_log_scaled:
            log_sum = value_sum
        else:
            log_sub = sub.copy()
            log_sub.data = np.log1p(log_sub.data)
            log_sum = np.asarray(log_sub.sum(axis=0)).ravel()
        for i, panel_row in present.iterrows():
            rows.append(
                {
                    "source_layer": source_layer,
                    "expression_scale": expression_scale,
                    "broad_class": broad,
                    "panel": panel_row["panel"],
                    "gene": panel_row["gene"],
                    "gene_key": panel_row["gene_key"],
                    "candidate_role": panel_row["candidate_role"],
                    "notes": panel_row["notes"],
                    "n_cells": n_cells,
                    "detection_fraction": float(nonzero[i] / n_cells),
                    "mean_value_all_cells": float(value_sum[i] / n_cells),
                    "mean_log1p_value_all_cells": float(log_sum[i] / n_cells),
                    "present_in_matrix": True,
                }
            )

    absent = panel_presence.loc[~panel_presence["present_in_matrix"]].copy()
    for _, row in absent.iterrows():
        rows.append(
            {
                "source_layer": source_layer,
                "expression_scale": expression_scale,
                "broad_class": "not_present",
                "panel": row["panel"],
                "gene": row["gene"],
                "gene_key": row["gene_key"],
                "candidate_role": row["candidate_role"],
                "notes": row["notes"],
                "n_cells": 0,
                "detection_fraction": 0.0,
                "mean_value_all_cells": 0.0,
                "mean_log1p_value_all_cells": 0.0,
                "present_in_matrix": False,
            }
        )
    return pd.DataFrame(rows)


def human_core_meta() -> pd.DataFrame:
    obs = pd.read_csv(HUMAN_CORE_OBS, sep="\t", usecols=["cell_id"])
    labels = pd.read_csv(
        HUMAN_CORE_LABELS,
        sep="\t",
        usecols=["cell_id", "dataset", "replicate_unit", "tuned_label"],
        low_memory=False,
    )
    meta = obs.merge(labels, on="cell_id", how="left", sort=False)
    meta["analysis_class"] = meta["tuned_label"].map(classify_human_core)
    meta["broad_class"] = meta["analysis_class"].map(broad_class)
    meta["analysis_include"] = meta["tuned_label"].notna()
    return meta


def gse325391_meta() -> pd.DataFrame:
    obs = pd.read_csv(GSE325391_OBS, sep="\t", usecols=["cell_id"])
    labels = pd.read_csv(
        GSE325391_LABELS,
        sep="\t",
        usecols=["cell_id", "sample", "group", "source_anchor_label", "analysis_include"],
        low_memory=False,
    )
    meta = obs.merge(labels, on="cell_id", how="left", sort=False)
    meta["analysis_class"] = meta["source_anchor_label"].map(classify_gse325391)
    meta["broad_class"] = meta["analysis_class"].map(broad_class)
    meta["analysis_include"] = meta["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    return meta


def gse268609_meta() -> pd.DataFrame:
    obs = pd.read_csv(GSE268609_OBS, sep="\t", usecols=["cell_id"])
    labels = pd.read_csv(
        GSE268609_LABELS,
        sep="\t",
        usecols=["cell_id", "sample_id", "diagnosis", "projected_label", "analysis_include"],
        low_memory=False,
    )
    meta = obs.merge(labels, on="cell_id", how="left", sort=False)
    meta["analysis_class"] = meta["projected_label"].map(classify_gse268609)
    meta["broad_class"] = meta["analysis_class"].map(broad_class)
    meta["analysis_include"] = meta["analysis_include"].astype(str).str.lower().isin(["true", "1", "yes"])
    return meta


def build_gene_summary(panels: pd.DataFrame) -> pd.DataFrame:
    chunks = [
        matrix_group_summary(
            source_layer="human_core_tuned",
            expression_scale="log1p_cp10k_selected_genes",
            X_path=HUMAN_CORE_X,
            var_path=HUMAN_CORE_VAR,
            meta=human_core_meta(),
            panels=panels,
            already_log_scaled=True,
        ),
        matrix_group_summary(
            source_layer="gse325391_adult_dg",
            expression_scale="raw_selected_counts_with_log1p_summary",
            X_path=GSE325391_X,
            var_path=GSE325391_VAR,
            meta=gse325391_meta(),
            panels=panels,
            already_log_scaled=False,
        ),
        matrix_group_summary(
            source_layer="gse268609_hippocampus_rna",
            expression_scale="raw_selected_counts_with_log1p_summary",
            X_path=GSE268609_X,
            var_path=GSE268609_VAR,
            meta=gse268609_meta(),
            panels=panels,
            already_log_scaled=False,
        ),
    ]
    out = pd.concat(chunks, ignore_index=True, sort=False)
    return out.sort_values(["candidate_role", "panel", "gene", "source_layer", "broad_class"])


def tier_for_gene(row: pd.Series) -> str:
    role = row["candidate_role"]
    detected = int(row["dentate_sources_detected_at_5pct"])
    median_det = float(row["dentate_candidate_median_detection"])
    background_det = float(row["background_or_ambiguous_median_detection"])
    if role == "shared_structural_executor":
        if detected >= 3 and median_det >= HIGH_DETECTION_SUPPORT:
            return "high_priority_shared_structural_executor"
        if detected >= 2 and median_det >= MIN_DETECTION_SUPPORT:
            return "moderate_priority_shared_structural_executor"
        return "low_current_human_bridge_support"
    if role == "regional_dentate_identity":
        if detected >= 2 and median_det >= MIN_DETECTION_SUPPORT and median_det >= background_det:
            return "dentate_identity_supported_in_human_bridge"
        return "dentate_identity_panel_gene_needs_context"
    if role == "regional_cerebellar_identity":
        if median_det < MIN_DETECTION_SUPPORT:
            return "clean_low_human_dentate_detection_cerebellar_identity"
        return "cerebellar_identity_gene_with_human_bridge_leakage_warning"
    if role == "shared_granule_neuronal_state":
        if detected >= 3 and median_det >= HIGH_DETECTION_SUPPORT:
            return "broad_shared_neuronal_state_high_detection"
        return "shared_neuronal_state_context_dependent"
    return "supporting_or_contextual"


def build_packet(gene_summary: pd.DataFrame, panels: pd.DataFrame) -> pd.DataFrame:
    present = gene_summary.loc[gene_summary["present_in_matrix"]].copy()
    dentate = present.loc[present["broad_class"].eq("dentate_candidate")]
    background = present.loc[present["broad_class"].isin(["non_dentate_background", "other_or_ambiguous"])]
    warning = present.loc[present["broad_class"].eq("broad_neuronal_structural_warning")]

    rows = []
    for _, panel_row in panels.drop_duplicates("gene_key").iterrows():
        gene = panel_row["gene"]
        key = panel_row["gene_key"]
        d = dentate.loc[dentate["gene_key"].eq(key)]
        b = background.loc[background["gene_key"].eq(key)]
        w = warning.loc[warning["gene_key"].eq(key)]
        source_detection = d.groupby("source_layer")["detection_fraction"].max() if len(d) else pd.Series(dtype=float)
        rows.append(
            {
                "gene": gene,
                "panel": panel_row["panel"],
                "candidate_role": panel_row["candidate_role"],
                "notes": panel_row["notes"],
                "human_bridge_sources_present": int(present.loc[present["gene_key"].eq(key), "source_layer"].nunique()),
                "dentate_sources_detected_at_5pct": int((source_detection >= MIN_DETECTION_SUPPORT).sum()),
                "dentate_sources_detected_at_20pct": int((source_detection >= HIGH_DETECTION_SUPPORT).sum()),
                "dentate_candidate_median_detection": float(d["detection_fraction"].median()) if len(d) else 0.0,
                "dentate_candidate_max_detection": float(d["detection_fraction"].max()) if len(d) else 0.0,
                "dentate_candidate_median_log1p_value": float(d["mean_log1p_value_all_cells"].median()) if len(d) else 0.0,
                "background_or_ambiguous_median_detection": float(b["detection_fraction"].median()) if len(b) else 0.0,
                "broad_warning_median_detection": float(w["detection_fraction"].median()) if len(w) else 0.0,
            }
        )

    out = pd.DataFrame(rows)
    out["support_tier"] = out.apply(tier_for_gene, axis=1)
    out["dentate_vs_background_detection_delta"] = (
        out["dentate_candidate_median_detection"] - out["background_or_ambiguous_median_detection"]
    )
    role_order = {
        "shared_structural_executor": 0,
        "regional_dentate_identity": 1,
        "regional_cerebellar_identity": 2,
        "shared_granule_neuronal_state": 3,
        "supporting_metabolic_validation": 4,
    }
    out["_role_order"] = out["candidate_role"].map(role_order).fillna(9)
    out = out.sort_values(
        [
            "_role_order",
            "dentate_sources_detected_at_20pct",
            "dentate_sources_detected_at_5pct",
            "dentate_candidate_median_detection",
            "gene",
        ],
        ascending=[True, False, False, False, True],
    ).drop(columns=["_role_order"])
    return out


def write_md(packet: pd.DataFrame, gene_summary: pd.DataFrame) -> None:
    stats = pd.read_csv(INTEGRATED_STATS, sep="\t") if INTEGRATED_STATS.exists() else pd.DataFrame()
    source = pd.read_csv(INTEGRATED_SUMMARY, sep="\t") if INTEGRATED_SUMMARY.exists() else pd.DataFrame()
    high_structural = packet.loc[packet["support_tier"].eq("high_priority_shared_structural_executor")].head(20)
    dentate_identity = packet.loc[packet["candidate_role"].eq("regional_dentate_identity")].head(12)
    cerebellar_identity = packet.loc[packet["candidate_role"].eq("regional_cerebellar_identity")].head(12)
    rows = [
        "# Human Bridge Candidate Gene Packet",
        "",
        "Date built: 2026-06-21",
        "",
        "## Scope",
        "",
        "This packet summarizes refined marker-panel genes across the constructed human dentate/hippocampal bridge objects: `human_core_tuned`, `GSE325391`, and `GSE268609`.",
        "",
        "It is a first candidate table for manuscript planning. It does not replace full differential-expression, ortholog-aware cross-species modeling, or source-taxonomy refinement for `GSE268609`.",
        "",
        "## Integrated Context",
        "",
    ]
    if not stats.empty:
        keep = stats.loc[stats["comparison"].isin(["dentate_candidate_vs_cerebellar_candidate", "dentate_candidate_vs_non_dentate_background"])]
        for _, row in keep.iterrows():
            rows.append(
                f"- `{row['comparison']}` / `{row['metric']}`: median delta {row['delta_a_minus_b']:.4f}, BH-adjusted p={row['p_adj_bh']:.3g}."
            )
    if not source.empty:
        structural = source.loc[source["broad_class"].isin(["dentate_candidate", "cerebellar_candidate"])]
        for _, row in structural.iterrows():
            rows.append(
                f"- `{row['source_layer']}` / `{row['broad_class']}`: structural-rank median {row['median_structural_rank']:.4f}."
            )
    rows.extend(
        [
            "",
            "## High-Priority Shared Structural Executors",
            "",
        ]
    )
    for _, row in high_structural.iterrows():
        rows.append(
            f"- `{row['gene']}` ({row['panel']}): dentate detection median {row['dentate_candidate_median_detection']:.3f}; "
            f"{int(row['dentate_sources_detected_at_20pct'])} human bridge sources at >=20% detection."
        )
    rows.extend(["", "## Regional Identity Examples", "", "Dentate identity panel:"])
    for _, row in dentate_identity.iterrows():
        rows.append(
            f"- `{row['gene']}`: tier `{row['support_tier']}`, dentate detection median {row['dentate_candidate_median_detection']:.3f}."
        )
    rows.append("")
    rows.append("Cerebellar identity panel:")
    for _, row in cerebellar_identity.iterrows():
        rows.append(
            f"- `{row['gene']}`: tier `{row['support_tier']}`, human dentate-bridge detection median {row['dentate_candidate_median_detection']:.3f}."
        )
    rows.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The strongest manuscript direction remains identity separation plus structural convergence, not a claim that dentate and cerebellar granule cells are transcriptionally identical.",
            "- High-priority structural/executor candidates are detected across multiple human dentate/hippocampal bridge sources and belong to morphogenesis, cytoskeletal, adhesion, guidance, or synaptic panels.",
            "- Cerebellar identity genes with high human dentate-bridge detection should be treated as specificity warnings, not as evidence of cerebellar identity in DG cells.",
            "",
            "## Outputs",
            "",
            f"- Gene-level human bridge summary: `{OUT_GENE_SUMMARY.relative_to(ROOT)}`",
            f"- Candidate gene packet: `{OUT_PACKET.relative_to(ROOT)}`",
            f"- Structural executor candidate plot: `{OUT_STRUCTURAL_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(rows))


def plot_structural_candidates(packet: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = packet.loc[packet["candidate_role"].eq("shared_structural_executor")].head(20).copy()
    top = top.sort_values("dentate_candidate_median_detection", ascending=True)
    y = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    ax.barh(
        y,
        top["dentate_candidate_median_detection"],
        color="#2a9d8f",
        alpha=0.86,
        label="dentate candidate",
    )
    ax.scatter(
        top["background_or_ambiguous_median_detection"],
        y,
        color="#6b7280",
        s=28,
        label="background/ambiguous median",
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlabel("Median detection fraction across human bridge sources")
    ax.set_title("Shared structural executor candidates")
    ax.set_xlim(0, min(1.0, max(0.1, float(top["dentate_candidate_median_detection"].max()) + 0.08)))
    ax.grid(axis="x", linewidth=0.4, color="#d9d9d9")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT_STRUCTURAL_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))
    panels = load_panels()
    gene_summary = build_gene_summary(panels)
    gene_summary.to_csv(OUT_GENE_SUMMARY, sep="\t", index=False, float_format="%.6g")
    packet = build_packet(gene_summary, panels)
    packet.to_csv(OUT_PACKET, sep="\t", index=False, float_format="%.6g")
    plot_structural_candidates(packet)
    write_md(packet, gene_summary)
    print(f"Wrote {OUT_GENE_SUMMARY}")
    print(f"Wrote {OUT_PACKET}")
    print(f"Wrote {OUT_STRUCTURAL_PLOT}")
    print(f"Wrote {OUT_MD}")
    print(f"genes={packet['gene'].nunique()}; gene_summary_rows={len(gene_summary)}")


if __name__ == "__main__":
    main()
