#!/usr/bin/env python3
"""Infer cerebellar conditioned-medium candidate factors from sequencing data.

The 2005 conditioned-medium experiment implicated peptide-like factors from
postnatal cerebellar cultures, with TGF-beta2 and BDNF accounting for a major
fraction of the anti-proliferative/pro-differentiating activity. This script
asks what additional secreted or ligand-like genes are supported by the current
primary-core cerebellar sequencing data.

This is a nomination screen. mRNA abundance can support source plausibility but
cannot prove secretion, extracellular concentration, protein processing, or
bioactivity.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
FULL_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_primary_core_aim2_niche_pathway_model import PATHWAY_MODULES  # noqa: E402


OUT_UNITS = RESULTS / "cerebellar_conditioned_medium_secretome_units.tsv"
OUT_CONTRASTS = RESULTS / "cerebellar_conditioned_medium_secretome_candidate_contrasts.tsv"
OUT_SUMMARY = RESULTS / "cerebellar_conditioned_medium_secretome_ranked_candidates.tsv"
OUT_STAGE = RESULTS / "cerebellar_conditioned_medium_secretome_gse122357_stage.tsv"
OUT_PLOT = RESULTS / "cerebellar_conditioned_medium_secretome_candidates.png"
OUT_MD = RESULTS / "cerebellar_conditioned_medium_secretome_candidates.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


CEREBELLAR_DATASETS = ["GSE122357", "GSE165657", "GSE312658"]
SAMPLE_STAGE = {
    "GSM3464549_P0": ("P0", 0.0),
    "GSM3464550_P8a": ("P8a", 8.1),
    "GSM3464551_P8b": ("P8b", 8.2),
    "Cerebellum_aggr": ("aggregate", 100.0),
    "Ctrl": ("adult_or_ctrl", 200.0),
    "cKO": ("perturbed_cKO", 201.0),
}


CURATED_SECRETED = {
    # Historical anchors.
    "TGFB2": ("TGF-beta superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "2005 anchor; SMAD-linked stop/maturation factor"),
    "BDNF": ("neurotrophin", "secreted_medium_compatible", "maturation_survival_context_dependent", "2005 anchor; TrkB/MAPK and SMAD-convergent factor"),
    "TGFB1": ("TGF-beta superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "TGF-beta family candidate"),
    "TGFB3": ("TGF-beta superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "TGF-beta family candidate"),
    # BMP/GDF/activin branch.
    "BMP2": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "BMP differentiation/neurogenic brake candidate"),
    "BMP4": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "BMP differentiation/neurogenic brake candidate"),
    "BMP6": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "BMP differentiation/neurogenic brake candidate"),
    "BMP7": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "BMP differentiation/neurogenic brake candidate"),
    "GDF10": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "TGF-beta family differentiation candidate"),
    "GDF11": ("BMP/GDF superfamily", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "TGF-beta family differentiation candidate"),
    "GDF15": ("BMP/GDF superfamily", "secreted_medium_compatible", "stress_maturation_context_dependent", "stress/differentiation-associated TGF-beta family candidate"),
    "INHBA": ("activin/inhibin", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "activin branch candidate"),
    "INHBB": ("activin/inhibin", "secreted_medium_compatible", "anti_proliferative_or_differentiation", "activin branch candidate"),
    "NODAL": ("TGF-beta superfamily", "secreted_medium_compatible", "developmental_context_dependent", "developmental TGF-beta family candidate"),
    # WNT antagonists and neurogenic brake candidates.
    "DKK1": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "WNT-antagonist candidate"),
    "DKK3": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "WNT-antagonist candidate"),
    "SFRP1": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "secreted frizzled-related candidate"),
    "SFRP2": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "secreted frizzled-related candidate"),
    "SFRP4": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "secreted frizzled-related candidate"),
    "WIF1": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "WNT inhibitory factor candidate"),
    "FRZB": ("WNT antagonist", "secreted_medium_compatible", "anti_neurogenic_or_differentiation", "secreted frizzled-related candidate"),
    # Migration-stop/guidance ligands.
    "RELN": ("migration-stop/guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "Reelin migration-stop/maturation candidate"),
    "SEMA3A": ("semaphorin guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted semaphorin guidance candidate"),
    "SEMA3C": ("semaphorin guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted semaphorin guidance candidate"),
    "SEMA3E": ("semaphorin guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted semaphorin guidance candidate"),
    "SLIT1": ("slit/robo guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted guidance candidate"),
    "SLIT2": ("slit/robo guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted guidance candidate"),
    "SLIT3": ("slit/robo guidance", "secreted_medium_compatible", "migration_stop_or_differentiation", "secreted guidance candidate"),
    "NTN1": ("netrin guidance", "secreted_medium_compatible", "migration_or_survival_context_dependent", "secreted guidance candidate"),
    "NTN3": ("netrin guidance", "secreted_medium_compatible", "migration_or_survival_context_dependent", "secreted guidance candidate"),
    # Neurotrophins/cytokines.
    "NGF": ("neurotrophin", "secreted_medium_compatible", "maturation_survival_context_dependent", "tested in 2005 but not main neutralized activity"),
    "NTF3": ("neurotrophin", "secreted_medium_compatible", "maturation_survival_context_dependent", "neurotrophin maturation candidate"),
    "NTF4": ("neurotrophin", "secreted_medium_compatible", "maturation_survival_context_dependent", "neurotrophin maturation candidate"),
    "CNTF": ("cytokine/neurotrophic", "secreted_medium_compatible", "maturation_survival_context_dependent", "neurotrophic cytokine candidate"),
    "LIF": ("cytokine/neurotrophic", "secreted_medium_compatible", "differentiation_context_dependent", "IL6-family differentiation candidate"),
    "IL6": ("cytokine/neuroimmune", "secreted_medium_compatible", "context_dependent", "cytokine candidate"),
    "CXCL12": ("chemokine", "secreted_medium_compatible", "migration_or_progenitor_context_dependent", "chemokine niche candidate"),
    # ECM and matricellular proteins likely to survive conditioned-medium concentration.
    "TGFBI": ("TGF-beta response/ECM", "secreted_medium_compatible", "maturation_ecm_context_dependent", "secreted TGF-beta-induced ECM candidate"),
    "THBS1": ("matricellular/synaptogenic", "secreted_medium_compatible", "maturation_synaptogenic_context_dependent", "secreted matricellular candidate"),
    "THBS2": ("matricellular/synaptogenic", "secreted_medium_compatible", "maturation_synaptogenic_context_dependent", "secreted matricellular candidate"),
    "SPARC": ("matricellular/ECM", "secreted_medium_compatible", "maturation_ecm_context_dependent", "secreted ECM candidate"),
    "SPARCL1": ("matricellular/ECM", "secreted_medium_compatible", "maturation_ecm_context_dependent", "secreted ECM candidate"),
    "CCN1": ("CCN matricellular", "secreted_medium_compatible", "context_dependent", "secreted matricellular candidate"),
    "CCN2": ("CCN matricellular", "secreted_medium_compatible", "context_dependent", "secreted matricellular candidate"),
    "CCN3": ("CCN matricellular", "secreted_medium_compatible", "context_dependent", "secreted matricellular candidate"),
    # Counter-signals and soluble factors to keep visible as alternatives.
    "SHH": ("SHH pathway", "secreted_medium_compatible", "proliferative_counter_signal", "known cerebellar granule precursor expansion cue"),
    "FGF2": ("FGF pathway", "secreted_medium_compatible", "proliferative_or_survival_counter_signal", "growth/permissive counter-signal"),
    "FGF8": ("FGF pathway", "secreted_medium_compatible", "proliferative_or_survival_counter_signal", "growth/permissive counter-signal"),
    "FGF9": ("FGF pathway", "secreted_medium_compatible", "proliferative_or_survival_counter_signal", "growth/permissive counter-signal"),
    "FGF10": ("FGF pathway", "secreted_medium_compatible", "proliferative_or_survival_counter_signal", "growth/permissive counter-signal"),
    "FGF17": ("FGF pathway", "secreted_medium_compatible", "proliferative_or_survival_counter_signal", "growth/permissive counter-signal"),
    "WNT3A": ("WNT pathway", "secreted_medium_compatible", "neurogenic_permissive_counter_signal", "neurogenic/permissive counter-signal"),
    "WNT5A": ("WNT pathway", "secreted_medium_compatible", "context_dependent", "WNT pathway candidate"),
    "WNT7A": ("WNT pathway", "secreted_medium_compatible", "neurogenic_permissive_counter_signal", "neurogenic/permissive counter-signal"),
    "WNT7B": ("WNT pathway", "secreted_medium_compatible", "neurogenic_permissive_counter_signal", "neurogenic/permissive counter-signal"),
    "IGF1": ("IGF pathway", "secreted_medium_compatible", "survival_or_proliferative_counter_signal", "growth/survival counter-signal"),
    "IGF2": ("IGF pathway", "secreted_medium_compatible", "survival_or_proliferative_counter_signal", "growth/survival counter-signal"),
    "MDK": ("growth factor", "secreted_medium_compatible", "survival_or_proliferative_counter_signal", "growth/survival counter-signal"),
    "PTN": ("growth factor", "secreted_medium_compatible", "survival_or_proliferative_counter_signal", "growth/survival counter-signal"),
    # Cell-bound ligands: biologically relevant, weaker conditioned-medium candidates.
    "JAG1": ("Notch ligand", "cell_bound_or_shed", "differentiation_context_dependent", "cell-bound Notch ligand; weaker CM candidate"),
    "JAG2": ("Notch ligand", "cell_bound_or_shed", "differentiation_context_dependent", "cell-bound Notch ligand; weaker CM candidate"),
    "DLL1": ("Notch ligand", "cell_bound_or_shed", "differentiation_context_dependent", "cell-bound Notch ligand; weaker CM candidate"),
    "DLL3": ("Notch ligand", "cell_bound_or_shed", "differentiation_context_dependent", "cell-bound Notch ligand; weaker CM candidate"),
    "DLL4": ("Notch ligand", "cell_bound_or_shed", "differentiation_context_dependent", "cell-bound Notch ligand; weaker CM candidate"),
    "SEMA5A": ("semaphorin guidance", "membrane_or_shed", "migration_stop_or_differentiation", "membrane semaphorin; weaker CM candidate"),
    "SEMA6A": ("semaphorin guidance", "membrane_or_shed", "migration_stop_or_differentiation", "membrane semaphorin; weaker CM candidate"),
    "SEMA6D": ("semaphorin guidance", "membrane_or_shed", "migration_stop_or_differentiation", "membrane semaphorin; weaker CM candidate"),
    "SEMA7A": ("semaphorin guidance", "membrane_or_shed", "migration_stop_or_differentiation", "membrane/GPI semaphorin; weaker CM candidate"),
}


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def finite_median(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def finite_mean(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.mean(arr))


def annotations() -> pd.DataFrame:
    rows = []
    pathway_roles: dict[str, set[str]] = {}
    for module in PATHWAY_MODULES:
        for gene, role in module["genes"].items():
            pathway_roles.setdefault(canon(gene), set()).add(f"{module['pathway_id']}:{role}")
    for gene, (mechanism, secretome_class, expected_effect, rationale) in CURATED_SECRETED.items():
        rows.append(
            {
                "canonical_gene": gene,
                "mechanism_class": mechanism,
                "secretome_class": secretome_class,
                "medium_compatibility": "high" if secretome_class == "secreted_medium_compatible" else "lower",
                "expected_effect": expected_effect,
                "rationale": rationale,
                "pathway_roles": ";".join(sorted(pathway_roles.get(gene, []))),
            }
        )
    return pd.DataFrame(rows)


def load_units(genes: set[str]) -> pd.DataFrame:
    use_cols = [
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "n_cells",
        "gene",
        "canonical_gene",
        "source_gene_symbol",
        "detection_fraction",
        "mean_expression",
        "mean_log1p_expression",
        "eligible_class",
        "mean_log1p_rank_within_sample_gene",
        "source_path",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(FULL_EXPR, sep="\t", usecols=use_cols, chunksize=100_000, low_memory=False):
        chunk["canonical_gene"] = chunk["canonical_gene"].map(canon)
        sub = chunk.loc[
            chunk["dataset"].isin(CEREBELLAR_DATASETS)
            & chunk["canonical_gene"].isin(genes)
            & chunk["eligible_class"].astype(str).str.lower().isin({"true", "1", "yes"})
        ].copy()
        if not sub.empty:
            pieces.append(sub)
    if not pieces:
        raise RuntimeError("No secretome candidate genes found in cerebellar expression table.")
    units = pd.concat(pieces, ignore_index=True)
    units["stage_label"] = units["sample"].map(lambda x: SAMPLE_STAGE.get(str(x), (str(x), np.nan))[0])
    units["stage_order"] = units["sample"].map(lambda x: SAMPLE_STAGE.get(str(x), (str(x), np.nan))[1])
    units = units.sort_values(["canonical_gene", "dataset", "sample", "broad_class"]).reset_index(drop=True)
    return units


def build_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (dataset, sample, gene), sub in units.groupby(["dataset", "sample", "canonical_gene"], sort=False):
        candidate = sub.loc[sub["broad_class"].eq("cerebellar_candidate")]
        background = sub.loc[~sub["broad_class"].eq("cerebellar_candidate")]
        if candidate.empty:
            continue
        cand = candidate.iloc[0]
        bg_rank = finite_median(background["mean_log1p_rank_within_sample_gene"])
        bg_det = finite_median(background["detection_fraction"])
        rows.append(
            {
                "dataset": dataset,
                "sample": sample,
                "stage_label": cand["stage_label"],
                "stage_order": cand["stage_order"],
                "canonical_gene": gene,
                "source_gene_symbol": cand["source_gene_symbol"],
                "candidate_n_cells": int(cand["n_cells"]),
                "candidate_detection_fraction": float(cand["detection_fraction"]),
                "candidate_mean_log1p_expression": float(cand["mean_log1p_expression"]),
                "candidate_rank": float(cand["mean_log1p_rank_within_sample_gene"]),
                "background_median_rank": bg_rank,
                "candidate_minus_background_rank": float(cand["mean_log1p_rank_within_sample_gene"] - bg_rank)
                if np.isfinite(bg_rank)
                else np.nan,
                "background_median_detection": bg_det,
                "candidate_minus_background_detection": float(cand["detection_fraction"] - bg_det)
                if np.isfinite(bg_det)
                else np.nan,
                "background_classes": ",".join(sorted(set(background["broad_class"].astype(str)))),
            }
        )
    return pd.DataFrame(rows).sort_values(["canonical_gene", "dataset", "stage_order"])


def effect_weight(effect: str, medium: str) -> float:
    effect = str(effect)
    medium = str(medium)
    weight = 0.0
    if "anti_proliferative" in effect or "anti_neurogenic" in effect:
        weight += 2.0
    elif "migration_stop" in effect:
        weight += 1.5
    elif "maturation" in effect or "differentiation" in effect:
        weight += 1.0
    elif "counter_signal" in effect or "proliferative" in effect:
        weight -= 1.0
    if medium == "high":
        weight += 0.5
    else:
        weight -= 0.5
    return weight


def build_summary(contrasts: pd.DataFrame, annot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage_rows: list[dict[str, object]] = []
    for gene, sub in contrasts.loc[contrasts["dataset"].eq("GSE122357")].groupby("canonical_gene", sort=False):
        p0 = sub.loc[sub["stage_label"].eq("P0")]
        p8 = sub.loc[sub["stage_label"].isin(["P8a", "P8b"])]
        if p0.empty or p8.empty:
            continue
        stage_rows.append(
            {
                "canonical_gene": gene,
                "p0_rank": finite_median(p0["candidate_rank"]),
                "p8_median_rank": finite_median(p8["candidate_rank"]),
                "p8_minus_p0_rank": finite_median(p8["candidate_rank"]) - finite_median(p0["candidate_rank"]),
                "p0_detection": finite_median(p0["candidate_detection_fraction"]),
                "p8_median_detection": finite_median(p8["candidate_detection_fraction"]),
                "p8_minus_p0_detection": finite_median(p8["candidate_detection_fraction"])
                - finite_median(p0["candidate_detection_fraction"]),
            }
        )
    stage = pd.DataFrame(stage_rows)

    rows: list[dict[str, object]] = []
    for gene, sub in contrasts.groupby("canonical_gene", sort=False):
        cereb = sub.copy()
        positive = pd.to_numeric(cereb["candidate_minus_background_rank"], errors="coerce") > 0
        detected = pd.to_numeric(cereb["candidate_detection_fraction"], errors="coerce") > 0
        row = {
            "canonical_gene": gene,
            "n_candidate_units": int(len(cereb)),
            "n_detected_candidate_units": int(detected.sum()),
            "detected_candidate_fraction": float(detected.mean()) if len(detected) else np.nan,
            "median_candidate_detection": finite_median(cereb["candidate_detection_fraction"]),
            "median_candidate_rank": finite_median(cereb["candidate_rank"]),
            "median_candidate_minus_background_rank": finite_median(cereb["candidate_minus_background_rank"]),
            "positive_candidate_background_units": int(positive.sum()),
            "positive_candidate_background_fraction": float(positive.mean()) if len(positive) else np.nan,
        }
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary = summary.merge(stage, on="canonical_gene", how="left").merge(annot, on="canonical_gene", how="left")
    summary["effect_weight"] = summary.apply(
        lambda row: effect_weight(row["expected_effect"], row["medium_compatibility"]), axis=1
    )
    summary["priority_score"] = (
        summary["effect_weight"]
        + 2.0 * np.maximum(summary["median_candidate_rank"].fillna(0) - 0.5, 0)
        + 2.0 * np.maximum(summary["median_candidate_minus_background_rank"].fillna(0), 0)
        + np.maximum(summary["p8_minus_p0_rank"].fillna(0), 0)
        + summary["median_candidate_detection"].fillna(0)
    )
    high_effect = summary["expected_effect"].astype(str).str.contains(
        "anti_proliferative|anti_neurogenic|migration_stop|differentiation", regex=True
    )
    medium_ok = summary["medium_compatibility"].eq("high")
    expression_ok = (summary["median_candidate_detection"] >= 0.02) | (summary["p8_median_detection"] >= 0.02)
    rank_ok = (summary["median_candidate_rank"] >= 0.5) | (summary["p8_minus_p0_rank"] > 0)
    summary["conditioned_medium_inference_tier"] = np.select(
        [
            medium_ok & high_effect & expression_ok & rank_ok & (summary["priority_score"] >= 3.0),
            medium_ok & high_effect & expression_ok & rank_ok,
            medium_ok & expression_ok,
        ],
        ["high_priority_cm_candidate", "supported_cm_candidate", "detected_secreted_context_candidate"],
        default="low_or_counter_signal",
    )
    summary = summary.sort_values(
        ["conditioned_medium_inference_tier", "priority_score", "median_candidate_detection"],
        ascending=[True, False, False],
    )
    tier_order = {
        "high_priority_cm_candidate": 0,
        "supported_cm_candidate": 1,
        "detected_secreted_context_candidate": 2,
        "low_or_counter_signal": 3,
    }
    summary["tier_order"] = summary["conditioned_medium_inference_tier"].map(tier_order)
    summary = summary.sort_values(["tier_order", "priority_score"], ascending=[True, False]).drop(columns="tier_order")
    return summary, stage


def plot_top(summary: pd.DataFrame) -> None:
    keep = summary.loc[summary["conditioned_medium_inference_tier"].isin(["high_priority_cm_candidate", "supported_cm_candidate"])].copy()
    if keep.empty:
        keep = summary.head(20).copy()
    keep = keep.sort_values("priority_score", ascending=False).head(20)
    keep = keep.sort_values("priority_score")
    colors = keep["conditioned_medium_inference_tier"].map(
        {
            "high_priority_cm_candidate": "#8c2d04",
            "supported_cm_candidate": "#d95f0e",
            "detected_secreted_context_candidate": "#3182bd",
            "low_or_counter_signal": "#737373",
        }
    ).fillna("#737373")

    fig, ax = plt.subplots(figsize=(9, max(5, 0.32 * len(keep))))
    ax.barh(keep["canonical_gene"], keep["priority_score"], color=colors)
    ax.set_xlabel("Conditioned-medium candidate priority score")
    ax.set_ylabel("")
    ax.set_title("Inferred cerebellar conditioned-medium candidate factors")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def fmt(value: object, digits: int = 3) -> str:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(val):
        return "NA"
    return f"{val:.{digits}f}"


def write_markdown(summary: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    high = summary.loc[summary["conditioned_medium_inference_tier"].eq("high_priority_cm_candidate")].copy()
    supported = summary.loc[summary["conditioned_medium_inference_tier"].eq("supported_cm_candidate")].copy()
    detected = summary.loc[summary["conditioned_medium_inference_tier"].eq("detected_secreted_context_candidate")].copy()
    counter = summary.loc[summary["expected_effect"].astype(str).str.contains("counter_signal", na=False)].copy()
    anchors = summary.loc[summary["canonical_gene"].isin(["TGFB2", "BDNF", "NGF"])].copy()

    def bullet_rows(df: pd.DataFrame, n: int = 12) -> list[str]:
        rows: list[str] = []
        for _, row in df.head(n).iterrows():
            rows.append(
                f"- `{row['canonical_gene']}` ({row['mechanism_class']}): "
                f"score {fmt(row['priority_score'])}, candidate detection {fmt(row['median_candidate_detection'])}, "
                f"candidate-background rank delta {fmt(row['median_candidate_minus_background_rank'])}, "
                f"P8-P0 rank delta {fmt(row['p8_minus_p0_rank'])}. {row['rationale']}."
            )
        return rows

    lines = [
        "# Cerebellar Conditioned-Medium Secretome Candidate Inference",
        "",
        "## Question",
        "",
        "Can current sequencing data nominate additional cerebellar conditioned-medium factors that might inhibit proliferation or promote differentiation besides TGF-beta2 and BDNF?",
        "",
        "## Answer",
        "",
        "Yes, partially. The data can nominate secreted or ligand-like genes with cerebellar granule-cell source plausibility and developmental timing, but it cannot prove that the proteins are secreted into medium, properly processed, abundant, or bioactive. This screen should therefore be treated as a prioritized validation list for proteomics, ELISA/Luminex, neutralization, or recombinant-factor rescue experiments.",
        "",
        "## Historical Anchors",
        "",
        *bullet_rows(anchors.sort_values("priority_score", ascending=False), 5),
        "",
        "`TGFB2` is well supported by the present cerebellar transcriptomic data, especially by the P0-to-P8 rise in `GSE122357`. `BDNF` remains experimentally important from the 2005 paper, but its mRNA is sparse in these sequencing tables; that means transcriptomics alone would under-prioritize it relative to antibody/functional evidence.",
        "",
        "## Highest-Priority Candidates",
        "",
        *bullet_rows(high, 15),
        "",
        "## Supported Candidates",
        "",
        *bullet_rows(supported, 15),
        "",
        "## Detected Context or Counter-Signal Candidates",
        "",
        "Some secreted factors are detected but are more likely to be permissive, survival, proliferative, or context-dependent rather than direct stop factors. These matter because conditioned medium can contain mixed activities.",
        "",
        *bullet_rows(pd.concat([detected, counter]).drop_duplicates("canonical_gene").sort_values("priority_score", ascending=False), 12),
        "",
        "## Interpretation",
        "",
        "The most manuscript-useful conclusion is that TGF-beta2 and BDNF remain supported historical anchors, but the modern data make a broader peptide/secreted-factor model testable. The strongest new classes to consider are BMP/GDF/activin-family factors, secreted migration-stop/guidance cues such as Reelin/Semaphorin/Slit, WNT antagonists, and matricellular proteins that could shift proliferating hippocampal progenitors toward differentiation or altered niche adhesion.",
        "",
        "Because the 2005 assay concentrated factors above 6 kDa and found multiple chromatographic activities, a multi-factor medium model is more plausible than a single-factor replacement model. Sequencing should be used to prioritize neutralization/proteomics, not as final proof.",
        "",
        "## Outputs",
        "",
        f"- Unit table: `{OUT_UNITS.relative_to(ROOT)}` ({len(contrasts):,} candidate contrasts after unit extraction).",
        f"- Ranked candidates: `{OUT_SUMMARY.relative_to(ROOT)}` ({len(summary):,} genes).",
        f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    annot = annotations()
    genes = set(annot["canonical_gene"])
    units = load_units(genes).merge(annot, on="canonical_gene", how="left")
    contrasts = build_contrasts(units).merge(annot, on="canonical_gene", how="left")
    summary, stage = build_summary(contrasts, annot)

    units.to_csv(OUT_UNITS, sep="\t", index=False)
    contrasts.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    stage.to_csv(OUT_STAGE, sep="\t", index=False)
    plot_top(summary)
    write_markdown(summary, contrasts)

    print(f"Wrote {OUT_MD}")
    print(f"Secretome units: {len(units):,}")
    print(f"Candidate contrasts: {len(contrasts):,}")
    print(f"Ranked candidates: {len(summary):,}")
    print("Top candidates:")
    cols = [
        "canonical_gene",
        "conditioned_medium_inference_tier",
        "priority_score",
        "mechanism_class",
        "median_candidate_detection",
        "median_candidate_minus_background_rank",
        "p8_minus_p0_rank",
    ]
    print(summary[cols].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
