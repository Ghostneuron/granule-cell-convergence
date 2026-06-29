#!/usr/bin/env python3
"""Triage expanded pseudobulk hits into manuscript-useful mechanism classes."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

STATS = RESULTS / "primary_core_expanded_gene_pseudobulk_statistics.tsv"
SHARED_HITS = RESULTS / "primary_core_expanded_gene_pseudobulk_shared_hits.tsv"

OUT_TRIAGE = RESULTS / "primary_core_expanded_gene_mechanism_triage.tsv"
OUT_MD = RESULTS / "primary_core_expanded_gene_mechanism_triage.md"


AXON_GUIDANCE_PREFIXES = (
    "ROBO",
    "SLIT",
    "SEMA",
    "PLXN",
    "EPH",
    "NRP",
    "CNTN",
    "NCAM",
    "L1CAM",
    "NRXN",
    "NLGN",
    "CADM",
    "PCDH",
    "DSCAM",
    "NTRK",
    "DCC",
)

STRUCTURAL_PREFIXES = (
    "MAP",
    "TUB",
    "ACT",
    "STMN",
    "DPYSL",
    "CFL",
    "CAPZ",
    "CORO",
    "ADD",
    "ABLIM",
    "KIF",
    "DYN",
    "DNM",
    "RTN",
    "BASP",
    "MARCK",
    "SPT",
    "MYO",
    "NEFL",
    "NEFM",
    "INA",
    "GAP43",
    "GPM6",
)

SYNAPTIC_PREFIXES = (
    "SYN",
    "SYT",
    "SNAP",
    "STX",
    "VAMP",
    "SLC17A",
    "GRIN",
    "GABR",
    "CAMK",
    "CACN",
    "KCN",
    "SCN",
    "DLG",
    "SHANK",
    "HOMER",
    "PPP3C",
    "CALM",
)

REGULATORY_PREFIXES = (
    "TCF",
    "NFIA",
    "NFIB",
    "NEUROD",
    "PAX",
    "ZIC",
    "EOMES",
    "SOX",
    "MEF",
    "CREB",
    "CAMTA",
    "BCL",
    "FOX",
    "LHX",
    "EMX",
    "DLX",
    "EGR",
    "RFX",
)

METABOLIC_PREFIXES = (
    "NDUF",
    "COX",
    "ATP5",
    "RPS",
    "RPL",
    "EIF",
    "EEF",
    "HSP",
    "PRDX",
    "SOD",
    "UBC",
)

RNA_PROCESSING_PREFIXES = (
    "DDX",
    "HNRNP",
    "SRSF",
    "CELF",
    "ELAVL",
    "RBFOX",
    "LUC7L",
    "NOVA",
    "KHDRBS",
)

BROAD_NEURONAL_GENES = {
    "MAP2",
    "SNAP25",
    "SYT1",
    "SYN1",
    "TUBB3",
    "SLC17A7",
    "NEUROD1",
    "NEUROD2",
    "DCX",
    "TBR1",
    "ST18",
    "EOMES",
    "PAX6",
}

SPECIFICITY_WARNING_ROLES = {
    "regional_cerebellar_identity",
    "regional_dentate_identity",
}


def startswith_any(gene: str, prefixes: tuple[str, ...]) -> bool:
    return any(gene.startswith(prefix) for prefix in prefixes)


def mechanism_class(row: pd.Series) -> str:
    gene = str(row["canonical_gene"]).upper()
    role = str(row.get("candidate_role", ""))
    if role in SPECIFICITY_WARNING_ROLES:
        return "regional_identity_or_specificity_warning"
    if role == "shared_structural_executor":
        return "curated_shared_structural_executor"
    if gene in BROAD_NEURONAL_GENES or role == "shared_granule_neuronal_state":
        return "broad_neuronal_state"
    if startswith_any(gene, AXON_GUIDANCE_PREFIXES):
        return "axon_guidance_adhesion"
    if startswith_any(gene, STRUCTURAL_PREFIXES):
        return "cytoskeleton_morphogenesis"
    if startswith_any(gene, SYNAPTIC_PREFIXES):
        return "synaptic_wiring"
    if startswith_any(gene, REGULATORY_PREFIXES):
        return "regulatory_morphogenesis_candidate"
    if startswith_any(gene, METABOLIC_PREFIXES) or role == "supporting_metabolic_validation":
        return "metabolic_or_housekeeping_support"
    if startswith_any(gene, RNA_PROCESSING_PREFIXES):
        return "rna_processing_neuronal_state"
    return "unclassified_shared_neuronal_or_context"


def manuscript_use(row: pd.Series) -> str:
    cls = row["mechanism_class"]
    if not bool(row["shared_positive_rank_delta"]):
        if row["dentate_rank_delta_vs_background"] > 0 and row["cerebellar_rank_delta_vs_background"] <= 0:
            return "branch_specific_dentate_context"
        if row["cerebellar_rank_delta_vs_background"] > 0 and row["dentate_rank_delta_vs_background"] <= 0:
            return "branch_specific_cerebellar_context"
        return "not_prioritized"
    if cls in {
        "curated_shared_structural_executor",
        "axon_guidance_adhesion",
        "cytoskeleton_morphogenesis",
        "synaptic_wiring",
        "regulatory_morphogenesis_candidate",
    }:
        if bool(row["shared_strict_bh_0_10"]):
            return "mechanism_figure_candidate"
        return "mechanism_followup_candidate"
    if cls == "regional_identity_or_specificity_warning":
        return "specificity_warning_or_branch_context"
    if cls in {"broad_neuronal_state", "metabolic_or_housekeeping_support", "rna_processing_neuronal_state"}:
        return "supporting_context_not_core_executor"
    return "exploratory_context"


def priority_score(row: pd.Series) -> float:
    score = float(row["combined_rank_delta"])
    score += 0.25 * float(row["minimum_branch_detection"] if pd.notna(row["minimum_branch_detection"]) else 0)
    if bool(row["shared_strict_bh_0_10"]):
        score += 1.0
    elif bool(row["shared_strict_bh_0_20"]):
        score += 0.5
    if row["mechanism_class"] in {
        "curated_shared_structural_executor",
        "axon_guidance_adhesion",
        "cytoskeleton_morphogenesis",
        "synaptic_wiring",
        "regulatory_morphogenesis_candidate",
    }:
        score += 0.5
    if row["mechanism_class"] in {"metabolic_or_housekeeping_support", "regional_identity_or_specificity_warning"}:
        score -= 0.5
    if bool(row["is_original_candidate_gene"]):
        score += 0.25
    return score


def build_triage() -> pd.DataFrame:
    stats = pd.read_csv(STATS, sep="\t")
    stats["mechanism_class"] = stats.apply(mechanism_class, axis=1)
    stats["manuscript_use"] = stats.apply(manuscript_use, axis=1)
    stats["mechanism_priority_score"] = stats.apply(priority_score, axis=1)
    stats = stats.sort_values(
        ["manuscript_use", "mechanism_priority_score", "combined_rank_delta", "minimum_branch_detection"],
        ascending=[True, False, False, False],
    )
    return stats


def write_report(triage: pd.DataFrame) -> None:
    shared = triage.loc[triage["shared_positive_rank_delta"]].copy()
    mechanism = triage.loc[triage["manuscript_use"].eq("mechanism_figure_candidate")].copy()
    followup = triage.loc[triage["manuscript_use"].eq("mechanism_followup_candidate")].copy()
    context = triage.loc[triage["manuscript_use"].eq("supporting_context_not_core_executor")].copy()
    warnings = triage.loc[triage["manuscript_use"].eq("specificity_warning_or_branch_context")].copy()
    class_counts = (
        shared.groupby(["mechanism_class", "manuscript_use"], dropna=False)
        .size()
        .reset_index(name="n_genes")
        .sort_values(["manuscript_use", "n_genes"], ascending=[True, False])
    )

    lines = [
        "# Expanded Gene Mechanism Triage",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "The expanded 2,169-gene screen is intentionally broad. This triage separates candidate morphology/wiring mechanisms from broad neuronal state, metabolic/supporting, RNA-processing, and regional-identity warning genes.",
        "",
        "## Summary",
        "",
        f"- Shared-positive genes in expanded screen: {len(shared):,}.",
        f"- Mechanism figure candidates: {len(mechanism):,}.",
        f"- Mechanism follow-up candidates: {len(followup):,}.",
        f"- Supporting-context genes: {len(context):,}.",
        f"- Specificity-warning or branch-context genes: {len(warnings):,}.",
        "",
        "## Shared-Positive Class Counts",
        "",
    ]
    for _, row in class_counts.iterrows():
        lines.append(f"- `{row['mechanism_class']}` / `{row['manuscript_use']}`: {int(row['n_genes'])} genes.")

    lines.extend(["", "## Top Mechanism Figure Candidates", ""])
    top = mechanism.head(40)
    if top.empty:
        lines.append("- No mechanism figure candidates passed the triage criteria.")
    else:
        for _, row in top.iterrows():
            lines.append(
                f"- `{row['gene']}` ({row['mechanism_class']}): dentate delta {row['dentate_rank_delta_vs_background']:.3f}, "
                f"cerebellar delta {row['cerebellar_rank_delta_vs_background']:.3f}, "
                f"original_packet={bool(row['is_original_candidate_gene'])}."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The best mechanism candidates are not simply the highest shared-positive genes; they are shared-positive genes that plausibly affect neurite growth, cytoskeletal remodeling, axon guidance, adhesion, synaptic wiring, or upstream morphogenesis regulation.",
            "- Broad neuronal and metabolic genes are still useful as controls/context, but they should not be the central mechanism claim.",
            "- Regional identity genes that appear shared-positive are warning markers: they help refine specificity rather than prove shared morphology.",
            "",
            "## Output",
            "",
            f"- Triage table: `{OUT_TRIAGE.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    triage = build_triage()
    triage.to_csv(OUT_TRIAGE, sep="\t", index=False)
    write_report(triage)
    print(f"Wrote {len(triage):,} triaged genes")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
