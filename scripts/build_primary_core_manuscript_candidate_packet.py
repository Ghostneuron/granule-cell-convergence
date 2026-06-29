#!/usr/bin/env python3
"""Build a compact manuscript candidate-tier packet from formal rank results."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import build_primary_core_mgi_ortholog_meta_model as meta


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

FORMAL_MECHANISM = RESULTS / "primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv"
FORMAL_SHARED = RESULTS / "primary_core_mgi_ortholog_formal_rank_shared_hits.tsv"

OUT_TSV = RESULTS / "primary_core_manuscript_candidate_tiers.tsv"
OUT_MD = RESULTS / "primary_core_manuscript_candidate_tiers.md"

CORE_SEED = ["GPM6A", "NFIB", "NFIA", "KCNK1", "RFX3", "GABRA2"]
SECOND_TIER = ["PPP3CA", "CACNA2D1", "KCNJ6", "GABRB3", "GRIN2B", "KCNJ3", "KCND2", "STXBP5L", "ROBO2"]

GENE_ROLE = {
    "GPM6A": "membrane/neurite outgrowth structural executor",
    "NFIB": "developmental transcriptional regulator",
    "NFIA": "developmental transcriptional regulator",
    "RFX3": "ciliogenesis/transcriptional regulatory candidate",
    "KCNK1": "ion-channel/excitability tuning candidate",
    "GABRA2": "GABA receptor and synaptic maturation candidate",
    "PPP3CA": "calcineurin/synaptic plasticity signaling candidate",
    "CACNA2D1": "calcium-channel auxiliary subunit and wiring candidate",
    "KCNJ6": "inward-rectifier potassium-channel candidate",
    "GABRB3": "GABA receptor and synaptic maturation candidate",
    "GRIN2B": "glutamatergic synapse maturation candidate",
    "KCNJ3": "inward-rectifier potassium-channel candidate",
    "KCND2": "voltage-gated potassium-channel candidate",
    "STXBP5L": "synaptic vesicle/exocytosis regulatory candidate",
    "ROBO2": "axon guidance and neurite-patterning executor",
}


def rel(path: Path) -> str:
    return meta.rel(path)


def bool_value(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def first_present(row: pd.Series, *cols: str) -> object:
    for col in cols:
        if col in row.index and pd.notna(row[col]) and str(row[col]) != "":
            return row[col]
    return np.nan


def mechanism_class(row: pd.Series) -> str:
    value = first_present(
        row,
        "genome_mechanism_class",
        "selected_mechanism_class",
        "genome_mechanism_class_formal",
        "selected_mechanism_class_formal",
    )
    return "unclassified" if pd.isna(value) else str(value)


def suggested_role(row: pd.Series) -> str:
    gene = str(row["gene"])
    if gene in GENE_ROLE:
        return GENE_ROLE[gene]
    cls = mechanism_class(row)
    if "regulatory" in cls:
        return "regulatory morphogenesis candidate"
    if "synaptic" in cls:
        return "synaptic wiring or excitability candidate"
    if "cytoskeleton" in cls:
        return "cytoskeleton/neurite morphogenesis candidate"
    if "axon_guidance" in cls or "adhesion" in cls:
        return "axon-guidance or adhesion candidate"
    if "structural" in cls:
        return "shared structural-executor candidate"
    return "exploratory shared ortholog candidate"


def branch_summary(row: pd.Series) -> str:
    cols = [
        ("selected_dentate", "branch_nominal_support_selected_dentate"),
        ("selected_cerebellar", "branch_nominal_support_selected_cerebellar"),
        ("full_dentate", "branch_nominal_support_full_matrix_dentate"),
        ("full_cerebellar", "branch_nominal_support_full_matrix_cerebellar"),
    ]
    supported = [label for label, col in cols if bool_value(row.get(col, False))]
    return ",".join(supported)


def assign_mechanism_tier(row: pd.Series) -> tuple[int, str, str]:
    gene = str(row["gene"])
    formal_tier = str(row.get("formal_rank_tier", ""))
    if gene in CORE_SEED:
        return (
            1,
            "Tier 1 core convergent program",
            "Primary manuscript seed gene; use in model schematic and main heatmap.",
        )
    if gene in SECOND_TIER and "both_screens" in formal_tier:
        return (
            2,
            "Tier 2 high-confidence wiring/synaptic executor",
            "Use as pathway-level support around the Tier 1 seed set.",
        )
    if "both_screens" in formal_tier:
        return (
            3,
            "Tier 3 broad both-screen mechanism support",
            "Use as supporting biology; screen for broad neuronal or housekeeping interpretation.",
        )
    return (
        4,
        "Tier 4 screen-specific mechanism support",
        "Use as follow-up or pathway context rather than as a central claim.",
    )


def build_mechanism_packet() -> pd.DataFrame:
    df = pd.read_csv(FORMAL_MECHANISM, sep="\t", low_memory=False)
    records: list[dict[str, object]] = []
    for _, row in df.iterrows():
        tier_rank, tier, manuscript_use = assign_mechanism_tier(row)
        records.append(
            {
                "tier_rank": tier_rank,
                "manuscript_tier": tier,
                "gene": row["gene"],
                "human_symbol": row.get("human_symbol", row["gene"]),
                "mouse_symbol": row.get("mouse_symbol", ""),
                "ortholog_symbol_class": row.get("ortholog_symbol_class", row.get("ortholog_symbol_class_formal", "")),
                "mechanism_class": mechanism_class(row),
                "suggested_role": suggested_role(row),
                "formal_rank_tier": row.get("formal_rank_tier", ""),
                "ortholog_meta_tier": row.get("ortholog_meta_tier", ""),
                "mechanism_hit_tier": row.get("mechanism_hit_tier", ""),
                "formal_n_nominal_branches": row.get("formal_n_nominal_branches", np.nan),
                "formal_n_available_branches": row.get("formal_n_available_branches", np.nan),
                "formal_model_n_nominal_branches": row.get("formal_model_n_nominal_branches", np.nan),
                "formal_rank_priority_score": row.get("formal_rank_priority_score", np.nan),
                "formal_fisher_q_bh_all_available_branches": row.get(
                    "formal_fisher_q_bh_all_available_branches", np.nan
                ),
                "nominal_branch_support": branch_summary(row),
                "median_delta_selected_dentate": row.get("median_dataset_rank_delta_selected_dentate", np.nan),
                "median_delta_selected_cerebellar": row.get("median_dataset_rank_delta_selected_cerebellar", np.nan),
                "median_delta_full_matrix_dentate": row.get("median_dataset_rank_delta_full_matrix_dentate", np.nan),
                "median_delta_full_matrix_cerebellar": row.get(
                    "median_dataset_rank_delta_full_matrix_cerebellar", np.nan
                ),
                "manuscript_use": manuscript_use,
            }
        )
    return pd.DataFrame(records)


def build_nonidentical_packet(existing_genes: set[str], limit: int = 30) -> pd.DataFrame:
    shared = pd.read_csv(FORMAL_SHARED, sep="\t", low_memory=False)
    sub = shared.loc[
        shared["ortholog_symbol_class"].eq("nonidentical_symbol")
        & ~shared["gene"].astype(str).isin(existing_genes)
    ].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = sub.sort_values(
        ["formal_n_nominal_branches", "formal_rank_priority_score"],
        ascending=[False, False],
    ).head(limit)
    records: list[dict[str, object]] = []
    for _, row in sub.iterrows():
        records.append(
            {
                "tier_rank": 5,
                "manuscript_tier": "Tier 5 exploratory non-identical ortholog",
                "gene": row["gene"],
                "human_symbol": row.get("human_symbol", row["gene"]),
                "mouse_symbol": row.get("mouse_symbol", ""),
                "ortholog_symbol_class": row.get("ortholog_symbol_class", ""),
                "mechanism_class": mechanism_class(row),
                "suggested_role": suggested_role(row),
                "formal_rank_tier": row.get("formal_rank_tier", ""),
                "ortholog_meta_tier": row.get("ortholog_meta_tier", ""),
                "mechanism_hit_tier": "nonidentical_ortholog_follow_up",
                "formal_n_nominal_branches": row.get("formal_n_nominal_branches", np.nan),
                "formal_n_available_branches": row.get("formal_n_available_branches", np.nan),
                "formal_model_n_nominal_branches": np.nan,
                "formal_rank_priority_score": row.get("formal_rank_priority_score", np.nan),
                "formal_fisher_q_bh_all_available_branches": row.get(
                    "formal_fisher_q_bh_all_available_branches", np.nan
                ),
                "nominal_branch_support": branch_summary(row),
                "median_delta_selected_dentate": row.get("median_dataset_rank_delta_selected_dentate", np.nan),
                "median_delta_selected_cerebellar": row.get("median_dataset_rank_delta_selected_cerebellar", np.nan),
                "median_delta_full_matrix_dentate": row.get("median_dataset_rank_delta_full_matrix_dentate", np.nan),
                "median_delta_full_matrix_cerebellar": row.get(
                    "median_dataset_rank_delta_full_matrix_cerebellar", np.nan
                ),
                "manuscript_use": "Keep as exploratory ortholog-completeness evidence; validate before central claims.",
            }
        )
    return pd.DataFrame(records)


def write_report(packet: pd.DataFrame) -> None:
    counts = packet["manuscript_tier"].value_counts().sort_index()
    tier1 = packet.loc[packet["tier_rank"].eq(1)]
    tier2 = packet.loc[packet["tier_rank"].eq(2)]
    tier3 = packet.loc[packet["tier_rank"].eq(3)]
    tier4 = packet.loc[packet["tier_rank"].eq(4)]
    tier5 = packet.loc[packet["tier_rank"].eq(5)]

    lines = [
        "# Manuscript Candidate Tiers",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This packet distills the formal MGI ortholog rank-meta validation into manuscript-facing candidate tiers. It does not add a new statistical test; it organizes the formal results into seed genes, pathway-support genes, and exploratory ortholog-completeness candidates.",
        "",
        "## Tier Counts",
        "",
    ]
    for tier, count in counts.items():
        lines.append(f"- {tier}: {int(count)} genes.")

    lines.extend(["", "## Tier 1 Core Convergent Program", ""])
    for _, row in tier1.iterrows():
        lines.append(
            f"- `{row['gene']}`: {row['suggested_role']}; "
            f"{int(row['formal_n_nominal_branches'])}/{int(row['formal_n_available_branches'])} formal branches."
        )

    lines.extend(["", "## Tier 2 High-Confidence Support", ""])
    for _, row in tier2.iterrows():
        lines.append(
            f"- `{row['gene']}`: {row['suggested_role']}; "
            f"{int(row['formal_n_nominal_branches'])}/{int(row['formal_n_available_branches'])} formal branches."
        )

    lines.extend(["", "## Tier 3 Broad Both-Screen Support", ""])
    for _, row in tier3.head(20).iterrows():
        lines.append(f"- `{row['gene']}`: {row['suggested_role']}.")

    lines.extend(["", "## Tier 4 Screen-Specific Support", ""])
    for _, row in tier4.head(20).iterrows():
        lines.append(f"- `{row['gene']}`: {row['suggested_role']}.")

    lines.extend(["", "## Tier 5 Exploratory Ortholog Completeness", ""])
    if tier5.empty:
        lines.append("- No non-identical-symbol exploratory genes were added.")
    else:
        for _, row in tier5.head(20).iterrows():
            lines.append(
                f"- `{row['gene']}` / mouse `{row['mouse_symbol']}`: "
                f"{row['formal_rank_tier']}; {int(row['formal_n_nominal_branches'])}/"
                f"{int(row['formal_n_available_branches'])} formal branches."
            )

    lines.extend(
        [
            "",
            "## Recommended Manuscript Use",
            "",
            "- Build the main model around Tier 1: shared downstream morphology/excitability/regulatory programs, not shared regional identity.",
            "- Use Tier 2 as pathway support for synaptic wiring, calcium signaling, potassium/GABA/glutamate receptor maturation, and axon guidance.",
            "- Use Tier 3 and Tier 4 as supportive context after checking broad neuronal or housekeeping interpretations.",
            "- Keep Tier 5 outside central claims until raw-count/object-level or external validation supports those non-identical ortholog mappings.",
            "",
            "## Outputs",
            "",
            f"- Candidate tier table: `{rel(OUT_TSV)}`",
            f"- Source formal rank model: `{rel(RESULTS / 'primary_core_mgi_ortholog_formal_rank_model.md')}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    mechanism = build_mechanism_packet()
    nonidentical = build_nonidentical_packet(set(mechanism["gene"].astype(str)))
    packet = pd.concat([mechanism, nonidentical], ignore_index=True, sort=False)
    packet = packet.sort_values(
        [
            "tier_rank",
            "formal_n_nominal_branches",
            "formal_model_n_nominal_branches",
            "formal_rank_priority_score",
        ],
        ascending=[True, False, False, False],
    )
    packet.to_csv(OUT_TSV, sep="\t", index=False)
    write_report(packet)
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Candidate tier rows: {len(packet):,}")
    print(packet["manuscript_tier"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
