#!/usr/bin/env python3
"""Triage genome-wide same-symbol pseudobulk hits into mechanism classes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import triage_expanded_pseudobulk_hits as triage


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

STATS = RESULTS / "primary_core_genomewide_symbol_pseudobulk_statistics.tsv"
CANDIDATE_PACKET = RESULTS / "human_bridge_candidate_gene_packet.tsv"

OUT_TRIAGE = RESULTS / "primary_core_genomewide_symbol_mechanism_triage.tsv"
OUT_MD = RESULTS / "primary_core_genomewide_symbol_mechanism_triage.md"


def canon_gene(gene: object) -> str:
    if pd.isna(gene):
        return ""
    return str(gene).strip().upper()


def build_triage() -> pd.DataFrame:
    stats = pd.read_csv(STATS, sep="\t")
    packet = pd.read_csv(CANDIDATE_PACKET, sep="\t")
    packet["canonical_gene"] = packet["gene"].map(canon_gene)
    packet_meta = packet.set_index("canonical_gene")[["panel", "candidate_role", "support_tier"]]
    stats = stats.merge(packet_meta, on="canonical_gene", how="left", suffixes=("", "_packet"))
    stats["is_original_candidate_gene"] = stats["candidate_role_packet"].notna()
    stats["candidate_role"] = stats["candidate_role_packet"].fillna(stats["candidate_role"])
    stats["panel"] = stats["panel"].fillna("genomewide_same_symbol")
    stats["support_tier"] = stats["support_tier"].fillna("full_matrix_symbol_universe")
    stats["shared_strict_bh_0_20"] = (
        stats["shared_positive_rank_delta"]
        & stats["dentate_rank_p_adj_bh"].lt(0.20)
        & stats["cerebellar_rank_p_adj_bh"].lt(0.20)
    )
    stats["mechanism_class"] = stats.apply(triage.mechanism_class, axis=1)
    stats["manuscript_use"] = stats.apply(triage.manuscript_use, axis=1)
    stats["mechanism_priority_score"] = stats.apply(triage.priority_score, axis=1)
    stats = stats.sort_values(
        ["manuscript_use", "mechanism_priority_score", "combined_rank_delta", "minimum_branch_detection"],
        ascending=[True, False, False, False],
    )
    return stats


def write_report(df: pd.DataFrame) -> None:
    shared = df.loc[df["shared_positive_rank_delta"]].copy()
    mechanism = df.loc[df["manuscript_use"].eq("mechanism_figure_candidate")].copy()
    followup = df.loc[df["manuscript_use"].eq("mechanism_followup_candidate")].copy()
    original = shared.loc[shared["is_original_candidate_gene"]].copy()
    class_counts = (
        shared.groupby(["mechanism_class", "manuscript_use"], dropna=False)
        .size()
        .reset_index(name="n_genes")
        .sort_values(["manuscript_use", "n_genes"], ascending=[True, False])
    )

    lines = [
        "# Genome-Wide Same-Symbol Mechanism Triage",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This triage classifies the full-matrix same-symbol pseudobulk hits into likely morphology/wiring mechanisms, broad neuronal state, metabolic/supporting genes, RNA-processing state, and regional-identity warnings.",
        "",
        "## Summary",
        "",
        f"- Shared-positive same-symbol hits: {len(shared):,}.",
        f"- Mechanism figure candidates: {len(mechanism):,}.",
        f"- Mechanism follow-up candidates: {len(followup):,}.",
        f"- Original 67-gene packet genes recovered among shared-positive hits: {len(original):,}.",
        "",
        "## Shared-Positive Class Counts",
        "",
    ]
    for _, row in class_counts.iterrows():
        lines.append(f"- `{row['mechanism_class']}` / `{row['manuscript_use']}`: {int(row['n_genes'])} genes.")

    lines.extend(["", "## Top Mechanism Figure Candidates", ""])
    for _, row in mechanism.head(50).iterrows():
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
            "- The full-matrix screen recovers broad neuronal and morphogenesis-related programs across a much larger symbol universe.",
            "- The mechanism figure list should be treated as a prioritization table. It still needs curated ortholog mapping and a model with dataset/sample effects before final manuscript claims.",
            "- Genes that also appeared in the selected-gene triage are the best near-term candidates because they survive two independent feature-universe definitions.",
            "",
            "## Output",
            "",
            f"- Triage table: `{OUT_TRIAGE.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    df = build_triage()
    df.to_csv(OUT_TRIAGE, sep="\t", index=False)
    write_report(df)
    print(f"Wrote {len(df):,} genome-wide triaged genes")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
