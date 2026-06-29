#!/usr/bin/env python3
"""Build consensus mechanism candidates across selected-gene and full-matrix screens."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SELECTED_TRIAGE = RESULTS / "primary_core_expanded_gene_mechanism_triage.tsv"
GENOME_TRIAGE = RESULTS / "primary_core_genomewide_symbol_mechanism_triage.tsv"

OUT_TSV = RESULTS / "primary_core_cross_screen_mechanism_consensus.tsv"
OUT_MD = RESULTS / "primary_core_cross_screen_mechanism_consensus.md"


def load_screen(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    keep = [
        "gene",
        "canonical_gene",
        "mechanism_class",
        "manuscript_use",
        "is_original_candidate_gene",
        "dentate_rank_delta_vs_background",
        "cerebellar_rank_delta_vs_background",
        "shared_positive_rank_delta",
        "mechanism_priority_score",
        "minimum_branch_detection",
    ]
    present = [col for col in keep if col in df.columns]
    out = df[present].copy()
    return out.rename(columns={col: f"{prefix}_{col}" for col in present if col not in {"gene", "canonical_gene"}})


def consensus_tier(row: pd.Series) -> str:
    selected_use = row.get("selected_manuscript_use", "")
    genome_use = row.get("genome_manuscript_use", "")
    if selected_use == "mechanism_figure_candidate" and genome_use == "mechanism_figure_candidate":
        return "consensus_figure_candidate"
    if selected_use == "mechanism_figure_candidate" and genome_use == "mechanism_followup_candidate":
        return "selected_figure_genome_followup"
    if selected_use == "mechanism_followup_candidate" and genome_use == "mechanism_figure_candidate":
        return "genome_figure_selected_followup"
    if bool(row.get("selected_shared_positive_rank_delta", False)) and bool(row.get("genome_shared_positive_rank_delta", False)):
        return "shared_positive_both_screens"
    return "screen_specific_or_low_priority"


def main() -> None:
    selected = load_screen(SELECTED_TRIAGE, "selected")
    genome = load_screen(GENOME_TRIAGE, "genome")
    merged = selected.merge(genome, on=["gene", "canonical_gene"], how="outer")
    merged["consensus_tier"] = merged.apply(consensus_tier, axis=1)
    merged["combined_priority_score"] = merged[["selected_mechanism_priority_score", "genome_mechanism_priority_score"]].fillna(0).sum(axis=1)
    merged["selected_in_screen"] = merged["selected_manuscript_use"].notna()
    merged["genome_in_screen"] = merged["genome_manuscript_use"].notna()
    merged = merged.sort_values(["consensus_tier", "combined_priority_score"], ascending=[True, False])
    merged.to_csv(OUT_TSV, sep="\t", index=False)

    counts = merged.groupby("consensus_tier").size().reset_index(name="n_genes")
    top = merged.loc[merged["consensus_tier"].eq("consensus_figure_candidate")].head(40)
    lines = [
        "# Cross-Screen Mechanism Consensus",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This table compares the 2,169-gene selected-feature mechanism triage with the full-matrix same-symbol mechanism triage. Genes that survive both screens are the best near-term mechanism candidates.",
        "",
        "## Consensus Counts",
        "",
    ]
    for _, row in counts.iterrows():
        lines.append(f"- `{row['consensus_tier']}`: {int(row['n_genes'])} genes.")
    lines.extend(["", "## Consensus Figure Candidates", ""])
    for _, row in top.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['selected_mechanism_class']} / {row['genome_mechanism_class']}): "
            f"selected delta {row['selected_dentate_rank_delta_vs_background']:.3f}/{row['selected_cerebellar_rank_delta_vs_background']:.3f}, "
            f"full-matrix delta {row['genome_dentate_rank_delta_vs_background']:.3f}/{row['genome_cerebellar_rank_delta_vs_background']:.3f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Consensus figure candidates are stronger than either screen alone because they survive both the selected-feature bridge and the independent full-matrix same-symbol pass.",
            "- This table should drive the next mechanism shortlist, while curated ortholog mapping and mixed-effect DE remain the final validation step.",
            "",
            "## Output",
            "",
            f"- Consensus table: `{OUT_TSV.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote {len(merged):,} consensus rows")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
