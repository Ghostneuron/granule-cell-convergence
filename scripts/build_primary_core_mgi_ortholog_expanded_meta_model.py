#!/usr/bin/env python3
"""Dataset-aware meta-model using the expanded MGI ortholog full-matrix screen."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import build_primary_core_mgi_ortholog_meta_model as meta


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

GENOME_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"

OUT_DELTAS = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_unit_deltas.tsv.gz"
OUT_BRANCH = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_branch_summary.tsv"
OUT_GENE = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv"
OUT_HITS = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_shared_hits.tsv"
OUT_MECHANISM_HITS = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv"
OUT_PLOT = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_top_hits.png"
OUT_MD = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model.md"


def rel(path: Path) -> str:
    return meta.rel(path)


def add_ortholog_symbol_class(df: pd.DataFrame, ortho: pd.DataFrame) -> pd.DataFrame:
    keep = ortho[["canonical_gene", "same_canonical_symbol"]].drop_duplicates("canonical_gene").copy()
    keep["ortholog_symbol_class"] = np.where(keep["same_canonical_symbol"], "same_symbol", "nonidentical_symbol")
    return df.drop(columns=["ortholog_symbol_class"], errors="ignore").merge(keep, on="canonical_gene", how="left")


def write_report(
    *,
    ortho: pd.DataFrame,
    mgi_summary: dict[str, int],
    selected_expr: pd.DataFrame,
    genome_expr: pd.DataFrame,
    deltas: pd.DataFrame,
    branch_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    hits: pd.DataFrame,
    mechanism_hits: pd.DataFrame,
) -> None:
    strict_both = hits.loc[hits["ortholog_meta_tier"].eq("strict_shared_both_screens")]
    supported_both = hits.loc[hits["ortholog_meta_tier"].eq("supported_shared_both_screens")]
    full_only = hits.loc[hits["ortholog_meta_tier"].str.contains("full_matrix_only", na=False)]
    selected_only = hits.loc[hits["ortholog_meta_tier"].str.contains("selected_only", na=False)]
    nonidentical_hits = hits.loc[hits["ortholog_symbol_class"].eq("nonidentical_symbol")]
    nonidentical_strict = nonidentical_hits.loc[nonidentical_hits["ortholog_meta_tier"].str.contains("strict", na=False)]
    robust_consensus = gene_summary.loc[meta.bool_col(gene_summary, "consensus_candidate_dataset_robust_all_available")]

    lines = [
        "# Expanded MGI Ortholog Dataset-Aware Meta-Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This model re-runs the dataset-aware candidate-versus-background meta-analysis using the MGI one-to-one full-matrix expression layer. Unlike the previous conservative model, this version includes non-identical human/mouse ortholog symbols resolved through MGI.",
        "",
        "Selected-feature expression remains limited to the selected human-core feature universe. Full-matrix evidence now comes from `primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`.",
        "",
        "## Ortholog Scope",
        "",
        f"- MGI report rows: {mgi_summary['mgi_rows']:,}.",
        f"- MGI human-mouse homology classes: {mgi_summary['mgi_classes_with_human_mouse']:,}.",
        f"- One-to-one human-mouse pairs: {mgi_summary['one_to_one_pairs']:,}.",
        f"- Same-symbol one-to-one pairs: {mgi_summary['strict_same_symbol_one_to_one_pairs']:,}.",
        f"- Non-identical one-to-one pairs: {mgi_summary['one_to_one_pairs'] - mgi_summary['strict_same_symbol_one_to_one_pairs']:,}.",
        f"- One-to-one pairs represented in selected-feature expression rows: {selected_expr['canonical_gene'].nunique():,}.",
        f"- One-to-one pairs represented in MGI full-matrix expression rows: {genome_expr['canonical_gene'].nunique():,}.",
        "",
        "## Meta-Model Summary",
        "",
        f"- Unit delta rows: {len(deltas):,}.",
        f"- Branch summary rows: {len(branch_summary):,}.",
        f"- Gene summary rows: {len(gene_summary):,}.",
        f"- Shared strict both-screen hits: {len(strict_both):,}.",
        f"- Shared supported both-screen hits: {len(supported_both):,}.",
        f"- Shared full-matrix-only hits: {len(full_only):,}.",
        f"- Shared selected-only hits: {len(selected_only):,}.",
        f"- Non-identical-symbol shared hits: {len(nonidentical_hits):,}.",
        f"- Non-identical-symbol strict shared hits: {len(nonidentical_strict):,}.",
        f"- Mechanism-prioritized shared hits: {len(mechanism_hits):,}.",
        "",
        "A branch is supported when at least two datasets contribute, at least 75% of datasets have positive candidate-versus-background deltas, and the median dataset delta is positive. A branch is strict when it also passes the dataset-level sign-test threshold p<=0.25.",
        "",
        "## Dataset-Robust Consensus Genes",
        "",
    ]
    if robust_consensus.empty:
        lines.append("- No dataset-robust consensus genes were retained.")
    else:
        genes = ", ".join(f"`{g}`" for g in robust_consensus["gene"])
        lines.append(f"- Retained in the expanded MGI model: {genes}.")

    lines.extend(["", "## Mechanism-Prioritized Hits", ""])
    for _, row in mechanism_hits.head(30).iterrows():
        mechanism_class = row.get("genome_mechanism_class") or row.get("selected_mechanism_class") or "unclassified"
        lines.append(
            f"- `{row['gene']}` ({row['mechanism_hit_tier']}; {mechanism_class}; "
            f"{row['ortholog_symbol_class']}): {int(row['n_supported_screen_branches'])}/"
            f"{int(row['n_available_screen_branches'])} supported screen/branches."
        )

    lines.extend(["", "## Top Non-Identical Ortholog Hits", ""])
    if nonidentical_hits.empty:
        lines.append("- No non-identical-symbol shared hits were retained.")
    else:
        for _, row in nonidentical_hits.head(30).iterrows():
            lines.append(
                f"- `{row['gene']}` / mouse `{row['mouse_symbol']}` ({row['ortholog_meta_tier']}): "
                f"{int(row['n_supported_screen_branches'])}/{int(row['n_available_screen_branches'])} supported screen/branches."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This is now the best ortholog-aware ranking layer because it includes same-symbol and non-identical one-to-one MGI pairs.",
            "- The top manuscript mechanism claims should still prioritize genes that are both biologically interpretable and robust across datasets, not merely broad neuronal/context hits.",
            "- This expanded MGI layer now feeds the formal rank-meta validation; raw-count/object-level DE remains a later optional strengthening step.",
            "",
            "## Outputs",
            "",
            f"- Unit deltas: `{rel(OUT_DELTAS)}`",
            f"- Branch summary: `{rel(OUT_BRANCH)}`",
            f"- Gene summary: `{rel(OUT_GENE)}`",
            f"- Shared hits: `{rel(OUT_HITS)}`",
            f"- Mechanism-prioritized hits: `{rel(OUT_MECHANISM_HITS)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    ortho, mgi_summary = meta.load_mgi_ortholog_map()
    all_genes = set(ortho["canonical_gene"])

    meta.OUT_PLOT = OUT_PLOT
    selected_expr = meta.load_expression(meta.SELECTED_EXPR, "selected", all_genes)
    genome_expr = meta.load_expression(GENOME_EXPR, "full_matrix", all_genes)
    expr = pd.concat([selected_expr, genome_expr], ignore_index=True, sort=False)

    deltas = meta.build_unit_deltas(expr, ortho)
    branch_summary = meta.summarize_branches(deltas)
    gene_summary = meta.build_gene_summary(branch_summary, ortho)
    gene_summary = add_ortholog_symbol_class(gene_summary, ortho)
    hits = gene_summary.loc[gene_summary["ortholog_meta_tier"].ne("not_shared_or_incomplete")].copy()
    hits = hits.sort_values(
        [
            "shared_meta_strict_both_screens",
            "shared_meta_support_both_screens",
            "n_strict_screen_branches",
            "n_supported_screen_branches",
            "ortholog_meta_priority_score",
        ],
        ascending=[False, False, False, False, False],
    )
    mechanism_hits = meta.mechanism_hits_from_gene_summary(gene_summary)
    mechanism_hits = add_ortholog_symbol_class(mechanism_hits, ortho)

    deltas.to_csv(OUT_DELTAS, sep="\t", index=False, compression="gzip")
    branch_summary.to_csv(OUT_BRANCH, sep="\t", index=False)
    gene_summary.to_csv(OUT_GENE, sep="\t", index=False)
    hits.to_csv(OUT_HITS, sep="\t", index=False)
    mechanism_hits.to_csv(OUT_MECHANISM_HITS, sep="\t", index=False)
    meta.plot_top_hits(gene_summary)
    write_report(
        ortho=ortho,
        mgi_summary=mgi_summary,
        selected_expr=selected_expr,
        genome_expr=genome_expr,
        deltas=deltas,
        branch_summary=branch_summary,
        gene_summary=gene_summary,
        hits=hits,
        mechanism_hits=mechanism_hits,
    )

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Unit deltas: {len(deltas):,}")
    print(f"Shared hits: {len(hits):,}")
    print(f"Mechanism hits: {len(mechanism_hits):,}")


if __name__ == "__main__":
    main()
