#!/usr/bin/env python3
"""Build the compact supplementary-table packet for the DGD analysis."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
OUT = ROOT / "Project/manuscript/Revised_supplementary_tables_20260817"
OUT.mkdir(parents=True, exist_ok=True)

ARCHIVE = Path(
    os.environ.get(
        "GRANULE_ARCHIVE",
        "/Volumes/VV 2021 backup drive 01/Codex_Project_Archive/Hippocanpus&Cerebellum",
    )
)
ORIGINAL_TABLES = (
    ARCHIVE
    / "Project/manuscript/Development_upload_files/direct_tables_tsv_short_names"
)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    tables: list[tuple[str, str, pd.DataFrame, str, str]] = []

    tables.append(
        (
            "Table_S1_primary_core_datasets.tsv",
            "S1_dataset_frame",
            read_tsv(ORIGINAL_TABLES / "Table_S1_primary_core_dataset_frame.tsv"),
            "Primary-core transcriptomic datasets and their analytical roles.",
            "Discovery resources are heterogeneous; within-dataset ranks are used rather than pooled expression.",
        )
    )
    tables.append(
        (
            "Table_S2_candidate_tiers.tsv",
            "S2_candidate_tiers",
            read_tsv(RESULTS / "primary_core_manuscript_candidate_tiers.tsv"),
            "Candidate tier assignments, ortholog status, branch support and rank-priority fields.",
            "Tiers summarize discovery evidence and are not independent validation.",
        )
    )
    tables.append(
        (
            "Table_S3_species_screen_support.tsv",
            "S3_species_screen",
            read_tsv(RESULTS / "dgd_species_stratified_candidates.tsv"),
            "Mouse- and human-bridge candidate support separated by screen and branch.",
            "Human rows are bridge evidence; the main discovery interpretation is mouse-first.",
        )
    )
    tables.append(
        (
            "Table_S4_leave_one_dataset_out.tsv",
            "S4_leave_one_out",
            read_tsv(RESULTS / "dgd_candidate_leave_one_dataset_out.tsv"),
            "Candidate leave-one-dataset-out results with dataset-level medians.",
            "This is an selection-conditioned robustness analysis of a preselected candidate set.",
        )
    )

    independent = read_tsv(RESULTS / "dgd_dataset_level_configuration.tsv")
    tables.append(
        (
            "Table_S5_independent_dataset_scores.tsv",
            "S5_dataset_scores",
            independent,
            "One configuration-delta summary per independent dataset.",
            "Datasets, not nested contrasts, are the inferential units.",
        )
    )
    tables.append(
        (
            "Table_S6_module_inference.tsv",
            "S6_module_tests",
            read_tsv(RESULTS / "dgd_module_level_inference.tsv"),
            "Five-module summaries and the exact downstream-versus-upstream comparison.",
            "The module comparison is directional and descriptive (exact one-sided p=0.10).",
        )
    )

    allen_counts = read_tsv(RESULTS / "dgd_allen_consensus_population_counts.tsv")
    allen_counts = (
        allen_counts[allen_counts["n_cells"].ge(50)]
        .groupby("population", as_index=False)
        .agg(
            source_matrices=("matrix", lambda values: ",".join(sorted(set(values)))),
            n_libraries=("library", "nunique"),
            n_cells=("n_cells", "sum"),
        )
    )
    tables.append(
        (
            "Table_S7_Allen_population_coverage.tsv",
            "S7_Allen_coverage",
            allen_counts,
            "Allen common-matrix population coverage after the minimum-cell filter.",
            "Libraries are independent units; cell counts describe coverage only.",
        )
    )

    allen_gene = read_tsv(RESULTS / "dgd_allen_consensus_local_gene_contrasts.tsv")
    tiers = read_tsv(RESULTS / "primary_core_manuscript_candidate_tiers.tsv")
    tiers = tiers[tiers["manuscript_tier"].isin([
        "Tier 1 core convergent program",
        "Tier 2 high-confidence wiring/synaptic executor",
    ])]
    tier_cols = ["mouse_symbol", "gene", "manuscript_tier", "mechanism_class"]
    allen_gene = allen_gene[allen_gene["gene_symbol"].isin(tiers["mouse_symbol"])]
    allen_gene = tiers[tier_cols].merge(
        allen_gene, left_on="mouse_symbol", right_on="gene_symbol", how="inner"
    )
    tables.append(
        (
            "Table_S8_Allen_candidate_contrasts.tsv",
            "S8_Allen_candidates",
            allen_gene,
            "Allen branch-local candidate contrasts against Purkinje and CA1/CA3 comparators.",
            "Positive values indicate target-minus-local-comparator expression, not strict cell-type specificity.",
        )
    )
    tables.append(
        (
            "Table_S9_Allen_module_contrasts.tsv",
            "S9_Allen_modules",
            read_tsv(RESULTS / "dgd_allen_consensus_local_module_contrasts.tsv"),
            "Allen branch-local contrasts for the five curated modules.",
            "Module effects are heterogeneous and do not establish uniform adult convergence.",
        )
    )
    tables.append(
        (
            "Table_S10_Allen_pair_similarity.tsv",
            "S10_Allen_similarity",
            read_tsv(RESULTS / "dgd_allen_consensus_pair_similarity.tsv"),
            "Direct pairwise module similarities with library-and-gene bootstrap intervals.",
            "The adult cerebellar-versus-dentate target pair does not show downstream-over-upstream similarity.",
        )
    )
    tables.append(
        (
            "Table_S11_Allen_matched_null.tsv",
            "S11_Allen_null",
            read_tsv(RESULTS / "dgd_allen_consensus_candidate_matched_null.tsv"),
            "Expression- and detection-matched null analysis for the preselected candidate sets.",
            "This is external sensitivity evidence for recurrence, not proof of developmental causality.",
        )
    )

    manifest_rows = []
    for filename, sheet, frame, purpose, caveat in tables:
        frame.to_csv(OUT / filename, sep="\t", index=False)
        manifest_rows.append(
            {
                "table": filename.removesuffix(".tsv").replace("_", " "),
                "file": filename,
                "workbook_sheet": sheet,
                "rows": len(frame),
                "columns": len(frame.columns),
                "purpose": purpose,
                "interpretive_caveat": caveat,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_name = "Table_S12_table_manifest.tsv"
    manifest.to_csv(OUT / manifest_name, sep="\t", index=False)

    workbook = OUT / "Revised_Supplementary_Tables_S1-S12.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        for _, sheet, frame, _, _ in tables:
            frame.to_excel(writer, sheet_name=sheet[:31], index=False)
        manifest.to_excel(writer, sheet_name="S12_manifest", index=False)
        for worksheet in writer.book.worksheets:
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions

    print(f"Wrote {len(tables) + 1} supplementary tables and {workbook}")


if __name__ == "__main__":
    main()
