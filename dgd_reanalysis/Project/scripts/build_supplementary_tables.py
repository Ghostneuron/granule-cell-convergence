#!/usr/bin/env python3
"""Build machine-readable TSV tables for the compact DGD supplement."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
OUT = ROOT / "Project/manuscript/source_tables"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    tables: list[tuple[str, str, pd.DataFrame, str, str]] = []

    tables.append(
        (
            "Table_S1_primary_core_datasets.tsv",
            "S1_dataset_frame",
            read_tsv(OUT / "Table_S1_primary_core_datasets.tsv"),
            "Primary-core transcriptomic datasets and their analytical roles.",
            "Discovery resources are heterogeneous; within-dataset ranks are used rather than pooled expression.",
        )
    )
    candidate_tiers = read_tsv(RESULTS / "primary_core_manuscript_candidate_tiers.tsv").rename(
        columns={
            "manuscript_tier": "candidate_tier",
            "manuscript_use": "interpretive_role",
        }
    )
    candidate_tiers["interpretive_role"] = (
        candidate_tiers["interpretive_role"]
        .str.replace("Primary manuscript seed gene; use in model schematic and main heatmap.",
                     "Primary seed gene for the model schematic and candidate heatmap.", regex=False)
        .str.replace("Use in extended candidate heatmap and supporting interpretation.",
                     "Extended candidate for the heatmap and supporting interpretation.", regex=False)
        .str.replace("Keep as context-specific or secondary candidate.",
                     "Context-dependent secondary candidate.", regex=False)
    )
    tables.append(
        (
            "Table_S2_candidate_tiers.tsv",
            "S2_candidate_tiers",
            candidate_tiers,
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
            "This sensitivity analysis uses a preselected candidate set and is not independent validation.",
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
    specificity = read_tsv(
        RESULTS / "dgd_allen_consensus_candidate_specificity.tsv"
    )
    profile_cols = [
        column
        for column in specificity.columns
        if column.startswith("mean_log2_")
        or column
        in {
            "gene_symbol",
            "minimum_target_mean_log2",
            "maximum_comparator_mean_log2",
            "target_min_minus_comparator_max",
            "strict_target_pair_specific",
        }
    ]
    allen_gene = allen_gene.merge(
        specificity[profile_cols], on="gene_symbol", how="left"
    )
    tables.append(
        (
            "Table_S8_Allen_candidate_contrasts.tsv",
            "S8_Allen_candidates",
            allen_gene,
            "Allen candidate contrasts against Purkinje, CA1 and CA3, with expression across all nine populations.",
            "Positive values indicate target-minus-comparator expression. CA1-only and CA3-only fields test reference-weighting sensitivity; the population profile quantifies non-exclusivity.",
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
            "This tests generalization of a direction selected in discovery, not blind validation or developmental causality.",
        )
    )

    transfer = read_tsv(RESULTS / "dgd_allen_cross_region_specificity_summary.tsv")
    matched_panels = read_tsv(
        RESULTS / "dgd_allen_cross_region_matched_gene_null.tsv.gz"
    )
    matched_iqr = (
        matched_panels.groupby("feature_set")["minimum_bidirectional_auc"]
        .agg(
            matched_null_q25_minimum_auc=lambda values: values.quantile(0.25),
            matched_null_q75_minimum_auc=lambda values: values.quantile(0.75),
        )
        .reset_index()
    )
    transfer = transfer.merge(matched_iqr, on="feature_set", how="left")
    concordance = read_tsv(
        RESULTS / "dgd_allen_cross_region_contrast_concordance.tsv"
    )
    stage = read_tsv(RESULTS / "dgd_allen_cross_region_stage_sensitivity.tsv")
    transfer_cols = [
        "feature_set",
        "feature_set_label",
        "n_features",
        "dentate_to_cerebellum_auc",
        "cerebellum_to_dentate_auc",
        "minimum_bidirectional_auc",
        "intersection_union_label_permutation_p",
        "intersection_union_label_permutation_q_bh",
        "matched_null_median_minimum_auc",
        "matched_null_q25_minimum_auc",
        "matched_null_q75_minimum_auc",
        "matched_null_95ci_low_minimum_auc",
        "matched_null_95ci_high_minimum_auc",
        "matched_gene_panel_empirical_p",
        "matched_gene_panel_q_bh",
        "transfer_interpretation",
    ]
    concordance_cols = [
        "feature_set",
        "cosine_concordance",
        "spearman_contrast_correlation",
        "same_sign_fraction",
        "both_positive_fraction",
        "matched_null_p_greater_cosine_concordance",
        "matched_null_p_greater_spearman_contrast_correlation",
        "matched_null_p_greater_same_sign_fraction",
        "matched_null_p_greater_both_positive_fraction",
    ]
    stage_cols = [
        "feature_set",
        "n_shared_libraries",
        "median_mature_margin",
        "median_immature_margin",
        "median_immature_minus_mature_margin",
        "paired_wilcoxon_p_two_sided",
        "paired_wilcoxon_q_bh",
        "scope",
    ]
    transfer = transfer[transfer_cols].merge(
        concordance[concordance_cols], on="feature_set", how="left"
    ).merge(stage[stage_cols], on="feature_set", how="left")
    tables.append(
        (
            "Table_S12_Allen_cross_region.tsv",
            "S12_cross_region",
            transfer,
            "Allen cross-region transfer, contrast concordance and adult-state sensitivity.",
            "Positive ranking did not exceed both configuration nulls; the immature group is an adult cell state.",
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
    manifest_name = "Table_S13_table_manifest.tsv"
    manifest.to_csv(OUT / manifest_name, sep="\t", index=False)

    print(f"Wrote {len(tables) + 1} supplementary TSV tables to {OUT}")


if __name__ == "__main__":
    main()
