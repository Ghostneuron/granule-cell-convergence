#!/usr/bin/env python3
"""Build cluster-supported GSE322785 marker groups and epigenomic scores."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SELECTED_SUMMARY = RESULTS / "gse322785_human_h5_selected_matrix_summary.tsv"
CLUSTER_ASSIGNMENTS = RESULTS / "gse322785_human_h5_cluster_validation_barcode_assignments.tsv.gz"
CLUSTER_ENRICHMENT = RESULTS / "gse322785_human_h5_cluster_validation_marker_call_enrichment.tsv"

OUT_RULES = RESULTS / "gse322785_human_h5_cluster_supported_marker_rules.tsv"
OUT_BARCODES = RESULTS / "gse322785_human_h5_cluster_supported_marker_barcodes.tsv.gz"
OUT_FEATURE = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_feature_scores.tsv.gz"
OUT_MODULE = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_module_scores.tsv"
OUT_CONTRASTS = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_contrasts.tsv"
OUT_MD = RESULTS / "gse322785_human_h5_cluster_supported_epigenomic_scoring.md"

MIN_CLUSTER_MARKER_BARCODES = 20
MIN_FRACTION_OF_CLUSTER = 0.02
MIN_FOLD_ENRICHMENT = 2.0

EXCLUDED_CALLS = {
    "ambiguous_neuronal",
    "ambiguous_non_neuronal_or_niche",
    "low_information_or_low_qc",
}

COMPARATORS = [
    "purkinje_candidate",
    "astrocyte_bergmann_candidate",
    "oligodendrocyte_candidate",
    "opc_candidate",
    "microglia_candidate",
    "vascular_candidate",
    "inhibitory_interneuron_candidate",
]


def split_semicolon(value: object) -> list[str]:
    text = str(value)
    if not text or text == "nan":
        return []
    return [part for part in text.split(";") if part]


def build_rules(assignments: pd.DataFrame, enrichment: pd.DataFrame) -> pd.DataFrame:
    donor_totals = assignments.groupby("donor_id").size().rename("n_donor_basic_qc_barcodes")
    call_totals = (
        assignments.groupby(["donor_id", "marker_call"])
        .size()
        .rename("n_marker_call_barcodes_in_donor")
        .reset_index()
    )
    rules = enrichment.merge(call_totals, on=["donor_id", "marker_call"], how="left")
    rules = rules.merge(donor_totals.reset_index(), on="donor_id", how="left")
    rules["donor_marker_call_fraction"] = (
        rules["n_marker_call_barcodes_in_donor"] / rules["n_donor_basic_qc_barcodes"]
    )
    rules["fold_enrichment_over_donor"] = rules["fraction_of_cluster"] / rules[
        "donor_marker_call_fraction"
    ].replace(0, np.nan)
    rules["cluster_support_tier"] = "not_supported"
    eligible = (
        ~rules["marker_call"].isin(EXCLUDED_CALLS)
        & (rules["n_marker_call_barcodes_in_cluster"] >= MIN_CLUSTER_MARKER_BARCODES)
        & (rules["fraction_of_cluster"] >= MIN_FRACTION_OF_CLUSTER)
        & (rules["fold_enrichment_over_donor"] >= MIN_FOLD_ENRICHMENT)
    )
    high = eligible & (
        (rules["fold_enrichment_over_donor"] >= 3.0)
        | (rules["fraction_of_cluster"] >= 0.40)
    )
    rules.loc[eligible, "cluster_support_tier"] = "cluster_enriched_moderate"
    rules.loc[high, "cluster_support_tier"] = "cluster_enriched_high"
    rules["cluster_support_rule"] = (
        f"n>={MIN_CLUSTER_MARKER_BARCODES}; cluster_fraction>={MIN_FRACTION_OF_CLUSTER}; "
        f"fold_enrichment>={MIN_FOLD_ENRICHMENT}; ambiguous/low-info calls excluded"
    )
    return rules.sort_values(
        ["donor_id", "cluster_id", "cluster_support_tier", "marker_call"],
        ascending=[True, True, True, True],
    )


def build_supported_barcodes(assignments: pd.DataFrame, rules: pd.DataFrame) -> pd.DataFrame:
    supported_rules = rules.loc[~rules["cluster_support_tier"].eq("not_supported")].copy()
    key_cols = ["sample_accession", "donor_id", "cluster_id", "marker_call"]
    keep_rules = supported_rules[
        key_cols
        + [
            "cluster_support_tier",
            "n_marker_call_barcodes_in_cluster",
            "n_cluster_barcodes",
            "fraction_of_cluster",
            "n_marker_call_barcodes_in_donor",
            "donor_marker_call_fraction",
            "fold_enrichment_over_donor",
        ]
    ]
    supported = assignments.merge(keep_rules, on=key_cols, how="inner")
    supported = supported.loc[supported["marker_confidence"].isin(["high", "medium"])].copy()
    supported["cluster_supported_group"] = supported["marker_call"]
    return supported


def feature_scores_for_group(
    sample_accession: str,
    donor_id: str,
    marker_call: str,
    support_scope: str,
    X_group: sparse.csr_matrix,
    var: pd.DataFrame,
) -> pd.DataFrame:
    n = X_group.shape[0]
    sums = np.asarray(X_group.sum(axis=0)).ravel()
    detected = X_group.getnnz(axis=0)
    mean_count = sums / max(n, 1)
    detection_fraction = detected / max(n, 1)
    score = np.log1p(mean_count) + detection_fraction
    out = var[
        [
            "selected_feature_index",
            "source_feature_index_0based",
            "feature_id",
            "feature_name",
            "feature_type",
            "linked_genes",
            "best_priorities",
            "model_terms_supported",
            "target_sets",
            "peak_categories",
            "selection_reasons",
            "selected_for_epigenomic_target",
            "selected_for_marker_panel",
        ]
    ].copy()
    out.insert(0, "sample_accession", sample_accession)
    out.insert(1, "donor_id", donor_id)
    out.insert(2, "marker_call", marker_call)
    out.insert(3, "cluster_support_scope", support_scope)
    out.insert(4, "n_barcodes", n)
    out["mean_count"] = mean_count
    out["detected_barcodes"] = detected
    out["detection_fraction"] = detection_fraction
    out["feature_score_log1p_mean_plus_detection"] = score
    return out


def explode_module_scores(feature_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rec in feature_scores.to_dict("records"):
        model_terms = split_semicolon(rec["model_terms_supported"]) or ["marker_or_unassigned"]
        target_sets = split_semicolon(rec["target_sets"]) or ["marker_or_unassigned"]
        peak_categories = split_semicolon(rec["peak_categories"]) or ["gene_expression_or_marker"]
        linked_genes = split_semicolon(rec["linked_genes"]) or [str(rec["feature_name"])]
        for term in model_terms:
            for target_set in target_sets:
                for peak_category in peak_categories:
                    rows.append(
                        {
                            "sample_accession": rec["sample_accession"],
                            "donor_id": rec["donor_id"],
                            "marker_call": rec["marker_call"],
                            "cluster_support_scope": rec["cluster_support_scope"],
                            "feature_type": rec["feature_type"],
                            "model_term_supported": term,
                            "target_set": target_set,
                            "peak_category": peak_category,
                            "n_barcodes": rec["n_barcodes"],
                            "feature_name": rec["feature_name"],
                            "linked_genes": ";".join(linked_genes),
                            "mean_count": rec["mean_count"],
                            "detection_fraction": rec["detection_fraction"],
                            "feature_score": rec["feature_score_log1p_mean_plus_detection"],
                        }
                    )
    exploded = pd.DataFrame(rows)
    if exploded.empty:
        return exploded
    return (
        exploded.groupby(
            [
                "sample_accession",
                "donor_id",
                "marker_call",
                "cluster_support_scope",
                "feature_type",
                "model_term_supported",
                "target_set",
                "peak_category",
            ],
            dropna=False,
        )
        .agg(
            n_barcodes=("n_barcodes", "first"),
            n_features=("feature_name", "nunique"),
            n_linked_gene_strings=("linked_genes", "nunique"),
            mean_count=("mean_count", "mean"),
            mean_detection_fraction=("detection_fraction", "mean"),
            mean_feature_score=("feature_score", "mean"),
            median_feature_score=("feature_score", "median"),
        )
        .reset_index()
    )


def score_supported_groups(supported: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(SELECTED_SUMMARY, sep="\t")
    feature_tables = []
    for row in summary.to_dict("records"):
        donor_id = str(row["donor_id"])
        donor_supported = supported.loc[supported["donor_id"].eq(donor_id)].copy()
        if donor_supported.empty:
            continue
        X = sparse.load_npz(ROOT / row["matrix_path"]).tocsr()
        var = pd.read_csv(ROOT / row["var_path"], sep="\t")
        sample_accession = str(row["sample_accession"])
        for marker_call, sub in donor_supported.groupby("marker_call", sort=True):
            idx = sorted(sub["barcode_index_0based"].astype(int).unique())
            if not idx:
                continue
            feature_tables.append(
                feature_scores_for_group(
                    sample_accession,
                    donor_id,
                    marker_call,
                    "cluster_enriched_all",
                    X[idx, :],
                    var,
                )
            )
    feature_scores = pd.concat(feature_tables, ignore_index=True) if feature_tables else pd.DataFrame()
    module_scores = explode_module_scores(feature_scores) if not feature_scores.empty else pd.DataFrame()
    return feature_scores, module_scores


def build_contrasts(module: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if module.empty:
        return pd.DataFrame()
    key_cols = ["feature_type", "model_term_supported", "target_set", "peak_category"]
    for donor_id, sub in module.groupby("donor_id"):
        granule = sub.loc[sub["marker_call"].eq("cerebellar_granule_candidate")]
        if granule.empty:
            continue
        for comparator in COMPARATORS:
            comp = sub.loc[sub["marker_call"].eq(comparator)]
            if comp.empty:
                continue
            merged = granule.merge(comp, on=key_cols, suffixes=("_granule", "_comparator"))
            for rec in merged.to_dict("records"):
                rows.append(
                    {
                        "donor_id": donor_id,
                        "comparator": comparator,
                        "feature_type": rec["feature_type"],
                        "model_term_supported": rec["model_term_supported"],
                        "target_set": rec["target_set"],
                        "peak_category": rec["peak_category"],
                        "n_barcodes_granule": rec["n_barcodes_granule"],
                        "n_barcodes_comparator": rec["n_barcodes_comparator"],
                        "mean_feature_score_granule": rec["mean_feature_score_granule"],
                        "mean_feature_score_comparator": rec["mean_feature_score_comparator"],
                        "delta_granule_minus_comparator": rec["mean_feature_score_granule"]
                        - rec["mean_feature_score_comparator"],
                        "mean_detection_granule": rec["mean_detection_fraction_granule"],
                        "mean_detection_comparator": rec["mean_detection_fraction_comparator"],
                        "delta_detection_granule_minus_comparator": rec["mean_detection_fraction_granule"]
                        - rec["mean_detection_fraction_comparator"],
                    }
                )
    contrasts = pd.DataFrame(rows)
    if contrasts.empty:
        return contrasts
    pooled = (
        contrasts.groupby(["comparator", "feature_type", "model_term_supported", "target_set", "peak_category"], dropna=False)
        .agg(
            n_donors=("donor_id", "nunique"),
            mean_delta_score=("delta_granule_minus_comparator", "mean"),
            median_delta_score=("delta_granule_minus_comparator", "median"),
            mean_delta_detection=("delta_detection_granule_minus_comparator", "mean"),
            median_delta_detection=("delta_detection_granule_minus_comparator", "median"),
        )
        .reset_index()
    )
    pooled.insert(0, "donor_id", "pooled")
    pooled["n_barcodes_granule"] = ""
    pooled["n_barcodes_comparator"] = ""
    pooled["mean_feature_score_granule"] = ""
    pooled["mean_feature_score_comparator"] = ""
    pooled["delta_granule_minus_comparator"] = pooled["mean_delta_score"]
    pooled["mean_detection_granule"] = ""
    pooled["mean_detection_comparator"] = ""
    pooled["delta_detection_granule_minus_comparator"] = pooled["mean_delta_detection"]
    pooled = pooled[
        [
            "donor_id",
            "comparator",
            "feature_type",
            "model_term_supported",
            "target_set",
            "peak_category",
            "n_barcodes_granule",
            "n_barcodes_comparator",
            "mean_feature_score_granule",
            "mean_feature_score_comparator",
            "delta_granule_minus_comparator",
            "mean_detection_granule",
            "mean_detection_comparator",
            "delta_detection_granule_minus_comparator",
            "n_donors",
            "mean_delta_score",
            "median_delta_score",
            "mean_delta_detection",
            "median_delta_detection",
        ]
    ]
    contrasts["n_donors"] = ""
    contrasts["mean_delta_score"] = ""
    contrasts["median_delta_score"] = ""
    contrasts["mean_delta_detection"] = ""
    contrasts["median_delta_detection"] = ""
    return pd.concat([contrasts, pooled], ignore_index=True, sort=False)


def write_markdown(rules: pd.DataFrame, supported: pd.DataFrame, module: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    supported_rules = rules.loc[~rules["cluster_support_tier"].eq("not_supported")]
    group_counts = (
        supported.groupby("marker_call").size().sort_values(ascending=False)
        if not supported.empty
        else pd.Series(dtype=int)
    )
    group_lines = [f"- `{idx}`: {int(value)} supported barcodes." for idx, value in group_counts.items()]
    pooled = contrasts.loc[contrasts["donor_id"].eq("pooled")].copy() if not contrasts.empty else pd.DataFrame()
    top = pooled.sort_values("mean_delta_score", ascending=False).head(6) if not pooled.empty else pd.DataFrame()
    top_lines = [
        "- `{}` / `{}` / `{}` / `{}`: delta {:.4g}.".format(
            row.comparator,
            row.feature_type,
            row.model_term_supported,
            row.target_set,
            float(row.mean_delta_score),
        )
        for row in top.itertuples()
    ]

    lines = [
        "# GSE322785 Cluster-Supported Epigenomic Scoring",
        "",
        "Date built: 2026-06-26",
        "",
        "## Support Rule",
        "",
        f"Supported marker clusters required at least {MIN_CLUSTER_MARKER_BARCODES} barcodes of the marker call, cluster fraction at least {MIN_FRACTION_OF_CLUSTER}, and at least {MIN_FOLD_ENRICHMENT}-fold enrichment over the donor-level marker-call frequency. Ambiguous and low-information calls were excluded.",
        "",
        "## Scope",
        "",
        f"- Supported cluster-marker rules: {len(supported_rules)}.",
        f"- Supported barcodes: {len(supported)}.",
        f"- Supported marker groups: {supported['marker_call'].nunique() if not supported.empty else 0}.",
        f"- Module-score rows: {len(module)}.",
        f"- Granule-versus-supported-comparator contrast rows: {len(contrasts)}.",
        "",
        "## Supported Barcode Counts",
        "",
        *(group_lines or ["- No supported barcodes passed the rule."]),
        "",
        "## Top Pooled Granule-Positive Supported Contrasts",
        "",
        *(top_lines or ["- No pooled granule/comparator contrasts available."]),
        "",
        "## Interpretation",
        "",
        "This stricter layer reduces false confidence by using only marker calls located in donor-specific clusters enriched for that same call. It is better suited for sensitivity analysis than the broader provisional marker-score layer, but it still does not replace source-author taxonomy or full multimodal clustering.",
        "",
        "## Outputs",
        "",
        f"- Cluster support rules: `{OUT_RULES.relative_to(ROOT)}`",
        f"- Supported barcodes: `{OUT_BARCODES.relative_to(ROOT)}`",
        f"- Supported feature scores: `{OUT_FEATURE.relative_to(ROOT)}`",
        f"- Supported module scores: `{OUT_MODULE.relative_to(ROOT)}`",
        f"- Supported granule/comparator contrasts: `{OUT_CONTRASTS.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    assignments = pd.read_csv(CLUSTER_ASSIGNMENTS, sep="\t")
    enrichment = pd.read_csv(CLUSTER_ENRICHMENT, sep="\t")
    rules = build_rules(assignments, enrichment)
    supported = build_supported_barcodes(assignments, rules)
    feature_scores, module_scores = score_supported_groups(supported)
    contrasts = build_contrasts(module_scores)

    rules.to_csv(OUT_RULES, sep="\t", index=False)
    supported.to_csv(OUT_BARCODES, sep="\t", index=False, compression="gzip")
    feature_scores.to_csv(OUT_FEATURE, sep="\t", index=False, compression="gzip")
    module_scores.to_csv(OUT_MODULE, sep="\t", index=False)
    contrasts.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    write_markdown(rules, supported, module_scores, contrasts)


if __name__ == "__main__":
    main()
