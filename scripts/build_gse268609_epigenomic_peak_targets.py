#!/usr/bin/env python3
"""Map GSE268609 ATAC peak rows to epigenomic-extension target genes.

This script does not require the large multiome count matrix. It uses the
already-downloaded GSE268609 feature table to identify which ATAC peak rows sit
near curated target genes, so a later matrix extraction can be selective.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

FEATURES = EXTERNAL / "GEO/GSE268609/GSE268609_features.tsv.gz"
TARGETS = RESULTS / "epigenomic_extension_regulatory_targets.tsv"

OUT_PEAK_TARGETS = RESULTS / "gse268609_epigenomic_peak_targets.tsv"
OUT_GENE_SUMMARY = RESULTS / "gse268609_epigenomic_peak_gene_summary.tsv"
OUT_MANIFEST = RESULTS / "gse268609_epigenomic_selective_extraction_manifest.tsv"
OUT_MD = RESULTS / "gse268609_epigenomic_peak_targeting.md"

WINDOWS = {
    "gene_edge_2kb": 2_000,
    "proximal_10kb": 10_000,
    "distal_100kb": 100_000,
}

PRIORITY_ORDER = {
    "highest": 0,
    "high": 1,
    "medium_high": 2,
    "medium": 3,
    "low": 4,
}


def parse_features() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    with gzip.open(FEATURES, "rt") as handle:
        for row_idx, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            feature_id, name, feature_type, chrom, start, end = parts[:6]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            rows.append(
                {
                    "feature_row_1based": row_idx,
                    "feature_id": feature_id,
                    "feature_name": name,
                    "feature_type": feature_type,
                    "chrom": chrom,
                    "start": start_i,
                    "end": end_i,
                }
            )
    features = pd.DataFrame(rows)
    genes = features.loc[features["feature_type"].eq("Gene Expression")].copy()
    genes["gene"] = genes["feature_name"].str.upper()
    peaks = features.loc[features["feature_type"].eq("Peaks")].copy()
    return genes, peaks


def collapse_targets(targets: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    present = targets.loc[targets["present_as_gse268609_gene_expression_feature"].astype(bool)].copy()
    groups = []
    for gene, sub in present.groupby("gene", sort=True):
        gene_upper = str(gene).upper()
        gene_rows = genes.loc[genes["gene"].eq(gene_upper)]
        if gene_rows.empty:
            continue
        gene_row = gene_rows.iloc[0]
        priority = min(sub["priority"].astype(str), key=lambda x: PRIORITY_ORDER.get(x, 99))
        groups.append(
            {
                "gene": gene_upper,
                "gene_feature_row_1based": int(gene_row["feature_row_1based"]),
                "gene_feature_id": gene_row["feature_id"],
                "gene_feature_name": gene_row["feature_name"],
                "chrom": gene_row["chrom"],
                "gene_start": int(gene_row["start"]),
                "gene_end": int(gene_row["end"]),
                "target_sources": ";".join(sorted(sub["target_source"].astype(str).unique())),
                "target_sets": ";".join(sorted(sub["target_set"].astype(str).unique())),
                "model_terms_supported": ";".join(sorted(sub["model_term_supported"].astype(str).unique())),
                "best_priority": priority,
                "is_tf_or_regulatory_candidate": bool(sub["is_tf_or_regulatory_candidate"].astype(bool).any()),
                "n_target_rows": len(sub),
            }
        )
    return pd.DataFrame(groups)


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def distance_to_interval(start: int, end: int, ref_start: int, ref_end: int) -> int:
    if overlap(start, end, ref_start, ref_end):
        return 0
    if end < ref_start:
        return ref_start - end
    return start - ref_end


def classify_peak(peak_start: int, peak_end: int, gene_start: int, gene_end: int) -> str:
    if overlap(peak_start, peak_end, gene_start, gene_end):
        return "gene_body_overlap"
    for label, window in WINDOWS.items():
        if overlap(peak_start, peak_end, gene_start - window, gene_end + window):
            return label
    return "outside_window"


def map_peaks(target_genes: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    peak_rows = []
    peaks_by_chrom = {chrom: sub.copy() for chrom, sub in peaks.groupby("chrom", sort=False)}
    for gene in target_genes.to_dict("records"):
        chrom_peaks = peaks_by_chrom.get(gene["chrom"])
        if chrom_peaks is None or chrom_peaks.empty:
            continue
        query_start = int(gene["gene_start"]) - WINDOWS["distal_100kb"]
        query_end = int(gene["gene_end"]) + WINDOWS["distal_100kb"]
        nearby = chrom_peaks.loc[
            (chrom_peaks["start"].le(query_end)) & (chrom_peaks["end"].ge(query_start))
        ]
        for peak in nearby.to_dict("records"):
            category = classify_peak(
                int(peak["start"]),
                int(peak["end"]),
                int(gene["gene_start"]),
                int(gene["gene_end"]),
            )
            if category == "outside_window":
                continue
            distance = distance_to_interval(
                int(peak["start"]),
                int(peak["end"]),
                int(gene["gene_start"]),
                int(gene["gene_end"]),
            )
            peak_rows.append(
                {
                    **gene,
                    "peak_feature_row_1based": int(peak["feature_row_1based"]),
                    "peak_feature_id": peak["feature_id"],
                    "peak_chrom": peak["chrom"],
                    "peak_start": int(peak["start"]),
                    "peak_end": int(peak["end"]),
                    "peak_category": category,
                    "distance_to_gene_body_bp": int(distance),
                    "peak_width_bp": int(peak["end"]) - int(peak["start"]) + 1,
                    "extraction_priority": extraction_priority(gene["best_priority"], category),
                }
            )
    return pd.DataFrame(peak_rows)


def extraction_priority(gene_priority: str, peak_category: str) -> str:
    if gene_priority == "highest" and peak_category in {"gene_body_overlap", "gene_edge_2kb"}:
        return "tier1_core_peak"
    if gene_priority in {"highest", "high"} and peak_category in {"gene_body_overlap", "gene_edge_2kb"}:
        return "high_priority_promoter_gene_peak"
    if peak_category == "gene_body_overlap":
        return "gene_body_peak"
    if peak_category == "gene_edge_2kb":
        return "gene_edge_peak"
    if peak_category == "proximal_10kb":
        return "proximal_peak"
    return "distal_context_peak"


def summarize_gene_peaks(target_genes: pd.DataFrame, peak_targets: pd.DataFrame) -> pd.DataFrame:
    if peak_targets.empty:
        summary = target_genes.copy()
        for col in [
            "n_peaks_total_100kb",
            "n_gene_body_overlap",
            "n_gene_edge_2kb",
            "n_proximal_10kb",
            "n_distal_100kb",
        ]:
            summary[col] = 0
        return summary

    counts = (
        peak_targets.pivot_table(
            index="gene",
            columns="peak_category",
            values="peak_feature_row_1based",
            aggfunc="count",
            fill_value=0,
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for col in ["gene_body_overlap", "gene_edge_2kb", "proximal_10kb", "distal_100kb"]:
        if col not in counts.columns:
            counts[col] = 0
    counts["n_peaks_total_100kb"] = counts[
        ["gene_body_overlap", "gene_edge_2kb", "proximal_10kb", "distal_100kb"]
    ].sum(axis=1)
    counts = counts.rename(
        columns={
            "gene_body_overlap": "n_gene_body_overlap",
            "gene_edge_2kb": "n_gene_edge_2kb",
            "proximal_10kb": "n_proximal_10kb",
            "distal_100kb": "n_distal_100kb",
        }
    )
    summary = target_genes.merge(counts, on="gene", how="left")
    count_cols = [
        "n_peaks_total_100kb",
        "n_gene_body_overlap",
        "n_gene_edge_2kb",
        "n_proximal_10kb",
        "n_distal_100kb",
    ]
    summary[count_cols] = summary[count_cols].fillna(0).astype(int)
    return summary.sort_values(["best_priority", "gene"], key=priority_sort)


def priority_sort(series: pd.Series) -> pd.Series:
    if series.name == "best_priority":
        return series.map(lambda x: PRIORITY_ORDER.get(str(x), 99))
    return series


def build_manifest(target_genes: pd.DataFrame, peak_targets: pd.DataFrame) -> pd.DataFrame:
    gene_manifest = target_genes[
        [
            "gene_feature_row_1based",
            "gene_feature_id",
            "gene_feature_name",
            "gene",
            "chrom",
            "gene_start",
            "gene_end",
            "best_priority",
            "model_terms_supported",
            "target_sets",
        ]
    ].rename(
        columns={
            "gene_feature_row_1based": "feature_row_1based",
            "gene_feature_id": "feature_id",
            "gene_feature_name": "feature_name",
            "chrom": "feature_chrom",
            "gene_start": "feature_start",
            "gene_end": "feature_end",
        }
    )
    gene_manifest["feature_type"] = "Gene Expression"
    gene_manifest["selection_reason"] = "target_gene_expression_row"
    gene_manifest["linked_gene"] = gene_manifest["gene"]
    gene_manifest["peak_category"] = ""

    peak_manifest = peak_targets[
        [
            "peak_feature_row_1based",
            "peak_feature_id",
            "gene",
            "peak_chrom",
            "peak_start",
            "peak_end",
            "best_priority",
            "model_terms_supported",
            "target_sets",
            "peak_category",
            "extraction_priority",
        ]
    ].rename(
        columns={
            "peak_feature_row_1based": "feature_row_1based",
            "peak_feature_id": "feature_id",
            "peak_chrom": "feature_chrom",
            "peak_start": "feature_start",
            "peak_end": "feature_end",
            "gene": "linked_gene",
            "extraction_priority": "selection_reason",
        }
    )
    peak_manifest["feature_type"] = "Peaks"
    peak_manifest["feature_name"] = peak_manifest["feature_id"]
    peak_manifest["gene"] = peak_manifest["linked_gene"]

    manifest = pd.concat([gene_manifest, peak_manifest], ignore_index=True, sort=False)
    manifest = manifest.drop_duplicates(["feature_row_1based", "linked_gene", "selection_reason"])
    cols = [
        "feature_row_1based",
        "feature_id",
        "feature_name",
        "feature_type",
        "feature_chrom",
        "feature_start",
        "feature_end",
        "linked_gene",
        "best_priority",
        "model_terms_supported",
        "target_sets",
        "peak_category",
        "selection_reason",
    ]
    return manifest[cols].sort_values(["feature_row_1based", "linked_gene"])


def main() -> None:
    genes, peaks = parse_features()
    targets = pd.read_csv(TARGETS, sep="\t")
    target_genes = collapse_targets(targets, genes)
    peak_targets = map_peaks(target_genes, peaks)
    gene_summary = summarize_gene_peaks(target_genes, peak_targets)
    manifest = build_manifest(target_genes, peak_targets)

    peak_targets.to_csv(OUT_PEAK_TARGETS, sep="\t", index=False)
    gene_summary.to_csv(OUT_GENE_SUMMARY, sep="\t", index=False)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)

    n_genes = len(target_genes)
    n_peak_rows = len(peak_targets)
    n_unique_peak_features = peak_targets["peak_feature_row_1based"].nunique() if n_peak_rows else 0
    n_gene_rows = len(target_genes)
    n_manifest_features = manifest["feature_row_1based"].nunique()
    n_full_features = len(genes) + len(peaks)
    manifest_fraction = n_manifest_features / n_full_features if n_full_features else 0.0
    top = gene_summary.sort_values("n_peaks_total_100kb", ascending=False).head(8)
    top_lines = [
        f"- `{r.gene}`: {r.n_peaks_total_100kb} nearby peaks "
        f"({r.n_gene_body_overlap} gene-body, {r.n_gene_edge_2kb} edge-2kb, "
        f"{r.n_proximal_10kb} proximal, {r.n_distal_100kb} distal)."
        for r in top.itertuples()
    ]
    lines = [
        "# GSE268609 Epigenomic Peak Targeting",
        "",
        "Date built: 2026-06-26",
        "",
        "## Purpose",
        "",
        "This no-heavy-download step maps curated epigenomic target genes to nearby `GSE268609` ATAC peak feature rows. It prepares a selective extraction manifest for future matrix processing.",
        "",
        "## Scope",
        "",
        f"- Target genes with gene-expression feature rows: {n_genes}.",
        f"- Peak-target rows within gene body plus 100 kb: {n_peak_rows}.",
        f"- Unique nearby ATAC peak features: {n_unique_peak_features}.",
        f"- Target gene-expression rows: {n_gene_rows}.",
        f"- Unique feature rows in selective manifest: {n_manifest_features} ({manifest_fraction:.2%} of all `GSE268609` feature rows).",
        "",
        "## Highest Peak-Count Genes",
        "",
        *top_lines,
        "",
        "## Interpretation",
        "",
        "The manifest makes later peak-count processing substantially narrower than the full 177,931-row GSE268609 matrix. It still does not fit chromatin accessibility because the large count matrix is not currently local.",
        "",
        "## Outputs",
        "",
        f"- Peak targets: `{OUT_PEAK_TARGETS.relative_to(ROOT)}`",
        f"- Gene summary: `{OUT_GENE_SUMMARY.relative_to(ROOT)}`",
        f"- Selective extraction manifest: `{OUT_MANIFEST.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
