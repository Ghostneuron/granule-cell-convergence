#!/usr/bin/env python3
"""Map human GSE322785 cerebellar H5 features to epigenomic target genes."""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE322785"
RESULTS = ROOT / "Project/results"

TARGETS = RESULTS / "epigenomic_extension_regulatory_targets.tsv"
PLAN = RESULTS / "gse322785_cerebellar_multiome_download_plan.tsv"

OUT_H5_INVENTORY = RESULTS / "gse322785_human_h5_feature_inventory.tsv"
OUT_GENE_SUMMARY = RESULTS / "gse322785_human_h5_epigenomic_gene_summary.tsv"
OUT_PEAK_TARGETS = RESULTS / "gse322785_human_h5_epigenomic_peak_targets.tsv"
OUT_MANIFEST = RESULTS / "gse322785_human_h5_epigenomic_selective_manifest.tsv"
OUT_MD = RESULTS / "gse322785_human_h5_epigenomic_targeting.md"

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


def parse_interval(interval: str) -> tuple[str, int, int] | None:
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", str(interval))
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def read_targets() -> pd.DataFrame:
    targets = pd.read_csv(TARGETS, sep="\t")
    rows = []
    for gene, sub in targets.groupby("gene", sort=True):
        priority = min(sub["priority"].astype(str), key=lambda x: PRIORITY_ORDER.get(x, 99))
        rows.append(
            {
                "gene": str(gene).upper(),
                "target_sources": ";".join(sorted(sub["target_source"].astype(str).unique())),
                "target_sets": ";".join(sorted(sub["target_set"].astype(str).unique())),
                "model_terms_supported": ";".join(sorted(sub["model_term_supported"].astype(str).unique())),
                "best_priority": priority,
                "is_tf_or_regulatory_candidate": bool(sub["is_tf_or_regulatory_candidate"].astype(bool).any()),
                "n_target_rows": len(sub),
            }
        )
    return pd.DataFrame(rows)


def human_h5_files() -> pd.DataFrame:
    plan = pd.read_csv(PLAN, sep="\t")
    human = plan.loc[plan["download_tier"].eq("tier1_human_cerebellar_h5")].copy()
    human["local_path"] = human["file_name"].map(lambda x: str(BASE / x))
    human["local_exists"] = human["local_path"].map(lambda x: Path(x).exists())
    return human


def read_h5_features(path: Path) -> tuple[pd.DataFrame, tuple[int, int], int]:
    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][:])
        nnz = int(matrix["data"].shape[0])
        features = pd.DataFrame(
            {
                "feature_index_0based": range(len(matrix["features/name"])),
                "feature_id": [x.decode() for x in matrix["features/id"][:]],
                "feature_name": [x.decode() for x in matrix["features/name"][:]],
                "feature_type": [x.decode() for x in matrix["features/feature_type"][:]],
                "genome": [x.decode() for x in matrix["features/genome"][:]],
                "interval": [x.decode() for x in matrix["features/interval"][:]],
            }
        )
    parsed = features["interval"].map(parse_interval)
    features["chrom"] = parsed.map(lambda x: x[0] if x else "")
    features["start"] = parsed.map(lambda x: x[1] if x else pd.NA)
    features["end"] = parsed.map(lambda x: x[2] if x else pd.NA)
    features["gene"] = features["feature_name"].str.upper()
    return features, shape, nnz


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


def map_one_file(sample: pd.Series, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    path = Path(sample["local_path"])
    features, shape, nnz = read_h5_features(path)
    genes = features.loc[features["feature_type"].eq("Gene Expression")].copy()
    peaks = features.loc[features["feature_type"].eq("Peaks")].copy()
    peaks_by_chrom = {chrom: sub for chrom, sub in peaks.groupby("chrom", sort=False)}

    inventory = {
        "sample_accession": sample["sample_accession"],
        "donor_id": sample["donor_id"],
        "organism": sample["organism"],
        "file_name": sample["file_name"],
        "n_matrix_rows": shape[0],
        "n_barcodes": shape[1],
        "n_nonzero": nnz,
        "n_gene_expression_features": int(genes.shape[0]),
        "n_peak_features": int(peaks.shape[0]),
        "local_size_bytes": path.stat().st_size,
    }

    gene_rows = []
    peak_rows = []
    for rec in targets.to_dict("records"):
        gene = rec["gene"]
        matching = genes.loc[genes["gene"].eq(gene)]
        if matching.empty:
            gene_rows.append(
                {
                    **rec,
                    "sample_accession": sample["sample_accession"],
                    "donor_id": sample["donor_id"],
                    "present_in_h5": False,
                    "gene_feature_index_0based": "",
                    "gene_feature_id": "",
                    "gene_feature_name": "",
                    "chrom": "",
                    "gene_start": "",
                    "gene_end": "",
                    "n_peaks_total_100kb": 0,
                    "n_gene_body_overlap": 0,
                    "n_gene_edge_2kb": 0,
                    "n_proximal_10kb": 0,
                    "n_distal_100kb": 0,
                }
            )
            continue
        g = matching.iloc[0]
        chrom = g["chrom"]
        gene_start = int(g["start"])
        gene_end = int(g["end"])
        nearby = peaks_by_chrom.get(chrom, pd.DataFrame())
        if not nearby.empty:
            nearby = nearby.loc[
                nearby["start"].le(gene_end + WINDOWS["distal_100kb"])
                & nearby["end"].ge(gene_start - WINDOWS["distal_100kb"])
            ]
        counts = {
            "gene_body_overlap": 0,
            "gene_edge_2kb": 0,
            "proximal_10kb": 0,
            "distal_100kb": 0,
        }
        for peak in nearby.to_dict("records"):
            category = classify_peak(int(peak["start"]), int(peak["end"]), gene_start, gene_end)
            if category == "outside_window":
                continue
            counts[category] += 1
            dist = distance_to_interval(int(peak["start"]), int(peak["end"]), gene_start, gene_end)
            peak_rows.append(
                {
                    **rec,
                    "sample_accession": sample["sample_accession"],
                    "donor_id": sample["donor_id"],
                    "gene_feature_index_0based": int(g["feature_index_0based"]),
                    "gene_feature_id": g["feature_id"],
                    "gene_feature_name": g["feature_name"],
                    "chrom": chrom,
                    "gene_start": gene_start,
                    "gene_end": gene_end,
                    "peak_feature_index_0based": int(peak["feature_index_0based"]),
                    "peak_feature_id": peak["feature_id"],
                    "peak_chrom": peak["chrom"],
                    "peak_start": int(peak["start"]),
                    "peak_end": int(peak["end"]),
                    "peak_category": category,
                    "distance_to_gene_body_bp": int(dist),
                    "peak_width_bp": int(peak["end"]) - int(peak["start"]) + 1,
                    "extraction_priority": extraction_priority(rec["best_priority"], category),
                }
            )
        gene_rows.append(
            {
                **rec,
                "sample_accession": sample["sample_accession"],
                "donor_id": sample["donor_id"],
                "present_in_h5": True,
                "gene_feature_index_0based": int(g["feature_index_0based"]),
                "gene_feature_id": g["feature_id"],
                "gene_feature_name": g["feature_name"],
                "chrom": chrom,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "n_peaks_total_100kb": sum(counts.values()),
                "n_gene_body_overlap": counts["gene_body_overlap"],
                "n_gene_edge_2kb": counts["gene_edge_2kb"],
                "n_proximal_10kb": counts["proximal_10kb"],
                "n_distal_100kb": counts["distal_100kb"],
            }
        )

    gene_summary = pd.DataFrame(gene_rows)
    peak_targets = pd.DataFrame(peak_rows)
    manifest = build_manifest(gene_summary, peak_targets)
    return gene_summary, peak_targets, manifest, inventory


def build_manifest(gene_summary: pd.DataFrame, peak_targets: pd.DataFrame) -> pd.DataFrame:
    genes = gene_summary.loc[gene_summary["present_in_h5"].astype(bool)].copy()
    gene_manifest = genes[
        [
            "sample_accession",
            "donor_id",
            "gene_feature_index_0based",
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
            "gene_feature_index_0based": "feature_index_0based",
            "gene_feature_id": "feature_id",
            "gene_feature_name": "feature_name",
            "chrom": "feature_chrom",
            "gene_start": "feature_start",
            "gene_end": "feature_end",
        }
    )
    gene_manifest["feature_type"] = "Gene Expression"
    gene_manifest["linked_gene"] = gene_manifest["gene"]
    gene_manifest["peak_category"] = ""
    gene_manifest["selection_reason"] = "target_gene_expression_row"

    if peak_targets.empty:
        return gene_manifest

    peak_manifest = peak_targets[
        [
            "sample_accession",
            "donor_id",
            "peak_feature_index_0based",
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
            "peak_feature_index_0based": "feature_index_0based",
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
    return manifest.drop_duplicates(["sample_accession", "feature_index_0based", "linked_gene", "selection_reason"])


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    targets = read_targets()
    samples = human_h5_files()
    gene_tables = []
    peak_tables = []
    manifest_tables = []
    inventory_rows = []

    for sample in samples.loc[samples["local_exists"]].to_dict("records"):
        gene_summary, peak_targets, manifest, inventory = map_one_file(pd.Series(sample), targets)
        gene_tables.append(gene_summary)
        peak_tables.append(peak_targets)
        manifest_tables.append(manifest)
        inventory_rows.append(inventory)

    inventory = pd.DataFrame(inventory_rows)
    gene_summary = pd.concat(gene_tables, ignore_index=True) if gene_tables else pd.DataFrame()
    peak_targets = pd.concat(peak_tables, ignore_index=True) if peak_tables else pd.DataFrame()
    manifest = pd.concat(manifest_tables, ignore_index=True) if manifest_tables else pd.DataFrame()

    inventory.to_csv(OUT_H5_INVENTORY, sep="\t", index=False)
    gene_summary.to_csv(OUT_GENE_SUMMARY, sep="\t", index=False)
    peak_targets.to_csv(OUT_PEAK_TARGETS, sep="\t", index=False)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)

    n_samples = inventory.shape[0]
    n_present_rows = int(gene_summary["present_in_h5"].sum()) if not gene_summary.empty else 0
    n_gene_tests = len(gene_summary)
    n_unique_genes_all = targets["gene"].nunique()
    n_unique_genes_any = gene_summary.loc[gene_summary["present_in_h5"].astype(bool), "gene"].nunique()
    n_peak_rows = len(peak_targets)
    n_unique_manifest = manifest[["sample_accession", "feature_index_0based"]].drop_duplicates().shape[0] if not manifest.empty else 0
    n_total_rows = int(inventory["n_matrix_rows"].sum()) if not inventory.empty else 0
    manifest_fraction = n_unique_manifest / n_total_rows if n_total_rows else 0
    top = (
        gene_summary.groupby("gene", as_index=False)["n_peaks_total_100kb"]
        .sum()
        .sort_values("n_peaks_total_100kb", ascending=False)
        .head(8)
    )
    top_lines = [f"- `{r.gene}`: {int(r.n_peaks_total_100kb)} nearby peaks across human H5 files." for r in top.itertuples()]
    lines = [
        "# GSE322785 Human H5 Epigenomic Targeting",
        "",
        "Date built: 2026-06-26",
        "",
        "## Purpose",
        "",
        "This analysis maps the manuscript epigenomic target genes to gene-expression and nearby ATAC peak rows in the downloaded human adult cerebellar multiome H5 files.",
        "",
        "## Scope",
        "",
        f"- Human H5 files analyzed: {n_samples}.",
        f"- Target genes represented in any human H5: {n_unique_genes_any}/{n_unique_genes_all}.",
        f"- Per-sample target gene rows present: {n_present_rows}/{n_gene_tests}.",
        f"- Peak-target rows within gene body plus 100 kb: {n_peak_rows}.",
        f"- Selective manifest sample-feature rows: {n_unique_manifest} ({manifest_fraction:.2%} of summed H5 matrix rows across analyzed files).",
        "",
        "## Highest Peak-Count Genes",
        "",
        *top_lines,
        "",
        "## Interpretation",
        "",
        "The downloaded human cerebellar H5 files are usable for the epigenomic compatibility extension. They still lack manuscript-ready granule/Purkinje labels until cell-type annotation or label transfer is performed.",
        "",
        "## Outputs",
        "",
        f"- H5 inventory: `{OUT_H5_INVENTORY.relative_to(ROOT)}`",
        f"- Gene summary: `{OUT_GENE_SUMMARY.relative_to(ROOT)}`",
        f"- Peak targets: `{OUT_PEAK_TARGETS.relative_to(ROOT)}`",
        f"- Selective manifest: `{OUT_MANIFEST.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
