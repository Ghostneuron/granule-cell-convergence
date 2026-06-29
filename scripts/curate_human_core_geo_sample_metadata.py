#!/usr/bin/env python3
"""Curate GEO sample metadata for the built human core components."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data/GEO"
RESULTS = ROOT / "Project/results"

SOFT_FILES = {
    "GSE185277": EXTERNAL / "GSE185277/GSE185277_family.soft.gz",
    "GSE185553": EXTERNAL / "GSE185553/GSE185553_family.soft.gz",
    "GSE186538": EXTERNAL / "GSE186538/GSE186538_family.soft.gz",
}

QC_CELLS = RESULTS / "human_core_harmonized_cell_qc.tsv.gz"
SEED_SUMMARY = RESULTS / "human_seed_sparse_object_summary.tsv"
GSE186538_CELL_META = ROOT / "Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/cell_metadata.tsv.gz"
MARKER_CALL_SUMMARY = RESULTS / "human_core_marker_validation_call_summary.tsv"

OUT_SAMPLE_METADATA = RESULTS / "human_core_geo_sample_metadata.tsv"
OUT_COMPONENT_METADATA = RESULTS / "human_core_component_metadata_curated.tsv"
OUT_DONOR_SUMMARY = RESULTS / "human_core_gse186538_dg_donor_summary.tsv"
OUT_ENRICHED_CELLS = RESULTS / "human_core_enriched_cell_metadata.tsv.gz"
OUT_MD = RESULTS / "human_core_geo_metadata_curation_summary.md"


def norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def norm_file(value: str) -> str:
    name = Path(urlparse(value).path).name if "://" in value else Path(value).name
    name = re.sub(r"\.(txt|dge|deg|barcode)(\.gz)?$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[^A-Za-z0-9]+", "", name).lower()
    return name


def basename(value: str) -> str:
    return Path(urlparse(value).path).name if "://" in value else Path(value).name


def first_nonempty(values: list[str]) -> str:
    for value in values:
        if str(value).strip():
            return str(value).strip()
    return ""


def split_multi(value: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in str(value).split("/") if part.strip()]


def split_repeated(value: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def parse_age_years(value: str) -> float:
    text = str(value).strip().lower()
    if not text:
        return np.nan
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(yrs?|years?|y)\b", text)
    if match:
        return float(match.group(1))
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(mos?|months?)\b", text)
    if match:
        return float(match.group(1)) / 12.0
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(days?|d)\b", text)
    if match:
        return float(match.group(1)) / 365.25
    return np.nan


def parse_soft(dataset: str, path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    series: dict[str, str] = {"dataset": dataset, "soft_path": str(path.relative_to(ROOT))}
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    repeated: dict[str, list[str]] = defaultdict(list)

    def flush_current() -> None:
        nonlocal current, repeated
        if current is None:
            return
        for key, values in repeated.items():
            current[key] = "; ".join(values)
        samples.append(current)
        current = None
        repeated = defaultdict(list)

    with gzip.open(path, "rt", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                flush_current()
                current = {"dataset": dataset, "sample_accession": line.split(" = ", 1)[1]}
                continue
            if line.startswith("^") and not line.startswith("^SAMPLE = "):
                flush_current()
                continue
            if line.startswith("!Series_"):
                key, value = line[1:].split(" = ", 1)
                series[norm_key(key.replace("Series_", ""))] = value
                continue
            if current is None or not line.startswith("!Sample_") or " = " not in line:
                continue

            key, value = line[1:].split(" = ", 1)
            sample_key = norm_key(key.replace("Sample_", "").replace("_ch1", ""))
            if sample_key.startswith("characteristics"):
                if ":" in value:
                    char_key, char_value = value.split(":", 1)
                    current[norm_key(char_key)] = char_value.strip()
                else:
                    repeated["characteristics"].append(value)
            elif sample_key.startswith("supplementary_file"):
                repeated["supplementary_files"].append(value)
            elif sample_key in {"extract_protocol", "data_processing", "relation"}:
                repeated[sample_key].append(value)
            else:
                current[sample_key] = value
    flush_current()
    return series, samples


def sample_number_from_component(component_id: str) -> str:
    match = re.search(r"[Ss]ample[-_]?(\d+)|sample(\d+)", component_id)
    if not match:
        return ""
    return first_nonempty([item for item in match.groups() if item])


def donor_from_title(title: str) -> str:
    match = re.match(r"(HSB\d+)", str(title))
    return match.group(1) if match else ""


def specimen_metadata(sample: pd.Series, sample_number: str) -> dict[str, object]:
    title = str(sample.get("title", ""))
    specimen_ids = []
    match = re.search(r"Specimen[_ ](.+)$", title, flags=re.IGNORECASE)
    if match:
        specimen_ids = split_multi(match.group(1))

    stages = split_multi(str(sample.get("development_stage", "")))
    ages = split_multi(str(sample.get("age", "")))
    idx = 0
    match_method = "gsm_only"
    if sample_number and specimen_ids and sample_number in specimen_ids:
        idx = specimen_ids.index(sample_number)
        match_method = "gsm_and_specimen_number"
    elif len(specimen_ids) == 1:
        match_method = "single_specimen_gsm"
    elif sample_number:
        match_method = "gsm_matched_specimen_number_unresolved"

    age = ages[idx] if idx < len(ages) else first_nonempty(ages)
    stage = stages[idx] if idx < len(stages) else first_nonempty(stages)
    specimen_id = specimen_ids[idx] if idx < len(specimen_ids) else sample_number
    return {
        "specimen_id": specimen_id,
        "development_stage_curated": stage,
        "age_curated": age,
        "age_years_curated": parse_age_years(age),
        "sample_match_method": match_method,
    }


def build_sample_metadata() -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    all_samples: list[dict[str, object]] = []
    series_by_dataset: dict[str, dict[str, str]] = {}
    for dataset, path in SOFT_FILES.items():
        series, samples = parse_soft(dataset, path)
        series_by_dataset[dataset] = series
        for sample in samples:
            sample["series_title"] = series.get("title", "")
            sample["series_pubmed_id"] = series.get("pubmed_id", "")
            sample["soft_path"] = str(path.relative_to(ROOT))
            supplementary_files = split_repeated(str(sample.get("supplementary_files", "")))
            sample["supplementary_file_basenames"] = ";".join([basename(item) for item in supplementary_files])
            sample["supplementary_file_norms"] = ";".join([norm_file(item) for item in supplementary_files])
            sample["age_years"] = parse_age_years(str(sample.get("age", "")))
            sample["donor_hint"] = donor_from_title(str(sample.get("title", "")))
        all_samples.extend(samples)
    df = pd.DataFrame(all_samples)
    preferred_cols = [
        "dataset",
        "sample_accession",
        "title",
        "source_name",
        "tissue",
        "development_stage",
        "age",
        "age_years",
        "sex",
        "molecule_subtype",
        "donor_hint",
        "supplementary_file_basenames",
        "series_title",
        "series_pubmed_id",
        "soft_path",
    ]
    for col in preferred_cols:
        if col not in df:
            df[col] = ""
    extra_cols = [col for col in sorted(df.columns) if col not in preferred_cols + ["supplementary_file_norms"]]
    return df[preferred_cols + extra_cols + ["supplementary_file_norms"]], series_by_dataset


def call_counts_by_component() -> pd.DataFrame:
    if not MARKER_CALL_SUMMARY.exists():
        return pd.DataFrame()
    calls = pd.read_csv(MARKER_CALL_SUMMARY, sep="\t")
    pivot = calls.pivot_table(
        index=["dataset", "component_id"],
        columns="marker_call",
        values="n_cells",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    pivot.columns.name = None
    return pivot


def curate_seed_components(sample_df: pd.DataFrame, marker_counts: pd.DataFrame) -> pd.DataFrame:
    seed = pd.read_csv(SEED_SUMMARY, sep="\t")
    sample_lookup = sample_df.set_index(["dataset", "sample_accession"], drop=False)
    rows: list[dict[str, object]] = []

    for _, row in seed.iterrows():
        gsm_match = re.search(r"(GSM\d+)", row["library_id"])
        gsm = gsm_match.group(1) if gsm_match else ""
        sample = sample_lookup.loc[(row["dataset"], gsm)] if (row["dataset"], gsm) in sample_lookup.index else pd.Series(dtype=object)
        sample_number = sample_number_from_component(row["library_id"])
        spec = specimen_metadata(sample, sample_number) if len(sample) else {}
        supp_norms = split_repeated(str(sample.get("supplementary_file_norms", ""))) if len(sample) else []
        source_norm = norm_file(row["source_member"])
        matched_supp = ""
        for supp in supp_norms:
            if supp and (source_norm in supp or supp in source_norm):
                matched_supp = supp
                break

        rows.append(
            {
                "dataset": row["dataset"],
                "component_id": row["library_id"],
                "component_type": "library",
                "source_member": row["source_member"],
                "gsm": gsm,
                "geo_sample_title": sample.get("title", ""),
                "specimen_id": spec.get("specimen_id", sample_number),
                "geo_tissue": sample.get("tissue", sample.get("source_name", "")),
                "geo_development_stage": spec.get("development_stage_curated", sample.get("development_stage", "")),
                "geo_age": spec.get("age_curated", sample.get("age", "")),
                "geo_age_years": spec.get("age_years_curated", sample.get("age_years", np.nan)),
                "geo_sex": sample.get("sex", ""),
                "sample_match_method": spec.get("sample_match_method", "unmatched"),
                "supplementary_file_match": "matched" if matched_supp else "gsm_only_or_unmatched_supplement",
                "n_cells": int(row["n_cells"]),
                "nnz": int(row["nnz"]),
                "total_counts": int(row["total_counts"]),
            }
        )
    component = pd.DataFrame(rows)
    if not marker_counts.empty:
        component = component.merge(marker_counts, on=["dataset", "component_id"], how="left")
    return component


def curate_gse186538_donors(sample_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cell_meta = pd.read_csv(GSE186538_CELL_META, sep="\t")
    human_dg_samples = sample_df[
        (sample_df["dataset"] == "GSE186538")
        & (sample_df["donor_hint"] != "")
        & (sample_df["tissue"].astype(str).str.lower() == "dentate gyrus")
    ].copy()

    donor_rows: list[dict[str, object]] = []
    for donor, group in cell_meta.groupby("samplename", dropna=False):
        donor_samples = human_dg_samples[human_dg_samples["donor_hint"] == donor]
        cluster_counts = group["cluster"].value_counts().to_dict()
        donor_rows.append(
            {
                "dataset": "GSE186538",
                "donor_id": donor,
                "n_dg_gc_cells": len(group),
                "dg_gc_prox1_sgcz_cells": int(cluster_counts.get("DG GC PROX1 SGCZ", 0)),
                "dg_gc_prox1_pdlim5_cells": int(cluster_counts.get("DG GC PROX1 PDLIM5", 0)),
                "geo_gsm": ";".join(donor_samples["sample_accession"].astype(str)),
                "geo_sample_titles": ";".join(donor_samples["title"].astype(str)),
                "geo_tissue": first_nonempty(donor_samples["tissue"].astype(str).tolist()),
                "geo_age": first_nonempty(donor_samples["age"].astype(str).tolist()),
                "geo_age_years": first_nonempty([str(item) for item in donor_samples["age_years"].dropna().unique()]),
                "geo_sex": first_nonempty(donor_samples["sex"].astype(str).tolist()),
                "sample_match_method": "donor_id_to_human_dg_soft_records" if len(donor_samples) else "unmatched",
            }
        )
    donor_summary = pd.DataFrame(donor_rows).sort_values("donor_id")

    component = pd.DataFrame(
        [
            {
                "dataset": "GSE186538",
                "component_id": "DG_GC_candidates",
                "component_type": "curated_subset",
                "source_member": "GSE186538_Human_counts.mtx.gz",
                "gsm": ";".join(sorted(set(";".join(donor_summary["geo_gsm"]).split(";")) - {""})),
                "geo_sample_title": "DG GC candidate subset aggregated across human donors",
                "specimen_id": ";".join(donor_summary["donor_id"].astype(str)),
                "geo_tissue": "Dentate gyrus",
                "geo_development_stage": "Adult/Aging",
                "geo_age": ";".join(donor_summary["geo_age"].astype(str)),
                "geo_age_years": ";".join(donor_summary["geo_age_years"].astype(str)),
                "geo_sex": ";".join(sorted(set(donor_summary["geo_sex"].astype(str)) - {""})),
                "sample_match_method": "donor_id_to_human_dg_soft_records",
                "supplementary_file_match": "matrix_and_metadata_downloaded",
                "n_cells": int(donor_summary["n_dg_gc_cells"].sum()),
                "nnz": "",
                "total_counts": "",
            }
        ]
    )
    marker_counts = call_counts_by_component()
    if not marker_counts.empty:
        component = component.merge(marker_counts, on=["dataset", "component_id"], how="left")
    return donor_summary, component


def write_enriched_cells(component_metadata: pd.DataFrame, donor_summary: pd.DataFrame) -> pd.DataFrame:
    qc = pd.read_csv(QC_CELLS, sep="\t", low_memory=False)
    component_cols = [
        "dataset",
        "component_id",
        "gsm",
        "geo_sample_title",
        "specimen_id",
        "geo_tissue",
        "geo_development_stage",
        "geo_age",
        "geo_age_years",
        "geo_sex",
        "sample_match_method",
    ]
    enriched = qc.merge(component_metadata[component_cols], on=["dataset", "component_id"], how="left")

    donor_cols = ["donor_id", "geo_age", "geo_age_years", "geo_sex", "geo_gsm", "geo_sample_titles", "sample_match_method"]
    donor = donor_summary[donor_cols].rename(
        columns={
            "donor_id": "sample_hint",
            "geo_age": "donor_geo_age",
            "geo_age_years": "donor_geo_age_years",
            "geo_sex": "donor_geo_sex",
            "geo_gsm": "donor_geo_gsm",
            "geo_sample_titles": "donor_geo_sample_titles",
            "sample_match_method": "donor_sample_match_method",
        }
    )
    enriched = enriched.merge(donor, on="sample_hint", how="left")

    is_gse186538 = enriched["dataset"] == "GSE186538"
    enriched.loc[is_gse186538, "geo_age"] = enriched.loc[is_gse186538, "donor_geo_age"]
    enriched.loc[is_gse186538, "geo_age_years"] = enriched.loc[is_gse186538, "donor_geo_age_years"]
    enriched.loc[is_gse186538, "geo_sex"] = enriched.loc[is_gse186538, "donor_geo_sex"]
    enriched.loc[is_gse186538, "gsm"] = enriched.loc[is_gse186538, "donor_geo_gsm"]
    enriched.loc[is_gse186538, "geo_sample_title"] = enriched.loc[is_gse186538, "donor_geo_sample_titles"]
    enriched.loc[is_gse186538, "sample_match_method"] = enriched.loc[is_gse186538, "donor_sample_match_method"]

    drop_cols = [col for col in enriched.columns if col.startswith("donor_")]
    enriched = enriched.drop(columns=drop_cols)
    enriched.to_csv(OUT_ENRICHED_CELLS, sep="\t", index=False, compression="gzip")
    return enriched


def write_summary(sample_df: pd.DataFrame, component_df: pd.DataFrame, donor_df: pd.DataFrame, enriched: pd.DataFrame) -> None:
    matched_components = int((component_df["sample_match_method"] != "unmatched").sum())
    gse186538_donor_matches = int((donor_df["sample_match_method"] != "unmatched").sum())
    lines = [
        "# Human Core GEO Metadata Curation",
        "",
        "Date built: 2026-06-21",
        "",
        "## Inputs",
        "",
        "- `GSE185277_family.soft.gz`, `GSE185553_family.soft.gz`, and `GSE186538_family.soft.gz` were downloaded from NCBI GEO FTP and parsed locally.",
        "- Built sparse-object metadata were taken from the human seed summaries, `GSE186538` DG GC cell metadata, and the harmonized QC table.",
        "",
        "## Curation Result",
        "",
        f"- Parsed {len(sample_df)} GEO sample records across the three human series.",
        f"- Curated {len(component_df)} built human core components; {matched_components} have a non-empty GEO/sample match method.",
        f"- Matched {gse186538_donor_matches} `GSE186538` human DG donors to donor-level SOFT records.",
        f"- Wrote enriched per-cell metadata for {len(enriched)} cells/nuclei.",
        "",
        "## Notes",
        "",
        "- `GSE185277` sample records map cleanly to individual specimens and ages.",
        "- `GSE185553` includes pooled GSM records with multiple specimen IDs; component-level specimen numbers were parsed from supplementary file names where possible.",
        "- `GSE186538` is represented locally as one DG GC subset, so GEO age/sex annotations are donor-level rather than component-level.",
        "",
        "## Outputs",
        "",
        f"- GEO sample metadata: `{OUT_SAMPLE_METADATA.relative_to(ROOT)}`",
        f"- Built component metadata: `{OUT_COMPONENT_METADATA.relative_to(ROOT)}`",
        f"- `GSE186538` donor summary: `{OUT_DONOR_SUMMARY.relative_to(ROOT)}`",
        f"- Enriched cell metadata: `{OUT_ENRICHED_CELLS.relative_to(ROOT)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    sample_df, _ = build_sample_metadata()
    marker_counts = call_counts_by_component()
    seed_component_df = curate_seed_components(sample_df, marker_counts)
    donor_summary, gse186538_component_df = curate_gse186538_donors(sample_df)
    component_df = pd.concat([seed_component_df, gse186538_component_df], ignore_index=True, sort=False)

    sample_df.to_csv(OUT_SAMPLE_METADATA, sep="\t", index=False)
    component_df.to_csv(OUT_COMPONENT_METADATA, sep="\t", index=False)
    donor_summary.to_csv(OUT_DONOR_SUMMARY, sep="\t", index=False)
    enriched = write_enriched_cells(component_df, donor_summary)
    write_summary(sample_df, component_df, donor_summary, enriched)

    print(f"Wrote {OUT_SAMPLE_METADATA}")
    print(f"Wrote {OUT_COMPONENT_METADATA}")
    print(f"Wrote {OUT_DONOR_SUMMARY}")
    print(f"Wrote {OUT_ENRICHED_CELLS}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
