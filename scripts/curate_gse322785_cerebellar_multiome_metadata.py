#!/usr/bin/env python3
"""Curate lightweight metadata and download plan for GSE322785."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE322785"
RESULTS = ROOT / "Project/results"

SOFT = BASE / "GSE322785_family.soft.gz"
FILELIST = BASE / "filelist.txt"

OUT_SAMPLES = RESULTS / "gse322785_cerebellar_multiome_sample_metadata.tsv"
OUT_FILES = RESULTS / "gse322785_cerebellar_multiome_file_inventory.tsv"
OUT_DONORS = RESULTS / "gse322785_cerebellar_multiome_donor_summary.tsv"
OUT_PLAN = RESULTS / "gse322785_cerebellar_multiome_download_plan.tsv"
OUT_MD = RESULTS / "gse322785_cerebellar_multiome_metadata.md"

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9558nnn"


SPECIES_BY_PREFIX = {
    "C": "Callithrix jacchus",
    "H": "Homo sapiens",
    "P": "Pan troglodytes",
    "R": "Macaca mulatta",
}


def parse_filelist() -> pd.DataFrame:
    rows = []
    with FILELIST.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 5 or parts[0] != "File":
                continue
            _, name, time, size, file_type = parts
            sample_accession = name.split("_", 1)[0]
            donor_match = re.match(r"^GSM\d+_([^_]+)", name)
            donor_id = donor_match.group(1) if donor_match else ""
            if donor_id.endswith("CBm"):
                donor_id = donor_id[:-3]
            library_type = "ATAC" if "atac_fragments" in name else "RNA_multiome_h5"
            rows.append(
                {
                    "sample_accession": sample_accession,
                    "file_name": name,
                    "donor_id": donor_id,
                    "species_prefix": donor_id[:1],
                    "species_inferred_from_prefix": SPECIES_BY_PREFIX.get(donor_id[:1], ""),
                    "library_file_class": library_type,
                    "file_type": file_type,
                    "size_bytes": int(size),
                    "size_mb": int(size) / 1_000_000,
                    "file_time": time,
                    "sample_file_url": f"{FTP_BASE}/{sample_accession}/suppl/{name}",
                    "raw_tar_member": f"GSE322785_RAW/{name}",
                }
            )
    return pd.DataFrame(rows)


def parse_soft() -> pd.DataFrame:
    rows = []
    current: dict[str, object] | None = None
    with gzip.open(SOFT, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                current = {"sample_accession": line.split("=", 1)[1].strip()}
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["sample_title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_geo_accession = "):
                current["sample_accession"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_organism_ch1 = "):
                current["organism"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source_name"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, val = value.split(":", 1)
                    key = key.strip().lower().replace(" ", "_")
                    current[key] = val.strip()
            elif line.startswith("!Sample_relation = "):
                value = line.split("=", 1)[1].strip()
                if value.startswith("BioSample:"):
                    current["biosample"] = value.split(":", 1)[1].strip()
                elif value.startswith("SRA:"):
                    current["sra"] = value.split(":", 1)[1].strip()
            elif line.startswith("!Sample_supplementary_file"):
                current["supplementary_file"] = line.split("=", 1)[1].strip()
        if current:
            rows.append(current)
    samples = pd.DataFrame(rows)
    if not samples.empty:
        samples["donor_id"] = samples["sample_title"].str.replace(r"_(ATAC|RNA)$", "", regex=True)
        samples["donor_id"] = samples["donor_id"].str.replace("CBm", "", regex=False)
        samples["species_prefix"] = samples["donor_id"].str[:1]
        samples["species_inferred_from_prefix"] = samples["species_prefix"].map(SPECIES_BY_PREFIX).fillna("")
    return samples


def build_download_plan(files: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    merged = files.merge(
        samples[
            [
                "sample_accession",
                "sample_title",
                "organism",
                "tissue",
                "developmental_stage",
                "library_type",
                "donor_id",
            ]
        ],
        on=["sample_accession", "donor_id"],
        how="left",
    )
    rows = []
    for rec in merged.to_dict("records"):
        if rec["library_file_class"] == "RNA_multiome_h5" and rec["species_inferred_from_prefix"] == "Homo sapiens":
            tier = "tier1_human_cerebellar_h5"
            reason = "Human cerebellar multiome H5 is the closest counterpart to the human dentate/hippocampal extension."
        elif rec["library_file_class"] == "RNA_multiome_h5":
            tier = "tier2_cross_primate_h5"
            reason = "Cross-primate H5 files can support orthology and conservation checks after human H5 feasibility."
        elif rec["species_inferred_from_prefix"] == "Homo sapiens":
            tier = "tier3_human_atac_fragments"
            reason = "Human ATAC fragments are large and should wait until H5 label/gene-activity feasibility is established."
        else:
            tier = "tier4_cross_primate_atac_fragments"
            reason = "Large fragment files are not first-line for this manuscript extension."
        rows.append(
            {
                **rec,
                "download_tier": tier,
                "download_reason": reason,
                "recommended_now": tier in {"tier1_human_cerebellar_h5"},
            }
        )
    order = {
        "tier1_human_cerebellar_h5": 1,
        "tier2_cross_primate_h5": 2,
        "tier3_human_atac_fragments": 3,
        "tier4_cross_primate_atac_fragments": 4,
    }
    plan = pd.DataFrame(rows)
    plan["download_order"] = plan["download_tier"].map(order)
    return plan.sort_values(["download_order", "size_bytes", "sample_accession"])


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    files = parse_filelist()
    samples = parse_soft()
    donors = (
        samples.groupby(["donor_id", "species_inferred_from_prefix"], dropna=False)
        .agg(
            n_geo_samples=("sample_accession", "nunique"),
            library_types=("library_type", lambda x: ";".join(sorted(set(map(str, x))))),
            sample_titles=("sample_title", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )
    plan = build_download_plan(files, samples)

    files.to_csv(OUT_FILES, sep="\t", index=False)
    samples.to_csv(OUT_SAMPLES, sep="\t", index=False)
    donors.to_csv(OUT_DONORS, sep="\t", index=False)
    plan.to_csv(OUT_PLAN, sep="\t", index=False)

    h5 = plan.loc[plan["library_file_class"].eq("RNA_multiome_h5")]
    frag = plan.loc[plan["library_file_class"].eq("ATAC")]
    recommended = plan.loc[plan["recommended_now"]]
    species_counts = donors["species_inferred_from_prefix"].value_counts().to_dict()
    lines = [
        "# GSE322785 Cerebellar Multiome Metadata",
        "",
        "Date built: 2026-06-26",
        "",
        "## Purpose",
        "",
        "This lightweight curation prepares the adult primate cerebellar multiome dataset as the cerebellar epigenomic counterpart to the current dentate/hippocampal extension.",
        "",
        "## Inventory",
        "",
        f"- Donor/sample prefixes by species: {species_counts}.",
        f"- GEO sample records: {len(samples)}.",
        f"- Filelist records: {len(files)}.",
        f"- H5 feature-barcode files: {len(h5)}, total {h5['size_bytes'].sum() / 1e9:.2f} GB.",
        f"- ATAC fragment files: {len(frag)}, total {frag['size_bytes'].sum() / 1e9:.2f} GB.",
        f"- Raw tar size: 19.61 GB by GEO filelist.",
        "",
        "## Recommended First Download",
        "",
        f"- Human H5 files recommended now: {len(recommended)} files, {recommended['size_bytes'].sum() / 1e6:.1f} MB total.",
        "- Fragment files should be deferred until cell-label/gene-activity feasibility is established from H5 files.",
        "",
        "## Outputs",
        "",
        f"- Sample metadata: `{OUT_SAMPLES.relative_to(ROOT)}`",
        f"- File inventory: `{OUT_FILES.relative_to(ROOT)}`",
        f"- Donor summary: `{OUT_DONORS.relative_to(ROOT)}`",
        f"- Download plan: `{OUT_PLAN.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
