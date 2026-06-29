#!/usr/bin/env python3
"""Curate GEO sample metadata for the GSE268609 human hippocampal multiome series."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE268609"
RESULTS = ROOT / "Project/results"

SOFT = BASE / "GSE268609_family.soft.gz"

OUT_SAMPLE_METADATA = RESULTS / "gse268609_geo_sample_metadata.tsv"
OUT_RNA_SUMMARY = RESULTS / "gse268609_geo_rna_sample_summary.tsv"
OUT_DOWNLOAD_STATUS = RESULTS / "gse268609_download_status.tsv"
OUT_MD = RESULTS / "gse268609_geo_metadata_summary.md"

EXPECTED_BYTES = {
    "GSE268609_barcodes.tsv.gz": 1449649,
    "GSE268609_features.tsv.gz": 2946511,
    "GSE268609_matrix.mtx.gz": 4573344069,
    "GSE268609_family.soft.gz": 14120,
}


def norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def basename(value: str) -> str:
    return Path(urlparse(value).path).name if "://" in value else Path(value).name


def parse_float(value: object) -> float | None:
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"[-+]?[0-9]+(?:\.[0-9]+)?", text)
    return float(match.group(0)) if match else None


def library_name_from_description(value: object) -> str:
    for part in str(value).split("; "):
        match = re.search(r"Library name:\s*(.+?)\s*$", part)
        if match:
            return re.sub(r"\s+", " ", match.group(1)).strip()
    return ""


def sample_id_from_title(title: object) -> str:
    text = str(title).strip()
    match = re.match(r"^(\d+)(?:_|\s+)", text)
    return match.group(1) if match else ""


def library_type_from_title(title: object) -> str:
    text = str(title).strip().lower()
    if "rna" in text or "gex" in text:
        return "RNA"
    if "atac" in text:
        return "ATAC"
    return ""


def parse_soft() -> tuple[dict[str, str], pd.DataFrame]:
    series: dict[str, str] = {}
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    repeated: dict[str, list[str]] = defaultdict(list)

    def flush() -> None:
        nonlocal current, repeated
        if current is None:
            return
        for key, values in repeated.items():
            current[key] = "; ".join(values)
        samples.append(current)
        current = None
        repeated = defaultdict(list)

    with gzip.open(SOFT, "rt", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                flush()
                current = {"dataset": "GSE268609", "sample_accession": line.split(" = ", 1)[1]}
                continue
            if line.startswith("^") and not line.startswith("^SAMPLE = "):
                flush()
                continue
            if line.startswith("!Series_") and " = " in line:
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
            elif sample_key in {"description", "extract_protocol", "data_processing", "relation"}:
                repeated[sample_key].append(value)
            elif sample_key.startswith("supplementary_file"):
                repeated["supplementary_files"].append(value)
            else:
                current[sample_key] = value
    flush()

    df = pd.DataFrame(samples)
    if df.empty:
        return series, df

    df["library_name"] = df.get("description", "").map(library_name_from_description)
    df["sample_id_from_title"] = df.get("title", "").map(sample_id_from_title)
    df["sample_id"] = df.get("sample_id", "").astype(str).str.strip()
    missing_sample_id = df["sample_id"].eq("") | df["sample_id"].eq("nan")
    df.loc[missing_sample_id, "sample_id"] = df.loc[missing_sample_id, "sample_id_from_title"]
    df["library_type"] = df.get("library_type", "").astype(str).str.upper().replace({"NAN": ""})
    missing_library_type = df["library_type"].eq("") | df["library_type"].eq("NAN")
    df.loc[missing_library_type, "library_type"] = df.loc[missing_library_type, "title"].map(library_type_from_title)
    df["sample_id_int"] = pd.to_numeric(df["sample_id"], errors="coerce").astype("Int64")
    df["age_at_death_years"] = df.get("age_at_death", "").map(parse_float)
    df["pmi_hours"] = df.get("pmi", "").map(parse_float)
    df["diagnosis"] = df.get("diagnosis", "").astype(str).str.strip()
    df["tissue"] = df.get("tissue", df.get("source_name", "")).astype(str)
    df["supplementary_file_basenames"] = df.get("supplementary_files", "").map(
        lambda x: ";".join(basename(part.strip()) for part in str(x).split(";") if part.strip())
    )

    preferred = [
        "dataset",
        "sample_accession",
        "title",
        "sample_id",
        "sample_id_int",
        "sample_id_from_title",
        "library_type",
        "library_name",
        "diagnosis",
        "age_at_death",
        "age_at_death_years",
        "pmi",
        "pmi_hours",
        "tissue",
        "source_name",
        "organism",
        "instrument_model",
        "library_strategy",
        "library_source",
        "library_selection",
        "molecule",
        "relation",
        "supplementary_file_basenames",
    ]
    for col in preferred:
        if col not in df:
            df[col] = ""
    extras = [col for col in sorted(df.columns) if col not in preferred]
    return series, df[preferred + extras].sort_values(["sample_id_int", "library_type", "sample_accession"])


def summarize_rna(samples: pd.DataFrame) -> pd.DataFrame:
    if samples.empty:
        return pd.DataFrame()
    rna = samples.loc[samples["library_type"].astype(str).str.upper().isin(["RNA", "GEX"])].copy()
    if rna.empty:
        return pd.DataFrame()
    out = (
        rna.groupby("diagnosis", dropna=False)
        .agg(
            n_rna_samples=("sample_accession", "nunique"),
            sample_ids=("sample_id_int", lambda x: ";".join(str(int(v)) for v in sorted(x.dropna().unique()))),
            age_min=("age_at_death_years", "min"),
            age_median=("age_at_death_years", "median"),
            age_max=("age_at_death_years", "max"),
            pmi_median=("pmi_hours", "median"),
            tissues=("tissue", lambda x: ";".join(sorted(set(str(v) for v in x if str(v))))),
        )
        .reset_index()
        .sort_values("diagnosis")
    )
    return out


def write_download_status() -> None:
    rows = []
    for name, expected in EXPECTED_BYTES.items():
        path = BASE / name
        actual = path.stat().st_size if path.exists() else 0
        rows.append(
            {
                "dataset": "GSE268609",
                "file": str(path.relative_to(ROOT)),
                "expected_bytes": expected,
                "actual_bytes": actual,
                "download_complete_by_size": actual == expected,
                "status": "downloaded" if actual == expected else "partial_or_missing",
            }
        )
    pd.DataFrame(rows).to_csv(OUT_DOWNLOAD_STATUS, sep="\t", index=False)


def write_md(series: dict[str, str], samples: pd.DataFrame, rna_summary: pd.DataFrame) -> None:
    status = pd.read_csv(OUT_DOWNLOAD_STATUS, sep="\t")
    n_rna = int((samples["library_type"].astype(str).str.upper() == "RNA").sum()) if len(samples) else 0
    n_atac = int((samples["library_type"].astype(str).str.upper() == "ATAC").sum()) if len(samples) else 0
    complete = int(status["download_complete_by_size"].sum())
    lines = [
        "# GSE268609 GEO Metadata Summary",
        "",
        "Date curated: 2026-06-21",
        "",
        "## Series",
        "",
        f"- Title: {series.get('title', '')}",
        f"- GEO status: {series.get('status', '')}",
        f"- Samples parsed: {len(samples)} ({n_rna} RNA, {n_atac} ATAC).",
        f"- Downloaded expected small/matrix files: {complete} / {len(status)} complete by byte size.",
        "",
        "## RNA Diagnosis Summary",
        "",
    ]
    if rna_summary.empty:
        lines.append("No RNA samples parsed.")
    else:
        for _, row in rna_summary.iterrows():
            lines.append(
                f"- `{row['diagnosis']}`: {int(row['n_rna_samples'])} RNA samples, "
                f"sample IDs {row['sample_ids']}, age median {row['age_median']:.1f} years "
                f"(range {row['age_min']:.1f}-{row['age_max']:.1f}), median PMI {row['pmi_median']:.1f} h."
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This is a primary human dentate/hippocampal candidate because the GEO design explicitly isolates dentate gyrus/hippocampal nuclei and includes multiome RNA with barcode suffixes matching sample IDs.",
            "- The bundled sparse matrix is mixed gene-expression plus ATAC-peak features; RNA extraction must restrict to `Gene Expression` rows before cross-dataset projection.",
            "",
            "## Outputs",
            "",
            f"- Sample metadata: `{OUT_SAMPLE_METADATA.relative_to(ROOT)}`",
            f"- RNA summary: `{OUT_RNA_SUMMARY.relative_to(ROOT)}`",
            f"- Download status: `{OUT_DOWNLOAD_STATUS.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    series, samples = parse_soft()
    rna_summary = summarize_rna(samples)
    samples.to_csv(OUT_SAMPLE_METADATA, sep="\t", index=False)
    rna_summary.to_csv(OUT_RNA_SUMMARY, sep="\t", index=False)
    write_download_status()
    write_md(series, samples, rna_summary)

    print(f"Wrote {OUT_SAMPLE_METADATA}")
    print(f"Wrote {OUT_RNA_SUMMARY}")
    print(f"Wrote {OUT_DOWNLOAD_STATUS}")
    print(f"Wrote {OUT_MD}")
    print(f"samples={len(samples)}; rna={(samples['library_type'].astype(str).str.upper() == 'RNA').sum()}; atac={(samples['library_type'].astype(str).str.upper() == 'ATAC').sum()}")


if __name__ == "__main__":
    main()
