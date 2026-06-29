#!/usr/bin/env python3
"""Curate GEO sample metadata for the GSE325391 adult dentate object."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "External_Data/GEO/GSE325391"
RESULTS = ROOT / "Project/results"

SOFT = BASE / "GSE325391_family.soft.gz"
RDS = BASE / "GSE325391_adultgc_filtered.RDS"

OUT_SAMPLE_METADATA = RESULTS / "gse325391_geo_sample_metadata.tsv"
OUT_GROUP_SUMMARY = RESULTS / "gse325391_geo_group_summary.tsv"
OUT_DOWNLOAD_STATUS = RESULTS / "gse325391_download_status.tsv"
OUT_MD = RESULTS / "gse325391_geo_metadata_summary.md"


def norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def basename(value: str) -> str:
    return Path(urlparse(value).path).name if "://" in value else Path(value).name


def parse_age(value: object) -> float | None:
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)", text)
    return float(match.group(1)) if match else None


def split_title(title: object) -> dict[str, str]:
    text = str(title)
    match = re.match(r"Sample(\d+)\s+of\s+Run(\d+),\s*(.*?),\s*Dentate gyrus punch", text, flags=re.I)
    if not match:
        return {
            "sample_number_from_title": "",
            "run_from_title": "",
            "condition_from_title": "",
        }
    return {
        "sample_number_from_title": match.group(1),
        "run_from_title": match.group(2),
        "condition_from_title": re.sub(r"\s+", " ", match.group(3)).strip(),
    }


def library_name_from_description(value: object) -> str:
    for part in str(value).split("; "):
        match = re.search(r"Library name:\s*(\S+)", part)
        if match:
            return match.group(1)
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
                current = {"dataset": "GSE325391", "sample_accession": line.split(" = ", 1)[1]}
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

    title_parts = df["title"].map(split_title).apply(pd.Series)
    df = pd.concat([df, title_parts], axis=1)
    df["library_name"] = df.get("description", "").map(library_name_from_description)
    df["age_years"] = df.get("age", "").map(parse_age)
    df["geo_group"] = df.get("group", "").astype(str)
    df["batch_or_run"] = df.get("batch", "").astype(str).str.replace("Run ", "", regex=False)
    df["supplementary_file_basenames"] = df.get("supplementary_files", "").map(
        lambda x: ";".join(basename(part.strip()) for part in str(x).split(";") if part.strip())
    )

    preferred = [
        "dataset",
        "sample_accession",
        "title",
        "sample_number_from_title",
        "run_from_title",
        "condition_from_title",
        "library_name",
        "source_name",
        "tissue",
        "sex",
        "age",
        "age_years",
        "geo_group",
        "batch",
        "batch_or_run",
        "organism",
        "instrument_model",
        "library_strategy",
        "library_source",
        "molecule",
        "relation",
        "supplementary_file_basenames",
    ]
    for col in preferred:
        if col not in df:
            df[col] = ""
    extras = [col for col in sorted(df.columns) if col not in preferred]
    return series, df[preferred + extras]


def summarize_groups(samples: pd.DataFrame) -> pd.DataFrame:
    if samples.empty:
        return pd.DataFrame()
    out = (
        samples.groupby(["geo_group", "condition_from_title"], dropna=False)
        .agg(
            n_samples=("sample_accession", "nunique"),
            n_male=("sex", lambda x: int((x.astype(str).str.upper() == "M").sum())),
            n_female=("sex", lambda x: int((x.astype(str).str.upper() == "F").sum())),
            age_min=("age_years", "min"),
            age_median=("age_years", "median"),
            age_max=("age_years", "max"),
            runs=("batch_or_run", lambda x: ";".join(sorted(set(str(v) for v in x if str(v))))),
            libraries=("library_name", lambda x: ";".join(sorted(set(str(v) for v in x if str(v))))),
        )
        .reset_index()
    )
    return out


def write_download_status() -> None:
    expected_bytes = 3515050479
    actual_bytes = RDS.stat().st_size if RDS.exists() else 0
    rows = [
        {
            "dataset": "GSE325391",
            "file": str(RDS.relative_to(ROOT)),
            "expected_bytes": expected_bytes,
            "actual_bytes": actual_bytes,
            "download_complete_by_size": actual_bytes == expected_bytes,
            "status": "downloaded" if actual_bytes == expected_bytes else "partial_or_missing",
        }
    ]
    pd.DataFrame(rows).to_csv(OUT_DOWNLOAD_STATUS, sep="\t", index=False)


def write_md(series: dict[str, str], samples: pd.DataFrame, groups: pd.DataFrame) -> None:
    status = pd.read_csv(OUT_DOWNLOAD_STATUS, sep="\t")
    status_row = status.iloc[0].to_dict()
    lines = [
        "# GSE325391 GEO Metadata Summary",
        "",
        "Date curated: 2026-06-21",
        "",
        "## Series",
        "",
        f"- Title: {series.get('title', '')}",
        f"- GEO status: {series.get('status', '')}",
        f"- Samples parsed: {len(samples)}",
        f"- Adult RDS bytes: {status_row['actual_bytes']} / {status_row['expected_bytes']} ({status_row['status']})",
        "",
        "## Group Summary",
        "",
    ]
    if groups.empty:
        lines.append("No groups parsed.")
    else:
        for _, row in groups.iterrows():
            lines.append(
                f"- `{row['geo_group']}` / {row['condition_from_title']}: "
                f"{int(row['n_samples'])} samples, age median {row['age_median']:.1f} years "
                f"(range {row['age_min']:.1f}-{row['age_max']:.1f}), runs {row['runs']}."
            )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Sample metadata: `{OUT_SAMPLE_METADATA.relative_to(ROOT)}`",
            f"- Group summary: `{OUT_GROUP_SUMMARY.relative_to(ROOT)}`",
            f"- Download status: `{OUT_DOWNLOAD_STATUS.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    series, samples = parse_soft()
    groups = summarize_groups(samples)
    samples.to_csv(OUT_SAMPLE_METADATA, sep="\t", index=False)
    groups.to_csv(OUT_GROUP_SUMMARY, sep="\t", index=False)
    write_download_status()
    write_md(series, samples, groups)
    print(f"Wrote {OUT_SAMPLE_METADATA}")
    print(f"Wrote {OUT_GROUP_SUMMARY}")
    print(f"Wrote {OUT_DOWNLOAD_STATUS}")
    print(f"Wrote {OUT_MD}")
    print(f"samples={len(samples)} groups={len(groups)}")


if __name__ == "__main__":
    main()
