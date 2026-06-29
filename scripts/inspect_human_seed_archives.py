#!/usr/bin/env python3
"""Inspect first-pass human dentate/hippocampal seed archives."""

from __future__ import annotations

import csv
import gzip
import io
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"

ARCHIVES = [
    ("GSE185277", ROOT / "External_Data/GEO/GSE185277/GSE185277_RAW.tar"),
    ("GSE185553", ROOT / "External_Data/GEO/GSE185553/GSE185553_RAW.tar"),
    ("GSE185553", ROOT / "External_Data/GEO/GSE185553/GSE185553.barcodes.tar.gz"),
]

OUT = RESULTS / "human_seed_archive_inventory.tsv"
SUMMARY_OUT = RESULTS / "human_seed_archive_summary.md"


def open_member_text(tar: tarfile.TarFile, member: tarfile.TarInfo) -> io.TextIOWrapper:
    raw = tar.extractfile(member)
    if raw is None:
        raise ValueError(f"Could not read {member.name}")
    if member.name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace")
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")


def classify_member(name: str) -> str:
    lower = name.lower()
    if "barcode" in lower:
        return "barcode_list"
    if lower.endswith(".txt.gz") or lower.endswith(".dge.txt.gz") or lower.endswith(".deg.txt.gz"):
        return "expression_table"
    return "other"


def infer_orientation(first_fields: list[str], second_fields: list[str], third_fields: list[str]) -> str:
    if len(first_fields) > 10 and len(second_fields) > 10:
        return "wide_matrix_gene_by_barcode"
    if len(first_fields) in (2, 3) and len(second_fields) in (2, 3) and len(third_fields) in (2, 3):
        return "long_sparse_or_pair_table"
    if len(first_fields) == 1 and len(second_fields) == 1:
        return "single_column_list"
    return "unknown_tabular"


def inspect_member(dataset: str, archive: Path, tar: tarfile.TarFile, member: tarfile.TarInfo) -> dict[str, str | int]:
    kind = classify_member(member.name)
    first_lines: list[str] = []
    line_count = 0
    max_tab_fields = 0
    max_whitespace_fields = 0
    first_tab_fields: list[str] = []
    second_tab_fields: list[str] = []
    third_tab_fields: list[str] = []
    first_ws_fields: list[str] = []
    second_ws_fields: list[str] = []
    third_ws_fields: list[str] = []

    with open_member_text(tar, member) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            line_count += 1
            if len(first_lines) < 3:
                first_lines.append(line[:500])
            tab_fields = line.split("\t")
            ws_fields = line.split()
            max_tab_fields = max(max_tab_fields, len(tab_fields))
            max_whitespace_fields = max(max_whitespace_fields, len(ws_fields))
            if line_count == 1:
                first_tab_fields = tab_fields
                first_ws_fields = ws_fields
            elif line_count == 2:
                second_tab_fields = tab_fields
                second_ws_fields = ws_fields
            elif line_count == 3:
                third_tab_fields = tab_fields
                third_ws_fields = ws_fields

    parse_fields = first_tab_fields
    parse_second_fields = second_tab_fields
    parse_third_fields = third_tab_fields
    delimiter_hint = "tab"
    if max_tab_fields == 1 and max_whitespace_fields > 1:
        parse_fields = first_ws_fields
        parse_second_fields = second_ws_fields
        parse_third_fields = third_ws_fields
        delimiter_hint = "whitespace"
    elif max_tab_fields == 1:
        delimiter_hint = "single_column"

    orientation = infer_orientation(parse_fields, parse_second_fields, parse_third_fields)
    if kind == "barcode_list":
        inferred_observations = line_count
        inferred_features = ""
    elif orientation == "wide_matrix_gene_by_barcode":
        inferred_features = max(line_count - 1, 0)
        first_token = parse_fields[0].strip('"') if parse_fields else ""
        if first_token.upper() == "GENE":
            inferred_observations = max(len(parse_fields) - 1, 0)
        else:
            inferred_observations = len(parse_fields)
    else:
        inferred_features = ""
        inferred_observations = ""

    return {
        "dataset": dataset,
        "archive": str(archive.relative_to(ROOT)),
        "member": member.name,
        "member_size_bytes": member.size,
        "kind": kind,
        "line_count": line_count,
        "max_tab_fields": max_tab_fields,
        "max_whitespace_fields": max_whitespace_fields,
        "delimiter_hint": delimiter_hint,
        "inferred_orientation": orientation,
        "inferred_features": inferred_features,
        "inferred_observations_or_barcodes": inferred_observations,
        "preview": " | ".join(first_lines),
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int]] = []
    for dataset, archive in ARCHIVES:
        mode = "r:gz" if archive.name.endswith(".tar.gz") else "r:"
        with tarfile.open(archive, mode) as tar:
            for member in tar.getmembers():
                if member.isfile():
                    rows.append(inspect_member(dataset, archive, tar, member))

    fieldnames = [
        "dataset",
        "archive",
        "member",
        "member_size_bytes",
        "kind",
        "line_count",
        "max_tab_fields",
        "max_whitespace_fields",
        "delimiter_hint",
        "inferred_orientation",
        "inferred_features",
        "inferred_observations_or_barcodes",
        "preview",
    ]
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    expression_rows = [row for row in rows if row["kind"] == "expression_table"]
    barcode_rows = [row for row in rows if row["kind"] == "barcode_list"]
    dataset_counts = {}
    for row in rows:
        dataset_counts.setdefault(row["dataset"], {"expression": 0, "barcode": 0})
        if row["kind"] == "expression_table":
            dataset_counts[row["dataset"]]["expression"] += 1
        if row["kind"] == "barcode_list":
            dataset_counts[row["dataset"]]["barcode"] += 1

    lines = [
        "# Human seed archive inspection",
        "",
        "Date inspected: 2026-06-21",
        "",
        "## Summary",
        "",
    ]
    for dataset, counts in sorted(dataset_counts.items()):
        lines.append(
            f"- `{dataset}`: {counts['expression']} expression tables and {counts['barcode']} barcode files inspected."
        )
    lines.extend(
        [
            "",
            "## Format inference",
            "",
            "The expression files are wide count matrices, with genes/features as rows and barcodes/cells as columns. The archive mixes tab-delimited matrices and whitespace-separated quoted matrices, so the parser must infer the delimiter per file. The barcode files are single-column barcode lists. This means the first build should stream each gzipped text file into a sparse matrix rather than fully materializing dense data frames.",
            "",
            "## Build implications",
            "",
            "- Construct one object per expression table/library first, then merge within dataset after metadata labels are stable.",
            "- Keep `GSE185277` and `GSE185553` separate at first, because they differ in sample grouping and some `GSE185553` libraries have explicit barcode sidecar files.",
            "- Use the archive inventory TSV as the source of truth for member names, inferred dimensions, and parsing mode.",
            "",
            f"Detailed table: `{OUT.relative_to(ROOT)}`",
            "",
        ]
    )
    SUMMARY_OUT.write_text("\n".join(lines))
    print(f"Wrote {OUT}")
    print(f"Wrote {SUMMARY_OUT}")
    print(f"Expression tables: {len(expression_rows)}; barcode files: {len(barcode_rows)}")


if __name__ == "__main__":
    main()
