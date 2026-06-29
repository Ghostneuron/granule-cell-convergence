# Human seed archive inspection

Date inspected: 2026-06-21

## Summary

- `GSE185277`: 7 expression tables and 0 barcode files inspected.
- `GSE185553`: 27 expression tables and 25 barcode files inspected.

## Format inference

The expression files are wide count matrices, with genes/features as rows and barcodes/cells as columns. The archive mixes tab-delimited matrices and whitespace-separated quoted matrices, so the parser must infer the delimiter per file. The barcode files are single-column barcode lists. This means the first build should stream each gzipped text file into a sparse matrix rather than fully materializing dense data frames.

## Build implications

- Construct one object per expression table/library first, then merge within dataset after metadata labels are stable.
- Keep `GSE185277` and `GSE185553` separate at first, because they differ in sample grouping and some `GSE185553` libraries have explicit barcode sidecar files.
- Use the archive inventory TSV as the source of truth for member names, inferred dimensions, and parsing mode.

Detailed table: `Project/results/human_seed_archive_inventory.tsv`
