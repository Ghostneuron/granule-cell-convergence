# Human seed sparse object build

Date built: 2026-06-21

## Summary

- `GSE185277`: 7 libraries, 89526 cells/barcodes, 60286599 non-zero counts, 99592073 total counts.
- `GSE185553`: 27 libraries, 97729 cells/barcodes, 103757887 non-zero counts, 205573021 total counts.

## Output Format

Each library has a compressed SciPy sparse matrix in cells-by-genes orientation, plus gzipped barcodes, genes, cell metadata, and gene metadata files. This is intentionally neutral: it can be imported into AnnData, Seurat, or SingleCellExperiment later.

## Immediate Use

- Use per-library QC first; do not merge libraries until sample labels and donor/age metadata are curated.
- Prefix cell IDs with dataset and library ID to avoid barcode collisions.
- Treat this as the first constructed human dentate/hippocampal branch for marker/QC work, not yet the final harmonized object.

Detailed table: `Project/results/human_seed_sparse_object_summary.tsv`
