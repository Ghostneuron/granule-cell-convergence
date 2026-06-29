# GSE186538 human DG GC sparse subset

Date built: 2026-06-21

## Summary

- Cells: 32067
- Genes: 33939
- Non-zero count entries: 113642650
- Total counts: 331331616

## Cluster Counts

- `DG GC PROX1 SGCZ`: 27821
- `DG GC PROX1 PDLIM5`: 4246

## Sample Counts

- `HSB179`: 9568
- `HSB181`: 7319
- `HSB628`: 4877
- `HSB237`: 4135
- `HSB231`: 3718
- `HSB282`: 2450

## Output Format

The object is saved as a compressed SciPy sparse matrix in cells-by-genes orientation, with matching cell, barcode, gene, and gene-metadata tables.

- Matrix: `Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/matrix_cells_by_genes.npz`
- Cell metadata: `Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/cell_metadata.tsv.gz`
- Gene metadata: `Project/processed/human_dg_taxonomy_sparse_objects/GSE186538/DG_GC_candidates/gene_metadata.tsv.gz`
- Summary TSV: `Project/results/gse186538_human_dg_gc_sparse_subset_summary.tsv`
