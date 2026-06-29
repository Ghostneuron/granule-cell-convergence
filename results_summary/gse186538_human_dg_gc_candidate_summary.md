# GSE186538 human DG granule-cell candidates

Date prepared: 2026-06-21

## Selection Rule

Cells were selected from the human metadata using `region == DG` and `cluster` names beginning with `DG GC`.

## Counts

- Candidate DG granule cells: 32067
- Source human metadata cells: 219058

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

## Outputs

- Candidate cell table: `Project/results/gse186538_human_dg_gc_candidate_cells.tsv`
- Candidate summary table: `Project/results/gse186538_human_dg_gc_candidate_summary.tsv`

The count matrix is already local as Matrix Market format and is oriented genes-by-cells. The next compute-heavy step is to stream-extract these 32,067 candidate columns from the full matrix.
