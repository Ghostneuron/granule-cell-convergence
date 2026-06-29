# Granule-Special Gene Named-Comparator Screen

This screen searches for genes that are higher in dentate granule-lineage groups than pyramidal comparators in `GSE104323` and higher in cerebellar granule-lineage groups than Purkinje cells in `GSE122357`.

A gene is called shared-positive when dentate delta is positive, median cerebellar delta is positive, at least two cerebellar samples are positive, and median granule detection is at least 5% in both branches. A strong shared call additionally requires deltas >= 0.05 and all three cerebellar samples positive.

## Screen Counts

- Granule lineage: tested 18,081 genes; shared-positive 10; strong shared 5.
- Postmitotic/granule-cell state: tested 18,081 genes; shared-positive 3; strong shared 2.

## Top Shared-Positive Candidates

### Granule lineage

`NFIA`, `HMGN2`, `RBFOX3`, `GM8730`, `NEUROD1`, `H3F3B`, `RPS24`, `RPL8`, `RPS4X`, `RPL9`

### Postmitotic/granule-cell state

`RBFOX3`, `NFIA`, `HMGN2`

## NFIA Position

- Granule lineage: rank 1 among shared-positive genes; dentate delta 0.490, cerebellar median delta 1.498.
- Postmitotic/granule-cell state: rank 2 among shared-positive genes; dentate delta 0.144, cerebellar median delta 1.386.

## Interpretation

- The screen identifies local named-comparator granule-enriched genes, not universal granule-cell-specific markers across the whole brain.
- The granule-lineage contrast is broader and includes precursor/neuroblast biology.
- The postmitotic/granule-cell-state contrast is more conservative for mature or differentiating granule-cell identity.
- Candidates should be cross-checked against broader cell atlases before being described as granule-cell-specific.

## Outputs

- Unit table: `Project/results/primary_core_granule_special_gene_named_comparator_units.tsv.gz`
- Summary table: `Project/results/primary_core_granule_special_gene_named_comparator_summary.tsv`
- Top candidates: `Project/results/primary_core_granule_special_gene_named_comparator_top_candidates.tsv`
- Plot: `Project/results/primary_core_granule_special_gene_named_comparator_top_candidates.png`
