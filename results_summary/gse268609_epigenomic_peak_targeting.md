# GSE268609 Epigenomic Peak Targeting

Date built: 2026-06-26

## Purpose

This no-heavy-download step maps curated epigenomic target genes to nearby `GSE268609` ATAC peak feature rows. It prepares a selective extraction manifest for future matrix processing.

## Scope

- Target genes with gene-expression feature rows: 153.
- Peak-target rows within gene body plus 100 kb: 2623.
- Unique nearby ATAC peak features: 2591.
- Target gene-expression rows: 153.
- Unique feature rows in selective manifest: 2744 (1.54% of all `GSE268609` feature rows).

## Highest Peak-Count Genes

- `ZBTB20`: 78 nearby peaks (65 gene-body, 0 edge-2kb, 1 proximal, 12 distal).
- `NFIA`: 61 nearby peaks (54 gene-body, 0 edge-2kb, 3 proximal, 4 distal).
- `NRXN1`: 56 nearby peaks (51 gene-body, 1 edge-2kb, 1 proximal, 3 distal).
- `TCF4`: 49 nearby peaks (39 gene-body, 0 edge-2kb, 0 proximal, 10 distal).
- `ROBO2`: 49 nearby peaks (47 gene-body, 0 edge-2kb, 0 proximal, 2 distal).
- `NFIB`: 39 nearby peaks (28 gene-body, 1 edge-2kb, 1 proximal, 9 distal).
- `DPYSL2`: 36 nearby peaks (6 gene-body, 2 edge-2kb, 2 proximal, 26 distal).
- `DCC`: 35 nearby peaks (34 gene-body, 0 edge-2kb, 0 proximal, 1 distal).

## Interpretation

The manifest makes later peak-count processing substantially narrower than the full 177,931-row GSE268609 matrix. It still does not fit chromatin accessibility because the large count matrix is not currently local.

## Outputs

- Peak targets: `Project/results/gse268609_epigenomic_peak_targets.tsv`
- Gene summary: `Project/results/gse268609_epigenomic_peak_gene_summary.tsv`
- Selective extraction manifest: `Project/results/gse268609_epigenomic_selective_extraction_manifest.tsv`
