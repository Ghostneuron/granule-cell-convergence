# GSE268609 GEO Metadata Summary

Date curated: 2026-06-21

## Series

- Title: A roadmap to human hippocampal neurogenesis in adulthood, aging and AD
- GEO status: Public on Dec 15 2025
- Samples parsed: 78 (39 RNA, 39 ATAC).
- Downloaded expected small/matrix files: 4 / 4 complete by byte size.

## RNA Diagnosis Summary

- `AD`: 8 RNA samples, sample IDs 1;2;3;4;26;28;30;34, age median 81.0 years (range 70.0-92.0), median PMI 5.7 h.
- `HA`: 9 RNA samples, sample IDs 5;6;7;8;23;27;29;31;36, age median 87.0 years (range 60.0-93.0), median PMI 7.6 h.
- `MCI`: 8 RNA samples, sample IDs 9;10;11;12;21;22;33;39, age median 87.5 years (range 80.0-94.0), median PMI 6.3 h.
- `SA`: 6 RNA samples, sample IDs 15;16;17;18;19;20, age median 91.0 years (range 86.0-100.0), median PMI 7.0 h.
- `YA`: 8 RNA samples, sample IDs 13;14;24;25;32;35;37;38, age median 32.0 years (range 21.0-38.0), median PMI 5.8 h.

## Interpretation

- This is a primary human dentate/hippocampal candidate because the GEO design explicitly isolates dentate gyrus/hippocampal nuclei and includes multiome RNA with barcode suffixes matching sample IDs.
- The bundled sparse matrix is mixed gene-expression plus ATAC-peak features; RNA extraction must restrict to `Gene Expression` rows before cross-dataset projection.

## Outputs

- Sample metadata: `Project/results/gse268609_geo_sample_metadata.tsv`
- RNA summary: `Project/results/gse268609_geo_rna_sample_summary.tsv`
- Download status: `Project/results/gse268609_download_status.tsv`
