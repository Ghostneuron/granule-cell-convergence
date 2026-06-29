# Human dentate core construction plan

Date revised: 2026-06-21

## Rationale

The earlier local core had a real asymmetry: it included a human cerebellar dataset (`GSE165657`) but no constructed human dentate or hippocampal granule-cell comparator. A staged human dentate/hippocampal core branch has now been constructed through sparse-object building, QC, marker triage, GEO metadata curation, a normalized reduced object, tuned labels, first-pass dataset-aware module tests, adult dentate anchoring, broader hippocampal RNA projection, and rank-level integration with the refined dentate/cerebellar backbone.

## Build Order

1. `GSE185277`: built first because the RAW archive was small enough for quick object construction and gives a human immature granule-cell marker scaffold; now sparse-built, QC-harmonized, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
2. `GSE185553`: built with `GSE185277` as the broader human hippocampal companion reference; now sparse-built, QC-harmonized, marker-scored, pooled specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
3. `GSE186538`: added as the stronger DG taxonomy and annotation-transfer reference; now a 32,067-cell sparse DG GC subset with QC, marker validation, donor-level metadata, normalized-object inclusion, and tuned-label use as the DG anchor.
4. `GSE325391`: adult dentate gyrus RDS added as the primary modern human dentate dataset; 59,075 adult DG nuclei are now projected into the tuned human-core convention.
5. `GSE268609`: RNA matrix/barcodes/features are now built as a broader human hippocampal aging/AD multiome RNA branch, with projected labels and first-pass diagnosis module tests; defer ATAC and the full Seurat RDS unless source taxonomy or regulatory analysis becomes necessary.
6. Integrated rank layer: the human bridge objects and refined backbone are now merged into `Project/results/human_bridge_backbone_rank_*` outputs for manuscript-scale comparison.
7. `GSE233363`: add later as mouse DG aging/neurogenic-lineage validation, not as the immediate fix for the missing human comparison.

## Practical Rule

The main comparative claim should include the constructed human dentate/hippocampal branch. The smallest useful first pass, `GSE185277` plus `GSE185553`, is now built, marker-scored, label-tuned, and module-tested; the stronger human core includes the `GSE186538` DG GC anchor, `GSE325391` as the direct adult human dentate anchor, `GSE268609` as a broader projected hippocampal aging/AD RNA branch, and the integrated rank bridge to the dentate/cerebellar backbone.

The remaining optional `GSE268609` work is source-taxonomy/regulatory: download the full Seurat RDS only if source cell annotations are needed, and keep the 76G ATAC fragments deferred until a regulatory aim is active.
