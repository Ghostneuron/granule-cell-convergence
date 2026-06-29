# Secondary dataset promotion and human dentate options

## Answer

The secondary label was not one thing. For `GSE292261` and `GSE214309`, it mostly reflected unfinished annotation/gene-identifier curation, so these two can be promoted to primary validation datasets. For `GSE214905` and `GSE242688`, the secondary label is intrinsic to the experimental design: `GSE214905` is small targeted patch-seq, and `GSE242688` is spatial/proteomics-linked spot-level data rather than cell-level single-cell/nucleus data.

Promoted local datasets: GSE292261, GSE214309.

## Local promotion decisions

- `GSE292261`: promote to primary dentate developmental validation. The count matrix already uses gene symbols, and the metadata now records postnatal stage, Leiden/louvain group, QC metrics, and matrix presence.
- `GSE214309`: promote to primary dentate maturation/activity validation. The sample metadata are strong and now split maturation state, activity state, timepoint, sex, and mouse ID. The remaining technical requirement is a full Ensembl-to-symbol map for whole-transcriptome object-level analysis.
- `GSE214905`: keep as supporting targeted validation. Gene symbols and metadata are usable, but the sample size and patch-seq/projection design make it a validation layer rather than a discovery backbone.
- `GSE242688`: keep as supporting context. Spatial spots should not be combined with cell-level observations in primary cell-level statistics.

## Human dentate/hippocampal core construction candidates

- `GSE185277`: first sparse human hippocampal/dentate imGC scaffold is built, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
- `GSE185553`: broader human hippocampal companion sparse reference is built, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
- `GSE186538`: human DG GC sparse subset is built from 32,067 candidate cells, marker-validated, donor metadata curated, and included as the tuned DG anchor in the normalized reduced object.
- `GSE325391`: primary modern adult human dentate acquisition is now downloaded, inspected, converted to a selected sparse bridge, and projected into the tuned human-core label convention.
- `GSE268609`: RNA matrix branch is now downloaded, selected-gene bridged, and projected into human-core labels as the broader human hippocampal aging/AD multiome expansion; defer ATAC/full Seurat object initially.
- `GSE198323`: keep as disease context after the healthy/reference human imGC branch is built.

`GSE216877` and `GSE317381` should be kept as disease/spatial context rather than the first human primary datasets.

## Practical next step

Build primary object-level analyses in two tiers:

1. Human dentate/hippocampal construction: use `GSE325391` as the primary adult human DG anchor alongside the tuned `GSE185277`, `GSE185553`, and `GSE186538` scaffold, with `GSE268609` as the broader aging/AD hippocampal RNA expansion.
2. Main comparative rerun: combine the new human dentate/hippocampal branch with the local mouse dentate and cerebellar backbone (`GSE104323`, `GSE95752`, `GSE122357`, `GSE165657`, `GSE312658`, plus promoted `GSE292261` and `GSE214309`).
