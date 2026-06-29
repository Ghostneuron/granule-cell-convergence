# User-suggested resource assessment

Date checked: 2026-06-21

## Summary

All three resources are useful and have now been integrated into the candidate-resource inventory. The Cell Reports 2025 human dentate gyrus paper is already present in the local `Literature` folder, but its spatial data are not incorporated into the current primary expression-analysis table. The Nature 2022 immature-neuron paper now contributes the first built human hippocampal/dentate seed branch through `GSE185277` and `GSE185553`.

## Resource-by-resource decision

| Resource | Title | Data/accession found | Current project status | Recommended role |
|---|---|---:|---|---|
| Nature Neuroscience 2025, DOI `10.1038/s41593-024-01848-4` | Multimodal transcriptomics reveal neurogenic aging trajectories and age-related regional inflammation in the dentate gyrus | `GSE233363` | Integrated into candidate resources; not downloaded or analyzed | Add as mouse dentate aging/neurogenic-lineage validation. Strong for maturation, neurogenic aging, spatial inflammation, and separating granule/neurogenic programs from niche effects. |
| Cell Reports 2025, DOI `10.1016/j.celrep.2025.115300` | Spatiotemporal analysis of gene expression in the human dentate gyrus reveals age-associated changes in cellular maturation and neuroinflammation | Zenodo `10.5281/zenodo.10126687`; raw data `10.5281/zenodo.10126688`; code `10.5281/zenodo.10126715`; spatialLIBD Shiny portal | PDF is present locally in `Literature`; integrated into candidate resources; data not downloaded or analyzed | Add as human DG spatial/lifespan context. Important for human dentate maturation and neuroinflammation, but spatial/tissue-level rather than single-cell primary. |
| Nature 2022, DOI `10.1038/s41586-022-04912-w`, PMCID `PMC9316413` | Molecular landscapes of human hippocampal immature neurons across lifespan | `GSE185277`, `GSE185553`, `GSE198323` | `GSE185277` and `GSE185553` raw archives downloaded and sparse objects built; `GSE198323` remains disease context, not downloaded | Use the built seed objects for human imGC marker scaffolding and cross-species caution alongside the built `GSE186538` human DG taxonomy subset. |

## Practical recommendation

1. Add `GSE233363` to the next acquisition list as a mouse dentate aging validation dataset.
2. Use the local Cell Reports 2025 human DG spatial PDF as a literature source, and acquire the Zenodo SpatialExperiment objects if we want a human DG spatial validation layer.
3. Use the built `GSE185277`, `GSE185553`, and `GSE186538` sparse objects as the first human dentate/hippocampal reference branch; keep `GSE198323` as later AD context.
