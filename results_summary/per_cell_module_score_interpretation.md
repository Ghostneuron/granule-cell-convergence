# Per-cell module-score interpretation

## What this adds

This pass moves beyond dataset-level marker presence and creates a per-cell/per-spot module-score layer for the local hippocampal, cerebellar, organoid, and spatial datasets. The score is a lightweight mean log1p marker-panel score computed by streaming only marker genes from each matrix. It is not a replacement for full Seurat/SCE normalization, but it is useful for testing whether the core project hypothesis has the right shape.

Generated outputs:

- `Project/results/per_cell_marker_module_scores.tsv.gz`
- `Project/results/per_cell_marker_module_score_summary.tsv`
- `Project/results/module_score_identity_structural_contrasts.tsv`
- `Project/results/module_score_identity_structural_scatter.png`

## Key result

The per-cell contrast supports the working model: dentate and cerebellar granule-cell datasets separate by region-identity modules, while both carry detectable structural/wiring modules.

Examples:

- `GSE104323` curated dentate granule lineage groups are dentate-high and cerebellar-low. `Immature-GC` has dentate identity 0.473 versus cerebellar identity 0.005, with structural-program score 0.552. `GC-adult` has dentate identity 0.227 versus cerebellar identity 0.0035, with structural-program score 0.226.
- `GSE292261` postnatal DG Smart-seq cells are strongly dentate-high. `DG_P5` has dentate identity 2.614 versus cerebellar identity 0.021, with structural-program score 3.496. `DG_P28` remains dentate-high but has lower structural score, consistent with a developmental signal that needs proper normalization before interpretation.
- `GSE214309` dentate granule nuclei are strongly dentate-high across immature, mature, active, and inactive groups. The 4 hr groups show especially high structural-program scores, suggesting this dataset may be useful for linking activity state to morphogenesis/wiring programs.
- `GSE122357`, `GSE312658`, and `GSE165657` cerebellar datasets are generally cerebellar-identity higher, while still showing structural-program signal. This is the crucial comparison for the paper: different regional identity, shared structural machinery.

## Important caveats

- These are targeted module scores, not fully normalized cross-dataset statistics. Platform effects are still strong.
- `GSE122357`, `GSE165657`, `GSE312658`, and `GSE242688` are currently treated as mixed cerebellum or spatial matrices unless stronger per-cell annotations are added.
- `GSE242688` appears dentate-identity higher despite being cerebellum/spatial. This should be treated as a warning that the current dentate/cerebellar panels need refinement before final claims.
- `GSE150153` organoid metadata rows do not match barcode counts, so the current grouping is deliberately `metadata_unmatched`.

## Project strategy from here

1. Refine marker panels before formal statistics. Remove ambiguous aliases such as `Tau`, check `Atp5f1a`/`Atp5a1`, and split region-identity genes from structural executor genes.
2. Add stronger cerebellar granule-cell annotation. The most important next task is to isolate true cerebellar granule cells from mixed cerebellum matrices using marker-based clustering or external annotations.
3. Build full Seurat/SCE objects for the main datasets and compute dataset-normalized module scores with control genes or rank-based scoring.
4. Add ortholog mapping for mouse, rat, and human before cross-species comparisons.
5. Test the manuscript model statistically: region-identity contrast should distinguish dentate versus cerebellar granule cells, while the structural-program score should converge across both.
6. Use developmental datasets to ask whether convergence is strongest during maturation, when axon extension, dendritic simplification, and synapse formation are active.
7. Use regulatory-network analysis to separate upstream identity regulators from downstream morphology executors.

## Current conclusion

The project remains a strong "yes." The analysis is now pointing toward a sharper and more publishable claim: similar morphology is unlikely to mean shared cell identity; it more likely reflects convergence onto a partially shared morphogenesis and wiring program downstream of distinct hippocampal and cerebellar developmental identities.
