# Primary-Core MGI Ortholog Meta-Model

Date built: 2026-06-22

## Purpose

This analysis adds a conservative ortholog-aware layer to the selected-feature and full-matrix pseudobulk screens. It uses the official MGI human-mouse homology report, keeps one_to_one human-mouse classes, and further restricts the strict model to classes where the human and mouse symbols have the same canonical symbol.

Because the current local matrix extraction was done in a same-symbol frame, non-identical one_to_one orthologs are intentionally deferred until a mouse-symbol-aware extraction is built.

## Ortholog Scope

- MGI report rows: 46,522.
- MGI human-mouse homology classes: 20,181.
- One_to_one human-mouse pairs: 17,611.
- Strict same-symbol one_to_one pairs: 16,245.
- Strict pairs represented in selected-feature expression rows: 1,903.
- Strict pairs represented in full-matrix expression rows: 15,417.

## Meta-Model Summary

- Unit delta rows: 241,690.
- Branch summary rows: 33,829.
- Gene summary rows: 15,345.
- Shared strict both-screen hits: 68.
- Shared supported both-screen hits: 90.
- Shared full-matrix-only hits: 1,140.
- Shared selected-only hits: 6.
- Mechanism-prioritized shared hits: 36.

A branch is supported when at least two datasets contribute, at least 75% of datasets have positive candidate-versus-background deltas, and the median dataset delta is positive. A branch is strict when it also passes the dataset-level sign-test threshold p<=0.25, which is intentionally permissive because some branches have only two or three independent datasets.

## Top Shared Ortholog Hits

- `CELF2` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.98.
- `FTH1` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.97.
- `GPM6A` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.93.
- `PPP3CA` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.89.
- `NFIB` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.89.
- `NFIA` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.82.
- `ANK3` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.76.
- `RBFOX3` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.74.
- `NBEA` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.72.
- `NRN1` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.66.
- `ATP2B1` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.62.
- `TNIK` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.59.
- `PPP2R2C` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.58.
- `CELF4` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.58.
- `ZBTB18` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.57.
- `PPM1H` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.47.
- `CACNA2D1` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.45.
- `ST18` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.44.
- `SSBP3` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.41.
- `BAIAP2` (strict_shared_both_screens): 4/4 supported screen/branches, score 14.41.

## Mechanism-Prioritized Hits

- `GPM6A` (dataset_robust_consensus_figure; cytoskeleton_morphogenesis): 4/4 supported screen/branches.
- `NFIB` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): 4/4 supported screen/branches.
- `NFIA` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): 4/4 supported screen/branches.
- `KCNK1` (dataset_robust_consensus_figure; synaptic_wiring): 4/4 supported screen/branches.
- `RFX3` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): 4/4 supported screen/branches.
- `GABRA2` (dataset_robust_consensus_figure; synaptic_wiring): 4/4 supported screen/branches.
- `PPP3CA` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `CACNA2D1` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `KCNJ6` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `GABRB3` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `GRIN2B` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `KCNJ3` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `KCND2` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `STXBP5L` (consensus_figure_candidate; synaptic_wiring): 4/4 supported screen/branches.
- `ROBO2` (dual_screen_mechanism_figure; curated_shared_structural_executor): 4/4 supported screen/branches.
- `DYNLL1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 4/4 supported screen/branches.
- `BASP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 4/4 supported screen/branches.
- `TUBA1A` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 4/4 supported screen/branches.
- `RFX7` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate): 4/4 supported screen/branches.
- `MAPKAP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 4/4 supported screen/branches.
- `MAP1B` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 3/4 supported screen/branches.
- `ACTB` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 3/4 supported screen/branches.
- `TCF4` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate): 3/4 supported screen/branches.
- `RTN3` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): 3/4 supported screen/branches.

## Consensus Candidate Check

The dataset-robust six-gene consensus shortlist remains inside this strict MGI one_to_one same-symbol frame: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`.

## Interpretation

- This is stronger than the same-symbol screen alone because it removes many-to-many and non-MGI-supported symbol matches.
- It is still not the final ortholog DE model because non-identical one_to_one orthologs are absent from the current same-symbol matrix extraction.
- The strongest manuscript-ready tier should come from genes that are supported in both selected-feature and full-matrix screens, and especially from genes that also passed the 24-candidate dataset-aware validation.
- Rat is not represented in the MGI report used here; a rat extension should be added only when rat primary datasets enter the core.

## Outputs

- Ortholog map: `Project/results/primary_core_mgi_ortholog_meta_model_map.tsv`
- Unit deltas: `Project/results/primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz`
- Branch summary: `Project/results/primary_core_mgi_ortholog_meta_model_branch_summary.tsv`
- Gene summary: `Project/results/primary_core_mgi_ortholog_meta_model_gene_summary.tsv`
- Shared hits: `Project/results/primary_core_mgi_ortholog_meta_model_shared_hits.tsv`
- Mechanism-prioritized hits: `Project/results/primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv`
- Plot: `Project/results/primary_core_mgi_ortholog_meta_model_top_hits.png`
