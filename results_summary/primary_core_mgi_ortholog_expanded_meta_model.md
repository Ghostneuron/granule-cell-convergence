# Expanded MGI Ortholog Dataset-Aware Meta-Model

Date built: 2026-06-22

## Purpose

This model re-runs the dataset-aware candidate-versus-background meta-analysis using the MGI one-to-one full-matrix expression layer. Unlike the previous conservative model, this version includes non-identical human/mouse ortholog symbols resolved through MGI.

Selected-feature expression remains limited to the selected human-core feature universe. Full-matrix evidence now comes from `primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`.

## Ortholog Scope

- MGI report rows: 46,522.
- MGI human-mouse homology classes: 20,181.
- One-to-one human-mouse pairs: 17,611.
- Same-symbol one-to-one pairs: 16,245.
- Non-identical one-to-one pairs: 1,366.
- One-to-one pairs represented in selected-feature expression rows: 1,957.
- One-to-one pairs represented in MGI full-matrix expression rows: 16,799.

## Meta-Model Summary

- Unit delta rows: 252,469.
- Branch summary rows: 36,303.
- Gene summary rows: 16,708.
- Shared strict both-screen hits: 68.
- Shared supported both-screen hits: 90.
- Shared full-matrix-only hits: 1,206.
- Shared selected-only hits: 6.
- Non-identical-symbol shared hits: 64.
- Non-identical-symbol strict shared hits: 11.
- Mechanism-prioritized shared hits: 36.

A branch is supported when at least two datasets contribute, at least 75% of datasets have positive candidate-versus-background deltas, and the median dataset delta is positive. A branch is strict when it also passes the dataset-level sign-test threshold p<=0.25.

## Dataset-Robust Consensus Genes

- Retained in the expanded MGI model: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`.

## Mechanism-Prioritized Hits

- `GPM6A` (dataset_robust_consensus_figure; cytoskeleton_morphogenesis; same_symbol): 4/4 supported screen/branches.
- `NFIB` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate; same_symbol): 4/4 supported screen/branches.
- `NFIA` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate; same_symbol): 4/4 supported screen/branches.
- `KCNK1` (dataset_robust_consensus_figure; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `RFX3` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate; same_symbol): 4/4 supported screen/branches.
- `GABRA2` (dataset_robust_consensus_figure; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `PPP3CA` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `CACNA2D1` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `KCNJ6` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `GABRB3` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `GRIN2B` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `KCNJ3` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `KCND2` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `STXBP5L` (consensus_figure_candidate; synaptic_wiring; same_symbol): 4/4 supported screen/branches.
- `ROBO2` (dual_screen_mechanism_figure; curated_shared_structural_executor; same_symbol): 4/4 supported screen/branches.
- `DYNLL1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 4/4 supported screen/branches.
- `BASP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 4/4 supported screen/branches.
- `TUBA1A` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 4/4 supported screen/branches.
- `RFX7` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate; same_symbol): 4/4 supported screen/branches.
- `MAPKAP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 4/4 supported screen/branches.
- `MAP1B` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 3/4 supported screen/branches.
- `ACTB` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 3/4 supported screen/branches.
- `TCF4` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate; same_symbol): 3/4 supported screen/branches.
- `RTN3` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 3/4 supported screen/branches.
- `ACTG1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 3/4 supported screen/branches.
- `TUBA1B` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis; same_symbol): 3/4 supported screen/branches.
- `STMN2` (selected_screen_mechanism_figure; curated_shared_structural_executor; same_symbol): 3/4 supported screen/branches.
- `CADM3` (selected_screen_mechanism_figure; axon_guidance_adhesion; same_symbol): 3/4 supported screen/branches.
- `KCNMB4` (selected_screen_mechanism_figure; synaptic_wiring; same_symbol): 3/4 supported screen/branches.
- `DCC` (selected_screen_mechanism_figure; axon_guidance_adhesion; same_symbol): 3/4 supported screen/branches.

## Top Non-Identical Ortholog Hits

- `ZNF148` / mouse `Zfp148` (strict_shared_full_matrix_only): 3/4 supported screen/branches.
- `C1orf21` / mouse `1700025G04Rik` (strict_shared_full_matrix_only): 3/4 supported screen/branches.
- `ZNF706` / mouse `Zfp706` (supported_shared_full_matrix_only): 3/4 supported screen/branches.
- `RAB7A` / mouse `Rab7` (supported_shared_full_matrix_only): 3/4 supported screen/branches.
- `ZNF292` / mouse `Zfp292` (supported_shared_full_matrix_only): 3/4 supported screen/branches.
- `ZNF827` / mouse `Zfp827` (supported_shared_full_matrix_only): 3/4 supported screen/branches.
- `KIAA1328` / mouse `AW554918` (supported_shared_full_matrix_only): 3/4 supported screen/branches.
- `ZNF536` / mouse `Zfp536` (strict_shared_full_matrix_only): 2/4 supported screen/branches.
- `TMEM178A` / mouse `Tmem178` (strict_shared_full_matrix_only): 2/4 supported screen/branches.
- `MIR124-1HG` / mouse `Mir124a-1hg` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF667` / mouse `Zfp667` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF226` / mouse `Zfp61` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `EGFEM1P` / mouse `Egfem1` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `RBM12B` / mouse `Rbm12b2` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF189` / mouse `Zfp189` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `CETN4P` / mouse `Cetn4` (strict_shared_full_matrix_only): 2/2 supported screen/branches.
- `TUBB` / mouse `Tubb5` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `C5orf34` / mouse `4833420G17Rik` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF280D` / mouse `Zfp280d` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF32` / mouse `Zfp637` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `TP53I11` / mouse `Trp53i11` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF821` / mouse `Zfp821` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `C16orf87` / mouse `4921524J17Rik` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF281` / mouse `Zfp281` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `C14orf119` / mouse `1700123O20Rik` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF410` / mouse `Zfp410` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `C9orf85` / mouse `1110059E24Rik` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `C17orf75` / mouse `5730455P16Rik` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF780A` / mouse `Zfp60` (supported_shared_full_matrix_only): 2/2 supported screen/branches.
- `ZNF260` / mouse `Zfp260` (supported_shared_full_matrix_only): 2/2 supported screen/branches.

## Interpretation

- This is now the best ortholog-aware ranking layer because it includes same-symbol and non-identical one-to-one MGI pairs.
- The top manuscript mechanism claims should still prioritize genes that are both biologically interpretable and robust across datasets, not merely broad neuronal/context hits.
- This expanded MGI layer now feeds the formal rank-meta validation; raw-count/object-level DE remains a later optional strengthening step.

## Outputs

- Unit deltas: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_unit_deltas.tsv.gz`
- Branch summary: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_branch_summary.tsv`
- Gene summary: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv`
- Shared hits: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_shared_hits.tsv`
- Mechanism-prioritized hits: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv`
- Plot: `Project/results/primary_core_mgi_ortholog_expanded_meta_model_top_hits.png`
