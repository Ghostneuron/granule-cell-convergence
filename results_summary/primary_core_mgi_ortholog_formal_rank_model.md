# Formal MGI Ortholog Rank-Meta Validation

Date built: 2026-06-22

## Purpose

This analysis adds a stricter statistical validation layer to the expanded MGI ortholog meta-model. It uses the existing pseudobulk rank deltas and tests whether granule-cell candidate classes are consistently above branch-specific background classes across independent datasets.

This remains a rank-meta pseudobulk model, not raw-count DESeq2/edgeR differential expression. Its strength is cross-dataset replication; its limitation is small independent dataset count in several branches.

## Input Scope

- Unit delta rows: 252,469.
- Dataset-level delta rows: 116,013.
- Branch tests: 36,303.
- Gene summaries: 16,708.
- Mechanism-prioritized genes modeled: 36.

## Test Definitions

- Replication support: at least 2 datasets, median dataset rank delta > 0, and positive dataset fraction >= 0.75.
- Nominal support: replication support plus best one-sided dataset-level p <= 0.25 from t, sign, or Wilcoxon tests.
- FDR10 support: replication support plus best within-screen/branch BH q <= 0.10.
- Shared support requires both dentate and cerebellar branches to pass in the selected screen, full-matrix screen, or both.

## Summary

- Branches with replication support: 11,337.
- Branches with nominal support: 11,337.
- Branches with FDR10 support: 1,102.
- Formally shared hits: 1,370.
- Replication-shared hits: 1,370.
- Nominal-shared hits: 1,370.
- FDR10-shared hits: 0.
- Both-screen shared hits: 158.
- Mechanism model methods: mixedlm_dataset_random_intercept: 120, cluster_robust_intercept_fallback: 16.

## Mechanism-Prioritized Result

- `GPM6A` (dataset_robust_consensus_figure; cytoskeleton_morphogenesis): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `NFIB` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `KCNK1` (dataset_robust_consensus_figure; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `RFX3` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `GABRA2` (dataset_robust_consensus_figure; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `NFIA` (dataset_robust_consensus_figure; regulatory_morphogenesis_candidate): formal_nominal_shared_both_screens, 4/4 nominal branches, 0 model-supported branches.
- `GABRB3` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 3 model-supported branches.
- `KCND2` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 3 model-supported branches.
- `PPP3CA` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `CACNA2D1` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `KCNJ6` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `GRIN2B` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `KCNJ3` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `STXBP5L` (consensus_figure_candidate; synaptic_wiring): formal_nominal_shared_both_screens, 4/4 nominal branches, 1 model-supported branches.
- `ROBO2` (dual_screen_mechanism_figure; curated_shared_structural_executor): formal_nominal_shared_both_screens, 4/4 nominal branches, 3 model-supported branches.
- `DYNLL1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_both_screens, 4/4 nominal branches, 4 model-supported branches.
- `BASP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_both_screens, 4/4 nominal branches, 4 model-supported branches.
- `RFX7` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate): formal_nominal_shared_both_screens, 4/4 nominal branches, 4 model-supported branches.
- `MAPKAP1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_both_screens, 4/4 nominal branches, 4 model-supported branches.
- `TUBA1A` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_both_screens, 4/4 nominal branches, 2 model-supported branches.
- `ACTB` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_full_matrix, 3/4 nominal branches, 3 model-supported branches.
- `RTN3` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_full_matrix, 3/4 nominal branches, 3 model-supported branches.
- `ACTG1` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_full_matrix, 3/4 nominal branches, 3 model-supported branches.
- `TUBA1B` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_full_matrix, 3/4 nominal branches, 3 model-supported branches.
- `MAP1B` (selected_screen_mechanism_figure; cytoskeleton_morphogenesis): formal_nominal_shared_selected, 3/4 nominal branches, 1 model-supported branches.
- `TCF4` (selected_screen_mechanism_figure; regulatory_morphogenesis_candidate): formal_nominal_shared_full_matrix, 3/4 nominal branches, 1 model-supported branches.
- `STMN2` (selected_screen_mechanism_figure; curated_shared_structural_executor): formal_nominal_shared_full_matrix, 3/4 nominal branches, 1 model-supported branches.
- `CADM3` (selected_screen_mechanism_figure; axon_guidance_adhesion): formal_nominal_shared_full_matrix, 3/4 nominal branches, 1 model-supported branches.
- `KCNMB4` (selected_screen_mechanism_figure; synaptic_wiring): formal_nominal_shared_full_matrix, 3/4 nominal branches, 1 model-supported branches.
- `DCC` (selected_screen_mechanism_figure; axon_guidance_adhesion): formal_nominal_shared_full_matrix, 3/4 nominal branches, 1 model-supported branches.

## Consensus Six-Gene Check

- The dataset-robust consensus genes retained after formal rank-meta validation are: `GPM6A`, `NFIB`, `KCNK1`, `RFX3`, `GABRA2`, `NFIA`.

## Top Formal Shared Hits

- `CELF2` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.98.
- `FTH1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.97.
- `GPM6A` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.93.
- `PPP3CA` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.89.
- `NFIB` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.89.
- `DDX5` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.88.
- `DYNLL1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.85.
- `YWHAE` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.83.
- `PTMA` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.83.
- `NFIA` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.82.
- `HMGB1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.81.
- `NAP1L1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.80.
- `BASP1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.79.
- `GSK3B` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.78.
- `MORF4L1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.77.
- `NCOR1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.77.
- `ANK3` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.76.
- `RBFOX3` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.74.
- `SRRM2` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.72.
- `SFPQ` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.72.
- `NBEA` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.72.
- `SF3B1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.70.
- `NRN1` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.66.
- `TUBA1A` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.65.
- `NDUFB8` (formal_nominal_shared_both_screens; same_symbol): 4/4 nominal branches, score 17.64.

## Interpretation

- This formal layer supports the project strategy: the strongest genes are not just morphology-plausible, they replicate as dentate-plus-cerebellar granule-cell enriched across independent datasets.
- The six-gene consensus (`GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`) remains the safest manuscript seed set.
- Additional synaptic and wiring genes are useful as pathway/context support, while non-identical-symbol ortholog hits should remain secondary until raw-count or external validation is added.

## Outputs

- Dataset deltas: `Project/results/primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz`
- Branch tests: `Project/results/primary_core_mgi_ortholog_formal_rank_branch_tests.tsv`
- Gene summary: `Project/results/primary_core_mgi_ortholog_formal_rank_gene_summary.tsv`
- Formal shared hits: `Project/results/primary_core_mgi_ortholog_formal_rank_shared_hits.tsv`
- Mechanism model long table: `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_model_long.tsv`
- Mechanism summary: `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv`
- Plot: `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.png`
