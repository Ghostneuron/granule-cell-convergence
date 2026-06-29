# Refined panel comparison interpretation

## What was refined

I created a refined marker-panel file that keeps strict identity genes separate from broad neurogenic, immature-neuron, and excitatory-neuron genes. The broad genes were moved into `shared_granule_neuronal` rather than being used as strict dentate or cerebellar identity evidence.

New/refined files:

- `Project/config/granule_marker_panels_refined.tsv`
- `Project/results/refined_per_cell_marker_module_scores.tsv.gz`
- `Project/results/refined_per_cell_marker_module_score_summary.tsv`
- `Project/results/refined_candidate_granule_cell_calls.tsv.gz`
- `Project/results/refined_candidate_granule_cell_call_summary.tsv`
- `Project/results/refined_candidate_granule_program_class_summary.tsv`
- `Project/results/refined_candidate_granule_program_statistics.tsv`
- `Project/results/refined_candidate_granule_program_comparison.png`

## Effect of refinement

The refined panels improved the core identity separation.

Original calls:

- 30,797 candidate cerebellar granule-like cells/spots
- 11,259 candidate dentate granule cells
- 22,000 cerebellum dentate-panel warning calls

Refined calls:

- 36,982 candidate cerebellar granule-like cells/spots
- 11,200 candidate dentate granule cells
- 17,122 cerebellum dentate-panel warning calls

The warning category decreased by 4,878 calls, and cerebellar candidate calls increased by 6,185, which means moving broad markers out of the strict dentate identity panel made the cerebellar classification cleaner.

## Refined statistics

Identity separation remains very strong:

- Dentate candidate versus cerebellar candidate identity contrast: median 0.339 versus -0.209; Cliff's delta 0.999.
- Dentate candidate versus known non-dentate reference identity contrast: median 0.339 versus 0.0236; Cliff's delta 0.796.
- Cerebellar candidate versus ambiguous cells cerebellar identity: median 0.272 versus 0.0999; Cliff's delta 0.747.

Structural-program ranks remain high in candidate populations:

- Dentate candidates: median structural rank 0.524.
- Cerebellar candidates: median structural rank 0.728.
- Ambiguous cells/spots: median structural rank 0.314.

This is the cleanest evidence so far for the working model: identity modules distinguish the two populations, while structural/morphogenesis modules are active in both.

## Interpretation

The refined analysis strengthens the project. It shows that the shared morphology question should be framed as convergence of downstream structural and wiring programs, not equivalence of dentate and cerebellar granule-cell identity.

The persistent warning category is still useful. It shows that some cerebellar cells/spots activate genes historically assigned to dentate-like or broad neurogenic modules. This is not a reason to stop; it is a reason to make the final paper more precise about strict identity genes versus shared neuronal morphogenesis genes.

## Recommended next move

Use the refined outputs as the current working analysis layer. The next technical step should be full object-level normalization and annotation:

1. Build Seurat/SCE objects for the main datasets.
2. Add refined module scores plus control-gene or rank-based scoring.
3. Subset candidate dentate and cerebellar granule cells.
4. Test identity separation and structural-program convergence with dataset-aware models.
5. Add ortholog mapping before any mouse-human comparison.

The refined result is now strong enough to support a project proposal or the first Results figure concept.
