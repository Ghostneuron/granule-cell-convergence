# Candidate granule-program comparison interpretation

## What was added

I compared the working candidate dentate and cerebellar granule-cell subsets using the candidate-call table. This stage asks the central manuscript question more directly:

Can dentate and cerebellar granule candidates be separated by regional identity while still sharing structural/morphogenesis program activity?

New files:

- `Project/scripts/compare_candidate_granule_programs.py`
- `Project/results/candidate_granule_program_class_summary.tsv`
- `Project/results/candidate_granule_program_statistics.tsv`
- `Project/results/candidate_granule_program_comparison.png`

## Main result

The candidate classes separate strongly by identity. Dentate candidates have positive dentate-minus-cerebellar identity contrast, while cerebellar candidates have negative contrast.

Key statistics:

- Dentate candidate versus cerebellar candidate identity contrast: median 0.356 versus -0.149; Cliff's delta 0.997; Mann-Whitney p reported as 0 due numerical underflow.
- Dentate candidate versus known non-dentate reference identity contrast: median 0.356 versus 0.0578; Cliff's delta 0.607.
- Cerebellar candidate versus ambiguous cells cerebellar-identity score: median 0.375 versus 0.116; Cliff's delta 0.828.

This supports the first half of the model: the two granule-cell populations are not just one generic granule identity. They retain distinct regional identity signatures.

## Structural-program signal

Both candidate classes show structural-program activity when measured as within-sample rank:

- Dentate candidates: median structural rank 0.529.
- Cerebellar candidates: median structural rank 0.678.
- Ambiguous/non-granule cells: median structural rank 0.276.

Cerebellar candidates are strongly enriched over ambiguous cells for structural-program rank, while dentate candidates are only modestly enriched over the known non-dentate reference set. This is not a failure; it is an important refinement. The non-dentate reference set includes other neurons and immature excitatory cells, so cytoskeletal, axon-guidance, and synaptic genes are expected to be active there too.

The correct interpretation is therefore not "structural genes identify granule cells." It is: distinct dentate and cerebellar granule-cell lineages may converge on a shared downstream structural/wiring toolkit that is also partly reused by other developing neurons.

## Warning category

The `cerebellum_warning` class has high structural rank and positive dentate-panel contrast. This says the current dentate identity panel still contains genes that are too broad for final claims. The most likely culprits are general neurogenic/excitatory or immature-neuron genes that should move into shared or developmental panels rather than strict dentate identity.

Before final statistics, the project should refine gene panels into:

- strict dentate-region identity
- strict cerebellar granule identity
- shared granule/excitatory neuronal state
- structural/morphogenesis/wiring executors

## Manuscript-level implication

The strongest current claim is:

Dentate and cerebellar granule cells preserve different regional identity programs but show partial convergence in downstream neuronal structural and wiring modules. This convergence provides a plausible molecular explanation for similar compact granule-cell morphology despite distinct developmental origins and circuit contexts.

## Next step

The immediate next step should be a panel-refinement and re-scoring round. After that, the analysis can move to full normalization, ortholog-aware scoring, and cell-level statistical modeling in Seurat/SCE.
