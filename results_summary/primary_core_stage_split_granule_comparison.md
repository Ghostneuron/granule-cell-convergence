# Stage-Split Granule-Cell Comparison

## Question

Can dentate and cerebellar granule-cell candidates be compared separately at immature and mature stages?

## Short Answer

Yes. The comparison is feasible and biologically useful, but the evidence is asymmetric: dentate has several explicit immature/mature resources, while cerebellum has one explicit staged developmental resource (`GSE122357`: P0 versus P8a/P8b). Therefore this layer should be treated as a stage-aware support analysis rather than a fully powered stage-by-branch statistical test.

## Stage Definitions

- Dentate immature: `GSE104323` Neuroblast/Immature-GC, `GSE214309` immature/immatureactive states, `GSE292261` DG_P5/DG_P7/DG_P10, and `GSE325391` DiffN labels.
- Dentate mature: `GSE104323` GC-juv/GC-adult, `GSE214309` mature/matureactive states, `GSE292261` DG_P15/DG_P28, and `GSE325391` MatN labels.
- Cerebellar immature: `GSE122357` P0 cerebellar granule candidates.
- Cerebellar mature/maturing: `GSE122357` P8a/P8b cerebellar granule candidates.
- `GSE268609` projected labels are retained as supporting, not strict, stage calls. Unstaged atlas resources are not forced into the binary split.

## Coverage

- Strict dentate stage datasets: 4.
- Strict cerebellar stage datasets: 1.
- Stage-called candidate groups: 36.

## Main Result

The stage split supports the manuscript model: the strongest interpretable comparison is not stem/progenitor identity, but postmitotic assembly and maturation. Downstream neurite/morphology and synaptic/excitability modules can be compared separately in immature and mature windows, while regional fate modules remain branch-biased.

Selected similarity scores below use `1 - abs(cerebellar median - dentate median)`, so higher means more similar between branches.

| module | immature similarity | mature similarity | mature - immature |
|---|---:|---:|---:|
| Neurite/morphology | 0.948 | 0.961 | 0.013 |
| Synaptic/excitability | 0.852 | 0.927 | 0.076 |
| Neuronal maturation | 0.812 | 0.838 | 0.025 |
| Immature/progenitor | 0.994 | 0.844 | -0.150 |
| Shared neurogenic niche | 0.877 | 0.940 | 0.063 |

## Interpretation

- If a module is already similar in the immature stage, it likely reflects early postmitotic assembly machinery.
- If similarity increases in the mature stage, it likely reflects later circuit/synaptic maturation or final geometry constraints.
- If similarity decreases, the two branches may begin with shared immature neuronal machinery but diverge as region-specific circuit implementation becomes stronger.
- Because cerebellar stage evidence is currently concentrated in `GSE122357`, any stage-specific cerebellar conclusion should be phrased cautiously.

## Outputs

- Stage-called groups: `Project/results/primary_core_stage_split_granule_group_calls.tsv`.
- Branch-stage module summary: `Project/results/primary_core_stage_split_granule_module_branch_summary.tsv`.
- Cross-branch stage similarity: `Project/results/primary_core_stage_split_granule_stage_similarity.tsv`.
- Mature-minus-immature transition table: `Project/results/primary_core_stage_split_granule_stage_transitions.tsv`.
- Plot: `Project/results/primary_core_stage_split_granule_comparison.png`.
