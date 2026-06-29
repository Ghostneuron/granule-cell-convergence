# First-pass marker-panel interpretation

## What was run

I regenerated the local dataset audit and extracted focused marker-panel expression summaries from the downloaded hippocampal, cerebellar, organoid, and spatial matrices. The current pass is designed as a fast feasibility screen: it checks marker availability and expression detection for dentate identity, cerebellar identity, shared granule-neuronal programs, cytoskeletal/morphogenesis genes, axon-guidance/synapse genes, and a small metabolic/proteomic validation panel.

Key scripts:

- `Project/scripts/audit_dataset_shapes.py`
- `Project/scripts/extract_marker_panel_expression.py`

Key outputs:

- `Project/results/matrix_dimension_audit.tsv`
- `Project/results/marker_gene_expression_summary.tsv`
- `Project/results/marker_panel_expression_summary.tsv`

## Parsing and data sanity

- `GSE214309` is comma-delimited and uses Ensembl gene IDs. It is now parsed as 46,615 genes by 401 nuclei, with sample-title groups recovered from the GEO series matrix, including mature/immature and 1 hr/4 hr activity groups.
- `GSE292261` is cell-by-gene rather than gene-by-cell. It is now parsed as 13,101 genes by 509 cells, grouped by the provided `Sample` metadata such as `DG_P5`, `DG_P7`, `DG_P10`, `DG_P15`, and `DG_P28`.
- Ensembl-to-symbol aliasing is built from local 10x feature files and the Visium archive feature files, so marker genes are recovered across symbol and Ensembl formats.
- Allen Institute files currently contribute taxonomy/reference context, not expression matrices, in this first-pass extraction.

## Main biological signal

The project is reasonable and worth pursuing. The first-pass data support a useful central hypothesis: hippocampal dentate granule cells and cerebellar granule cells are not the same cell type, but they may converge on a shared structural and wiring program that produces similar compact soma, simple dendritic architecture, and long-range axonal output logic.

The dentate datasets behave as expected. In `GSE104323`, hippocampal `GC-adult`, `GC-juv`, `Immature-GC`, and `Neuroblast` groups show high dentate identity detection and very low cerebellar identity detection. For example, `GC-adult` has dentate detection 0.238 versus cerebellar detection 0.0049, while `Immature-GC` has dentate detection 0.401 versus cerebellar detection 0.0071. `GSE292261` is also strongly dentate: all DG postnatal groups recover 12/12 dentate markers, while only 3/12 cerebellar identity markers are found in that Smart-seq table.

The cerebellar datasets show the complementary pattern, although many are currently whole cerebellum or mixed-cell matrices rather than purified granule-cell-only tables. `GSE122357`, `GSE165657`, and `GSE312658` show broad cerebellar identity, strong shared neuronal modules, and strong cytoskeletal/morphogenesis modules. This supports moving to cell-level granule-cell selection rather than interpreting whole-sample averages.

The most important positive signal is the shared module. Across DG and cerebellar resources, `shared_granule_neuronal`, `morphogenesis_cytoskeleton`, and `axon_guidance_synapse` panels are consistently detected. This is exactly the right shape for a novel project: distinct upstream regional identity programs, but partially convergent downstream morphogenesis and wiring programs.

## Caveats

These values are not final cross-platform statistics. The current tables mix 10x, Smart-seq, snRNA-seq, organoid data, and Visium-like spatial matrices, so raw mean expression values should not be compared directly across datasets. Detection fractions and marker coverage are useful for triage, but the next analysis must normalize within dataset, annotate cell types, subset granule cells, and then compare standardized module scores.

The marker panels should be refined before final analysis. `Tau` is a protein/common-name alias and should probably be represented by `Mapt` only. `Atp5f1a` is absent in several feature lists and may need an older alias check such as `Atp5a1`, depending on the dataset annotation.

## Recommended strategy

1. Build harmonized Seurat or SingleCellExperiment objects for each dataset, preserving dataset, species, age, region, genotype, and cell-type metadata.
2. Annotate and subset dentate granule cells and cerebellar granule cells before cross-region comparison.
3. Convert mouse, rat, and human genes through ortholog mapping, then score gene modules per cell using rank-based methods such as AUCell/UCell or carefully controlled Seurat module scores.
4. Test the core model explicitly: region-identity modules should separate DG and cerebellar granule cells, while morphogenesis, cytoskeletal, synapse, and axon-guidance modules should partially converge.
5. Add trajectory analysis for developmental datasets to ask whether convergence increases during granule-cell maturation.
6. Add regulatory-network analysis, such as SCENIC or motif enrichment, to separate upstream specification factors from downstream shared morphology executors.
7. Connect the transcriptomic modules back to the morphology question using your original thesis/paper, public morphology resources, and literature curation.

## Go/no-go judgment

Go. The project is reasonable, novel if framed mechanistically, and important enough for either a computational paper or a review-plus-reanalysis. The strongest paper angle is not simply "they look similar"; it is "two distinct granule-cell lineages converge on a conserved structural/wiring gene program while retaining region-specific identity control."
