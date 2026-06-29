# Human dentate core construction status

Date updated: 2026-06-22

## Completed First Seed Build

`GSE185277` and `GSE185553` are now local raw downloads and have first-pass sparse per-library objects.

- `GSE185277`: 7 libraries, 89,526 cells/barcodes, 60,286,599 non-zero count entries, 99,592,073 total counts.
- `GSE185553`: 27 libraries, 97,729 cells/barcodes, 103,757,887 non-zero count entries, 205,573,021 total counts.

The objects are stored under `Project/processed/human_seed_sparse_objects/` in cells-by-genes SciPy sparse format, with matching barcodes, genes, cell metadata, and gene metadata files.

## QC Harmonization

The built human branch now has a harmonized cell QC table covering `GSE185277`, `GSE185553`, and the `GSE186538` DG GC subset:

- `GSE185277`: 89,526 cells, median 518 counts, median 411 genes, median 1.43% MT, 50.5% preliminary QC pass.
- `GSE185553`: 97,729 cells, median 997 counts, median 708 genes, median 0.11% MT, 80.4% preliminary QC pass.
- `GSE186538`: 32,067 cells, median 9,229 counts, median 3,583 genes, median 0.58% MT, 99.8% preliminary QC pass.

The preliminary pass rule is diagnostic only: `n_counts >= 500`, `n_genes >= 300`, and `percent_mt <= 20`.

## Marker Validation

The human branch now has a marker-validation layer using refined granule-cell identity panels, human DG maturation panels, and non-neuronal background panels:

- `GSE185277`: 11,044 marker-supported human DG-like cells and 7,846 immature-neuron/neurogenic candidates, alongside low-QC/background/ambiguous cells.
- `GSE185553`: 13,197 marker-supported human DG-like cells and 6,825 immature-neuron/neurogenic candidates, with a larger expected broader-hippocampal/background component.
- `GSE186538`: 32,067 cells remain the curated human DG GC reference because the source taxonomy labels are `DG GC PROX1 SGCZ` or `DG GC PROX1 PDLIM5`.

The key marker result is that `GSE186538` carries strong DG/shared-neuronal signal with essentially zero median cerebellar-identity score, making it the current human dentate anchor.

## GEO Metadata Curation

Lightweight GEO SOFT metadata were parsed for `GSE185277`, `GSE185553`, and `GSE186538`. The curated tables now map all 35 built human components to sample/specimen metadata, and all six `GSE186538` human DG donors to donor-level age/sex/tissue records.

Main new metadata outputs:

- `Project/results/human_core_component_metadata_curated.tsv`
- `Project/results/human_core_enriched_cell_metadata.tsv.gz`
- `Project/results/human_core_gse186538_dg_donor_summary.tsv`

## Normalized Reduced Object

A portable normalized reduced object has now been built for the current human core. It uses preliminary QC-pass cells, log1p(CP10K) normalization, and a selected feature set combining high-information genes with marker-panel genes.

- Shape: 155,828 cells/nuclei by 2,169 selected genes.
- Included cells: `GSE185277` 45,219 / 89,526, `GSE185553` 78,605 / 97,729, and `GSE186538` 32,004 / 32,067.
- Feature set: 112 marker-panel genes plus high-information genes.
- First SVD check: 30 components computed; first five components explain 22.8% of variance.

Main object outputs:

- `Project/processed/human_core_normalized_reduced_object/X_log1p_cp10k_selected_genes.npz`
- `Project/processed/human_core_normalized_reduced_object/obs.tsv.gz`
- `Project/processed/human_core_normalized_reduced_object/var.tsv`
- `Project/results/human_core_normalized_svd_pc12.png`

## Label Tuning And Dataset-Aware Module Tests

The normalized reduced object has now been used to define a first stable human dentate/hippocampal label convention and to run replicate-aware module checks.

- Tuned-label object: 155,828 unique QC-pass cells/nuclei.
- Curated human DG GC anchor: 32,004 cells from `GSE186538`.
- Human DG-like high-confidence labels: 22,043 cells.
- Immature/neurogenic candidate labels: 13,376 cells.
- Non-neuronal/background and lower-support labels remain explicit rather than silently removed.

Dataset-aware module testing used specimen-level replicate units for `GSE185277` and `GSE185553`. High-confidence DG-like labels show higher dentate/shared/morphogenesis or axon-guidance module signal and lower background signal than non-neuronal background; the larger `GSE185553` replicate set gives the strongest adjusted significance.

Main new outputs:

- `Project/results/human_core_tuned_labels.tsv.gz`
- `Project/results/human_core_tuned_label_summary.tsv`
- `Project/results/human_core_dataset_aware_module_tests.tsv`
- `Project/results/human_core_tuned_label_module_heatmap.png`
- `Project/results/human_core_dataset_aware_module_deltas.png`
- `Project/results/human_core_label_tuning_and_module_tests.md`

## GSE325391 Adult Dentate Primary Anchor

`GSE325391` is now local and promoted from acquisition candidate to built primary adult human dentate anchor.

- Downloaded adult RDS: 3,515,050,479 bytes, matching the GEO file size.
- Source object: Seurat, 59,075 nuclei by 36,588 features, with RNA/spliced/unspliced assays and harmony/UMAP reductions.
- Source scope: adult dentate granule-cell nuclei only in the adult RDS; GEO series metadata also contains three fetal records, but the fetal RDS remains deferred.
- Selected sparse bridge: 59,075 nuclei by the same 2,169 human-core selected genes; 2,083 selected genes present and 86 absent/zero-filled.
- Label projection: 44,239 mature adult DG anchors, 12,611 differentiating adult DG anchors, 1,712 doublet flags, and 513 background-warning cells.

Main new outputs:

- `External_Data/GEO/GSE325391/GSE325391_adultgc_filtered.RDS`
- `Project/processed/gse325391_adult_dg_selected/matrix_cells_by_selected_genes.npz`
- `Project/results/gse325391_adult_rds_inspection.md`
- `Project/results/gse325391_selected_npz_summary.md`
- `Project/results/gse325391_human_core_label_projection.md`
- `Project/results/gse325391_human_core_module_tests.tsv`
- `Project/results/gse325391_human_core_label_umap.png`

## GSE268609 Human Hippocampal Multiome RNA Expansion

`GSE268609` is now local as an RNA-row selected-gene bridge and promoted from "download next" to built broader human hippocampal aging/AD context.

- GEO sample metadata: 78 records parsed, consisting of 39 RNA and 39 ATAC records.
- RNA sample scope: 8 AD, 9 HA, 8 MCI, 6 SA, and 8 YA RNA samples.
- Downloaded RNA files: matrix, barcodes, features, and SOFT metadata all match expected byte sizes.
- Source mixed multiome matrix: 177,931 rows by 366,175 columns with 1,451,618,555 non-zero entries; 36,601 rows are gene-expression features and 141,330 rows are ATAC peaks.
- Selected sparse RNA bridge: 366,175 nuclei/barcodes by 2,169 human-core selected genes; 2,083 selected genes present, 232,722,561 non-zero selected entries, and 1,023,207,067 selected counts.
- Basic RNA QC: 202,431 cells pass the source-like `nCount_RNA >= 2000` and `nFeature_RNA >= 1000` rule.
- Human-core label projection: 94,430 cells retained in the broad analysis include after basic QC plus marker/module filtering; 16,051 human DG-like candidates, 14,692 immature/neurogenic candidates, and 2,132 high-confidence DG-like projections.

Main new outputs:

- `External_Data/GEO/GSE268609/GSE268609_matrix.mtx.gz`
- `Project/processed/gse268609_rna_selected/matrix_cells_by_selected_genes.npz`
- `Project/results/gse268609_geo_metadata_summary.md`
- `Project/results/gse268609_rna_file_inspection.md`
- `Project/results/gse268609_selected_npz_summary.md`
- `Project/results/gse268609_human_core_label_projection.md`
- `Project/results/gse268609_human_core_diagnosis_module_tests.tsv`

## Integrated Human Bridge And Backbone Rank Layer

The built human dentate/hippocampal bridge objects have now been integrated with the existing refined mouse/human cerebellar-dentate backbone using within-sample rank metrics. This is a bridge-level sanity-check analysis, not yet a replacement for a single fully harmonized cross-species expression object.

The strict manuscript-scale primary core is now captured as a 10-dataset table:

- `Project/results/integrated_primary_core_datasets.tsv`
- `Project/results/integrated_primary_core_datasets.md`

- Integrated units: 635 summarized dataset/sample/label units.
- Statistics/plotting units after the `>=20` cells/spots filter: 546 units.
- Dentate candidates separate from cerebellar candidates by dentate-minus-cerebellar identity rank contrast: median delta 0.6472, BH-adjusted p=1.57e-10.
- Dentate candidates separate from non-dentate/background units by identity rank contrast: median delta 0.3423, BH-adjusted p=9.93e-22.
- Structural rank behaves as the expected convergence axis rather than a strict identity separator: dentate and cerebellar candidate structural ranks are both above the within-sample median, while their direct structural-rank difference is not significant.
- `GSE325391` and `GSE186538` are strongest as source-aware human dentate anchors; because they are DG-enriched, they are more informative for structural/neurogenic state than for within-sample dentate-versus-cerebellar identity contrast.
- `GSE268609` adds broader human hippocampal aging/AD context with useful projected dentate signal, but should still be described as projection-labeled unless the full source taxonomy is added.

Main new outputs:

- `Project/results/human_bridge_backbone_rank_units.tsv`
- `Project/results/human_bridge_backbone_rank_statistics.tsv`
- `Project/results/human_bridge_backbone_rank_source_summary.tsv`
- `Project/results/human_bridge_backbone_rank_units.png`
- `Project/results/human_bridge_backbone_rank_integration.md`

## Human Bridge Candidate Gene Packet

A first manuscript-planning candidate gene packet has been built from the refined marker panels across the constructed human bridge objects (`human_core_tuned`, `GSE325391`, and `GSE268609`).

- Gene scope: 67 refined panel genes.
- Summary rows: 670 source-layer/class/gene rows.
- High-priority shared structural/executor genes are those detected in dentate candidates across all three human bridge layers at `>=20%` detection.
- Current top shared structural/executor candidates include `NCAM1`, `MAPT`, `PLXNA4`, `DPYSL2`, `CNTN5`, `ROBO2`, `NRP1`, `EPHA4`, `ROBO1`, `PLXNA2`, `SLIT2`, `STMN2`, and `GAP43`.
- Dentate identity genes supported in the human bridge include `PROX1`, `CALB1`, `C1QL3`, `GLIS3`, `EGR3`, and `ITPKA`.
- Cerebellar panel genes with unexpectedly detectable human dentate-bridge signal, such as `ROR1`, `ETV1`, and `GABRA6`, are specificity warnings rather than evidence that DG cells adopt cerebellar identity.

Main new outputs:

- `Project/results/human_bridge_marker_gene_summary.tsv`
- `Project/results/human_bridge_candidate_gene_packet.tsv`
- `Project/results/human_bridge_structural_executor_candidates.png`
- `Project/results/human_bridge_candidate_gene_packet.md`

## Primary Core Ortholog-Aware Module Analysis

The strict 10-dataset primary core now has a first formal module-level analysis. This freezes the core datasets, adds a marker-panel mouse-human ortholog map, filters the integrated units to the 10 primary datasets, and runs unit/sample-aware tests with leave-one-dataset-out checks.

- Primary core datasets: 10.
- Integrated primary-core units: 418.
- Eligible units after the `>=20` cells/spots filter: 361.
- Identity separation is stable: dentate candidates exceed cerebellar candidates in dentate-minus-cerebellar identity rank, median delta 0.5356, BH-adjusted p=2.95e-08.
- Leave-one-dataset-out identity separation remains positive and BH-significant for every held-out dataset.
- Structural programs are elevated in both candidate classes: dentate candidate structural rank is above median, BH-adjusted p=0.000117; cerebellar candidate structural rank is above median, BH-adjusted p=0.0156.
- Structural magnitude is not equal: cerebellar candidates are higher than dentate candidates in this strict-core rank layer, so the manuscript should frame this as a shared elevated executor axis, not identical structural state.

Main new outputs:

- `Project/results/primary_core_marker_panel_ortholog_map.tsv`
- `Project/results/primary_core_integrated_module_units.tsv`
- `Project/results/primary_core_ortholog_module_statistics.tsv`
- `Project/results/primary_core_ortholog_module_leave_one_dataset_out.tsv`
- `Project/results/primary_core_identity_structural_modules.png`
- `Project/results/primary_core_ortholog_module_analysis.md`

## Primary Core Candidate-Gene Pseudobulk Analysis

The strict 10-dataset primary core now also has a DE-adjacent candidate-gene pseudobulk layer. This analysis aggregates the 67 human-bridge candidate genes by refined broad class inside each primary dataset/sample and tests whether proposed structural-executor genes rank higher in candidate dentate and cerebellar granule populations than in local backgrounds.

- Primary core datasets represented: 10/10.
- Pseudobulk rows: 16,719.
- Candidate genes tested: 67.
- Structural-executor genes with positive dentate and cerebellar candidate deltas: 14.
- Structural-executor genes passing the stricter exploratory BH<0.2 rule in both branches: 12.
- Strongest shared structural-executor signals: `CFL1`, `GAP43`, `ROBO2`, `STMN2`, `STMN3`, `CDK5R1`, `DPYSL2`, `EPHB2`, `L1CAM`, `MAPT`, `DPYSL3`, and `ELAVL4`.
- `NCAM1` and `CNTN6` are positive in both branches but do not pass the stricter two-branch exploratory cutoff.
- Some axon-guidance genes remain dentate-biased or branch-skewed in this pass (`CNTN5`, `NRP1`, `EPHA4`, `EPHB1`, `NRP2`, `PLXNA4`), which is useful for separating universal morphology executors from region-specific wiring programs.

Main new outputs:

- `Project/results/primary_core_candidate_gene_pseudobulk_expression.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_hits.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_effects.png`
- `Project/results/primary_core_candidate_gene_pseudobulk_analysis.md`

## Expanded Selected-Gene Pseudobulk And Mechanism Triage

The strict 10-dataset primary core now has an expanded discovery screen over the 2,169-gene human-core selected-feature universe. This is not yet whole-transcriptome mixed-effect DE, but it is the first broad ortholog-ready discovery layer that keeps all 10 primary datasets in one target-gene frame.

- Selected-gene universe: 2,169 genes.
- Pseudobulk expression rows: 537,070.
- Primary datasets represented: 10/10.
- Shared-positive rank genes: 1,279.
- Shared-positive genes passing BH<0.10 in both branches: 579.
- Original 67-gene packet genes recovered among shared-positive hits: 42.
- Mechanism triage nominates 71 figure-level candidates and 76 follow-up candidates.
- Strongest recovered structural candidates: `STMN2`, `GAP43`, `ROBO2`, `CFL1`, and `STMN3`.
- Strong new mechanism candidates include `GPM6A`, `MAP1B`, `TCF4`, `NFIA`, `NFIB`, `PPP3CA`, `CALM1`, `RTN3`, `RTN1`, `NRXN1`, `BASP1`, `CAMTA1`, `CADM1`, and `DCC`.

Main new outputs:

- `Project/results/primary_core_expanded_gene_pseudobulk_expression.tsv.gz`
- `Project/results/primary_core_expanded_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_branch_specific.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.png`
- `Project/results/primary_core_expanded_gene_pseudobulk_analysis.md`
- `Project/results/primary_core_expanded_gene_mechanism_triage.tsv`
- `Project/results/primary_core_expanded_gene_mechanism_triage.md`

## Full-Matrix Same-Symbol Pseudobulk And Consensus

The strict primary core now has a full-matrix same-symbol pseudobulk screen for the datasets where local full matrices are usable without exporting large source objects. This analysis uses the GSE186538 human DG full gene list as the target symbol universe and maps mouse genes by same-root upper-case symbol. It is ortholog-ready, not yet a curated ortholog model.

- Full-symbol target universe: 33,939 genes.
- Full-matrix expression rows: 628,339.
- Full-matrix datasets represented: 8/10 primary datasets.
- Datasets contributing to rank contrasts: 7 (`GSE104323`, `GSE122357`, `GSE165657`, `GSE214309`, `GSE292261`, `GSE312658`, `GSE95752`).
- Genes tested in contrast statistics: 21,253.
- Shared-positive same-symbol genes: 6,440.
- Shared-positive genes passing BH<0.10 in both branches: 283.
- Genome-wide mechanism triage nominates 33 figure-level candidates and 261 follow-up candidates.
- Cross-screen consensus identifies 24 figure-level mechanism candidates that survive both the 2,169-gene selected-feature screen and the full-matrix same-symbol screen.
- Top cross-screen consensus candidates include `GPM6A`, `NFIA`, `NFIB`, `PPP3CA`, `CAMTA1`, `MAPK1`, `STXBP1`, `CALM2`, `CACNA2D1`, `SYNPR`, `KCNK1`, `GABRB3`, `ADD2`, `KCNJ6`, `RFX3`, `KCNJ3`, `GRIN2B`, `MAP3K4`, `KCND3`, `KCND2`, `GABRA2`, `STXBP5L`, `CACNA1E`, and `STXBP5`.

Main new outputs:

- `Project/results/primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_coverage.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_statistics.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_branch_specific.tsv`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.png`
- `Project/results/primary_core_genomewide_symbol_pseudobulk_analysis.md`
- `Project/results/primary_core_genomewide_symbol_mechanism_triage.tsv`
- `Project/results/primary_core_genomewide_symbol_mechanism_triage.md`
- `Project/results/primary_core_cross_screen_mechanism_consensus.tsv`
- `Project/results/primary_core_cross_screen_mechanism_consensus.md`

## Dataset-Aware Consensus Candidate Validation

The 24 cross-screen consensus figure candidates have now been tested for robustness across individual dataset/sample units rather than only pooled pseudobulk summaries.

- Consensus candidates tested: 24.
- Dataset/sample/gene branch-delta rows: 1,608.
- Robustness rule: a screen/branch is robust when it has at least 2 dataset/sample units, at least 75% positive candidate-versus-background deltas, and median rank delta greater than 0.
- Genes robust across all available selected-feature and full-matrix screen/branch tests: `GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, and `RFX3`.
- Genes robust in 3/4 screen/branch tests: `CACNA2D1`, `GABRB3`, `GRIN2B`, `KCND2`, `KCNJ3`, `KCNJ6`, `PPP3CA`, and `STXBP5L`.
- This makes the strongest current manuscript mechanism shortlist smaller and more dataset-aware than the 24-gene consensus list.

Main new outputs:

- `Project/results/primary_core_consensus_candidate_dataset_deltas.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation_heatmap.png`
- `Project/results/primary_core_consensus_candidate_dataset_validation.md`

## MGI Ortholog Meta-Model

The strict primary core now has a conservative ortholog-aware meta-model using the official MGI human-mouse homology report. This model keeps one-to-one human-mouse homology classes and, because the current extraction is same-symbol based, restricts strict gene-level claims to one-to-one classes where the human and mouse symbols collapse to the same canonical symbol.

- MGI report rows: 46,522.
- MGI human-mouse homology classes: 20,181.
- One-to-one human-mouse pairs: 17,611.
- Strict same-symbol one-to-one pairs: 16,245.
- Strict pairs represented in selected-feature expression rows: 1,903.
- Strict pairs represented in full-matrix expression rows: 15,417.
- Dataset/sample candidate-versus-background unit deltas: 241,690.
- Shared hits after ortholog filtering: 1,304.
- Mechanism-prioritized shared hits: 36.
- Strongest dataset-robust consensus mechanism genes retained in the strict MGI frame: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- Additional consensus figure candidates with strict both-screen ortholog meta-support include `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, and `STXBP5L`; `ROBO2` is retained as a dual-screen curated structural-executor hit.

Main new outputs:

- `Project/results/primary_core_mgi_ortholog_meta_model_map.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_meta_model_branch_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_top_hits.png`
- `Project/results/primary_core_mgi_ortholog_meta_model.md`

## Expanded MGI Ortholog Full-Matrix And Meta-Model

The MGI ortholog layer has now been expanded beyond same-symbol genes by resolving mouse matrix rows through the MGI mouse-symbol side of one-to-one human-mouse homology classes.

- MGI one-to-one target genes: 17,611.
- Non-identical human/mouse symbol targets: 1,366.
- Ortholog full-matrix pseudobulk rows: 508,917.
- Genes tested in ortholog full-matrix contrast statistics: 16,704.
- Shared-positive ortholog genes: 6,413.
- Shared-positive non-identical-symbol genes: 246.
- Shared-positive genes passing BH<0.10 in both branches: 286.
- Expanded meta-model unit delta rows: 252,469.
- Expanded meta-model shared hits: 1,370.
- Expanded meta-model non-identical-symbol shared hits: 64, including 11 strict full-matrix-only hits.
- Mechanism-prioritized shared hits remain 36.
- The six dataset-robust consensus genes remain the strongest mechanism tier: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- The strongest non-identical-symbol hits are mostly full-matrix-only or context candidates, led by `ZNF148`/mouse `Zfp148`, `C1orf21`/mouse `1700025G04Rik`, `ZNF706`/mouse `Zfp706`, `RAB7A`/mouse `Rab7`, `ZNF292`/mouse `Zfp292`, and `ZNF827`/mouse `Zfp827`.

Main new outputs:

- `Project/results/primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_full_matrix_coverage.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_statistics.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_analysis.md`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_unit_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model.md`

## Formal MGI Ortholog Rank-Meta Validation

The expanded MGI ortholog layer now has a formal dataset-level rank-meta validation. This tests whether each gene's granule-cell candidate class is consistently above branch-specific background classes across independent datasets, then adds mixed/intercept model checks for the 36 mechanism-prioritized genes.

- Unit delta rows: 252,469.
- Dataset-level delta rows: 116,013.
- Branch tests: 36,303.
- Gene summary rows: 16,708.
- Formal shared hits: 1,370.
- Both-screen formal shared hits: 158.
- Mechanism-prioritized genes modeled: 36.
- The six dataset-robust consensus genes remain the strongest mechanism seed set: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- Strong second-tier formal mechanism candidates include `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`, and `ROBO2`.

Main new outputs:

- `Project/results/primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_formal_rank_branch_tests.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_model_long.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.png`
- `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`

## Manuscript Candidate Tier Packet

The formal rank-meta results have also been distilled into a compact manuscript-facing candidate table.

- Candidate tier rows: 66.
- Tier 1 core convergent-program genes: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- Tier 2 high-confidence wiring/synaptic/executor genes: `ROBO2`, `GABRB3`, `KCND2`, `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GRIN2B`, `KCNJ3`, and `STXBP5L`.
- Tier 5 contains 30 exploratory non-identical-symbol orthologs, kept outside central claims until external or raw-count validation.

Main new outputs:

- `Project/results/primary_core_manuscript_candidate_tiers.tsv`
- `Project/results/primary_core_manuscript_candidate_tiers.md`

## Mechanism Axis Model

The manuscript candidate tiers have been organized into mechanism axes for figure planning and Results drafting.

- Developmental regulatory control: 7 genes, including Tier 1 `NFIB`, `NFIA`, and `RFX3`.
- Neurite/cytoskeleton morphogenesis: 13 genes, including Tier 1 `GPM6A`.
- Axon guidance and adhesion: 4 genes, including Tier 2 `ROBO2`.
- Synaptic/excitability maturation: 12 genes, including Tier 1 `KCNK1` and `GABRA2`, plus eight Tier 2 genes.
- Exploratory ortholog completeness: 30 non-identical-symbol ortholog hits kept outside central claims.

Main new outputs:

- `Project/results/primary_core_mechanism_axis_gene_table.tsv`
- `Project/results/primary_core_mechanism_axis_summary.tsv`
- `Project/results/primary_core_mechanism_axis_branch_summary.tsv`
- `Project/results/primary_core_mechanism_axis_model.png`
- `Project/results/primary_core_mechanism_axis_model.md`

## Named-Comparator Specificity Audit

The current mechanism axes have been tested against explicit named non-granule comparators in the local primary datasets.

- `GSE104323`: dentate granule-lineage groups versus `CA3-Pyr` and `Immature-Pyr`.
- `GSE122357`: cerebellar granule-lineage groups versus `Purkinje cells`.
- Units scored: 192 group-axis units across 4 axes.
- Specificity result: 3 axes are not granule-specific versus the named comparators; developmental regulatory control is cerebellar granule-enriched versus Purkinje cells but not dentate granule-enriched versus pyramidal comparators.

Main new outputs:

- `Project/results/primary_core_granule_specificity_named_comparator_units.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_axis_summary.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_gene_coverage.tsv`
- `Project/results/primary_core_granule_specificity_named_comparators.png`
- `Project/results/primary_core_granule_specificity_named_comparators.md`

## Niche/Fate Versus Circuit/Morphology Model

The niche/circuit question now has a formal module layer.

- Modules tested: cerebellar fate/rhombic-lip/SHH, dentate fate/WNT/PROX1, shared neurogenic niche/progenitor state, downstream neurite/morphology, and downstream synaptic/excitability.
- Genes scored in the formal model: 92/92.
- Median formal convergence delta: -0.500 for upstream/niche modules and 0.500 for downstream circuit/morphology modules.
- Named-comparator result: cerebellar fate is cerebellar-granule enriched only, dentate fate is dentate-granule enriched only, and downstream modules are not uniquely granule-specific versus pyramidal/Purkinje comparators.

Main new outputs:

- `Project/results/primary_core_niche_circuit_module_gene_sets.tsv`
- `Project/results/primary_core_niche_circuit_module_formal_gene_scores.tsv`
- `Project/results/primary_core_niche_circuit_module_formal_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_named_comparator_units.tsv`
- `Project/results/primary_core_niche_circuit_module_named_comparator_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_model.png`
- `Project/results/primary_core_niche_circuit_module_model.md`

## Transcriptomic Configuration Model

The module-balance assembly-plan test is now complete.

- Configuration units: 48 local source-group units.
- Named granule-versus-comparator contrasts: 4.
- Combined configuration score is positive in 4/4 contrasts.
- Median combined configuration scores: dentate granule 1.052, pyramidal comparator 0.938, cerebellar granule 0.406, Purkinje comparator -0.312.
- Small-n Wilcoxon p for combined configuration greater than comparator: 0.0625.

Main new outputs:

- `Project/results/primary_core_transcriptomic_configuration_units.tsv`
- `Project/results/primary_core_transcriptomic_configuration_role_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_model.png`
- `Project/results/primary_core_transcriptomic_configuration_model.md`

## Primary-Core Configuration Validation

The transcriptomic assembly-plan score has been broadened from named local contrasts to primary-core pseudobulk candidate-background contrasts.

- Configuration class units: 210 across 10 datasets.
- Candidate-versus-background contrasts: 63 across 7 datasets.
- Combined configuration positive contrasts: 52/63.
- Median candidate-background configuration delta: 0.417.
- Sign-test p: 8.37e-08.
- Wilcoxon p: 4.89e-08.

Main new outputs:

- `Project/results/primary_core_transcriptomic_configuration_primary_units.tsv.gz`
- `Project/results/primary_core_transcriptomic_configuration_primary_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_coverage.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.png`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`

## Configuration Driver Audit

The assembly-plan score has been audited to separate downstream construction balance from regional fate polarity.

- Total contrasts audited: 67.
- Configuration-positive contrasts: 56/67.
- Both components positive: 28.
- Fate-driven positives: 27.
- Construction-driven positives: 1.
- Primary-core layer: 52/63 configuration-positive, split into 26 both-component positives and 26 fate-driven positives.

Main new outputs:

- `Project/results/primary_core_configuration_driver_audit_contrasts.tsv`
- `Project/results/primary_core_configuration_driver_audit_module_deltas.tsv`
- `Project/results/primary_core_configuration_driver_audit_summary.tsv`
- `Project/results/primary_core_configuration_driver_audit_gene_priorities.tsv`
- `Project/results/primary_core_configuration_driver_audit.png`
- `Project/results/primary_core_configuration_driver_audit.md`

## Aim 2b Stage-Resolved TGF-beta/BDNF Audit

The historical TGF-beta/BDNF maturation/stop mechanism has now been tested as a stage/state-resolved signal rather than only as a pooled pathway-readiness contrast.

- Stage/state gene units: 3,305.
- Pathway-stage units: 261.
- Signature-stage units: 116.
- Transition rows: 112.
- Main interpretation: TGF-beta/BDNF behaves as a timing- and state-sensitive maturation-readiness axis, not as a universal region-wide granule-cell stop switch.
- Cerebellar candidate granule cells show a P0-to-P8 increase in TGF-beta/BDNF/stop balance, whereas dentate data show strong progenitor/postnatal and activity-state dependence with lower adult mature GC scores.

Main new outputs:

- `Project/results/primary_core_aim2b_stage_tgf_bdnf_gene_units.tsv`
- `Project/results/primary_core_aim2b_stage_tgf_bdnf_pathway_units.tsv`
- `Project/results/primary_core_aim2b_stage_tgf_bdnf_signature_units.tsv`
- `Project/results/primary_core_aim2b_stage_tgf_bdnf_transitions.tsv`
- `Project/results/primary_core_aim2b_stage_tgf_bdnf_plot.png`
- `Project/results/primary_core_aim2b_stage_tgf_bdnf.md`

## Cerebellar Conditioned-Medium Secretome Candidate Screen

The 2005 conditioned-medium question now has a sequencing-derived candidate-factor list. This screen scores curated secreted/ligand genes across primary-core cerebellar datasets for cerebellar candidate granule-cell source plausibility, candidate-background enrichment, and P0-to-P8 developmental timing.

- Secretome expression units: 1,113.
- Candidate-background contrasts: 371.
- Ranked candidates: 68.
- Highest-priority inferred factors besides the historical anchors: `SFRP1`, `BMP6`, `RELN`, and `SFRP2`.
- Supported additional classes: BMP/GDF/activin-family factors, secreted guidance/migration-stop cues (`SEMA3A`, `SLIT2`), WNT antagonists, and matricellular/context factors.
- Key caution: this nominates possible conditioned-medium components but does not prove protein secretion, processing, concentration, or anti-proliferative bioactivity.

Main new outputs:

- `Project/results/cerebellar_conditioned_medium_secretome_units.tsv`
- `Project/results/cerebellar_conditioned_medium_secretome_candidate_contrasts.tsv`
- `Project/results/cerebellar_conditioned_medium_secretome_ranked_candidates.tsv`
- `Project/results/cerebellar_conditioned_medium_secretome_gse122357_stage.tsv`
- `Project/results/cerebellar_conditioned_medium_secretome_candidates.png`
- `Project/results/cerebellar_conditioned_medium_secretome_candidates.md`

## 2005 Paper Endpoint Trajectory/Pseudotime Audit

The remaining 2005 paper readouts have now been mapped onto stage and pseudotime-style RNA modules.

- Stage/module gene units: 2,709.
- Stage/module units: 261.
- Stage trajectory rows: 29.
- Stage correlation rows: 60.
- Cell-level graph pseudotime cells: 19,192 across `GSE104323`, `GSE292261`, and `GSE122357`.
- Graph pseudotime module correlation rows: 27.
- Main interpretation: the paper endpoints are stage/trajectory dependent, but not all are monotonic with chronological age.
- Adult dentate lineage is the cleanest trajectory: neuronal maturation rises along graph pseudotime; proliferation is trajectory-windowed rather than simply early versus late.
- Cerebellar candidates show a P0-rooted graph trajectory with TGF/SMAD and BDNF/ERK increasing along graph pseudotime, while P8a/P8b differ enough that age should not be treated as pseudotime.

Main new outputs:

- `Project/results/primary_core_2005_endpoint_pseudotime_gene_units.tsv`
- `Project/results/primary_core_2005_endpoint_pseudotime_module_units.tsv`
- `Project/results/primary_core_2005_endpoint_pseudotime_trajectory_scores.tsv`
- `Project/results/primary_core_2005_endpoint_pseudotime_correlations.tsv`
- `Project/results/primary_core_2005_endpoint_pseudotime_audit.png`
- `Project/results/primary_core_2005_endpoint_pseudotime_audit.md`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime_cell_scores.tsv.gz`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime_group_summary.tsv`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime_bin_summary.tsv`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime_module_correlations.tsv`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime.png`
- `Project/results/primary_core_2005_endpoint_graph_pseudotime.md`

## Full Primary-Core Diffusion/Pseudotime Upgrade

The endpoint-gene trajectory audit has now been upgraded to a primary-core highly-variable-gene diffusion/pseudotime layer.

- Primary datasets processed: 10/10.
- Cells/nuclei in trajectory tables: 56,892 after deterministic stratified sampling.
- Strict local full-transcriptome or full-DG-subset datasets: 8.
- Selected-feature bridge datasets: 2 (`GSE325391`, `GSE268609`).
- Ordered-stage validation: `GSE104323` rho 0.797, `GSE292261` rho 0.751, `GSE122357` rho 0.750, and `GSE214309` state-axis rho 0.684.
- Source-pseudotime validation: `GSE325391` diffusion pseudotime versus source sling pseudotime rho 0.579.
- Main interpretation: the full trajectory layer does not overturn Fig1-5, but it makes the model stage-aware. Fig5 should be revised to include pseudotime as an explicit layer between regional fate/niche and terminal morphology.

Main new outputs:

- `Project/results/primary_core_full_transcriptome_diffusion_cell_scores.tsv.gz`
- `Project/results/primary_core_full_transcriptome_diffusion_dataset_summary.tsv`
- `Project/results/primary_core_full_transcriptome_diffusion_group_summary.tsv`
- `Project/results/primary_core_full_transcriptome_diffusion_module_correlations.tsv`
- `Project/results/primary_core_full_transcriptome_diffusion_fig1_5_impact.tsv`
- `Project/results/primary_core_full_transcriptome_diffusion_overview.png`
- `Project/results/primary_core_full_transcriptome_diffusion.md`

## Manuscript Planning Packet

The first manuscript scaffold is now complete.

- Main packet: `Project/results/manuscript_planning_packet.md`
- Claim/evidence/caveat guardrail table: `Project/results/manuscript_claim_evidence_caveat_table.tsv`
- Figure plan table: `Project/results/manuscript_figure_plan.tsv`
- Recommended central language: distinct regional fate programs converge on an identity-coupled transcriptomic assembly configuration.
- Main caution: do not frame the current result as a unique granule-specific pathway or pure morphology-only transcriptomic code.

## Interpretation

This corrects the first major gap in the previous core design: the project no longer has only a human cerebellum dataset without a human dentate/hippocampal construction branch. The human branch now includes the lighter scaffold (`GSE185277`, `GSE185553`, `GSE186538`), a direct adult human dentate anchor (`GSE325391`), and a broader aging/AD hippocampal RNA expansion (`GSE268609`), with QC/metadata, tuned or projected labels, selected-feature sparse objects, first-pass dataset-aware module tests, an integrated rank-level comparison against the existing dentate/cerebellar backbone, and a primary-core trajectory layer.

## Next Build Step

The revised manuscript figure assembly is now complete in `Project/results/manuscript_figures/`, with six composite figures including a stage-windowed Fig5 and an integrated Fig6 working model. Aim 2 now has pooled pathway-readiness, stage-resolved TGF-beta/BDNF, conditioned-medium secretome-prioritization, endpoint-gene pseudotime, and full primary-core diffusion layers. The next build step should move from figure assembly into manuscript text drafting and optional external validation. The named-comparator audit, niche/circuit model, transcriptomic configuration model, primary-core configuration validation, driver audit, Aim 2b stage audit, conditioned-medium secretome screen, endpoint pseudotime analyses, full diffusion layer, and Fig6 working model should be used to phrase the central model as distinct upstream fate programs converging on a partly shared downstream neuronal morphogenesis, excitability, and maturation toolkit whose identity-coupled configuration is stage-windowed. For `GSE268609` and `GSE325391`, source full-object export remains the key practical refinement if the manuscript requires raw-count or genome-wide human dentate/hippocampal trajectories beyond the current selected-gene bridge.
