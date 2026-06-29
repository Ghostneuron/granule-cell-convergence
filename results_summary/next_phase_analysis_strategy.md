# Next-phase analysis strategy

## Decision

The project should proceed as an original in silico analysis, with a review/synthesis paper as a fallback or companion format.

The current evidence supports a focused and novel thesis:

Dentate gyrus granule cells and cerebellar granule cells retain distinct regional identity programs, but partially converge on downstream structural, morphogenesis, axon-guidance, and synaptic programs. This convergence may help explain why two developmentally and circuit-distinct neuron classes share a compact granule-cell morphology.

## Current local core datasets

These are already local and remain valuable, but they should now be treated as the local starting core rather than the complete analysis core:

- `GSE104323`: mouse dentate gyrus 10x, annotated dentate granule lineage plus non-dentate controls. Best dentate reference and control dataset.
- `GSE95752`: mouse dentate gyrus C1 single-cell data. Useful for dentate maturation and small-cell validation.
- `GSE292261`: mouse postnatal dentate gyrus Smart-seq2. Promoted to primary dentate developmental validation after curation; gene symbols are already present and the metadata now preserve postnatal stage, cluster labels, and QC metrics.
- `GSE214309`: mouse adult/activity-state dentate snRNA-seq. Promoted to primary dentate maturation/activity validation after curation; the remaining technical requirement is full Ensembl-to-symbol mapping for whole-transcriptome object-level analysis.
- `GSE122357`: mouse cerebellar postnatal single-cell data. Best cerebellar developmental comparison currently in the local data.
- `GSE165657`: human cerebellum 10x aggregate. Strong large human cerebellar validation dataset.
- `GSE312658`: mouse cerebellum control and cKO 10x. Useful for cerebellar validation and perturbation sensitivity.

The strict 10-dataset integrated primary core is now defined separately in:

- `Project/results/integrated_primary_core_datasets.tsv`
- `Project/results/integrated_primary_core_datasets.md`

## Constructed human dentate/hippocampal core extension

The earlier local core had a human cerebellum anchor but no constructed human dentate or hippocampal counterpart. That asymmetry weakened the cross-species and human relevance of the project. The human dentate/hippocampal branch is now constructed through sparse objects, QC, marker triage, GEO metadata curation, a normalized reduced object, tuned labels, first-pass dataset-aware module tests for `GSE185277`, `GSE185553`, and `GSE186538`, and a direct adult DG anchor from `GSE325391`.

Recommended construction order and current status:

1. `GSE185277`: Nature 2022 human hippocampal/dentate immature-neuron resource. Its RAW archive has been downloaded, converted into 7 sparse per-library objects, QC-harmonized, marker-scored, mapped to GEO specimen/age metadata, included in the normalized reduced object, label-tuned, and module-tested.
2. `GSE185553`: broader human hippocampal companion reference for separating imGC/DG markers from general hippocampal neuronal states. Its RAW archive has been downloaded, converted into 27 sparse per-library objects, QC-harmonized, marker-scored, mapped to pooled GEO specimen/age metadata, included in the normalized reduced object, label-tuned, and module-tested.
3. `GSE186538`: human DG taxonomy and cross-species reference; human counts, metadata, and genes are downloaded, and 32,067 DG granule-cell candidates are extracted as a sparse, marker-validated, donor-curated subset and included as the tuned normalized DG anchor.
4. `GSE325391`: primary modern adult human dentate acquisition; adult DG RDS is downloaded, inspected, converted to a selected sparse bridge, and mapped into the tuned-label convention.
5. `GSE268609`: human hippocampal aging/AD multiome expansion; RNA matrix/barcodes/features are downloaded, converted into a selected sparse bridge, projected into the tuned-label convention, and first-pass diagnosis module tests are complete. ATAC fragments and the full Seurat object remain deferred.

`GSE233363` should become a strong mouse DG aging/neurogenic-lineage validation dataset after this human branch is present. It does not solve the missing human dentate comparator by itself.

## Supporting datasets

Use these as validation or context, not main discovery:

- `GSE242688`: cerebellum Visium/spatial or proteomics-linked matrix. Useful for spatial support, but spots are not single cells.
- `GSE214905`: dentate patch-seq. Useful because it links transcription to physiology, but the sample is small.
- `GSE216877`: human epilepsy hippocampus context; useful as disease validation, not as the only human primary dataset.
- `GSE317381`: human DG spatial/annotation-transfer context; useful as spatial support, not cell-level discovery.
- `spatial_DG_lifespan`: Cell Reports 2025 human dentate gyrus spatial transcriptomics, already present as a local PDF. Use as human spatial/lifespan context; processed SpatialExperiment objects are available through Zenodo `10.5281/zenodo.10126687`, with raw data at `10.5281/zenodo.10126688` and code at `10.5281/zenodo.10126715`.
- `GSE198323`: Nature 2022 human hippocampal immature-neuron Alzheimer's disease context. Useful after the healthy/reference human imGC branch is defined.
- `GSE233363`: mouse DG aging/neurogenic-lineage scRNA-seq plus spatial transcriptomics. Use after human DG construction as validation.
- `GSE150153`: human organoid granule-like comparison. Keep as context because metadata were previously flagged as unmatched.

## Atlas references

- Allen Institute WMB taxonomy: use as marker and annotation reference, especially the `DG Glut`, `DG-PIR Ex IMN`, and `CB Granule Glut` subclasses. This should guide marker refinement and figure labeling, not be treated as a local expression matrix unless expression-level files are added.

## Current integrated evidence layer

The human dentate/hippocampal bridge has now been merged with the refined dentate/cerebellar backbone at a rank-summary level. Because the older backbone and newer human selected-gene bridges were not produced from one shared full-transcriptome matrix, the integrated analysis uses within-sample rank metrics.

- Integrated table: 635 summarized units; 546 units pass the `>=20` cells/spots threshold for statistics and plotting.
- Regional identity separation is strong: dentate candidates have higher dentate-minus-cerebellar identity rank than cerebellar candidates (median delta 0.6472, BH-adjusted p=1.57e-10) and non-dentate/background units (median delta 0.3423, BH-adjusted p=9.93e-22).
- Structural-program rank behaves as a convergence axis: dentate and cerebellar candidates are both above the within-sample structural median, while their direct structural-rank difference is not significant.
- Human primary anchors have different roles: `GSE325391` and `GSE186538` are source-aware DG anchors, while `GSE268609` supplies a broader hippocampal aging/AD context with projected labels.

Main outputs:

- `Project/results/human_bridge_backbone_rank_units.tsv`
- `Project/results/human_bridge_backbone_rank_statistics.tsv`
- `Project/results/human_bridge_backbone_rank_source_summary.tsv`
- `Project/results/human_bridge_backbone_rank_units.png`
- `Project/results/human_bridge_backbone_rank_integration.md`

## First candidate-gene packet

A first manuscript-planning gene packet has been built from the refined marker panels across the constructed human bridge. This packet is useful for framing the mechanistic model, but it should be treated as a hypothesis shortlist until full ortholog-aware differential expression is run.

- Gene scope: 67 refined panel genes.
- High-priority shared structural/executor examples: `NCAM1`, `MAPT`, `PLXNA4`, `DPYSL2`, `CNTN5`, `ROBO2`, `NRP1`, `EPHA4`, `ROBO1`, `PLXNA2`, `SLIT2`, `STMN2`, and `GAP43`.
- Human-bridge dentate identity examples: `PROX1`, `CALB1`, `C1QL3`, `GLIS3`, `EGR3`, and `ITPKA`.
- Specificity warnings: `ROR1`, `ETV1`, and `GABRA6` show detectable human dentate-bridge signal and should not be treated as clean cerebellar-only markers without context.

Main outputs:

- `Project/results/human_bridge_marker_gene_summary.tsv`
- `Project/results/human_bridge_candidate_gene_packet.tsv`
- `Project/results/human_bridge_structural_executor_candidates.png`
- `Project/results/human_bridge_candidate_gene_packet.md`

## First formal primary-core module result

The strict 10-dataset primary core now has a marker-panel ortholog-aware module analysis.

- Core size: 10 datasets, 418 integrated units, 361 eligible units.
- Identity separation is robust: dentate candidates versus cerebellar candidates have identity-rank delta 0.5356, BH-adjusted p=2.95e-08.
- Leave-one-dataset-out identity separation remains positive and BH-significant.
- Both dentate and cerebellar candidates are above the structural median, but cerebellar candidates are higher in magnitude; phrase this as shared elevated structural/morphogenesis executor activity, not identical structural state.

Main outputs:

- `Project/results/primary_core_marker_panel_ortholog_map.tsv`
- `Project/results/primary_core_integrated_module_units.tsv`
- `Project/results/primary_core_ortholog_module_statistics.tsv`
- `Project/results/primary_core_ortholog_module_leave_one_dataset_out.tsv`
- `Project/results/primary_core_identity_structural_modules.png`
- `Project/results/primary_core_ortholog_module_analysis.md`

## Candidate-gene pseudobulk screen

The strict 10-dataset primary core now has a focused candidate-gene pseudobulk screen across the 67 human-bridge candidate genes.

- Coverage: 10/10 primary datasets represented by at least one candidate gene.
- Scale: 16,719 pseudobulk expression rows and 67 candidate-gene statistics.
- Structural-executor support: 14 structural-executor genes have positive dentate and cerebellar candidate rank deltas; 12 pass the exploratory BH<0.2 cutoff in both branches.
- Strongest shared structural-executor signals: `CFL1`, `GAP43`, `ROBO2`, `STMN2`, `STMN3`, `CDK5R1`, `DPYSL2`, `EPHB2`, `L1CAM`, `MAPT`, `DPYSL3`, and `ELAVL4`.
- Interpretation: the mechanism shortlist is now more than a human-bridge marker packet; it has primary-core pseudobulk support. It still needs genome-wide ortholog-aware DE before being treated as the final mechanism model.

Main outputs:

- `Project/results/primary_core_candidate_gene_pseudobulk_expression.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_hits.tsv`
- `Project/results/primary_core_candidate_gene_pseudobulk_effects.png`
- `Project/results/primary_core_candidate_gene_pseudobulk_analysis.md`

## Expanded selected-gene discovery screen

The strict 10-dataset primary core now also has an expanded pseudobulk screen over the 2,169-gene human-core selected-feature universe.

- Coverage: 10/10 primary datasets represented.
- Scale: 537,070 pseudobulk expression rows and 2,169 gene statistics.
- Shared-positive rank genes: 1,279.
- Shared-positive genes passing BH<0.10 in both branches: 579.
- Original 67-gene packet genes recovered among shared-positive hits: 42.
- Mechanism triage: 71 figure-level candidates and 76 follow-up candidates.
- Strongest recovered structural genes: `STMN2`, `GAP43`, `ROBO2`, `CFL1`, and `STMN3`.
- Strong new candidate mechanisms: `GPM6A`, `MAP1B`, `TCF4`, `NFIA`, `NFIB`, `PPP3CA`, `CALM1`, `RTN3`, `RTN1`, `NRXN1`, `BASP1`, `CAMTA1`, `CADM1`, and `DCC`.

Main outputs:

- `Project/results/primary_core_expanded_gene_pseudobulk_expression.tsv.gz`
- `Project/results/primary_core_expanded_gene_pseudobulk_coverage.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_statistics.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_branch_specific.tsv`
- `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.png`
- `Project/results/primary_core_expanded_gene_pseudobulk_analysis.md`
- `Project/results/primary_core_expanded_gene_mechanism_triage.tsv`
- `Project/results/primary_core_expanded_gene_mechanism_triage.md`

## Full-matrix same-symbol screen and consensus

The strict primary core now has a full-matrix same-symbol pseudobulk screen for the local full matrices, plus a cross-screen consensus table.

- Full-symbol universe: 33,939 genes from the GSE186538 human DG taxonomy gene list.
- Scale: 628,339 pseudobulk expression rows and 21,253 same-symbol gene statistics.
- Full-matrix expression coverage: 8/10 primary datasets.
- Rank-contrast datasets: 7 datasets (`GSE104323`, `GSE122357`, `GSE165657`, `GSE214309`, `GSE292261`, `GSE312658`, `GSE95752`).
- Shared-positive same-symbol genes: 6,440.
- Shared-positive genes passing BH<0.10 in both branches: 283.
- Genome-wide mechanism triage: 33 figure-level candidates and 261 follow-up candidates.
- Cross-screen consensus: 24 figure-level mechanism candidates survive both the selected-feature and full-matrix same-symbol screens.
- Top consensus candidates: `GPM6A`, `NFIA`, `NFIB`, `PPP3CA`, `CAMTA1`, `MAPK1`, `STXBP1`, `CALM2`, `CACNA2D1`, `SYNPR`, `KCNK1`, `GABRB3`, `ADD2`, `KCNJ6`, `RFX3`, `KCNJ3`, `GRIN2B`, `MAP3K4`, `KCND3`, `KCND2`, `GABRA2`, `STXBP5L`, `CACNA1E`, and `STXBP5`.

Main outputs:

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

## Dataset-aware consensus validation

The 24 cross-screen figure candidates have now been re-tested at the dataset/sample level to ask whether their candidate-versus-background deltas are consistently positive, rather than driven by pooled summaries.

- Dataset/sample/gene branch-delta rows: 1,608.
- Summary rows: 96, covering 24 genes across selected-feature/full-matrix screens and dentate/cerebellar branches.
- Robust across all available screen/branch tests: `GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, and `RFX3`.
- Robust in 3/4 screen/branch tests: `CACNA2D1`, `GABRB3`, `GRIN2B`, `KCND2`, `KCNJ3`, `KCNJ6`, `PPP3CA`, and `STXBP5L`.
- Interpretation: the 24-gene cross-screen consensus is useful for context, but the six all-available robust genes are the current strongest figure-level mechanism candidates.

Main outputs:

- `Project/results/primary_core_consensus_candidate_dataset_deltas.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation.tsv`
- `Project/results/primary_core_consensus_candidate_dataset_validation_heatmap.png`
- `Project/results/primary_core_consensus_candidate_dataset_validation.md`

## MGI-filtered ortholog meta-model

The strict primary core now has a conservative MGI-filtered ortholog meta-model layered on top of the selected-feature and full-matrix pseudobulk screens.

- Official MGI human-mouse homology report downloaded locally: `External_Data/Orthology/HOM_MouseHumanSequence.rpt`.
- Strict same-symbol one-to-one human-mouse pairs: 16,245.
- Unit delta rows: 241,690.
- Gene summary rows: 15,345.
- Shared hits after ortholog filtering: 1,304.
- Mechanism-prioritized shared hits: 36.
- Strongest dataset-robust consensus mechanism genes retained: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- Important second-tier consensus/structural genes include `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`, and `ROBO2`.

Main outputs:

- `Project/results/primary_core_mgi_ortholog_meta_model_map.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_meta_model_branch_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_meta_model_top_hits.png`
- `Project/results/primary_core_mgi_ortholog_meta_model.md`

## Expanded MGI ortholog full-matrix screen

The MGI model has now been expanded from same-symbol one-to-one pairs to all MGI one-to-one human-mouse pairs by resolving mouse source rows through the MGI mouse symbol.

- MGI one-to-one target genes: 17,611.
- Non-identical human/mouse symbol targets: 1,366.
- Full-matrix ortholog pseudobulk rows: 508,917.
- Genes tested in full-matrix ortholog statistics: 16,704.
- Shared-positive ortholog genes: 6,413.
- Shared-positive non-identical-symbol genes: 246.
- Expanded meta-model shared hits: 1,370.
- Expanded meta-model non-identical-symbol shared hits: 64.
- Mechanism-prioritized shared hits: 36.
- The six dataset-robust consensus genes remain unchanged: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.

Main outputs:

- `Project/results/primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_full_matrix_statistics.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_full_matrix_analysis.md`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_expanded_meta_model.md`

## Formal MGI ortholog rank-meta validation

The expanded MGI one-to-one ortholog model now has a formal dataset-level rank-meta validation layer. This is the current strongest statistical evidence layer for manuscript gene claims, while still remaining rank/pseudobulk meta-analysis rather than raw-count DESeq2/edgeR.

- Unit delta rows: 252,469.
- Dataset-level delta rows: 116,013.
- Branch tests: 36,303.
- Gene summaries: 16,708.
- Formal shared hits: 1,370.
- Both-screen formal shared hits: 158.
- FDR10-supported branches: 1,102, but no gene passes shared dentate-plus-cerebellar FDR10 in this small-dataset setting.
- Mechanism-prioritized genes modeled: 36.
- The six-gene consensus remains the safest seed set: `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`.
- Second-tier formal both-screen mechanism candidates include `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`, and `ROBO2`.

Main outputs:

- `Project/results/primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz`
- `Project/results/primary_core_mgi_ortholog_formal_rank_branch_tests.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_gene_summary.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_shared_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_model_long.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv`
- `Project/results/primary_core_mgi_ortholog_formal_rank_mechanism_hits.png`
- `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`
- `Project/results/primary_core_manuscript_candidate_tiers.tsv`
- `Project/results/primary_core_manuscript_candidate_tiers.md`

## Mechanism axis model

The manuscript candidate tiers have now been translated into biological mechanism axes:

- Developmental regulatory control: 7 genes, including Tier 1 `NFIB`, `NFIA`, and `RFX3`.
- Neurite/cytoskeleton morphogenesis: 13 genes, including Tier 1 `GPM6A`.
- Axon guidance and adhesion: 4 genes, including Tier 2 `ROBO2`.
- Synaptic/excitability maturation: 12 genes, including Tier 1 `KCNK1` and `GABRA2`, plus eight Tier 2 genes.
- Exploratory ortholog completeness: 30 non-identical-symbol ortholog hits kept outside central claims.

Main outputs:

- `Project/results/primary_core_mechanism_axis_gene_table.tsv`
- `Project/results/primary_core_mechanism_axis_summary.tsv`
- `Project/results/primary_core_mechanism_axis_branch_summary.tsv`
- `Project/results/primary_core_mechanism_axis_model.png`
- `Project/results/primary_core_mechanism_axis_model.md`

## Named-comparator specificity audit

The first direct named-comparator test is now complete. It scores the Tier 1-4 mechanism axes in `GSE104323` dentate granule-lineage groups versus pyramidal comparators (`CA3-Pyr`, `Immature-Pyr`) and in `GSE122357` cerebellar granule-lineage groups versus Purkinje cells.

- Units scored: 192 group-axis units across 4 mechanism axes.
- Developmental regulatory control is higher in cerebellar granule-lineage groups than Purkinje cells, but not higher in dentate granule-lineage groups than pyramidal comparators.
- Neurite/cytoskeleton morphogenesis, axon guidance/adhesion, and synaptic/excitability maturation do not pass the strict named-comparator specificity rule.
- Interpretation: the shared axes should be described as convergent neuronal morphology/maturation programs, not as pathways unique to granule cells. Specificity may reside in gene combinations, developmental timing, regional regulatory context, and cell-shape implementation.

Main outputs:

- `Project/results/primary_core_granule_specificity_named_comparator_units.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_axis_summary.tsv`
- `Project/results/primary_core_granule_specificity_named_comparator_gene_coverage.tsv`
- `Project/results/primary_core_granule_specificity_named_comparators.png`
- `Project/results/primary_core_granule_specificity_named_comparators.md`

## Niche/fate versus circuit/morphology model

The niche/circuit question has now been translated into five hypothesis modules: cerebellar fate/rhombic-lip/SHH, dentate fate/WNT/PROX1, shared neurogenic niche/progenitor state, downstream neurite/morphology, and downstream synaptic/excitability.

- Genes scored: 92/92 in the formal MGI rank model.
- Median formal convergence delta for upstream/niche modules: -0.500.
- Median formal convergence delta for downstream circuit/morphology modules: 0.500.
- Downstream greater than upstream/niche by Mann-Whitney test: p=1.56e-09.
- Named-comparator layer: cerebellar fate is cerebellar-granule enriched only; dentate fate is dentate-granule enriched only; downstream modules are convergent in the formal model but not uniquely granule-specific versus pyramidal/Purkinje comparators.

Main outputs:

- `Project/results/primary_core_niche_circuit_module_gene_sets.tsv`
- `Project/results/primary_core_niche_circuit_module_formal_gene_scores.tsv`
- `Project/results/primary_core_niche_circuit_module_formal_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_named_comparator_units.tsv`
- `Project/results/primary_core_niche_circuit_module_named_comparator_summary.tsv`
- `Project/results/primary_core_niche_circuit_module_model.png`
- `Project/results/primary_core_niche_circuit_module_model.md`

## Transcriptomic configuration model

The assembly-plan idea has now been tested directly as a module-balance score rather than as single-module specificity.

- Downstream construction balance = mean(neurite/morphology, synaptic/excitability) minus neurogenic niche/progenitor rank.
- Regional fate balance = branch-matched fate rank minus branch-opposed fate rank.
- Combined configuration score = construction balance plus regional fate balance.
- Named granule-versus-comparator contrasts tested: 4.
- Positive combined configuration score: 4/4 contrasts.
- Wilcoxon p for combined configuration greater than comparator: 0.0625; treat as promising small-n support, not final proof.

Main outputs:

- `Project/results/primary_core_transcriptomic_configuration_units.tsv`
- `Project/results/primary_core_transcriptomic_configuration_role_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_model.png`
- `Project/results/primary_core_transcriptomic_configuration_model.md`

## Primary-core configuration validation

The transcriptomic assembly-plan score has now been validated across broader primary-core pseudobulk expression layers.

- Configuration class units: 210 across 10 datasets.
- Candidate-versus-background contrasts: 63 across 7 datasets.
- Expression layers: full MGI one-to-one ortholog matrix and selected-feature matrix.
- Positive combined configuration contrasts: 52/63.
- Median candidate-background combined configuration delta: 0.417.
- Sign-test p: 8.37e-08.
- Wilcoxon p: 4.89e-08.
- Caveat: this broad layer uses local background classes, not explicit pyramidal/Purkinje labels in every dataset.

Main outputs:

- `Project/results/primary_core_transcriptomic_configuration_primary_units.tsv.gz`
- `Project/results/primary_core_transcriptomic_configuration_primary_contrasts.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_summary.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_coverage.tsv`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.png`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`

## Configuration driver audit

The configuration score has now been decomposed into its two components: downstream construction balance and regional fate polarity.

- Total contrasts audited: 67, including 63 primary candidate-background contrasts and 4 local named-comparator contrasts.
- Configuration-positive contrasts: 56/67.
- Both construction and fate components positive: 28.
- Fate-driven positives: 27.
- Construction-driven positives: 1.
- Primary-core layer: 52/63 configuration-positive, with 26 both-component positives and 26 fate-driven positives.
- Interpretation: the assembly-plan evidence is strong but identity-coupled; downstream construction balance is selective, not a pure morphology-only code.

Main outputs:

- `Project/results/primary_core_configuration_driver_audit_contrasts.tsv`
- `Project/results/primary_core_configuration_driver_audit_module_deltas.tsv`
- `Project/results/primary_core_configuration_driver_audit_summary.tsv`
- `Project/results/primary_core_configuration_driver_audit_gene_priorities.tsv`
- `Project/results/primary_core_configuration_driver_audit.png`
- `Project/results/primary_core_configuration_driver_audit.md`

## Analysis work packages

1. Object building and QC
   - Convert the core matrices into AnnData, Seurat, or SingleCellExperiment objects.
   - Preserve raw counts, source metadata, sample labels, species, platform, region, and candidate-call annotations.
   - Normalize within dataset first; avoid direct cross-dataset expression comparisons until ortholog mapping and batch handling are explicit.

2. Annotation curation
   - Curate dentate granule, immature dentate granule, non-dentate excitatory, inhibitory, glial, cerebellar granule, and ambiguous classes.
   - Use Allen WMB subclass markers as an external reference.
   - Treat the `cerebellum_warning` class as a diagnostic group for genes that are broad neurogenic or structural markers rather than strict dentate identity markers.

3. Ortholog and gene-symbol harmonization
   - Build a mouse-human-rat ortholog table before any cross-species score comparison.
   - Separate strict identity genes from shared developmental and structural genes.
   - Remove broad genes from strict identity panels if they also mark immature neurons or general excitatory neurons.

4. Formal module scoring
   - Compute normalized module scores within each dataset using control-gene matched scoring or rank-based scoring.
   - Retain these modules: strict dentate identity, strict cerebellar granule identity, shared granule/excitatory neuronal state, morphogenesis/cytoskeleton, axon guidance/synapse.
   - Use within-sample ranks for cross-dataset summary plots.

5. Statistical testing
   - Test identity separation with pseudobulk or mixed-effect models using dataset/sample as the replicate structure.
   - Test structural convergence by asking whether dentate and cerebellar granule candidates are both enriched for structural modules, while not requiring them to be identical.
   - Test developmental timing: convergence should be strongest during maturation, migration, axon extension, dendrite remodeling, and synapse formation.

6. Regulatory-network layer
   - Use the curated cells to separate upstream identity regulators from downstream morphology executors.
   - Expected identity regulators include regionally biased transcriptional programs, while expected shared executors include cytoskeletal, neurite, axon-guidance, and synaptic genes.

## Manuscript figure plan

1. Figure 1: Biological question and historical context from the thesis/paper. Similar morphology, distinct anatomical/circuit origins.
2. Figure 2: Dataset map and candidate-call workflow. Include the next-phase priority table and Allen reference classes.
3. Figure 3: Regional identity separation. Dentate candidates are dentate-high/cerebellar-low; cerebellar candidates are cerebellar-high/dentate-low.
4. Figure 4: Structural-program convergence. Both candidate classes show elevated structural/morphogenesis rank, but the module is not granule-specific.
5. Figure 5: Developmental timing and regulatory model. Similar morphology emerges from shared downstream executors downstream of distinct regional identity programs.

## Go/no-go criteria

Proceed as an original analysis paper if the object-level analysis confirms:

- Dentate and cerebellar candidates remain identity-separated after dataset/sample-aware statistics.
- Both groups show enriched structural or wiring modules within their own datasets.
- The convergence signal is not fully explained by generic immature-neuron status.
- At least one mouse-mouse comparison and one cross-species or atlas-validation comparison support the model.

If these criteria weaken, the project is still valuable as a review with a smaller computational synthesis. In that version, the central contribution would be a conceptual framework: "granule-cell morphology as convergent deployment of a shared neuronal morphogenesis toolkit."

## Immediate next step

The next computational step should move beyond rank-meta validation into manuscript packaging and external validation. The marker-panel module result, 67-gene screen, 2,169-gene selected-feature discovery/triage, full-matrix same-symbol screen, cross-screen consensus, consensus candidate robustness layer, expanded MGI one-to-one meta-model, formal rank-meta validation, named-comparator specificity audit, niche/circuit module model, transcriptomic configuration model, primary-core configuration validation, and configuration driver audit are now strong enough to justify manuscript figures.

The revised manuscript figure assembly is now available in `Project/results/manuscript_figures/`, with six composite PNG/PDF figures plus `manuscript_figure_manifest.tsv` and `manuscript_figure_assembly.md`. Fig5 includes the stage-windowed pseudotime refinement, and Fig6 integrates the Fig1/2/4/5 flow charts into one working model. The current manuscript scaffold is `Project/results/manuscript_planning_packet.md`, with claim guardrails in `Project/results/manuscript_claim_evidence_caveat_table.tsv` and figure planning in `Project/results/manuscript_figure_plan.tsv`. A manuscript candidate-tier packet is available in `Project/results/primary_core_manuscript_candidate_tiers.md`, the current biological interpretation layer is `Project/results/primary_core_mechanism_axis_model.md`, the specificity caveat is documented in `Project/results/primary_core_granule_specificity_named_comparators.md`, the niche/circuit model is documented in `Project/results/primary_core_niche_circuit_module_model.md`, Aim 2 is documented in `Project/results/primary_core_aim2_niche_pathway_model.md`, Aim 3 is documented in `Project/results/primary_core_aim3_sparse_coding_model.md`, and the transcriptomic assembly-plan evidence is documented in `Project/results/primary_core_transcriptomic_configuration_model.md`, `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`, and `Project/results/primary_core_configuration_driver_audit.md`. The practical next outputs should be manuscript text drafting and targeted validation against external spatial, morphology-linked, aging, and disease resources if needed for journal strength. Raw-count/object-level DE remains useful if the paper needs claim-strength beyond rank-meta evidence, especially for `GSE325391`, `GSE268609`, and any added spatial or Allen-derived expression layer.
