# Manuscript Planning Packet

Date updated: 2026-06-24

## Purpose

This packet turns the completed primary-core evidence chain into a manuscript scaffold. It is meant to guide either an original in silico analysis paper or a computationally anchored review with new analyses. It now includes the completed fitted Aim 2 stage-window model, Aim 3 empirical sparse-coding calibration, NeuroMorpho/DANDI validation layer, and the 2005 paper support/revision table.

## Recommended Direction

Proceed. The project is reasonable, novel, and important if the central claim is framed precisely:

Dentate gyrus and cerebellar granule cells do not appear to share one simple granule-cell fate program, and the current data do not support a unique granule-specific pathway signature. Instead, distinct regional fate programs appear to converge on an identity-coupled transcriptomic assembly configuration that combines branch-matched fate polarity with downstream neurite, synaptic, excitability, and maturation machinery.

The safest manuscript type is an original in silico synthesis paper with a strong conceptual model. If later external validation is thinner than expected, it remains publishable as a review or perspective-style computational synthesis.

## Title Options

1. Distinct Granule-Cell Lineages Converge on an Identity-Coupled Transcriptomic Assembly Program
2. Transcriptomic Assembly Configuration Links Dentate and Cerebellar Granule-Cell Morphology
3. Shared Granule-Cell Morphology Reflects Convergent Downstream Assembly Rather Than Shared Fate Identity
4. Convergent Neuronal Construction Programs in Dentate and Cerebellar Granule Cells

Recommended current title: Distinct Granule-Cell Lineages Converge on an Identity-Coupled Transcriptomic Assembly Program.

## Central Claim

The cerebellar and dentate granule-cell similarity is unlikely to be chance, but it is also not explained by a single shared granule-cell identity program. Across a strict 10-dataset primary core, the strongest defensible model is that branch-specific upstream fate programs place each population into its local identity context, while partially shared downstream construction modules contribute to compact granule-cell morphology and circuit maturation. The transcriptomic signal is best described as an identity-coupled assembly configuration, not a pure morphology-only code.

## Why This Is Reasonable

- The strict primary core now includes mouse dentate, cerebellar, and human dentate/hippocampal branches.
- Ortholog-aware rank-meta analysis avoids overclaiming from one species or one dataset.
- The evidence chain includes positive support and claim-safety audits: formal rank model, candidate tiers, mechanism axes, named comparators, niche/circuit decomposition, primary-core configuration validation, driver audit, Aim 2 pathway-readiness analysis, fitted stage-window analysis, Aim 3 sparse-coding model, and empirical calibration against public morphology/activity data.
- The project answers a real biological puzzle from the thesis and Development 2005 work: why two anatomically distinct granule-cell populations look unexpectedly similar.

## What Is Novel

- It treats morphology similarity as a configuration problem rather than a single marker-gene or pathway-overlap problem.
- It explicitly separates upstream fate/niche modules from downstream neurite/morphology and synaptic/excitability modules.
- It includes a negative control logic: named pyramidal and Purkinje comparators temper the claim and make the model more credible.
- It integrates human dentate/hippocampal data into a question that is often addressed mainly with mouse developmental data.

## What Is Important

- The project could clarify how similar neuronal geometries emerge in different brain regions.
- It may provide a general framework for studying convergent neuronal morphologies: distinct fate programs, shared construction machinery, and local circuit constraints.
- It produces a practical prioritized gene set for future wet-lab or external in silico validation.

## Abstract Draft

Dentate gyrus granule cells and cerebellar granule cells share a compact excitatory-neuron design despite arising in different brain regions, developmental contexts, and circuit architectures. Earlier work from this laboratory showed that cerebellar conditioned medium can suppress hippocampal granule-cell proliferation and promote differentiation through TGF-beta2, BDNF, SMAD, and MAPK signaling, raising a broader question: does similar granule-cell morphology reflect shared molecular identity, convergent construction programs, niche signaling, or circuit-level constraint? To address this, we built a strict 10-dataset primary core spanning mouse dentate gyrus, mouse and human cerebellum, and human dentate/hippocampal datasets. We performed ortholog-aware rank-meta modeling, candidate-gene tiering, named-comparator specificity tests, niche/fate versus downstream-construction module analysis, transcriptomic configuration scoring, full-primary-core diffusion/pseudotime analysis, fitted stage-window modeling, and sparse expansion-coding simulation with empirical morphology/activity calibration. Dentate and cerebellar granule-cell candidates remained strongly regionally identity-separated, but downstream neurite/morphology and synaptic/excitability modules showed stronger convergence than upstream fate/niche modules. A formal MGI ortholog rank-meta model identified 1,370 shared hits, with a conservative six-gene seed set comprising `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, and `GABRA2`. Named-comparator analysis argued against a simple uniquely granule-specific pathway signature. Instead, candidate granule populations showed broad elevation of an identity-coupled transcriptomic assembly-configuration score. TGF-beta/BDNF/SMAD/MAPK programs behaved as stage-windowed maturation/readiness overlays rather than a universal cerebellar stop signal. NeuroMorpho and DANDI analyses provided first external morphology and activity support, while sparse expansion-coding simulation showed that granule-like architectures are computationally plausible under wiring and activity constraints. Together, these results support a model in which similar granule-cell morphology emerges from convergent downstream assembly machinery deployed within distinct regional identities and circuit niches.

## Results Spine

1. A 10-dataset strict primary core supports cross-region and cross-species granule-cell analysis.
2. Full-primary-core diffusion trajectories justify a stage-windowed interpretation rather than a static stop-factor model.
3. Ortholog-aware rank-meta modeling identifies a formal shared candidate space while preserving dentate/cerebellar identity separation.
4. Manuscript candidate tiers nominate a conservative Tier 1 seed set and Tier 2 wiring/synaptic support genes.
5. Named pyramidal and Purkinje comparator tests reject the over-simple claim that the shared axes are uniquely granule-specific pathways.
6. Niche/fate versus circuit/morphology modeling supports distinct upstream fate programs with stronger downstream convergence.
7. Transcriptomic configuration scoring and driver decomposition show an identity-coupled assembly-configuration signal, not a pure morphology-only code.
8. Aim 2 pathway-readiness and fitted stage-window analyses reframe TGF-beta/BDNF/SMAD/MAPK as branch-specific maturation/readiness overlays.
9. NeuroMorpho and DANDI provide morphology/activity constraints for interpreting compact input-expansion design.
10. Aim 3 sparse expansion-coding simulation and empirical calibration support the computational plausibility of granule-like design under resource constraints.

## Figure Plan Summary

Use the structured figure table in `Project/results/manuscript_figure_plan.tsv`.

- Figure 1: biological question, historical morphology, dataset map, and developmental/circuit logic.
- Figure 2: ortholog-aware workflow, Tier 1/Tier 2 candidates, and mechanism-axis organization.
- Figure 3: named-comparator specificity audit plus niche/fate versus downstream convergence.
- Figure 4: transcriptomic configuration formula, primary-core validation, and driver audit.
- Figure 5: Aim 2 niche/pathway result, fitted TGF-beta/BDNF stage-window model, Aim 3 empirical sparse-coding calibration, and stage-windowed final model.
- Figure 6: integrated working model combining the Fig1 primary-core frame, Fig2 candidate tiers, Fig4 configuration logic, Fig5 stage-window/computation logic, and NeuroMorpho/DANDI validation layer.

## Claim Table

Use `Project/results/manuscript_claim_evidence_caveat_table.tsv` as the manuscript guardrail. The central writing rule is:

- Say "convergent downstream assembly machinery" when discussing shared neurite, synaptic, and excitability modules.
- Say "identity-coupled transcriptomic assembly configuration" when discussing the combined configuration score.
- Avoid saying "unique granule-specific pathway" unless future comparator analyses change the result.
- Avoid saying "morphology is fully written in the transcriptome"; the data support transcriptomic configuration, not direct geometry.

## Candidate Genes For Main Text

Tier 1 seed genes:

- `GPM6A`: membrane/neurite outgrowth structural executor.
- `NFIB`, `NFIA`: developmental transcriptional regulators.
- `RFX3`: ciliogenesis/transcriptional regulatory candidate.
- `KCNK1`: ion-channel/excitability tuning candidate.
- `GABRA2`: GABA receptor and synaptic maturation candidate.

Tier 2 support genes:

- `ROBO2`, `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`.

Downstream morphology/executor priorities from the driver audit:

- `GPM6A`, `ROBO2`, `DCC`, `CADM3`, `STMN2`, `STMN3`, `DPYSL2`, `MAP1B`, `BASP1`, `CFL1`, `RTN1`, `RTN3`.

## Current Caveats

- The strongest evidence is rank-meta and module-balance evidence, not direct perturbation.
- Named-comparator specificity is currently limited to explicit local comparator labels from `GSE104323` and `GSE122357`.
- The human dentate branch is now strong enough for the primary frame, but `GSE268609` labels are projected and should be described as such.
- The configuration score is useful, but the driver audit shows it is driven broadly by regional fate polarity and selectively by downstream construction balance.
- NeuroMorpho and DANDI are now included as first external validation layers, but they remain indirect with respect to causality.
- Direct secreted-protein bioactivity, proteomics, neutralization/add-back, and perturbation experiments remain future validation needs.

## Immediate Next Step

The next practical step is manuscript revision and figure polish, not broad dataset hunting:

1. Use Figure 1 through Figure 6 as the first complete manuscript figure packet.
2. Fold the 2005 support/revision table into Results 7, Discussion, or a supplement.
3. Decide whether Figure 6 remains the final Results figure or becomes the graphical abstract.
4. Treat any additional data download as reviewer-driven, not as a prerequisite for the first full draft.

## Key Source Outputs

- `Project/results/integrated_primary_core_datasets.md`
- `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`
- `Project/results/primary_core_manuscript_candidate_tiers.md`
- `Project/results/primary_core_granule_specificity_named_comparators.md`
- `Project/results/primary_core_niche_circuit_module_model.md`
- `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`
- `Project/results/primary_core_configuration_driver_audit.md`
- `Project/results/primary_core_aim2_niche_pathway_model.md`
- `Project/results/aim2_stage_window_model.md`
- `Project/results/primary_core_aim3_sparse_coding_model.md`
- `Project/results/aim3_empirical_calibration.md`
- `Project/results/neuromorpho_granule_morphometry_validation.md`
- `Project/results/dandi_000003_multisession_spatial_extension.md`
- `Project/results/paper_2005_support_revision_table.md`
- `Project/results/manuscript_claim_evidence_caveat_table.tsv`
- `Project/results/manuscript_figure_plan.tsv`
