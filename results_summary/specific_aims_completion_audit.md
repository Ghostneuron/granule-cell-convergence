# Specific Aims Completion Audit

Date updated: 2026-06-25

Source aims file: `Project/specific_aims_granule_cell_convergence.md`

## Overall Verdict

Yes, the original three aims now have first-pass computational completion, and the previously optional trajectory layer has now been built across the strict primary core. The project is ready to move from analysis building into manuscript drafting and figure revision, while keeping the following caveats explicit:

- Aim 1 is complete as a stage/module/rank-meta plus full primary-core diffusion/pseudotime analysis. The caveat is that `GSE325391` and `GSE268609` are selected-feature bridge trajectories, not strict full-transcriptome trajectories yet.
- Aim 2 is complete at the primary-core pathway-readiness, stage-resolved, conditioned-medium candidate, trajectory-overlay, fitted stage-window model, and focused raw-expression sender-receiver ligand-receptor levels, although the current analysis is not a spatial adjacency assay and does not directly measure secreted protein bioactivity.
- Aim 3 is complete as a sparse expansion-coding model linked to transcriptomic construction modules, with empirical calibration against NeuroMorpho and DANDI outputs. First direct NeuroMorpho morphometry, DANDI 000003 activity, and six-session position-linked spatial-coding pilots have also been completed; direct granule-labeled spatial evidence now comes from all six analyzed sessions, and task-specific pattern-separation analysis remains an optional upgrade.

## Aim-By-Aim Status

| Aim | Original goal | Current status | What is finished | What remains |
|---|---|---|---|---|
| Aim 1 | Define matched maturation trajectories for cerebellar and dentate granule cells. | Complete as a primary-core computational trajectory aim. | Built a strict 10-dataset primary core; included mouse dentate, cerebellum, and human dentate/hippocampus; built ortholog-aware rank-meta models; identified shared and divergent maturation/morphogenesis/synaptic candidate modules; added full-transcriptome or available-feature HVG diffusion/pseudotime across all 10 primary datasets. Ordered anchors validate pseudotime direction: `GSE104323` rho 0.797, `GSE292261` rho 0.751, `GSE122357` rho 0.750, `GSE214309` rho 0.684, and `GSE325391` source-pseudotime rho 0.579. | For the manuscript, phrase this as full primary-core diffusion/pseudotime with dataset-specific evidence grades. Do not overstate `GSE325391` and `GSE268609` as strict full-transcriptome trajectories until full sparse source objects are exported. |
| Aim 2 | Test whether niche-signaling programs recapitulate the TGF-beta2/BDNF/SMAD-MAPK mechanism. | Complete as a pathway-readiness, trajectory-overlay, fitted stage-window, and focused sender-receiver LR aim; not yet a spatial adjacency or secreted-protein assay. | Built a niche/fate versus circuit/morphology model; targeted pathway/ligand-receptor audit for TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP, Reelin, Semaphorin, SHH, WNT, FGF, and Notch; stage-resolved TGF-beta/BDNF audit; conditioned-medium secretome candidate screen; 2005 endpoint pseudotime audit; full primary-core trajectory overlay; formal fitted stage-window model; focused sender-receiver LR prediction using `GSE122357` cerebellar Purkinje/astroglial-proxy/microglia/endothelial senders and `GSE104323` SGZ astrocyte/endothelial/PVM-proxy senders toward granule-lineage receivers. Key result: TGF-beta/BDNF/SMAD/ERK behave as branch-specific, stage-windowed maturation/readiness overlays; the sender-receiver layer highlights cerebellar Purkinje SHH/IGF1/guidance interactions and dentate astrocyte/PVM/endothelial APOE/Semaphorin/Notch-like interactions. | Optional upgrade: add spatial adjacency/protein evidence if the manuscript needs direct niche-contact or conditioned-medium bioactivity claims; refine Bergmann-specific and microglia-specific labels with datasets that explicitly separate those classes; validate candidates by proteomics/ELISA/Luminex, neutralization, or add-back experiments. |
| Aim 3 | Link morphology-associated modules to sparse-coding and pattern-separation performance. | Complete as first-pass computational model, with direct morphology, activity, spatial-coding pilots, and empirical calibration. | Built morphology/executor candidate lists, named-comparator specificity audit, transcriptomic configuration score, primary-core validation, driver decomposition, sparse expansion-coding model, and empirical calibration against NeuroMorpho/DANDI outputs. Completed NeuroMorpho morphometry validation: 558 dentate granule reconstructions and 62 cerebellar granule reconstructions support similar median branch counts but different stem count and dendritic-field scale. Completed DANDI 000003 activity and six-session position-linked analysis: 124 total units, 26 labeled granule units across six granule-containing sessions, pooled granule median spatial information 0.7800 bits/spike, active spatial-bin fraction 0.5489, and granule population-vector far-minus-near Euclidean separation of 0.1842, 1.2293, 0.0116, 0.0615, 0.0280, and 0.2260. Empirical calibration shows raw dense expansion can score well, whereas resource/morphology-constrained calibration shifts the optimum toward intermediate or sparse granule-like designs. | Optional upgrade: refine NeuroMorpho filters for manuscript-grade statistics, split Aim 3 morphology into primary stems/claws and dendritic-field complexity, extend DANDI extraction beyond six local sessions only with yield/breadth prioritization or freed storage, and compute task/trajectory-specific population-vector separation. Cerebellar granule-cell activity remains a public-data gap unless a direct dataset is verified. |
| Expected product | Original in silico manuscript or strong hypothesis/review paper. | Manuscript scaffold, first draft text, revised figure packet, and fitted-model reports are complete. | Built manuscript planning packet, claim/evidence/caveat table, figure plan, six composite figures, trajectory impact audit for Fig1-5, integrated Fig6 working model, and fitted Aim 2/Aim 3 report layer. | Revise the manuscript text, then polish figure sizing and typography for the target journal. |

## Evidence Already Completed

- Strict 10-dataset primary core: `Project/results/integrated_primary_core_datasets.md`
- Formal MGI ortholog rank-meta model: `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`
- Manuscript candidate tiers: `Project/results/primary_core_manuscript_candidate_tiers.md`
- Mechanism-axis model: `Project/results/primary_core_mechanism_axis_model.md`
- Named-comparator specificity audit: `Project/results/primary_core_granule_specificity_named_comparators.md`
- Niche/fate versus circuit/morphology model: `Project/results/primary_core_niche_circuit_module_model.md`
- Aim 2 niche/pathway signaling audit: `Project/results/primary_core_aim2_niche_pathway_model.md`
- Aim 2 focused sender-receiver LR prediction: `Project/results/aim2_sender_receiver_lr.md`
- Aim 2 fitted stage-window model: `Project/results/aim2_stage_window_model.md`
- Aim 3 sparse expansion-coding model: `Project/results/primary_core_aim3_sparse_coding_model.md`
- Aim 3 empirical sparse-coding calibration: `Project/results/aim3_empirical_calibration.md`
- Transcriptomic configuration model: `Project/results/primary_core_transcriptomic_configuration_model.md`
- Primary-core configuration validation: `Project/results/primary_core_transcriptomic_configuration_primary_validation.md`
- Configuration driver audit: `Project/results/primary_core_configuration_driver_audit.md`
- Aim 2b stage-resolved TGF-beta/BDNF audit: `Project/results/primary_core_aim2b_stage_tgf_bdnf.md`
- Cerebellar conditioned-medium secretome candidate screen: `Project/results/cerebellar_conditioned_medium_secretome_candidates.md`
- 2005 endpoint stage/module pseudotime audit: `Project/results/primary_core_2005_endpoint_pseudotime_audit.md`
- 2005 endpoint cell-level graph pseudotime: `Project/results/primary_core_2005_endpoint_graph_pseudotime.md`
- Full primary-core diffusion/pseudotime trajectory layer: `Project/results/primary_core_full_transcriptome_diffusion.md`
- TGF-beta/BDNF old-versus-full trajectory comparison: `Project/results/primary_core_tgf_bdnf_old_vs_full_trajectory_comparison.tsv`
- Secondary physiology/morphology validation plan: `Project/results/secondary_phys_morph_validation_plan.md`
- NeuroMorpho granule morphometry validation: `Project/results/neuromorpho_granule_morphometry_validation.md`
- DANDI 000003 activity/sparsity pilot: `Project/results/dandi_000003_activity_sparsity_pilot.md`
- DANDI 000003 spatial pattern pilot: `Project/results/dandi_000003_spatial_pattern_pilot.md`
- DANDI 000003 multi-session spatial extension: `Project/results/dandi_000003_multisession_spatial_extension.md`
- DANDI 000003 targeted download priority: `Project/results/dandi_000003_targeted_download_priority.md`
- Manuscript figure assembly: `Project/results/manuscript_figures/manuscript_figure_assembly.md`
- Manuscript planning packet: `Project/results/manuscript_planning_packet.md`

## Refined Central Hypothesis

The original hypothesis remains viable but should be sharpened. The current evidence supports this version:

Distinct dentate and cerebellar developmental systems converge on a compact excitatory input-expansion neuron design through identity-coupled, stage-windowed transcriptomic assembly configuration. Upstream fate programs remain branch-specific, while downstream neurite, synaptic, excitability, and maturation machinery provides the strongest shared-convergence layer. TGF-beta/BDNF/SMAD/ERK and other niche signals are best framed as trajectory-windowed modulators of maturation/readiness rather than universal granule-cell stop signals. Sparse expansion-coding constraints are computationally plausible, especially when separation is evaluated under wiring/activity resource constraints and morphology-informed input-sampling priors.

## Current Figure Impact

- Fig1: no structural change; add trajectory support as a completed Aim 1/Aim 2 refinement.
- Fig2: candidate tiers remain unchanged until a trajectory-aware gene-level test is built.
- Fig3: refine language to say downstream construction modules are stage-windowed and not uniquely granule-specific.
- Fig4: strengthen the configuration model, but describe it as maturation-window dependent.
- Fig5: revise or expand. Add pseudotime/stage-window logic between regional fate/niche inputs and final granule-cell morphology.

## Recommended Next Work

1. Move manuscript work forward.
   - Use the revised six-figure packet as the current figure backbone; Fig5 now carries the stage-windowed model and Fig6 provides the integrated working model.
   - Update the abstract and result spine to state that the shared construction program is stage-windowed and identity-coupled.

2. Choose optional validation upgrades only if needed for target-journal strength.
   - Aim 1 upgrade: export full source objects for `GSE325391` and `GSE268609` if strict full-transcriptome human trajectories are needed.
   - Aim 2 upgrade: spatial adjacency/protein validation, or a human marker-projected sender extension if explicit human SGZ niche-cell labels are needed.
   - Aim 3 upgrade: stricter NeuroMorpho manuscript filters, extend DANDI 000003 beyond the current six-session local spatial pilot only if storage permits or after prioritizing likely granule-containing or subject-breadth files, add task/trajectory-specific population-vector analysis, Allen Cell Types comparator ephys calibration, and connectomics/literature priors for synaptic input count.

3. Keep the main trajectory caveat in the manuscript.
   - The full-transcriptome trajectory validates developmental ordering, but TGF-beta/BDNF modules are not identical to the previous endpoint-gene map. They should be interpreted as stage-windowed overlays, not as the global pseudotime driver.
