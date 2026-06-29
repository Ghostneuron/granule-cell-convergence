# 2005 Paper Endpoint Stage/Pseudotime Audit

## Question

Are the other experimental readouts in the 2005 Development paper stage-dependent, and do they support a trajectory/pseudotime interpretation?

## Answer

Yes. The paper's readouts are better interpreted as trajectory-dependent than as one static factor effect. RNA modules corresponding to BrdU/cell-cycle, p21/p27 cell-cycle exit, neuronal differentiation/maturation, TGF-beta/SMAD/PAI1 response, BDNF/ERK response, apoptosis/survival, and secreted stop-factor candidates vary across ordered dentate and cerebellar axes.

Important limitation: this analysis uses RNA proxies. It cannot directly measure BrdU incorporation, p21/p27 protein, MAP2/synapsin protein, ERK phosphorylation, SMAD nuclear translocation, PAI1 reporter activity, or apoptosis by TUNEL/Hoechst.

## Key Stage-Dependence Results

- Adult dentate lineage (`GSE104323`): RNA maturation pseudotime versus lineage order is rho 0.714, p 0.047. Proliferation delta is RGL_young to GC-adult: -0.625; differentiation delta is RGL_young to GC-adult: 0.625.
- Postnatal dentate candidates (`GSE292261`): RNA maturation pseudotime versus age is rho 0.100, p 0.873. Proliferation delta is P5 to P28: -0.600; TGF/SMAD delta is P5 to P28: -0.800.
- Postnatal cerebellar candidates (`GSE122357`): RNA maturation pseudotime versus age is rho -1.000, p 0.000. Proliferation delta is P0 to P8b: 0.333; TGF/SMAD delta is P0 to P8b: 0.333.
- Adult dentate activity/maturation (`GSE214309`): RNA maturation pseudotime versus ordered state is rho -0.381, p 0.352. This axis is activity/time-state ordered, not a clean developmental lineage.

## Interpretation For The 2005 Paper

- The anti-proliferation result maps to a trajectory shift: proliferative modules and immature/progenitor modules tend to separate from cell-cycle-exit, maturation, and stop-factor modules.
- The p21/p27 result should be treated as a stage-sensitive cell-cycle-exit module, not only as an acute response marker.
- The MAP2/TuJ1/synapsin differentiation result is strongly compatible with a maturation trajectory.
- The TGF-beta/SMAD and BDNF/ERK results are partly stage-dependent, but RNA cannot substitute for pSMAD nuclear translocation or pERK signaling assays.
- The apoptosis/neuroprotection fraction in the paper is the least safely inferred from RNA alone and should remain an experimental validation target.

## Why True Pseudotime Is Needed Next

This audit uses ordered stages and module-defined RNA pseudotime. A stronger next step is cell-level graph pseudotime, ideally using `GSE104323` dentate lineage cells and `GSE292261` postnatal dentate cells first, then `GSE122357` cerebellar P0/P8 cells. That would test whether TGF/BDNF/secretome modules rise at a specific trajectory window after proliferation and before mature neuronal markers.

## Outputs

- Gene units: `Project/results/primary_core_2005_endpoint_pseudotime_gene_units.tsv` (2,709 rows).
- Module units: `Project/results/primary_core_2005_endpoint_pseudotime_module_units.tsv` (261 rows).
- Trajectory scores: `Project/results/primary_core_2005_endpoint_pseudotime_trajectory_scores.tsv` (29 rows).
- Correlations: `Project/results/primary_core_2005_endpoint_pseudotime_correlations.tsv` (60 rows).
- Transitions: `Project/results/primary_core_2005_endpoint_pseudotime_transitions.tsv` (60 rows).
- Plot: `Project/results/primary_core_2005_endpoint_pseudotime_audit.png`.
