# Aim 2b Stage-Resolved TGF-beta/BDNF Audit

## Question

Does the historical TGF-beta/BDNF maturation or stop mechanism behave differently across granule-cell developmental stage, especially given that cerebellar granule cells develop mostly postnatally while dentate granule cells retain adult neurogenesis?

## Scope

This is a pathway-readiness audit, not a direct niche sender-receiver assay. Scores are percentile ranks within each dataset and axis for each gene, then summarized at pathway and signature level. Therefore values test whether a module is relatively higher at a stage/state inside the same dataset, not whether raw expression is larger across platforms.

Datasets included:

- `GSE104323`: adult mouse dentate lineage states from RGL/RGL_young through neuroblast, immature GC, juvenile GC, and adult GC.
- `GSE292261`: mouse postnatal dentate stages P5, P7, P10, P15, and P28, scored both as all DG cells and candidate dentate granule cells.
- `GSE122357`: mouse cerebellar candidate granule cells at P0 and P8 replicates.
- `GSE214309`: adult mouse dentate granule-cell maturation/activity states, including immature, mature, active immature, and active mature cells at 1 hr and 4 hr.

## Main Findings

- Dentate lineage (`GSE104323`) TGF-beta/BDNF peaks at nIPC-perin (0.844); RGL_young to GC-adult delta is -0.281, while neuroblast to immature-GC delta is -0.016.
- Postnatal dentate candidates (`GSE292261`) TGF-beta/BDNF peaks at P5 (0.800). The P15 to P28 candidate-only delta is -0.550, but P28 candidate-cell count is only 7, so P28 must be treated as a low-support endpoint.
- Cerebellar candidates (`GSE122357`) TGF-beta/BDNF peaks at P8a (0.833); P0 to P8a delta is 0.500 and P0 to P8b delta is 0.333.
- Adult dentate activity/maturation (`GSE214309`) TGF-beta/BDNF peaks at immatureactive_4hr (0.750). Immature-to-mature deltas are 0.125 at 1 hr and 0.094 at 4 hr.

## Interpretation

The 2005 TGF-beta/BDNF mechanism should not be framed as a simple cerebellum-versus-dentate regional effect. In the current primary-core data, it is better framed as a stage- and state-sensitive maturation/readiness axis. Dentate evidence is strongest because adult lineage and postnatal stages can be directly ordered, while cerebellar evidence is limited to P0 versus P8 candidate granule cells in this first pass.

This supports the project hypothesis in a more precise form: cerebellar and dentate granule cells can converge on similar morphology through shared downstream maturation and wiring modules, but the timing and upstream niche logic are region-specific. Dentate retains a lifelong progenitor-to-granule continuum, whereas cerebellum has a more developmentally bounded expansion-and-differentiation program.

## Manuscript Use

- Add Aim 2b as a stage-aware refinement under the niche pathway aim.
- Phrase the TGF-beta/BDNF result as `stage-dependent maturation/readiness`, not as a universal granule-cell stop switch.
- Use `GSE104323`, `GSE292261`, and `GSE214309` as dentate neurogenesis anchors; use `GSE122357` as the cerebellar postnatal comparator.
- Keep direct ligand-source claims for a future spatial or ligand-receptor sender analysis.

## Outputs

- Gene units: `Project/results/primary_core_aim2b_stage_tgf_bdnf_gene_units.tsv` (3,305 rows).
- Pathway units: `Project/results/primary_core_aim2b_stage_tgf_bdnf_pathway_units.tsv` (261 rows).
- Signature units: `Project/results/primary_core_aim2b_stage_tgf_bdnf_signature_units.tsv` (116 rows).
- Transitions: `Project/results/primary_core_aim2b_stage_tgf_bdnf_transitions.tsv` (112 rows).
- Summary: `Project/results/primary_core_aim2b_stage_tgf_bdnf_summary.tsv` (25 rows).
- Plot: `Project/results/primary_core_aim2b_stage_tgf_bdnf_plot.png`.
