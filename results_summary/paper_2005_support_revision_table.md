# 2005 Paper Support And Revision Table

Date updated: 2026-06-24

## Purpose

This table translates the 2005 conditioned-medium paper into current manuscript language: what modern sequencing supports, what it revises, and what still requires direct protein or perturbation validation.

## Summary

The sequencing data support the broad direction of the 2005 paper: cerebellar conditioned-medium effects are consistent with a shift away from proliferation/progenitor state and toward maturation/differentiation readiness. The revision is timing and mechanism. TGF-beta/BDNF/SMAD/MAPK should be framed as branch-specific, stage-windowed maturation/readiness signaling rather than a universal cerebellar stop signal. The modern data also nominate a broader multifactor secretome model.

## Manuscript-Ready Table

| 2005 element | Current support | Revision | Caveat |
|---|---|---|---|
| Conditioned medium suppresses proliferation | Proliferation modules fall along dentate maturation axes; adult dentate RGL-to-GC delta -0.625 and postnatal dentate P5-to-P28 delta -0.600. | A trajectory/state transition is better supported than a static stop-factor model. | RNA modules do not measure BrdU or conditioned-medium bioactivity. |
| Conditioned medium promotes differentiation | Differentiation/maturation modules rise along ordered axes; adult dentate differentiation delta +0.625. | The result fits maturation-readiness rather than only acute differentiation. | RNA cannot replace MAP2, TuJ1, synapsin, or morphology assays. |
| TGF-beta2/SMAD/PAI1 candidate pathway | `TGFB2` remains a supported anchor; fitted Aim 2 model supports branch-specific stage windows. | TGF-beta is not a universal cerebellar stop signal. | No direct pSMAD, PAI1 reporter, or ligand-bioactivity measurement. |
| BDNF/TrkB/MAPK candidate pathway | Combined TGF-beta/BDNF index is positive in 48/63 candidate contrasts; branch timing differs. | BDNF is historically functional but transcriptomically sparse. | RNA does not measure secreted BDNF, TrkB activation, or pERK. |
| p21/p27 cell-cycle exit | Cell-cycle-exit modules vary with ordered stage/pseudotime. | p21/p27 should be interpreted as stage-sensitive cell-cycle exit. | No protein-level p21/p27 measurement. |
| Multiple medium activities | Secretome inference nominates SFRP1, BMP6, RELN, TGFB2, SFRP2, GDF/activin, semaphorin/slit, and ECM candidates. | Multifactor medium activity is more plausible than one replacement factor. | Need proteomics, ELISA/Luminex, neutralization, and add-back tests. |
| Apoptosis/survival observations | RNA includes apoptosis/survival and counter-signal candidates such as PTN, MDK, and FGF-family factors. | Treat as context-dependent and secondary. | RNA is weak for direct apoptosis/survival inference. |
| Developmental timing | Dentate TGF-beta/BDNF index peaks early/intermediate; cerebellar peak is late in the fitted model. | Same pathway family, different branch-specific timing. | Cerebellar stage inference has limited stage points. |

## Files

- Full TSV: `Project/results/paper_2005_support_revision_table.tsv`
- Endpoint audit: `Project/results/primary_core_2005_endpoint_pseudotime_audit.md`
- Fitted Aim 2 model: `Project/results/aim2_stage_window_model.md`
- Secretome candidates: `Project/results/cerebellar_conditioned_medium_secretome_candidates.md`
