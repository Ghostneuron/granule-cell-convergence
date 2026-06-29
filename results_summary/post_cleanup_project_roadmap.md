# Post-Cleanup Project Roadmap

Date updated: 2026-06-24

## Executive Decision

Proceed toward an original in silico analysis manuscript. The project is now stronger as a manuscript than as an open-ended dataset hunt.

The central model should be:

Distinct dentate and cerebellar developmental lineages converge on a compact excitatory input-expansion neuron design through an identity-coupled, stage-windowed transcriptomic assembly configuration. Upstream fate and niche programs remain branch-specific, while downstream neurite, synaptic, excitability, and maturation machinery provides the strongest shared-convergence layer. Sparse expansion-coding constraints make the shared compact granule-like design computationally plausible.

The safest phrasing is not that cerebellar and dentate granule cells share one universal granule-cell fate program. The evidence supports a convergent assembly configuration, not a single shared identity.

## Data Freeze

The cleaned workspace should now be treated as a data-freeze point for manuscript drafting.

- Current free disk space after cleanup: about 124 GiB.
- Large raw downloads should now go to the external drive root `/Volumes/VV 2021 backup drive 01/Hippocampus_Cerebellum_downloads`, linked inside the workspace as `Project/external_data`.
- Broad downloading should stop.
- Existing derived tables and figures are sufficient for a full first manuscript draft.
- Deleted raw GEO/DANDI/cache files should not be regenerated unless a specific figure, claim, or reviewer-style concern requires them.
- Current DANDI six-session outputs are a derived snapshot combining five legacy local NWB files plus the externally stored `YutaMouse55-160907` file. All six analyzed sessions now contain source-labeled granule units. One previously downloaded comparator NWB file, `YutaMouse37-150617`, was removed after showing 0 source-labeled granule units and is excluded from renewed download priority.

## Manuscript Backbone

| Figure | Current role | Readiness | Next action |
|---|---|---|---|
| Figure 1 | Biological puzzle, primary-core dataset frame, developmental/circuit logic. | Draft assembled; still visually most improvable. | Polish the conceptual panels and make the dataset map clean enough for first submission. |
| Figure 2 | Ortholog-aware candidate tiers and mechanism axes. | Manuscript-ready as a draft. | Keep Tier 1/Tier 2 genes concise; avoid overfilling the panel. |
| Figure 3 | Specificity caveat plus upstream-versus-downstream convergence. | Manuscript-ready as a draft. | Use this to prevent overclaiming a unique granule-specific pathway. |
| Figure 4 | Transcriptomic configuration score, primary-core validation, driver decomposition. | Manuscript-ready as a draft. | Emphasize "identity-coupled assembly configuration" rather than "morphology-only code." |
| Figure 5 | Aim 2 fitted stage-window niche signaling, trajectory refinement, Aim 3 sparse-coding calibration, final stage-window/resource-constraint model. | Manuscript-ready as a draft. | This is the key synthesis figure; keep TGF-beta/BDNF as maturation/readiness overlays and sparse coding as a resource-constrained plausibility model. |
| Figure 6 | Integrated working model. | Useful as final synthesis or graphical abstract. | Use as the final conceptual model; decide later whether it stays as Figure 6 or becomes graphical abstract. |

Primary figure packet:

- `Project/results/manuscript_figures/manuscript_figure_assembly.md`
- `Project/results/manuscript_figures/manuscript_figure_manifest.tsv`
- `Project/results/manuscript_figure_plan.tsv`

## Current Claim Guardrails

Use these as manuscript writing rules:

- Strong claim: the strict primary core is complete enough for manuscript-scale analysis.
- Strong claim: dentate and cerebellar granule cells are regionally identity-separated.
- Strong claim: downstream neurite/morphology and synaptic/excitability modules converge more than upstream fate/niche modules.
- Strong claim: the transcriptomic configuration score is broadly positive across candidate granule populations.
- Refined claim: the signal is identity-coupled, not a pure morphology-only transcriptomic code.
- Moderate claim: prioritized genes such as `GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`, `ROBO2`, `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, and `GRIN2B` are candidate regulators/executors, not proven causal drivers.
- Moderate claim: sparse expansion coding is computationally plausible under wiring/activity constraints; raw dense expansion can score well, but the resource/morphology-constrained fit shifts the optimum toward intermediate or sparse granule-like designs.
- Avoid: "unique granule-specific pathway."
- Avoid: "TGF-beta/BDNF is the universal cerebellar stop signal."
- Avoid: "transcriptomics directly encodes morphology."

Primary guardrail table:

- `Project/results/manuscript_claim_evidence_caveat_table.tsv`

## Priority Work Packages

### P0. Revise The First Manuscript Draft

This is the highest-value next step. The first draft skeleton now exists, and more analysis should be driven by manuscript gaps.

Recommended draft structure:

1. Introduction: the historical puzzle from the thesis and 2005 Development paper.
2. Results 1: strict primary core and cross-region/cross-species dataset construction.
3. Results 2: ortholog-aware candidate tiers.
4. Results 3: named-comparator specificity and upstream/downstream module separation.
5. Results 4: transcriptomic configuration score and driver decomposition.
6. Results 5: stage-windowed TGF-beta/BDNF/niche signaling and conditioned-medium reinterpretation.
7. Results 6: sparse expansion-coding model plus DANDI/NeuroMorpho validation layer.
8. Discussion: convergent assembly configuration, circuit constraints, 2005-paper reinterpretation, and experimental predictions.

### P1. Polish Figure 1-6 For Draft Submission

Focus first on typography, panel balance, and simple labels. The science is mostly assembled; the figures need to read clearly without the reader knowing the project history.

Highest-priority visual fixes:

- Figure 1: make the biological question and dataset map immediately understandable.
- Figure 5: make the stage-window logic visible without overcrowding.
- Figure 6: decide whether this is a final model figure or graphical abstract.

### Completed No-Download Validation Upgrades

These were completed after cleanup and should now be cited directly in the manuscript draft.

1. Aim 2 fitted stage-window model completed: `Project/results/aim2_stage_window_model.md`.
2. Aim 3 empirical sparse-coding calibration completed: `Project/results/aim3_empirical_calibration.md`.
3. DANDI six-session spatial results and NeuroMorpho morphometry are now usable as manuscript-facing validation layers.
4. The compact 2005 paper support/revision table is complete: `Project/results/paper_2005_support_revision_table.md`.
5. The manuscript draft now includes trajectory-stage, morphology, physiology, and sparse-coding calibration language.

### P1. Polish Manuscript And Figures

The next value step is revision rather than further dataset construction.

1. Tighten the abstract and Results spine around the fitted Aim 2 and Aim 3 outputs.
2. Update Figure 5 labels so readers can distinguish raw separation, resource-constrained calibration, and morphology/activity validation.
3. Decide whether Figure 6 stays as the final Results figure or becomes the graphical abstract.
4. Fold the completed 2005 support/revision table into the Results/Discussion or supplement.

### P2. Optional Small Targeted Download

Only do this after the first manuscript draft shows a real need for more physiology evidence.

One targeted DANDI file has now been added on the external drive:

- `YutaMouse55-160907`, 8.78 GB, yield-first track; added 7 source-labeled granule units.

Rule: stop here for the current manuscript draft unless a specific figure or review concern requires more physiology. If a further DANDI expansion becomes necessary, the refreshed targeted table currently points to `YutaMouse55-160911` as the next yield-first candidate, while the size-ranked breadth option is `YutaMouse44-151128`.

### P2. Optional Raw-Object Rebuilds

Only consider this if a target journal or internal review insists that the selected-feature bridge trajectories are insufficient.

Possible rebuilds:

- `GSE325391`: full sparse source object for stricter human adult dentate trajectory.
- `GSE268609`: fuller RNA object for aging/AD hippocampal projection and trajectory support.
- `GSE186538`: full source re-download only if the DG candidate subset needs to be rebuilt from scratch.

These are not needed for the first manuscript draft.

### P3. Deferred Wet-Lab Or External Validation Ideas

These are good Discussion/future-work items, not blockers for drafting.

- Conditioned-medium proteomics or cytokine panel for cerebellar granule-cell secreted inhibitors.
- Neutralization/add-back tests for TGF-beta/BDNF plus the new secretome candidates.
- Perturbation of Tier 1/Tier 2 candidate genes in granule-cell differentiation or neurite morphology assays.
- Direct cerebellar granule-cell activity or connectomics dataset if a strong public source is later identified.

## Recommended Next Concrete Step

Revise the manuscript draft using the completed fitted-model outputs, then polish Figure 5 and Figure 6. Treat any remaining analysis as figure-specific support rather than exploration.

Current draft:

- `Project/manuscript/granule_cell_convergence_manuscript_draft.md`

## Immediate To-Do List

| Priority | Task | Output | Stop condition |
|---|---|---|---|
| P0 | Revise manuscript draft with fitted Aim 2/Aim 3 outputs. | `Project/manuscript/granule_cell_convergence_manuscript_draft.md` | Abstract, Results 7, Results 9, Discussion, and limitations all reflect the fitted model results. |
| P0 | Convert claim table into Results/Discussion guardrails. | Manuscript text. | No unsupported causal or uniqueness language remains. |
| P1 | Figure 1 polish. | Revised `figure1_primary_core_concept` files. | Biological question and dataset frame are understandable in one glance. |
| P1 | Figure 5/6 polish. | Revised final-model figure(s). | Stage-window, resource-constrained sparse coding, and empirical calibration logic are visible without dense explanatory text. |
| P1 | Integrate 2005 support/revision table. | Compact table/panel or supplement. | Manuscript distinguishes sequencing support from untested secreted-protein bioactivity. |
| P2 | One targeted DANDI download, if needed. | Updated DANDI validation table. | Stop after one file unless it materially improves the claim. |
