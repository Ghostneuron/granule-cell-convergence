# Aim 2 Stage-Window Model Fit

Date built: 2026-06-24

## Purpose

This file adds a formal fitted-model layer to Aim 2. The model uses existing stage-resolved pathway/signature scores and asks whether the 2005 TGF-beta/BDNF axis behaves as a linear switch or as a stage-windowed maturation/readiness signal.

## Model

`signature_score ~ stage_norm + stage_norm^2 + cerebellum + cerebellum:stage_norm + cerebellum:stage_norm^2`

Scores are percentile-like within-dataset signature scores. `stage_norm` is normalized from 0 to 1 within each dataset, axis type, comparison group, and signature. Fits use weighted least squares with HC3 robust standard errors; weights are proportional to available cell count when present and to signature pathway coverage. In the current signature-level table, cell counts are not carried forward, so the fitted weights are driven by pathway coverage.

## TGF-beta/BDNF Stage-Window Result

- `cerebellum`: median slope 0.418, median endpoint delta 0.333, median observed peak stage 0.99 (late_stage_or_activity_window).
- `dentate`: median slope -0.293, median endpoint delta -0.416, median observed peak stage 0.21 (early_or_intermediate_stage_window).

Branch interaction terms from the fitted quadratic model:

- `stage_x_cerebellum` beta 14.025, HC3 p=2.75e-22, q=1.44e-21.
- `stage2_x_cerebellum` beta -13.523, HC3 p=4.8e-22, q=1.44e-21.

## Interpretation

The fitted layer supports the same conclusion as the earlier audits, but makes it quantitative: TGF-beta/BDNF is better modeled as a stage-windowed maturation/readiness axis than as a simple monotonic cerebellar stop signal. Dentate datasets provide richer stage coverage and show early/intermediate or activity-linked windows; cerebellar inference is useful but limited by three postnatal candidate granule-cell stage points in `GSE122357`.

This is still not a spatial sender-receiver model or secreted-protein bioactivity assay.

## Diffusion/Pseudotime Support

- `GSE104323` `bdnf_erk_response` rho=-0.052 (weak_or_flat).
- `GSE122357` `bdnf_erk_response` rho=-0.025 (weak_or_flat).
- `GSE165657` `bdnf_erk_response` rho=-0.026 (weak_or_flat).
- `GSE186538` `bdnf_erk_response` rho=0.017 (weak_or_flat).
- `GSE214309` `bdnf_erk_response` rho=0.129 (increases_with_pseudotime).
- `GSE268609` `bdnf_erk_response` rho=-0.333 (decreases_with_pseudotime).
- `GSE292261` `bdnf_erk_response` rho=0.054 (weak_or_flat).
- `GSE312658` `bdnf_erk_response` rho=-0.008 (weak_or_flat).
- `GSE325391` `bdnf_erk_response` rho=0.208 (increases_with_pseudotime).
- `GSE95752` `bdnf_erk_response` rho=-0.009 (weak_or_flat).
- `GSE104323` `immature_progenitor_state` rho=-0.625 (decreases_with_pseudotime).
- `GSE122357` `immature_progenitor_state` rho=-0.055 (weak_or_flat).
- `GSE165657` `immature_progenitor_state` rho=-0.146 (decreases_with_pseudotime).
- `GSE186538` `immature_progenitor_state` rho=-0.136 (decreases_with_pseudotime).
- `GSE214309` `immature_progenitor_state` rho=-0.076 (weak_or_flat).
- `GSE268609` `immature_progenitor_state` rho=-0.247 (decreases_with_pseudotime).
- `GSE292261` `immature_progenitor_state` rho=-0.371 (decreases_with_pseudotime).
- `GSE312658` `immature_progenitor_state` rho=-0.111 (decreases_with_pseudotime).
- `GSE325391` `immature_progenitor_state` rho=0.107 (increases_with_pseudotime).
- `GSE95752` `immature_progenitor_state` rho=-0.370 (decreases_with_pseudotime).
- `GSE104323` `neuronal_differentiation_maturation` rho=0.602 (increases_with_pseudotime).
- `GSE122357` `neuronal_differentiation_maturation` rho=-0.090 (weak_or_flat).
- `GSE165657` `neuronal_differentiation_maturation` rho=0.368 (increases_with_pseudotime).
- `GSE186538` `neuronal_differentiation_maturation` rho=0.111 (increases_with_pseudotime).
- `GSE214309` `neuronal_differentiation_maturation` rho=-0.093 (weak_or_flat).
- `GSE268609` `neuronal_differentiation_maturation` rho=-0.433 (decreases_with_pseudotime).
- `GSE292261` `neuronal_differentiation_maturation` rho=0.284 (increases_with_pseudotime).
- `GSE312658` `neuronal_differentiation_maturation` rho=0.308 (increases_with_pseudotime).
- `GSE325391` `neuronal_differentiation_maturation` rho=0.038 (weak_or_flat).
- `GSE95752` `neuronal_differentiation_maturation` rho=0.284 (increases_with_pseudotime).
- `GSE104323` `secreted_stop_candidate_axis` rho=-0.095 (weak_or_flat).
- `GSE122357` `secreted_stop_candidate_axis` rho=-0.025 (weak_or_flat).
- `GSE165657` `secreted_stop_candidate_axis` rho=-0.057 (weak_or_flat).
- `GSE186538` `secreted_stop_candidate_axis` rho=-0.054 (weak_or_flat).
- `GSE214309` `secreted_stop_candidate_axis` rho=0.031 (weak_or_flat).
- `GSE268609` `secreted_stop_candidate_axis` rho=-0.339 (decreases_with_pseudotime).
- `GSE292261` `secreted_stop_candidate_axis` rho=-0.243 (decreases_with_pseudotime).
- `GSE312658` `secreted_stop_candidate_axis` rho=0.076 (weak_or_flat).
- `GSE325391` `secreted_stop_candidate_axis` rho=0.225 (increases_with_pseudotime).
- `GSE95752` `secreted_stop_candidate_axis` rho=0.117 (increases_with_pseudotime).
- `GSE104323` `tgf_smad_pai1_response` rho=-0.339 (decreases_with_pseudotime).
- `GSE122357` `tgf_smad_pai1_response` rho=0.109 (increases_with_pseudotime).
- `GSE165657` `tgf_smad_pai1_response` rho=-0.121 (decreases_with_pseudotime).
- `GSE186538` `tgf_smad_pai1_response` rho=-0.036 (weak_or_flat).
- `GSE214309` `tgf_smad_pai1_response` rho=0.010 (weak_or_flat).
- `GSE268609` `tgf_smad_pai1_response` rho=-0.291 (decreases_with_pseudotime).
- `GSE292261` `tgf_smad_pai1_response` rho=-0.138 (decreases_with_pseudotime).
- `GSE312658` `tgf_smad_pai1_response` rho=-0.271 (decreases_with_pseudotime).
- `GSE325391` `tgf_smad_pai1_response` rho=0.105 (increases_with_pseudotime).
- `GSE95752` `tgf_smad_pai1_response` rho=-0.079 (weak_or_flat).

## Outputs

- Coefficients: `Project/results/aim2_stage_window_model_coefficients.tsv`
- Dataset/group fits: `Project/results/aim2_stage_window_model_group_fits.tsv`
- Branch summary: `Project/results/aim2_stage_window_model_branch_summary.tsv`
- Diffusion support: `Project/results/aim2_stage_window_model_diffusion_support.tsv`
- Plot: `Project/results/aim2_stage_window_model.png`
