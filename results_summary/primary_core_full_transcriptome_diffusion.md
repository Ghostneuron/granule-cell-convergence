# Primary-Core Full-Transcriptome Diffusion Trajectories

## Scope

This layer upgrades the endpoint-gene graph audit to a highly-variable-gene diffusion/pseudotime workflow across the strict 10-dataset primary core.

For datasets with local full matrices, geometry uses all available genes before HVG selection. `GSE325391` and `GSE268609` are retained as selected-feature bridge trajectories because the local analysis objects contain 2,169 selected genes rather than full cell-by-gene matrices ready for in-memory trajectory analysis.

The method is dependency-light: log1p(CP10K) normalization, HVG selection, TruncatedSVD, kNN graph, diffusion components, and shortest-path pseudotime from curated or inferred roots.

## Coverage

- Primary datasets processed: 10/10.
- Cells/nuclei/spots in trajectory tables after deterministic sampling: 56,892.
- Strict local full-transcriptome or full-DG-subset datasets: 8.
- Selected-feature bridge datasets: 2.

## Dataset Notes

- `GSE104323`: 8,000 cells, 19,091 genes available, 2,000 HVGs, root `curated_RGL_RGL_young_root` (primary_lineage_pseudotime); rho vs known/stage order 0.797.
- `GSE95752`: 1,755 cells, 16,266 genes available, 2,000 HVGs, root `inferred_high_progenitor_low_mature_root` (supporting_intrinsic_diffusion).
- `GSE292261`: 505 cells, 13,101 genes available, 2,000 HVGs, root `earliest_postnatal_P5_root` (primary_stage_pseudotime); rho vs known/stage order 0.751.
- `GSE214309`: 384 cells, 30,209 genes available, 2,000 HVGs, root `immature_DGC_state_root` (state_axis_pseudotime); rho vs known/stage order 0.684.
- `GSE122357`: 8,000 cells, 17,823 genes available, 2,000 HVGs, root `earliest_postnatal_P0_root` (primary_stage_pseudotime); rho vs known/stage order 0.750.
- `GSE165657`: 8,000 cells, 22,108 genes available, 2,000 HVGs, root `inferred_high_progenitor_low_mature_root` (supporting_intrinsic_diffusion).
- `GSE312658`: 6,248 cells, 26,819 genes available, 2,000 HVGs, root `inferred_high_progenitor_low_mature_root` (perturbation_context_diffusion).
- `GSE186538`: 8,000 cells, 28,779 genes available, 2,000 HVGs, root `inferred_high_progenitor_low_mature_root` (human_DG_intrinsic_diffusion).
- `GSE325391`: 8,000 cells, 2,081 genes available, 2,000 HVGs, root `source_differentiating_or_low_sling_root` (selected_feature_source_pseudotime_validation); rho vs known/stage order 0.579; rho vs source pseudotime 0.579.
- `GSE268609`: 8,000 cells, 2,082 genes available, 2,000 HVGs, root `projected_immature_neurogenic_root` (selected_feature_human_context_diffusion).

## Main Interpretation

The diffusion layer does not overturn the existing Fig1-5 structure. It makes the mechanism more precise: the shared construction/configuration signal is stage-windowed, and age labels alone should not be treated as pseudotime.
The clearest manuscript-level change is Figure 5: pseudotime should become an explicit layer between regional fate/niche input and final granule-cell morphology. Figure 4 can keep the configuration model, but should describe it as maturation-window dependent.

## Figure Impact

- `Figure 1`: no structural change. Add note that trajectory support is now available as an Aim 2 refinement.
- `Figure 2`: candidate tiers unchanged for now. Do not rerank genes until trajectory-aware gene tests are built; use module overlays as supporting evidence.
- `Figure 3`: refines interpretation. Mention that downstream modules are stage-windowed and not uniquely granule-specific.
- `Figure 4`: strengthens but makes it stage-aware. Add a trajectory panel or caption sentence that configuration is identity-coupled and maturation-window dependent.
- `Figure 5`: should be revised/expanded. Add pseudotime as an explicit layer between regional fate/niche and terminal morphology.

## Caveats

- This is not Scanpy/Monocle, because those packages are not installed locally; it is an equivalent lightweight diffusion/kNN trajectory layer.
- Pseudotime direction is root-dependent. Curated developmental roots are stronger than marker-inferred roots.
- `GSE325391` and `GSE268609` should not be called strict full-transcriptome trajectories until their complete source matrices are converted into trajectory-ready sparse objects.
- RNA trajectories still do not directly measure BrdU, p21/p27 protein, pERK, pSMAD, or secreted protein bioactivity.

## Outputs

- Cell scores: `Project/results/primary_core_full_transcriptome_diffusion_cell_scores.tsv.gz` (56,892 rows).
- Dataset summary: `Project/results/primary_core_full_transcriptome_diffusion_dataset_summary.tsv` (10 rows).
- Group summary: `Project/results/primary_core_full_transcriptome_diffusion_group_summary.tsv`.
- Module correlations: `Project/results/primary_core_full_transcriptome_diffusion_module_correlations.tsv`.
- Fig1-5 impact table: `Project/results/primary_core_full_transcriptome_diffusion_fig1_5_impact.tsv`.
- Overview plot: `Project/results/primary_core_full_transcriptome_diffusion_overview.png`.
- Pseudotime scatter plot: `Project/results/primary_core_full_transcriptome_diffusion_pseudotime_scatter.png`.
