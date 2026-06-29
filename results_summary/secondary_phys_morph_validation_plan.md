# Secondary Physiology And Morphology Validation Plan

Date built: 2026-06-23

## Decision

Yes, public physiology and morphology data can strengthen Aim 3, but they should be added as a secondary validation layer rather than merged into the strict 10-dataset transcriptomic core. The clean division is:

- Use `NeuroMorpho.Org` for direct dendrite/stem/branch morphometry.
- Use `DANDI:000003` as the main public dentate granule/mossy-cell activity and behavior-linked physiology source.
- Use Allen Cell Types for intrinsic electrophysiology calibration and comparator-cell feature ranges, not as primary dentate/cerebellar granule evidence.
- Use curated literature, and possibly later connectomics, for synaptic input counts where public raw datasets remain incomplete.

## Why Allen Helps But Is Not A Primary Granule Dataset

Allen Cell Types is highly useful because it provides whole-cell current-clamp recordings, computed ephys features, morphology reconstructions, NWB files, and neuronal models. However, the public overview states that data generation has focused on selected cortical and thalamic neurons, and the local API probes found no direct dentate or cerebellar structure matches.

| Probe | Matched rows | Interpretation |
|---|---:|---|
| `all_celltypes_api_rows` | 2333 | Allen Cell Types has a large ephys/morphology/model table suitable for comparator calibration. |
| `dentate_structure_name` | 0 | No direct dentate Cell Types rows in this API probe. |
| `cerebell_structure_name` | 0 | No direct cerebellar Cell Types rows in this API probe. |
| `DG_structure_acronym` | 0 | No dentate-gyrus acronym matches in this API probe. |

Therefore, Allen should be framed as a comparator/calibration resource for `intrinsic_excitability`, not evidence that cerebellar and dentate granule cells have matched physiology.

## Direct Morphology Evidence

| Resource | Matched neurons | Human | Mouse | Rat | Example stems | Example branches | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---|
| `dentate_gyrus_granule_all_species` | 9672 | 77 | 7649 | 1532 | 1.0 | 15.0 | Strong public morphology source for dentate granule dendrite and branch metrics. |
| `cerebellum_granule_all_species` | 62 | 5 | 0 | 0 | 4.0 | 16.0 | Useful but smaller cerebellar granule morphology source; strict mouse/rat filters are sparse or absent in current query. |

These values are only probes, not final filtered estimates. The next quantitative step is to download matched NeuroMorpho metadata/morphometry tables, exclude disease/genotype or incomplete reconstructions when appropriate, and compare distributions of `n_stems`, `n_branch`, `n_bifs`, length, and branch order.

That first quantitative step has now been built in `Project/results/neuromorpho_granule_morphometry_validation.md`. The first-pass sample includes 558 dentate granule reconstructions and all 62 cerebellar granule reconstructions from the strict query. The key result is not identical geometry: dentate and cerebellar granule cells differ strongly in primary stem count and dendritic length, but have strikingly similar median branch counts. This supports a refined claim of convergent compact/limited-branch input-expansion architecture rather than literal dendritic sameness.

## Direct Activity Evidence

`DANDI:000003` is the strongest current public activity source for the dentate branch: 101 NWB files, 16 subjects, and measured variables `DecompositionSeries,LFP,Units,Position,ElectricalSeries`. It can support firing sparsity, active fraction, spatial coding, and behavior-linked population-separation analyses.

`DANDI:000165` is useful as supporting DG/CA3/CA1 network physiology, especially for LFP/unit state context, but it is not a pure dentate granule-cell dataset.

A first DANDI 000003 pilot has now been completed in `Project/results/dandi_000003_activity_sparsity_pilot.md`. The smallest NWB file, `YutaMouse41-150829`, was downloaded and read directly with `h5py`. The pilot contains 23 units, including 3 source-labeled granule cells. Across the full recording, granule-cell-labeled units have median firing rate 0.9131 Hz and median active 1 s-bin fraction 0.3943. This supports using DANDI for direct activity validation, but it also shows that "sparse" must be defined carefully: this pilot measures temporal firing across sleep/wake intervals, not task-specific spatial/pattern-separation sparsity.

A first position-linked DANDI pilot has also been completed in `Project/results/dandi_000003_spatial_pattern_pilot.md`. Using awake-moving position samples from the same NWB file, the 3 labeled granule units show median spatial information 0.3098 bits/spike, median spatial sparsity 0.6273, median active spatial-bin fraction 0.4634, and median awake-moving firing rate 0.6414 Hz. This confirmed that the dataset can support spatial coding analyses.

The DANDI spatial workflow has now been extended to six analyzed sessions in `Project/results/dandi_000003_multisession_spatial_extension.md`. The current set includes `YutaMouse41-150829`, `YutaMouse37-150609`, `YutaMouse42-151102`, `YutaMouse55-160908`, `YutaMouse55-160909`, and the external-drive `YutaMouse55-160907`, totaling 124 units. Direct granule-labeled evidence comes from all six sessions, with 26 labeled granule units. Across these granule-labeled units, the pooled median spatial information is 0.7800 bits/spike, median spatial sparsity is 0.4537, median active spatial-bin fraction is 0.5489, and median awake-moving firing rate is 0.0627 Hz. Session-level granule population-vector checks show far-minus-near Euclidean separation of 0.1842, 1.2293, 0.0116, 0.0615, 0.0280, and 0.2260. The previous `YutaMouse37-150617` comparator file was removed after showing 0 source-labeled granule units and is now excluded from renewed download priority. This remains a pilot rather than a final task-specific pattern-separation result because the population-vector evidence is heterogeneous.

## Parameter Map For The Sparse-Coding Model

| Model parameter | Best public resource | Evidence grade | Use in model |
|---|---|---|---|
| `expansion_ratio` | curated_anatomical_literature;Allen_WMB_taxonomy_for_cell_labels | partial_prior_only | Set plausible ranges for the sparse-expansion simulation. |
| `input_degree` | NeuroMorpho_DG_granule_morphometry;NeuroMorpho_cerebellar_granule_morphometry;curated_cerebellar_synapse_literature | direct_morphology_partial_synapse | Replace the arbitrary Aim 3 input-degree grid with empirically plausible DG/CB ranges. |
| `output_active_fraction` | DANDI_000003_DG_granule_mossy_activity;Allen_Cell_Types_ephys_as_calibration | direct_DG_activity_comparator_Allen | Empirically constrain the output sparsity term in the Aim 3 model. |
| `pattern_separation_behavior` | DANDI_000003_multisession_spatial_extension;DANDI_000003_DG_granule_mossy_activity | first_pass_six_session_local_five_granule_session_DG_spatial_pilot | Test whether fitted sparse parameters improve separation of similar behavioral states. |
| `intrinsic_excitability` | Allen_Cell_Types_ephys;GSE214905_patch_seq_DG | strong_comparator_weak_DG_direct | Connect synaptic/excitability transcriptomic modules to measured electrophysiology feature families. |
| `morphology_similarity` | NeuroMorpho.Org API | direct_but_unbalanced | Quantify the morphology part of the central hypothesis instead of relying on qualitative similarity. |
| `secreted_stop_or_maturation_inputs` | existing_project_secretome_screen;GSE242688_spatial_proteomics;future_validation_assays | RNA_candidate_not_bioactivity | Keep Aim 2 mechanistically connected to the 2005 conditioned-medium question. |

## Recommended Next Analysis

1. Refine the NeuroMorpho analysis with stricter archive, genotype/condition, species, and reconstruction-completeness filters before final manuscript statistics.
2. Refit the Aim 3 model using separate empirical priors for `primary_input_sampling_stems_or_claws` and `dendritic_field_complexity`, rather than a single morphology/input-degree knob.
3. Extend the DANDI:000003 pilot beyond the current six local sessions only after prioritizing likely granule-containing files, deliberately switching to new-subject breadth, or freeing storage; then add stronger task-epoch and trajectory-specific population-vector analyses for pattern separation.
4. Add Allen Cell Types ephys as an intrinsic-excitability calibration table for comparator neurons and for mapping transcriptomic excitability modules to measured ephys feature families.
5. Keep cerebellar granule-cell activity and synaptic input count as explicit evidence gaps unless a direct public cerebellar granule ephys/connectomics dataset is verified.

## Outputs

- Resource table: `Project/results/secondary_phys_morph_candidate_resources.tsv`
- Parameter map: `Project/results/secondary_phys_morph_parameter_map.tsv`
- Allen API probe: `Project/results/secondary_phys_morph_allen_celltypes_probe.tsv`
- NeuroMorpho API probe: `Project/results/secondary_phys_morph_neuromorpho_probe.tsv`
- NeuroMorpho morphometry validation: `Project/results/neuromorpho_granule_morphometry_validation.md`
- DANDI 000003 activity/sparsity pilot: `Project/results/dandi_000003_activity_sparsity_pilot.md`
- DANDI 000003 spatial pattern pilot: `Project/results/dandi_000003_spatial_pattern_pilot.md`
- DANDI 000003 multi-session spatial extension: `Project/results/dandi_000003_multisession_spatial_extension.md`
- DANDI 000003 targeted download priority: `Project/results/dandi_000003_targeted_download_priority.md`
