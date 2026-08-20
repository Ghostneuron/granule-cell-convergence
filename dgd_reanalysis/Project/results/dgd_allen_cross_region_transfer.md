# Cross-region granule-cell configuration transfer

## Question

Does a molecular pattern that distinguishes dentate granule cells from CA1/CA3 pyramidal cells also distinguish cerebellar granule cells from Purkinje cells, and vice versa?

## Design

- The unit of analysis is a library-population pseudobulk, not an individual cell.
- Training used only complete within-library sets: cerebellar granule plus Purkinje or dentate granule plus CA1 plus CA3.
- The fitted classifier was transferred to the other region without refitting or threshold tuning.
- Training-label permutations were blocked within library. Matched random panels were drawn from 2,500 genes present in all four Allen matrices and matched on overall expression and detection.
- The mature-state test used 63 cerebellar library pairs and 17 dentate library triplets. The immature dentate sensitivity test used nine triplets.

## Main result

- Tier 1 ranked the held-out populations in the expected direction (dentate-to-cerebellum AUC 1.000; cerebellum-to-dentate AUC 0.901; minimum AUC 0.901), but this was not exceptional under blocked training-label permutations (intersection-union P=0.2964) or expression/detection-matched panels (P=0.3598). The training-derived decision threshold also failed to transfer (balanced accuracies 0.508 and 0.500).
- Tier 1+2 showed the same descriptive pattern at lower rank performance (AUCs 1.000 and 0.784; minimum 0.784; matched-panel P=0.5302).
- The broad downstream panel did not transfer as a common configuration. Both directions were reversed relative to the trained label (AUCs 0.000 and 0.000).
- Direct contrast concordance separated recurrence from configuration: 0.833 of Tier 1 genes were positive in both local contrasts (matched-panel P=0.0004998), but their cross-region contrast correlation was only 0.143 (matched-panel P=0.5567). The broad downstream contrast vectors were anti-aligned (Spearman rho=-0.179).

## Stage sensitivity

- For the cerebellum-trained downstream model, the median paired margin changed from -6.495 in mature dentate granule cells to 18.042 in the adult immature-neuron state (immature-minus-mature 24.537; paired P=0.003906).
- Tier 1 remained positive in both dentate states; its immature-minus-mature paired-margin shift was 11.442 (paired P=0.003906).
- The immature group is an adult HPF transcriptomic state, not a developmental-age series. Cerebellar developmental-stage transfer cannot be estimated from this adult Allen matrix.

## Interpretation

The result separates same-direction candidate recurrence from a shared multigene configuration. Tier 1 is enriched for genes that are positive in both local granule-versus-comparator contrasts, consistent with the existing branch-local analysis. However, neither classifier transfer nor the relative cross-region contrast pattern exceeded the matched configuration null. The full downstream neurite/synaptic panel was anti-aligned across mature regional contrasts. The current data therefore support limited candidate recurrence but do not establish a granule-cell-specific molecular configuration or a causal explanation of morphology.

## Limits

- Both directions use one adult mouse platform; anatomical region and comparator identity remain linked.
- Candidate tiers were selected in the primary discovery datasets, so the Allen transfer is external validation of a prespecified panel, not de novo discovery.
- The analysis tests the populations represented here and cannot establish uniqueness relative to every neuronal class.
- Morphology is not measured in the expression matrix.

## Outputs

- `Project/results/dgd_allen_cross_region_transfer_metrics.tsv`
- `Project/results/dgd_allen_cross_region_specificity_summary.tsv`
- `Project/results/dgd_allen_cross_region_label_permutations.tsv.gz`
- `Project/results/dgd_allen_cross_region_matched_gene_null.tsv.gz`
- `Project/results/dgd_allen_cross_region_feature_weights.tsv`
- `Project/results/dgd_allen_cross_region_stage_sensitivity.tsv`
- `Project/results/dgd_allen_cross_region_contrast_concordance.tsv`
- `Project/results/dgd_allen_cross_region_gene_contrasts.tsv`
- `Project/results/dgd_allen_cross_region_transfer.png`
- `Project/results/dgd_allen_cross_region_transfer.pdf`
