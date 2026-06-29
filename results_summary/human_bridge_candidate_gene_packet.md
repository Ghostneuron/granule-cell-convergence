# Human Bridge Candidate Gene Packet

Date built: 2026-06-21

## Scope

This packet summarizes refined marker-panel genes across the constructed human dentate/hippocampal bridge objects: `human_core_tuned`, `GSE325391`, and `GSE268609`.

It is a first candidate table for manuscript planning. It does not replace full differential-expression, ortholog-aware cross-species modeling, or source-taxonomy refinement for `GSE268609`.

## Integrated Context

- `dentate_candidate_vs_cerebellar_candidate` / `identity_rank_contrast`: median delta 0.6472, BH-adjusted p=1.57e-10.
- `dentate_candidate_vs_non_dentate_background` / `identity_rank_contrast`: median delta 0.3423, BH-adjusted p=9.93e-22.
- `dentate_candidate_vs_cerebellar_candidate` / `structural_rank`: median delta -0.0965, BH-adjusted p=0.435.
- `backbone_refined` / `cerebellar_candidate`: structural-rank median 0.6831.
- `backbone_refined` / `dentate_candidate`: structural-rank median 0.6861.
- `gse268609_hippocampus_rna` / `dentate_candidate`: structural-rank median 0.5487.
- `gse325391_adult_dg` / `dentate_candidate`: structural-rank median 0.3991.
- `human_core_tuned` / `dentate_candidate`: structural-rank median 0.7818.

## High-Priority Shared Structural Executors

- `Ncam1` (morphogenesis_cytoskeleton): dentate detection median 0.904; 3 human bridge sources at >=20% detection.
- `Mapt` (morphogenesis_cytoskeleton): dentate detection median 0.827; 3 human bridge sources at >=20% detection.
- `Plxna4` (axon_guidance_synapse): dentate detection median 0.827; 3 human bridge sources at >=20% detection.
- `Dpysl2` (morphogenesis_cytoskeleton): dentate detection median 0.777; 3 human bridge sources at >=20% detection.
- `Cntn5` (axon_guidance_synapse): dentate detection median 0.773; 3 human bridge sources at >=20% detection.
- `Robo2` (axon_guidance_synapse): dentate detection median 0.768; 3 human bridge sources at >=20% detection.
- `Nrp1` (axon_guidance_synapse): dentate detection median 0.744; 3 human bridge sources at >=20% detection.
- `Epha4` (axon_guidance_synapse): dentate detection median 0.711; 3 human bridge sources at >=20% detection.
- `Robo1` (axon_guidance_synapse): dentate detection median 0.688; 3 human bridge sources at >=20% detection.
- `Plxna2` (axon_guidance_synapse): dentate detection median 0.563; 3 human bridge sources at >=20% detection.
- `Slit2` (axon_guidance_synapse): dentate detection median 0.540; 3 human bridge sources at >=20% detection.
- `Stmn2` (morphogenesis_cytoskeleton): dentate detection median 0.501; 3 human bridge sources at >=20% detection.
- `Gap43` (morphogenesis_cytoskeleton): dentate detection median 0.498; 3 human bridge sources at >=20% detection.
- `Nrp2` (axon_guidance_synapse): dentate detection median 0.478; 3 human bridge sources at >=20% detection.
- `Ephb1` (axon_guidance_synapse): dentate detection median 0.435; 3 human bridge sources at >=20% detection.
- `Cfl1` (morphogenesis_cytoskeleton): dentate detection median 0.364; 3 human bridge sources at >=20% detection.
- `Cdk5r1` (morphogenesis_cytoskeleton): dentate detection median 0.346; 3 human bridge sources at >=20% detection.
- `Elavl4` (morphogenesis_cytoskeleton): dentate detection median 0.333; 3 human bridge sources at >=20% detection.
- `Cntn4` (axon_guidance_synapse): dentate detection median 0.322; 3 human bridge sources at >=20% detection.
- `Ephb2` (axon_guidance_synapse): dentate detection median 0.292; 3 human bridge sources at >=20% detection.

## Regional Identity Examples

Dentate identity panel:
- `Prox1`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.706.
- `Calb1`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.401.
- `C1ql3`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.392.
- `Glis3`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.315.
- `Egr3`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.136.
- `Itpka`: tier `dentate_identity_supported_in_human_bridge`, dentate detection median 0.131.
- `Calb2`: tier `dentate_identity_panel_gene_needs_context`, dentate detection median 0.015.
- `Mex3a`: tier `dentate_identity_panel_gene_needs_context`, dentate detection median 0.006.

Cerebellar identity panel:
- `Ror1`: tier `cerebellar_identity_gene_with_human_bridge_leakage_warning`, human dentate-bridge detection median 0.283.
- `Etv1`: tier `cerebellar_identity_gene_with_human_bridge_leakage_warning`, human dentate-bridge detection median 0.238.
- `Gabra6`: tier `cerebellar_identity_gene_with_human_bridge_leakage_warning`, human dentate-bridge detection median 0.144.
- `Grin2c`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.006.
- `Zic1`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.004.
- `Zic2`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.001.
- `Cbln3`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.000.
- `Barhl2`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.000.
- `Barhl1`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.000.
- `En2`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.000.
- `Atoh1`: tier `clean_low_human_dentate_detection_cerebellar_identity`, human dentate-bridge detection median 0.000.

## Interpretation

- The strongest manuscript direction remains identity separation plus structural convergence, not a claim that dentate and cerebellar granule cells are transcriptionally identical.
- High-priority structural/executor candidates are detected across multiple human dentate/hippocampal bridge sources and belong to morphogenesis, cytoskeletal, adhesion, guidance, or synaptic panels.
- Cerebellar identity genes with high human dentate-bridge detection should be treated as specificity warnings, not as evidence of cerebellar identity in DG cells.

## Outputs

- Gene-level human bridge summary: `Project/results/human_bridge_marker_gene_summary.tsv`
- Candidate gene packet: `Project/results/human_bridge_candidate_gene_packet.tsv`
- Structural executor candidate plot: `Project/results/human_bridge_structural_executor_candidates.png`
