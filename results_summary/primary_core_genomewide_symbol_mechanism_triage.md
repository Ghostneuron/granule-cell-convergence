# Genome-Wide Same-Symbol Mechanism Triage

Date built: 2026-06-22

## Purpose

This triage classifies the full-matrix same-symbol pseudobulk hits into likely morphology/wiring mechanisms, broad neuronal state, metabolic/supporting genes, RNA-processing state, and regional-identity warnings.

## Summary

- Shared-positive same-symbol hits: 6,440.
- Mechanism figure candidates: 33.
- Mechanism follow-up candidates: 261.
- Original 67-gene packet genes recovered among shared-positive hits: 29.

## Shared-Positive Class Counts

- `unclassified_shared_neuronal_or_context` / `exploratory_context`: 5830 genes.
- `synaptic_wiring` / `mechanism_figure_candidate`: 17 genes.
- `cytoskeleton_morphogenesis` / `mechanism_figure_candidate`: 8 genes.
- `regulatory_morphogenesis_candidate` / `mechanism_figure_candidate`: 6 genes.
- `axon_guidance_adhesion` / `mechanism_figure_candidate`: 1 genes.
- `curated_shared_structural_executor` / `mechanism_figure_candidate`: 1 genes.
- `cytoskeleton_morphogenesis` / `mechanism_followup_candidate`: 123 genes.
- `synaptic_wiring` / `mechanism_followup_candidate`: 76 genes.
- `regulatory_morphogenesis_candidate` / `mechanism_followup_candidate`: 36 genes.
- `axon_guidance_adhesion` / `mechanism_followup_candidate`: 16 genes.
- `curated_shared_structural_executor` / `mechanism_followup_candidate`: 10 genes.
- `regional_identity_or_specificity_warning` / `specificity_warning_or_branch_context`: 2 genes.
- `metabolic_or_housekeeping_support` / `supporting_context_not_core_executor`: 227 genes.
- `rna_processing_neuronal_state` / `supporting_context_not_core_executor`: 78 genes.
- `broad_neuronal_state` / `supporting_context_not_core_executor`: 9 genes.

## Top Mechanism Figure Candidates

- `ROBO2` (curated_shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500, original_packet=True.
- `GPM6A` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `NFIB` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `NFIA` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `PPP3CA` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `CAMTA1` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `MAPK1` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `STXBP1` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `CACNA2D1` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `SYNPR` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCNC1` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCNK1` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `MAPK8IP2` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `GABRB3` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `BCL7A` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `ADD2` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCNJ6` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `RFX3` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCNJ3` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `TUBGCP2` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `MAPK14` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `FOXN2` (regulatory_morphogenesis_candidate): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `GRIN2B` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `MAP3K4` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCND3` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `CALM2` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `KCND2` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `GABRA2` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `STXBP5L` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `MAP2K7` (cytoskeleton_morphogenesis): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `CACNA1E` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `SEMA7A` (axon_guidance_adhesion): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.
- `STXBP5` (synaptic_wiring): dentate delta 0.500, cerebellar delta 0.500, original_packet=False.

## Interpretation

- The full-matrix screen recovers broad neuronal and morphogenesis-related programs across a much larger symbol universe.
- The mechanism figure list should be treated as a prioritization table. It still needs curated ortholog mapping and a model with dataset/sample effects before final manuscript claims.
- Genes that also appeared in the selected-gene triage are the best near-term candidates because they survive two independent feature-universe definitions.

## Output

- Triage table: `Project/results/primary_core_genomewide_symbol_mechanism_triage.tsv`
