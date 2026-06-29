# Mechanism Axis Model

Date built: 2026-06-22

## Purpose

This report translates the formal manuscript candidate tiers into biological mechanism axes. The goal is to support a paper argument that dentate and cerebellar granule cells do not share regional identity, but do converge on a downstream toolkit for compact neuronal morphology, neurite patterning, and excitability maturation.

## Axis Definitions

- Developmental regulatory control: transcriptional or regulatory candidates that may help coordinate granule-cell maturation programs.
- Neurite/cytoskeleton morphogenesis: membrane, cytoskeletal, neurite-growth, and structural-executor candidates.
- Axon guidance and adhesion: genes plausibly linking compact soma/neurite morphology to wiring and local circuit integration.
- Synaptic/excitability maturation: receptor, ion-channel, calcium/signaling, and synaptic-release candidates.
- Exploratory ortholog completeness: non-identical mouse/human ortholog hits retained for follow-up, not central claims.

## Central Candidate Model

- Developmental regulatory control: `NFIB`, `NFIA`, `RFX3`.
- Neurite/cytoskeleton morphogenesis: `GPM6A`.
- Axon guidance and adhesion: `ROBO2`.
- Synaptic/excitability maturation: `KCNK1`, `GABRA2`, `PPP3CA`, `CACNA2D1`, `KCNJ6`, `GABRB3`, `GRIN2B`, `KCNJ3`, `KCND2`, `STXBP5L`.

## Tier 1 Seed Interpretation

- `NFIB` (Developmental regulatory control): developmental transcriptional regulator; 4/4 formal branches.
- `NFIA` (Developmental regulatory control): developmental transcriptional regulator; 4/4 formal branches.
- `RFX3` (Developmental regulatory control): ciliogenesis/transcriptional regulatory candidate; 4/4 formal branches.
- `GPM6A` (Neurite/cytoskeleton morphogenesis): membrane/neurite outgrowth structural executor; 4/4 formal branches.
- `KCNK1` (Synaptic/excitability maturation): ion-channel/excitability tuning candidate; 4/4 formal branches.
- `GABRA2` (Synaptic/excitability maturation): GABA receptor and synaptic maturation candidate; 4/4 formal branches.

## Tier 2 Pathway Support

- `ROBO2` (Axon guidance and adhesion): axon guidance and neurite-patterning executor; 4/4 formal branches.
- `PPP3CA` (Synaptic/excitability maturation): calcineurin/synaptic plasticity signaling candidate; 4/4 formal branches.
- `CACNA2D1` (Synaptic/excitability maturation): calcium-channel auxiliary subunit and wiring candidate; 4/4 formal branches.
- `KCNJ6` (Synaptic/excitability maturation): inward-rectifier potassium-channel candidate; 4/4 formal branches.
- `GABRB3` (Synaptic/excitability maturation): GABA receptor and synaptic maturation candidate; 4/4 formal branches.
- `GRIN2B` (Synaptic/excitability maturation): glutamatergic synapse maturation candidate; 4/4 formal branches.
- `KCNJ3` (Synaptic/excitability maturation): inward-rectifier potassium-channel candidate; 4/4 formal branches.
- `KCND2` (Synaptic/excitability maturation): voltage-gated potassium-channel candidate; 4/4 formal branches.
- `STXBP5L` (Synaptic/excitability maturation): synaptic vesicle/exocytosis regulatory candidate; 4/4 formal branches.

## Axis Summary

- Developmental regulatory control: 7 genes (Tier 1=3, Tier 2=0, Tier 3/4=4, Tier 5=0).
- Neurite/cytoskeleton morphogenesis: 13 genes (Tier 1=1, Tier 2=0, Tier 3/4=12, Tier 5=0).
- Axon guidance and adhesion: 4 genes (Tier 1=0, Tier 2=1, Tier 3/4=3, Tier 5=0).
- Synaptic/excitability maturation: 12 genes (Tier 1=2, Tier 2=8, Tier 3/4=2, Tier 5=0).
- Exploratory ortholog completeness: 30 genes (Tier 1=0, Tier 2=0, Tier 3/4=0, Tier 5=30).

## Interpretation

- The central model is not a single pathway. It is a layered toolkit: regulatory coordination, neurite/cytoskeleton execution, axon-guidance/adhesion, and synaptic/excitability maturation.
- Tier 1 should anchor the manuscript model because it is both compact and robust across all available formal branches.
- Tier 2 provides the strongest pathway-level biological breadth, especially synaptic/excitability maturation and ROBO2-mediated guidance.
- Tier 5 should stay exploratory because non-identical ortholog recovery improves completeness but does not yet have the same manuscript-readiness as the curated mechanism tiers.

## Outputs

- Gene-axis table: `Project/results/primary_core_mechanism_axis_gene_table.tsv`
- Axis summary: `Project/results/primary_core_mechanism_axis_summary.tsv`
- Branch summary: `Project/results/primary_core_mechanism_axis_branch_summary.tsv`
- Plot: `Project/results/primary_core_mechanism_axis_model.png`
