# Aim 2 Niche/Pathway Signaling Audit

Date built: 2026-06-22

## Purpose

This analysis addresses Specific Aim 2 by scoring targeted niche, ligand-receptor, and pathway-readiness programs in the primary-core pseudobulk expression layers.

The analysis is intentionally conservative: the available primary-core matrices are broad-class expression summaries, not spatial sender-receiver assays. Therefore the results test pathway/module readiness and candidate-versus-background enrichment, not direct cell-cell communication.

## Inputs

- Full MGI one-to-one ortholog expression layer.
- 2,169-gene selected-feature expression layer.
- Curated pathway modules for TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP/SMAD, Reelin, Semaphorin, SHH, WNT, FGF, and Notch.
- Curated ligand-receptor pairs centered on the 2005 TGF-beta2/BDNF mechanism and related niche pathways.

## Scale

- Pathway class units: 1,878 across 10 datasets.
- Pathway candidate-background contrasts: 566.
- Signature class units: 840.
- Signature candidate-background contrasts: 252.
- Ligand-receptor readiness units: 1,709.
- Ligand-receptor candidate-background contrasts: 285.

## Main Results

- Pathway modules were candidate-enriched in 336/566 contrasts (median delta 0.500, Wilcoxon p=2.61e-57).
- Composite signatures were candidate-enriched in 172/252 contrasts (median delta 0.250, Wilcoxon p=2.11e-24).
- Ligand-receptor readiness pairs were candidate-enriched in 139/285 contrasts (median delta 0.000, Wilcoxon p=0.0106).
- The TGF-beta/BDNF 2005 mechanism index was positive in 48/63 contrasts (median delta 0.500, Wilcoxon p=1.87e-09).
- Layer caveat: the composite-signature signal is driven mainly by the selected-feature matrix (162/220 positive, median delta 0.300); the full MGI layer is not broadly positive (10/32 positive, median delta -0.014).

Branch-level signature results:
- Dentate differentiation/stop: 48/51 positive, median delta 0.500, Wilcoxon p=2.32e-09.
- Cerebellar differentiation/stop: 2/12 positive, median delta -0.067, Wilcoxon p=0.847.
- Dentate neurogenic/permissive: 42/51 positive, median delta 0.500, Wilcoxon p=1.22e-07.
- Cerebellar neurogenic/permissive: 3/12 positive, median delta 0.000, Wilcoxon p=0.902.
- Dentate stop-minus-permissive: 23/51 positive, median delta 0.000, Wilcoxon p=3.6e-05.
- Cerebellar stop-minus-permissive: 6/12 positive, median delta 0.033, Wilcoxon p=0.374.

## Aim 2 Hypothesis Test Outcome

The starting prediction was that cerebellar granule-cell context might show stronger differentiation/stop signaling, while the dentate context might retain a more persistent neurogenic/permissive signature.
- Cerebellar candidates do not show broad differentiation/stop enrichment in this pseudobulk audit: differentiation/stop 2/12 positive, median delta -0.067; TGF-beta/BDNF index 1/12 positive, median delta -0.125.
- The clearest cerebellar pathway signal is SHH/PTCH/GLI: 10/12 positive, median delta 0.500.
- Dentate candidates show broad pathway readiness: differentiation/stop 48/51 positive, neurogenic/permissive 42/51 positive, and TGF-beta/BDNF index 47/51 positive.
- The TGF-beta/SMAD pathway itself is dentate-enriched rather than cerebellar-enriched in candidate-background contrasts: dentate median delta 0.500, cerebellar median delta -0.417.
Conclusion: Aim 2 supports context-dependent niche/maturation signaling, but not the simple version of a cerebellar-biased TGF-beta/BDNF stop-signaling program.

Pathway-level highlights:
- bdnf_trkb_mapk / dentate: median delta 0.500; 37/51 positive.
- bmp_smad / dentate: median delta 0.500; 42/51 positive.
- tgf_beta_smad / dentate: median delta 0.500; 46/51 positive.
- semaphorin_guidance / dentate: median delta 0.500; 44/51 positive.
- reelin_migration_stop / dentate: median delta 0.500; 31/51 positive.
- notch_hes / dentate: median delta 0.500; 26/51 positive.
- fgf_mapk / dentate: median delta 0.500; 37/51 positive.
- wnt_beta_catenin / dentate: median delta 0.500; 41/51 positive.
- shh_granule_expansion / cerebellar: median delta 0.500; 10/12 positive.
- notch_hes / cerebellar: median delta 0.250; 6/12 positive.

Minimum gene coverage across expression layers:
- Notch/HES: 1 genes present.
- SHH/PTCH/GLI: 1 genes present.
- TGF-beta/SMAD: 1 genes present.
- WNT/beta-catenin: 1 genes present.
- BMP/SMAD: 3 genes present.
- FGF/MAPK: 3 genes present.
- Reelin/migration-stop: 3 genes present.
- BDNF/TrkB/MAPK: 7 genes present.
- Semaphorin/plexin guidance: 9 genes present.

## Interpretation

- Aim 2 is now computably tested at the pathway-readiness level.
- The historical TGF-beta2/BDNF/SMAD-MAPK mechanism can be discussed as a prioritized niche/maturation axis, but in this primary-core audit it behaves as context-dependent and stronger in the dentate branch than in cerebellar candidate-background contrasts.
- SHH/PTCH/GLI is the strongest cerebellar pathway signal, consistent with cerebellar granule precursor biology.
- The key manuscript-safe interpretation should distinguish pathway readiness from true spatial niche signaling.
- If stronger Aim 2 evidence is required, the next layer should use spatial datasets, ligand-sender/receiver cell labels, or raw object-level cell-type communication tools.

## Outputs

- Pathway gene sets: `Project/results/primary_core_aim2_pathway_gene_sets.tsv`
- Pathway units: `Project/results/primary_core_aim2_pathway_units.tsv.gz`
- Pathway contrasts: `Project/results/primary_core_aim2_pathway_contrasts.tsv`
- Pathway summary: `Project/results/primary_core_aim2_pathway_summary.tsv`
- Pathway coverage: `Project/results/primary_core_aim2_pathway_coverage.tsv`
- Signature units: `Project/results/primary_core_aim2_signature_units.tsv.gz`
- Signature contrasts: `Project/results/primary_core_aim2_signature_contrasts.tsv`
- Signature summary: `Project/results/primary_core_aim2_signature_summary.tsv`
- Ligand-receptor pairs: `Project/results/primary_core_aim2_ligand_receptor_pairs.tsv`
- Ligand-receptor units: `Project/results/primary_core_aim2_ligand_receptor_units.tsv.gz`
- Ligand-receptor contrasts: `Project/results/primary_core_aim2_ligand_receptor_contrasts.tsv`
- Ligand-receptor summary: `Project/results/primary_core_aim2_ligand_receptor_summary.tsv`
- Plot: `Project/results/primary_core_aim2_niche_pathway_model.png`
