# Causal Perturbation Module-Shift Layer

This analysis converts public perturbation resources into preliminary module-shift signatures. The scored modules are the same curated fate, niche, pathway, neurite, synaptic, and excitability modules used in the primary manuscript analyses.

## Processing Status

- Processed accessions: GSE84786, GSE71916, GSE242199, GSE81962
- Deferred accessions: GSE107252
- Gene-level module effects: 1687
- Contrast-module summaries: 154

## Strongest Matched-Context Node-Module Shifts

| Node | Module | Direction | Median signed score | Contrasts | Context | Strongest contrast |
|---|---|---|---:|---:|---|---|
| RBFOX3/RBFOX1 | Downstream neurite/morphology | up_shift | 0.111 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| RBFOX3/RBFOX1 | Downstream synaptic/excitability | down_shift | -0.107 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| SHH/PTCH/Norrin | Downstream synaptic/excitability | up_shift | 0.104 | 4 | matched_context | GSE81962_Ptch_het_GNP_vs_WT_GNP |
| BDNF/NTRK2 | TGF-beta/SMAD | up_shift | 0.100 | 5 | matched_context | GSE242199_NTRK2ko_0h_vs_WT_0h |
| RBFOX3 | Semaphorin/plexin guidance | near_zero | 0.062 | 1 | matched_context | GSE84786_Rbfox3_KO_vs_WT |
| RBFOX3 | FGF/MAPK | near_zero | -0.057 | 1 | matched_context | GSE84786_Rbfox3_KO_vs_WT |
| RBFOX3/RBFOX1 | Semaphorin/plexin guidance | near_zero | 0.054 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| RBFOX3/RBFOX1 | Reelin/migration-stop | near_zero | 0.051 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| BDNF/NTRK2 | Downstream synaptic/excitability | near_zero | 0.038 | 5 | matched_context | GSE242199_NTRK2ko_3h_BDNF_vs_vehicle |
| SHH/PTCH/Norrin | Shared neurogenic niche/progenitor state | near_zero | -0.027 | 4 | matched_context | GSE81962_Ptch_het_GNP_vs_WT_GNP |
| RBFOX3/RBFOX1 | BDNF/TrkB/MAPK | near_zero | -0.023 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| RBFOX3 | BDNF/TrkB/MAPK | near_zero | -0.021 | 1 | matched_context | GSE84786_Rbfox3_KO_vs_WT |
| RBFOX3/RBFOX1 | FGF/MAPK | near_zero | -0.019 | 1 | matched_context | GSE71916_siRbfox1_3_vs_siNT |
| RBFOX3 | Reelin/migration-stop | near_zero | 0.017 | 1 | matched_context | GSE84786_Rbfox3_KO_vs_WT |
| SHH/PTCH/Norrin | Downstream neurite/morphology | near_zero | 0.017 | 4 | matched_context | GSE81962_Ndp_KO_GNP_vs_WT_GNP |
| BDNF/NTRK2 | FGF/MAPK | near_zero | -0.016 | 5 | matched_context | GSE242199_NTRK2ko_0h_vs_WT_0h |

Off-context shifts are retained in the tables as exploratory signals, but the manuscript interpretation should emphasize matched-context shifts.

## Strongest Matched-Context Contrast-Level Shifts

| Contrast | Node | Module | Direction | Signed score | Available genes | Top shifted genes |
|---|---|---|---|---:|---:|---|
| GSE242199_NTRK2ko_0h_vs_WT_0h | BDNF/NTRK2 | TGF-beta/SMAD | up_shift | 0.566 | 3 | ID1;SMAD3;ID2 |
| GSE242199_NTRK2ko_0h_vs_WT_0h | BDNF/NTRK2 | Downstream neurite/morphology | down_shift | -0.202 | 21 | DCC;NRXN1;STMN2;CADM3;PLXNA4 |
| GSE242199_NTRK2ko_3h_BDNF_vs_vehicle | BDNF/NTRK2 | Downstream synaptic/excitability | up_shift | 0.173 | 20 | KCNJ6;SLC17A7;STXBP5L;GABRA2;SNAP25 |
| GSE242199_NTRK2ko_3h_BDNF_vs_vehicle | BDNF/NTRK2 | TGF-beta/SMAD | up_shift | 0.168 | 3 | ID1;SMAD3;ID2 |
| GSE242199_NTRK2ko_0h_vs_WT_0h | BDNF/NTRK2 | Downstream synaptic/excitability | up_shift | 0.159 | 20 | SLC17A6;CAMK2B;GRIN2B;GABRA2;SYNPR |
| GSE242199_NTRK2ko_0h_vs_WT_0h | BDNF/NTRK2 | FGF/MAPK | up_shift | 0.149 | 4 | FGF8;MAP2K1;FRS2;MAPK1 |
| GSE242199_NTRK2ko_3h_BDNF_vs_vehicle | BDNF/NTRK2 | BDNF/TrkB/MAPK | down_shift | -0.148 | 8 | BDNF;SOS1;GRB2;MAPK1;NTRK2 |
| GSE81962_Ptch_het_GNP_vs_WT_GNP | SHH/PTCH/Norrin | Downstream synaptic/excitability | up_shift | 0.142 | 19 | SLC17A7;STXBP1;KCND2;PPP3CA;CACNA2D1 |
| GSE81962_Ndp_KO_GNP_vs_WT_GNP | SHH/PTCH/Norrin | Downstream synaptic/excitability | up_shift | 0.114 | 19 | SLC17A7;KCND2;SYT1;STXBP1;PPP3CA |
| GSE242199_NTRK2ko_3h_BDNF_vs_vehicle | BDNF/NTRK2 | Downstream neurite/morphology | up_shift | 0.111 | 21 | GAP43;PLXNA4;MAP1B;DCC;NRXN1 |
| GSE71916_siRbfox1_3_vs_siNT | RBFOX3/RBFOX1 | Downstream neurite/morphology | up_shift | 0.111 | 20 | NRP1;DPYSL3;PLXNA4;RTN1;STMN2 |
| GSE71916_siRbfox1_3_vs_siNT | RBFOX3/RBFOX1 | Downstream synaptic/excitability | down_shift | -0.107 | 20 | KCNJ3;SNAP25;SLC17A6;CACNA2D1;PPP3CA |
| GSE242199_WT_3h_BDNF_vs_vehicle | BDNF/NTRK2 | TGF-beta/SMAD | up_shift | 0.100 | 3 | ID1;SMAD3;ID2 |
| GSE81962_Ndp_KO_Ptch_het_GNP_vs_WT_GNP | SHH/PTCH/Norrin | Downstream synaptic/excitability | near_zero | 0.095 | 19 | KCND2;PPP3CA;STXBP1;SYT1;SLC17A7 |
| GSE81962_Ndp_KO_GNP_vs_WT_GNP | SHH/PTCH/Norrin | Cerebellar fate/rhombic-lip/SHH | near_zero | 0.082 | 13 | GLI2;EN1;BARHL1;ZIC3;MEIS1 |
| GSE84786_Rbfox3_KO_vs_WT | RBFOX3 | Semaphorin/plexin guidance | near_zero | 0.062 | 17 | NRP2;PLXNA1;SEMA3A;PLXNA4;SEMA3E |

## Interpretation

The output should be read as perturbation-sensitive module evidence, not as a causal mixed-effects model. It is strongest when a contrast is lineage-relevant, has adequate module-gene coverage, and shifts biologically matched modules in the expected branch. The current layer is therefore useful for prioritizing follow-up perturbations and for supporting the hierarchical integrative model, while still requiring matched single-cell or perturb-seq data for formal causality.

## Output Files

- Module catalog: `causal_perturbation_module_catalog.tsv`
- Gene-level effects: `causal_perturbation_module_shift_gene_effects.tsv`
- Contrast-module summary: `causal_perturbation_module_shift_summary.tsv`
- Node-module summary: `causal_perturbation_module_shift_node_summary.tsv`
- Processing status: `causal_perturbation_processing_status.tsv`
- Heatmap: `causal_perturbation_module_shift_heatmap.png`
