# Cerebellar Conditioned-Medium Secretome Candidate Inference

## Question

Can current sequencing data nominate additional cerebellar conditioned-medium factors that might inhibit proliferation or promote differentiation besides TGF-beta2 and BDNF?

## Answer

Yes, partially. The data can nominate secreted or ligand-like genes with cerebellar granule-cell source plausibility and developmental timing, but it cannot prove that the proteins are secreted into medium, properly processed, abundant, or bioactive. This screen should therefore be treated as a prioritized validation list for proteomics, ELISA/Luminex, neutralization, or recombinant-factor rescue experiments.

## Historical Anchors

- `TGFB2` (TGF-beta superfamily): score 3.646, candidate detection 0.146, candidate-background rank delta 0.000, P8-P0 rank delta 0.667. 2005 anchor; SMAD-linked stop/maturation factor.
- `BDNF` (neurotrophin): score 2.173, candidate detection 0.006, candidate-background rank delta 0.000, P8-P0 rank delta 0.333. 2005 anchor; TrkB/MAPK and SMAD-convergent factor.
- `NGF` (neurotrophin): score 1.500, candidate detection 0.000, candidate-background rank delta -0.500, P8-P0 rank delta NA. tested in 2005 but not main neutralized activity.

`TGFB2` is well supported by the present cerebellar transcriptomic data, especially by the P0-to-P8 rise in `GSE122357`. `BDNF` remains experimentally important from the 2005 paper, but its mRNA is sparse in these sequencing tables; that means transcriptomics alone would under-prioritize it relative to antibody/functional evidence.

## Highest-Priority Candidates

- `SFRP1` (WNT antagonist): score 4.826, candidate detection 0.326, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. secreted frizzled-related candidate.
- `BMP6` (BMP/GDF superfamily): score 4.747, candidate detection 0.247, candidate-background rank delta 0.500, P8-P0 rank delta NA. BMP differentiation/neurogenic brake candidate.
- `RELN` (migration-stop/guidance): score 4.739, candidate detection 0.739, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. Reelin migration-stop/maturation candidate.
- `TGFB2` (TGF-beta superfamily): score 3.646, candidate detection 0.146, candidate-background rank delta 0.000, P8-P0 rank delta 0.667. 2005 anchor; SMAD-linked stop/maturation factor.
- `SFRP2` (WNT antagonist): score 3.248, candidate detection 0.081, candidate-background rank delta 0.000, P8-P0 rank delta 0.333. secreted frizzled-related candidate.

## Supported Candidates

- `GDF11` (BMP/GDF superfamily): score 2.902, candidate detection 0.068, candidate-background rank delta 0.000, P8-P0 rank delta -0.333. TGF-beta family differentiation candidate.
- `INHBB` (activin/inhibin): score 2.861, candidate detection 0.028, candidate-background rank delta 0.000, P8-P0 rank delta -0.667. activin branch candidate.
- `TGFB3` (TGF-beta superfamily): score 2.524, candidate detection 0.024, candidate-background rank delta -0.250, P8-P0 rank delta 0.000. TGF-beta family candidate.
- `SEMA3A` (semaphorin guidance): score 2.415, candidate detection 0.082, candidate-background rank delta 0.000, P8-P0 rank delta -0.333. secreted semaphorin guidance candidate.
- `SLIT2` (slit/robo guidance): score 2.402, candidate detection 0.069, candidate-background rank delta 0.000, P8-P0 rank delta -0.167. secreted guidance candidate.

## Detected Context or Counter-Signal Candidates

Some secreted factors are detected but are more likely to be permissive, survival, proliferative, or context-dependent rather than direct stop factors. These matter because conditioned medium can contain mixed activities.

- `NTF3` (neurotrophin): score 3.570, candidate detection 0.070, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. neurotrophin maturation candidate.
- `CXCL12` (chemokine): score 2.708, candidate detection 0.208, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. chemokine niche candidate.
- `NTN1` (netrin guidance): score 2.630, candidate detection 0.130, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. secreted guidance candidate.
- `GDF10` (BMP/GDF superfamily): score 2.548, candidate detection 0.048, candidate-background rank delta -0.500, P8-P0 rank delta 0.000. TGF-beta family differentiation candidate.
- `DKK3` (WNT antagonist): score 2.545, candidate detection 0.045, candidate-background rank delta -0.500, P8-P0 rank delta 0.000. WNT-antagonist candidate.
- `PTN` (growth factor): score 2.377, candidate detection 0.877, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. growth/survival counter-signal.
- `FGF9` (FGF pathway): score 1.782, candidate detection 0.282, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. growth/permissive counter-signal.
- `SPARC` (matricellular/ECM): score 1.602, candidate detection 0.102, candidate-background rank delta -0.500, P8-P0 rank delta 0.000. secreted ECM candidate.
- `SPARCL1` (matricellular/ECM): score 1.593, candidate detection 0.093, candidate-background rank delta -0.500, P8-P0 rank delta 0.000. secreted ECM candidate.
- `FGF8` (FGF pathway): score 1.527, candidate detection 0.027, candidate-background rank delta 0.500, P8-P0 rank delta 0.000. growth/permissive counter-signal.
- `MDK` (growth factor): score 1.244, candidate detection 0.577, candidate-background rank delta 0.250, P8-P0 rank delta -0.167. growth/survival counter-signal.
- `WNT5A` (WNT pathway): score 0.534, candidate detection 0.034, candidate-background rank delta -0.500, P8-P0 rank delta 0.000. WNT pathway candidate.

## Interpretation

The most manuscript-useful conclusion is that TGF-beta2 and BDNF remain supported historical anchors, but the modern data make a broader peptide/secreted-factor model testable. The strongest new classes to consider are BMP/GDF/activin-family factors, secreted migration-stop/guidance cues such as Reelin/Semaphorin/Slit, WNT antagonists, and matricellular proteins that could shift proliferating hippocampal progenitors toward differentiation or altered niche adhesion.

Because the 2005 assay concentrated factors above 6 kDa and found multiple chromatographic activities, a multi-factor medium model is more plausible than a single-factor replacement model. Sequencing should be used to prioritize neutralization/proteomics, not as final proof.

## Outputs

- Unit table: `Project/results/cerebellar_conditioned_medium_secretome_units.tsv` (371 candidate contrasts after unit extraction).
- Ranked candidates: `Project/results/cerebellar_conditioned_medium_secretome_ranked_candidates.tsv` (68 genes).
- Plot: `Project/results/cerebellar_conditioned_medium_secretome_candidates.png`.
