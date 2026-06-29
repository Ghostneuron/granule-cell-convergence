# Public Perturbation Dataset Triage

This curation asks whether public datasets can move the project from a hierarchical evidence model toward causal or quasi-causal mixed-effects tests. The answer is uneven: `SHH`, `RBFOX3`, and `NFIA` have the most useful perturbation resources; `BDNF` and TGF-beta/SMAD have usable pathway-response resources; `HMGN2` currently has weak public perturbation support.

## Node Summary

| Node | Candidate datasets | Top accession | Priority call | Main caveat |
|---|---:|---|---|---|
| NFIA | 3 | GSE146793 | strong regulatory follow-up candidate | direct normal dentate granule-cell NFIA perturbation remains missing |
| BDNF | 4 | GSE242199 | usable pathway perturbation candidate | available datasets are neural progenitor, hippocampal injury/plasticity, or indirect cerebellar chromatin contexts |
| TGFBeta | 2 | GSE307973 | supporting secreted-factor/niche candidate | few direct granule-cell TGF-beta perturbation transcriptomes were found |
| TGFBeta_SMAD | 1 | GSE317025 | supporting hindbrain-patterning candidate | dual-SMAD patterning is not equivalent to TGF-beta ligand response |
| SHH | 5 | GSE260623 | best causal extension candidate | many datasets are medulloblastoma or tumor-susceptibility models rather than normal granule morphology |
| SHH_Niche | 1 | GSE281480 | supporting niche-response candidate | not granule lineage |
| RBFOX3 | 3 | GSE84786 | strong postmitotic functional follow-up candidate | best Rbfox3 KO RNA-seq has low sample count; Rbfox1/3 knockdown datasets are not RBFOX3-specific |
| HMGN2 | 2 | GSE186384 | weak public perturbation support | no clean neural or granule-lineage HMGN2-focused perturbation dataset was found |

## Recommended Use

1. Use `SHH`/cerebellar GNP perturbation resources first for the clearest causal extension of the cerebellar branch.
2. Use `RBFOX3` datasets as postmitotic dentate/hippocampal functional perturbation support, especially synaptic and plasticity modules.
3. Use `NFIA` datasets as regulatory-target evidence, strongest in cerebellar GNPs and weaker but useful in hippocampal astrocyte/niche context.
4. Treat `BDNF` and TGF-beta resources as pathway-response signatures rather than direct granule morphology perturbations.
5. Keep `HMGN2` as a hypothesis-generating chromatin-competence candidate until focused neural perturbation data are found or generated.

Detailed dataset table: `causal_perturbation_dataset_triage.tsv`
Node summary table: `causal_perturbation_node_summary.tsv`
