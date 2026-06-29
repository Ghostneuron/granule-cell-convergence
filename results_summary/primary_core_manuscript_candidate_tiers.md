# Manuscript Candidate Tiers

Date built: 2026-06-22

## Purpose

This packet distills the formal MGI ortholog rank-meta validation into manuscript-facing candidate tiers. It does not add a new statistical test; it organizes the formal results into seed genes, pathway-support genes, and exploratory ortholog-completeness candidates.

## Tier Counts

- Tier 1 core convergent program: 6 genes.
- Tier 2 high-confidence wiring/synaptic executor: 9 genes.
- Tier 3 broad both-screen mechanism support: 5 genes.
- Tier 4 screen-specific mechanism support: 16 genes.
- Tier 5 exploratory non-identical ortholog: 30 genes.

## Tier 1 Core Convergent Program

- `GPM6A`: membrane/neurite outgrowth structural executor; 4/4 formal branches.
- `NFIB`: developmental transcriptional regulator; 4/4 formal branches.
- `KCNK1`: ion-channel/excitability tuning candidate; 4/4 formal branches.
- `RFX3`: ciliogenesis/transcriptional regulatory candidate; 4/4 formal branches.
- `GABRA2`: GABA receptor and synaptic maturation candidate; 4/4 formal branches.
- `NFIA`: developmental transcriptional regulator; 4/4 formal branches.

## Tier 2 High-Confidence Support

- `ROBO2`: axon guidance and neurite-patterning executor; 4/4 formal branches.
- `GABRB3`: GABA receptor and synaptic maturation candidate; 4/4 formal branches.
- `KCND2`: voltage-gated potassium-channel candidate; 4/4 formal branches.
- `PPP3CA`: calcineurin/synaptic plasticity signaling candidate; 4/4 formal branches.
- `CACNA2D1`: calcium-channel auxiliary subunit and wiring candidate; 4/4 formal branches.
- `KCNJ6`: inward-rectifier potassium-channel candidate; 4/4 formal branches.
- `GRIN2B`: glutamatergic synapse maturation candidate; 4/4 formal branches.
- `KCNJ3`: inward-rectifier potassium-channel candidate; 4/4 formal branches.
- `STXBP5L`: synaptic vesicle/exocytosis regulatory candidate; 4/4 formal branches.

## Tier 3 Broad Both-Screen Support

- `DYNLL1`: cytoskeleton/neurite morphogenesis candidate.
- `BASP1`: cytoskeleton/neurite morphogenesis candidate.
- `RFX7`: regulatory morphogenesis candidate.
- `MAPKAP1`: cytoskeleton/neurite morphogenesis candidate.
- `TUBA1A`: cytoskeleton/neurite morphogenesis candidate.

## Tier 4 Screen-Specific Support

- `ACTB`: cytoskeleton/neurite morphogenesis candidate.
- `RTN3`: cytoskeleton/neurite morphogenesis candidate.
- `ACTG1`: cytoskeleton/neurite morphogenesis candidate.
- `TUBA1B`: cytoskeleton/neurite morphogenesis candidate.
- `TCF4`: regulatory morphogenesis candidate.
- `MAP1B`: cytoskeleton/neurite morphogenesis candidate.
- `STMN2`: shared structural-executor candidate.
- `CADM3`: axon-guidance or adhesion candidate.
- `KCNMB4`: synaptic wiring or excitability candidate.
- `DCC`: axon-guidance or adhesion candidate.
- `MAP3K13`: cytoskeleton/neurite morphogenesis candidate.
- `SLC17A6`: synaptic wiring or excitability candidate.
- `MAPK14`: cytoskeleton/neurite morphogenesis candidate.
- `FOXN2`: regulatory morphogenesis candidate.
- `SEMA7A`: axon-guidance or adhesion candidate.
- `BCL7A`: regulatory morphogenesis candidate.

## Tier 5 Exploratory Ortholog Completeness

- `ZNF706` / mouse `Zfp706`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `RAB7A` / mouse `Rab7`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `ZNF292` / mouse `Zfp292`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `ZNF148` / mouse `Zfp148`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `C1orf21` / mouse `1700025G04Rik`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `ZNF827` / mouse `Zfp827`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `KIAA1328` / mouse `AW554918`: formal_nominal_shared_full_matrix; 3/4 formal branches.
- `ZNF536` / mouse `Zfp536`: formal_nominal_shared_full_matrix; 2/4 formal branches.
- `TMEM178A` / mouse `Tmem178`: formal_nominal_shared_full_matrix; 2/4 formal branches.
- `TUBB` / mouse `Tubb5`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `MIR124-1HG` / mouse `Mir124a-1hg`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `C5orf34` / mouse `4833420G17Rik`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `ZNF280D` / mouse `Zfp280d`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `ZNF32` / mouse `Zfp637`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `TP53I11` / mouse `Trp53i11`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `ZNF821` / mouse `Zfp821`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `C16orf87` / mouse `4921524J17Rik`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `ZNF281` / mouse `Zfp281`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `C14orf119` / mouse `1700123O20Rik`: formal_nominal_shared_full_matrix; 2/2 formal branches.
- `ZNF410` / mouse `Zfp410`: formal_nominal_shared_full_matrix; 2/2 formal branches.

## Recommended Manuscript Use

- Build the main model around Tier 1: shared downstream morphology/excitability/regulatory programs, not shared regional identity.
- Use Tier 2 as pathway support for synaptic wiring, calcium signaling, potassium/GABA/glutamate receptor maturation, and axon guidance.
- Use Tier 3 and Tier 4 as supportive context after checking broad neuronal or housekeeping interpretations.
- Keep Tier 5 outside central claims until raw-count/object-level or external validation supports those non-identical ortholog mappings.

## Outputs

- Candidate tier table: `Project/results/primary_core_manuscript_candidate_tiers.tsv`
- Source formal rank model: `Project/results/primary_core_mgi_ortholog_formal_rank_model.md`
