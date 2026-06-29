# Distinct Dentate and Cerebellar Granule-Cell Lineages Converge through Niche and Circuit Constraints

## Author information

Jie Lu

Affiliation: [to be added]

Correspondence: [email to be added]

ORCID: https://orcid.org/0000-0001-6843-9720

## Short Title

Integrated granule-cell convergence

## Summary Statement

Distinct dentate and cerebellar granule-cell lineages preserve regional identity while converging on a granule-like assembly state shaped by niche maturation signals and sparse-expansion circuit constraints.

## Abstract

Granule cells in the cerebellum and dentate gyrus share compact excitatory input-expansion morphologies despite arising from different developmental territories and maturing in different circuits. Similar morphology could reflect hidden shared fate, convergent postmitotic construction, niche signaling, circuit constraint, or their integration. We tested these alternatives by combining a strict 10-dataset primary core spanning mouse dentate, mouse/human cerebellar, and human dentate/hippocampal resources with ortholog-aware rank-meta modeling, named-comparator tests, pseudotime/stage analysis, ligand-receptor prediction, regulatory-compatibility scoring, morphology/activity validation, and sparse expansion-coding simulation. Dentate and cerebellar granule-cell candidates remained regionally identity-separated, arguing against one universal granule-cell fate program. Instead, downstream neurite/morphology and synaptic/excitability modules converged more strongly than upstream fate/niche modules, and an identity-coupled configuration score was broadly positive across primary-core contrasts. A conservative seed set comprising GPM6A, NFIB, NFIA, KCNK1, RFX3 and GABRA2 marked a reusable assembly and maturation toolkit rather than an exclusive granule-cell barcode. TGF-beta/BDNF/SMAD/MAPK programs behaved as stage-windowed maturation/readiness overlays, while sender-receiver analysis nominated branch-specific developmental niche cues including cerebellar Purkinje-to-granule SHH/IGF1 and dentate SGZ astrocyte/vascular/immune interactions. NeuroMorpho, DANDI and sparse expansion-coding analyses supported compact input sampling under activity and resource constraints. An evidence-weighted hypothesis matrix disfavored a hidden shared-fate explanation and supported an integrated convergence model combining identity-coupled assembly, niche maturation signals and circuit constraints. These analyses suggest that compact granule-cell morphology emerges when distinct regional lineages deploy related postmitotic assembly machinery within compatible developmental windows and sparse-expansion circuit architectures.

## Key words

granule cell; dentate gyrus; cerebellum; developmental convergence; postmitotic assembly; neuronal morphology; single-cell transcriptomics; TGF-beta; BDNF; sparse coding/pattern separation

## Introduction

Granule cells are among the most abundant neurons in the mammalian brain, but "granule cell" is not a single lineage identity. Cerebellar granule cells arise from hindbrain/rhombic-lip progenitors, expand in the external granule layer, and mature under cerebellar niche signals including SHH. Dentate gyrus granule cells arise in a telencephalic hippocampal neurogenic context and continue to pass through immature-to-mature states during postnatal and, in rodents, adult neurogenesis. Despite these distinct developmental histories, both populations adopt compact excitatory-neuron designs with small somata, dense packing, restricted input-sampling structures, and roles in transforming input spaces into sparse downstream representations (Leutgeb et al., 2007; Bird et al., 2024).

This resemblance poses a developmental mechanism question. One possibility is that cerebellar and dentate granule cells share a hidden molecular fate. A second is that distinct lineage programs independently recruit overlapping downstream machinery for neurite outgrowth, synaptogenesis, excitability tuning, adhesion/guidance and maturation. A third is that local niche signals and circuit-level input-expansion constraints favor similar compact architectures even when lineage identity remains different. These possibilities are not mutually exclusive, but they make different predictions: a shared-fate model should reveal broad granule-cell specificity, whereas a developmental-convergence model should preserve regional identity while revealing later reuse of assembly modules.

This question is sharpened by earlier work showing that cerebellar conditioned medium suppresses hippocampal granule-cell proliferation and promotes differentiation, with TGF-beta2, BDNF, SMAD and MAPK signaling implicated as candidate mediators (Lu et al., 2005). That study suggested that extrinsic cerebellar factors can push hippocampal granule-lineage cells toward differentiation. Here, we revisit the same biological question at genome-scale resolution, asking whether distinct granule-cell lineages enter related developmental assembly states.

We therefore integrated a strict 10-dataset primary transcriptomic core with ortholog-aware rank-meta modeling, pseudotime and stage analysis, named non-granule comparators, niche/pathway scoring, focused sender-receiver ligand-receptor prediction, provisional regulatory-compatibility scoring, external morphology and activity validation, hierarchical evidence synthesis, and sparse expansion-coding simulation. Cerebellar and dentate granule cells remain regionally identity-separated, but they partially converge on downstream postmitotic construction modules. Similar granule-cell morphology is therefore modeled as the output of identity-coupled assembly configuration under developmental timing, niche and circuit constraints, rather than as a single shared fate identity or a simple morphology gene list.

## Results

### A strict 10-dataset primary core frames a developmental convergence test

The developmental-convergence model predicts that cerebellar and dentate granule cells should preserve branch-specific origin while sharing later assembly logic. These cells present a useful test case because they are not the same cell type, but they share a compact excitatory input-expansion design. Cerebellar granule cells arise from a rhombic-lip lineage and mature under cerebellar niche and SHH-associated developmental constraints (Wechsler-Reya and Scott, 1999; Lewis et al., 2004; Wizeman et al., 2019; Peng et al., 2019). Dentate granule cells arise in a hippocampal neurogenic context, retain prolonged postnatal and adult neurogenic programs in rodents, and show age-structured or controversial neurogenic signatures in humans (Kuhn et al., 1996; Artegiani et al., 2017; Eriksson et al., 1998; Spalding et al., 2013; Sorrells et al., 2018; Boldrini et al., 2018; Ramnauth et al., 2025). The project therefore tests convergence without assuming identity equivalence (Fig. 1a,b).

To make the comparison balanced, we assembled a strict 10-dataset primary core. The core contains four mouse dentate datasets (GSE104323, GSE95752, GSE292261, GSE214309), three cerebellar datasets (GSE122357, GSE165657, GSE312658), and three human dentate/hippocampal datasets (GSE186538, GSE325391, GSE268609). The mouse dentate arm used GSE104323 and GSE95752 from the postnatal dentate single-cell atlas (Hochgerner et al., 2018), GSE292261 from a postnatal hippocampal granule-cell synaptic-development resource (Lorente-Echeverria et al., 2025), and GSE214309 from an immature-versus-mature activity-state dentate granule-cell resource (Parylak et al., 2023). The cerebellar arm used GSE122357 (Peng et al., 2019), GSE165657 (Zhong et al., 2023), and GSE312658 (Chen et al., 2026). GSE312658 was retained as a cerebellar validation resource after within-dataset rank normalization; perturbation biology from that resource was not used to define baseline granule-cell identity. The human dentate/hippocampal arm used GSE186538 (Franjic et al., 2022), GSE325391, and GSE268609 (Disouky et al., 2026). The cerebellar datasets were treated as one regional/cell-lineage branch after ortholog-aware within-dataset ranking because only three cerebellar resources met primary-core criteria, including one human cerebellar validation dataset; therefore, cerebellar effects are interpreted as branch-level support rather than as species-resolved human-versus-mouse conclusions. Human scaffolding datasets (GSE185277 and GSE185553) were retained for label construction and marker tuning but were not counted in the strict primary core (Zhou et al., 2022; Fig. 1c). This separation allowed us to define granule-cell candidates, test cross-branch molecular convergence, and preserve branch-specific regional identity as distinct analytical tasks (Tables S1 and S2).

### Ortholog-aware rank-meta modeling identifies a conservative shared assembly-candidate set

Having defined a primary-core comparison, we next asked whether the data contained a reproducible shared molecular layer after species, platform and dataset effects were reduced, and whether that layer was more consistent with fate identity or postmitotic assembly. The formal MGI ortholog rank-meta model tested 17,611 target genes, generated 116,013 dataset-level deltas and 36,303 branch tests, and identified 1,370 formally shared hits. Because this broad set could reflect general neuronal maturation as well as granule-cell biology, we used it as a discovery layer rather than as a direct morphology program.

We therefore converted the rank-meta output into a tiered candidate framework (Fig. 2a). Genes were ranked within local candidate-versus-background contrasts, harmonized through one-to-one MGI orthology, and evaluated across selected-feature and full-MGI evidence layers. Candidates gained priority when support recurred across dentate and cerebellar branches and across more than one evidence route. Tier 1 represents the most conservative seed set, Tier 2 captures high-confidence wiring and synaptic/excitability support, Tiers 3 and 4 retain broader or screen-specific mechanism genes, and Tier 5 is exploratory because it depends on less direct non-identical ortholog evidence.

This filtering produced a compact developmental-mechanism candidate structure rather than a simple long gene list: 6 Tier 1 genes, 9 Tier 2 genes, 5 Tier 3 genes, 16 Tier 4 genes, and 30 Tier 5 exploratory genes. The Tier 1 seed set consists of GPM6A, NFIB, NFIA, KCNK1, RFX3, and GABRA2. These genes were prioritized for reproducible cross-branch support and mechanism relevance, not for exclusive granule-cell specificity. They span neurite or membrane outgrowth, developmental transcriptional regulation, ciliogenesis/transcriptional regulation, synaptic maturation, and excitability tuning. Tier 2 added stronger pathway context around this seed set, including ROBO2, PPP3CA, CACNA2D1, KCNJ6, GABRB3, GRIN2B, KCNJ3, KCND2, and STXBP5L.

Tier 1 and Tier 2 candidates showed positive candidate-background rank deltas across selected-feature and full-MGI branch layers, indicating reproducible relative enrichment rather than dependence on one dataset or one expression-processing route (Fig. 2b).

The prioritized candidates organized into four biologically interpretable axes: developmental regulatory control, neurite/cytoskeleton morphogenesis, synaptic/excitability maturation, and axon guidance/adhesion (Fig. 2c; Tables S3-S8). This axis structure reframes the shared gene signal as an implementation toolkit for neuronal assembly and wiring rather than as a single granule-cell identity barcode.

### Named comparators separate shared assembly modules from granule-cell identity

We first asked whether the four candidate-gene mechanism axes from Fig. 2c behaved as granule-cell-specific signatures. These axes classify prioritized genes by molecular function: developmental regulatory control, neurite/cytoskeleton morphogenesis, synaptic/excitability maturation, and axon guidance/adhesion. Using explicit source cell-type labels from GSE104323 and GSE122357, dentate granule-lineage groups were compared with the GSE104323 pyramidal-comparator labels CA3-Pyr and Immature-Pyr, while cerebellar granule-lineage groups were compared with Purkinje cells. None of the four axes passed strict dentate-granule > pyramidal and cerebellar-granule > Purkinje specificity (Fig. 3a). This negative result argues against a simple uniquely granule-specific pathway signature.

A complementary gene-level comparator screen gave the same caution at higher resolution. A small set of interpretable shared candidates, including NFIA, RBFOX3, and HMGN2, showed granule-enriched behavior relative to named pyramidal and Purkinje comparators, but the pattern did not collapse into a broad exclusive granule-cell marker program (Fig. S2).

Second, we reorganized the biology into five broader layers to ask where convergence sits in the developmental hierarchy. Three layers represented upstream or contextual programs: cerebellar fate/rhombic-lip/SHH, dentate fate/WNT/PROX1, and shared neurogenic niche/progenitor state. Two layers represented downstream construction programs: neurite/morphology and synaptic/excitability. These five layers are broader than the four mechanism axes in Fig. 3a; for example, axon guidance/adhesion is folded into the downstream neurite/morphology layer, and the synaptic/excitability layer includes additional maturation and activity genes beyond the candidate-axis set.

The formal module result strongly supported the layered model. Downstream modules had median convergence delta 0.500, whereas upstream/niche modules had median convergence delta -0.500, with Mann-Whitney p=1.56e-09. The strongest shared signal therefore occurs at the implementation layer, while lineage and niche programs preserve regional separation (Fig. 3b).

Ordered-state timing separated early branch identity from later construction reuse (Fig. S3). A developmental-origin control analysis further distinguished deep neural progenitor competence, anterior/telencephalic and medial pallium/dentate lineage patterning, hindbrain/rhombic-lip patterning, and later postmitotic/construction programs. Dentate states moved from high neural progenitor rank toward higher dentate-origin polarity and construction rank, whereas cerebellar granule precursor states retained high rhombic-lip/hindbrain rank and high shared maturation/construction ranks. The pattern supports distinct regional routes with later toolkit reuse rather than migration of one recent common granule-cell progenitor to both sites (Fig. S4; Tables S9-S13 and S41-S53).

### An identity-coupled transcriptomic configuration captures postmitotic assembly

To formalize convergence as a developmental state rather than a marker list, we defined a transcriptomic configuration score combining two components: construction-over-niche balance and regional fate polarity (Fig. 4a). The local named-comparator test showed positive combined configuration scores in 4/4 named granule-versus-comparator contrasts, but the small-n Wilcoxon p was 0.0625 (Fig. 4b). We therefore treated this named-comparator result as a directional specificity check rather than as the main statistical anchor.

We then broadened this to primary-core candidate-versus-background contrasts (Fig. 4c). Across 63 contrasts from seven datasets and two expression layers, 52/63 contrasts were positive. The median candidate-background combined configuration delta was 0.417, with sign-test p=8.37e-08 and Wilcoxon p=4.89e-08. The signal was positive in 8/8 full MGI ortholog matrix contrasts and 44/55 selected-feature matrix contrasts. Branch-level results were also positive: cerebellar 12/12 and dentate 40/51. The combined score reflects both construction-over-niche balance and branch-matched fate polarity; construction-over-niche balance alone was weaker.

Driver decomposition refined the claim (Fig. 4d). Regional fate polarity and downstream construction balance both contributed, with many positive contrasts driven by fate plus construction rather than construction alone. Thus, Fig. 4 supports a coupled assembly-configuration model, while Fig. 3 shows that downstream construction modules converge more strongly than upstream fate/niche modules.

These results support an identity-coupled assembly-configuration model rather than a single morphology-only gene signature (Fig. 4; Tables S14-S18).

### Stage-windowed niche signaling and sparse-coding constraints refine the 2005 model

The 2005 conditioned-medium result motivated a developmental-window analysis of TGF-beta/BDNF/SMAD/MAPK and related niche pathways (Lu et al., 2005). Across the primary core, we scored TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP, Reelin, Semaphorin, SHH, WNT, FGF, Notch, and ligand-receptor readiness modules (Fig. 5a). Pathway modules were candidate-enriched in 336/566 contrasts, composite signatures in 172/252 contrasts, and ligand-receptor readiness pairs in 139/285 contrasts. The historical TGF-beta/BDNF index was positive in 48/63 contrasts, but the branch pattern did not fit a simple cerebellar-biased stop-signal model: dentate candidates were positive in 47/51 contrasts, whereas cerebellar candidates were positive in only 1/12 contrasts. The clearest cerebellar pathway signal was SHH/PTCH/GLI, positive in 10/12 contrasts, consistent with established SHH control of cerebellar granule precursor expansion (Wechsler-Reya and Scott, 1999; Lewis et al., 2004). Thus, TGF-beta/BDNF/SMAD/MAPK remains an important maturation/readiness axis, but the sequencing data reframe it as stage- and context-dependent rather than as a universal cerebellar stop signal.

Stage-resolved modeling supported this interpretation (Fig. 5b). We fit quadratic branch curves to ask whether dentate and cerebellar observations differed in starting level, direction of change across developmental stage, or evidence for a peak-like window. For the 2005 TGF-beta/BDNF index, dentate scores were highest in an early/intermediate window (estimated peak stage 0.21) and declined toward later stages (endpoint delta -0.416). Cerebellar scores instead increased toward the latest sampled postnatal states (endpoint delta 0.333; estimated peak stage 0.99). With uneven stage support and cerebellar timing based mainly on GSE122357 P0/P8a/P8b, the branch differences support a stage-windowed model rather than one monotonic stop signal.

A complementary immature-versus-mature split supported the same conclusion: module similarity changed between immature and mature/maturing granule-cell bins rather than remaining fixed across the lineage (Fig. S1).

External morphology and activity data then constrained how this signaling model could map onto granule-cell design (Fig. 5c). NeuroMorpho analysis included 558 dentate and 62 cerebellar granule-cell reconstructions (Halavi et al., 2012; Polavaram et al., 2014). Dentate and cerebellar granule cells were not geometrically identical: dentate cells had fewer primary stems but larger dendritic fields, whereas cerebellar cells had more short stems or claw-like input structures. Median branch count was nevertheless similar, with 21 branches in dentate and 20 in cerebellar reconstructions. DANDI 000003 added a dentate activity constraint (Subash et al., 2023): across six NWB sessions, 124 units were analyzed, including 26 source-labeled granule units, with median spatial information 0.7800 bits/spike, spatial sparsity 0.4537, active spatial-bin fraction 0.5489, and awake-moving firing rate 0.0627 Hz. These data provide calibration for spatially structured dentate granule activity and compact input-sampling architecture.

The sparse expansion-coding model tested whether such compact input-expansion designs are computationally plausible (Fig. 5c). Dense high-activity expansion achieved strong raw useful scores, but it was expensive in wiring and activity load. Under resource-adjusted scoring, sparse granule-like designs were favored, and excessive sparsity performed poorly because it lost useful information. Calibration against NeuroMorpho and DANDI summaries sharpened the same point: the highest raw empirical grid point was an intermediate design with expansion ratio 0.5, input degree 2, observed output active fraction 0.0625, useful score 0.094, resource-adjusted useful score 1.507, and empirical calibration score 4.103. Across architecture families, raw empirical calibration favored dense_expansion, but emphasizing resource and morphology constraints shifted the leading family to intermediate, followed by granule_like_sparse_expansion and then dense_expansion. Thus, dense expansion can win under raw scoring, but biologically constrained nontrivial expansion favors intermediate or sparse granule-like designs, while excessive sparsity loses useful information.

Together, the stage-window, morphology, activity, and sparse-coding analyses refine the convergence model: distinct dentate and cerebellar lineages can pass through different signaling windows while converging on compact input-expansion designs constrained by morphology, excitability, and circuit resources (Fig. 5d; Tables S19-S37).

### Sender-receiver ligand-receptor prediction nominates testable developmental niche cues

To move from pathway readiness toward experimentally testable developmental niche hypotheses, we built a focused sender-receiver ligand-receptor analysis for the two primary datasets with suitable niche and receiver annotations. The cerebellar analysis used GSE122357 P0, P8a, and P8b mouse cerebellum; sender classes included Purkinje cells, astrocytes as a Bergmann/astroglial proxy, microglia, and endothelial cells, while receiver classes included granule precursors and granule cells, consistent with known cerebellar granule precursor dependence on niche SHH signaling (Wechsler-Reya and Scott, 1999; Lewis et al., 2004). The dentate SGZ analysis used GSE104323; sender classes included astrocyte states, endothelial cells, PVM/macrophage as a microglia-like proxy, and vascular/support classes, while receivers included RGL, nIPC, neuroblast, immature GC, juvenile GC, and adult GC states (Fig. 6a).

Across 47 curated ligand-receptor pairs, the screen generated 3,008 core-focus sender-receiver predictions, of which 208 were moderate/high expression-supported. Cerebellar samples contributed 97 supported predictions among 1,128 tested rows, with median supported score 0.002 and maximum score 0.502. Dentate SGZ samples contributed 111 supported predictions among 1,880 tested rows, with median supported score 0.001 and maximum score 0.557. High-scoring cerebellar examples included Purkinje-to-granule precursor IGF1-IGF1R and SHH-PTCH1/SMO interactions. High-scoring dentate examples included immature astrocyte-to-granule SEMA6A-PLXNA2/4, together with APOE/LRP, APOE/LDLR, C1QA/LRP1, and JAG1/NOTCH1 pathway support (Fig. 6b,c).

This sender-receiver layer sharpens the pathway-readiness model by requiring ligand expression in the sender and receptor expression in the receiver. The resulting predictions define a focused experimental menu for spatial, protein-level, and functional testing.

### Evidence-weighted hypothesis comparison favors an integrated convergence model

We finally compared the three conceptual hypotheses against the integrated evidence. H1, the hidden shared-fate hypothesis, predicts broad shared fate identity with weak regional separation. H2 predicts preserved regional identity plus shared downstream construction. H3 asks why a compact granule-like configuration would be favored through stage, niche, and circuit constraints. The evidence-index score was below neutral for H1 (39.3), but above neutral for H2 (69.1), H3 (76.9), and the integrated H2+H3 synthesis (71.8). Thus, the data do not support H1 alone; they support transcriptomic assembly convergence and niche/circuit constraints as complementary explanations (Fig. 7; Fig. S5; Tables S89-S93 and S102-S103).

A public perturbation module-shift audit nominated SHH/PTCH/Norrin, RBFOX-family, and NTRK2/BDNF-related contrasts as tractable perturbation contexts for downstream construction or maturation modules (Fig. S6; Tables S94-S100).

## Discussion

### Developmental convergence rather than shared fate identity

The central conclusion is developmental rather than taxonomic: dentate and cerebellar granule cells are not the same cell type, but distinct regional lineages can converge on related downstream postmitotic assembly configurations. Regional identity remains strongly separated in the primary core, whereas neurite/morphology, synaptic/excitability, and maturation programs provide the strongest shared signal. The shared morphology therefore appears less like a single granule-cell fate identity and more like repeated use of overlapping construction machinery in two different developmental and circuit settings (Fig. 3b; Fig. 4; Fig. 5d; Fig. 7).

The candidate-gene results support this developmental interpretation. The conservative Tier 1 seed set (GPM6A, NFIB, NFIA, KCNK1, RFX3, and GABRA2) did not behave as a granule-cell-exclusive barcode when tested against named pyramidal and Purkinje comparators. NFIA was the clearest shared granule-enriched seed, whereas the other seeds showed branch-specific, tied, or comparator-enriched behavior. Thus, the seed set is best interpreted as a reusable assembly and maturation toolkit whose deployment depends on lineage, stage, and circuit context (Fig. 2; Fig. 3a; Fig. S2).

The developmental analyses argue against a recent common dentate/cerebellar granule progenitor. The data are most consistent with a deep shared neural origin before anterior-posterior regional patterning. After that point, hippocampal granule cells follow a telencephalic/medial-pallial/dentate route, whereas cerebellar granule cells follow a hindbrain/rhombic-lip/external-granule-layer route. The relevant similarity emerges later, when postmitotic maturation and neurite/synapse construction programs become active in each branch (Fig. S3; Fig. S4).

### The 2005 conditioned-medium result as a stage-window clue

The sequencing data support the broad relevance of TGF-beta, BDNF, SMAD, MAPK and related maturation pathways, but they convert the original conditioned-medium observation into a stage-window hypothesis. The 2005 study measured secreted protein-level bioactivity from mixed postnatal cultures after days in vitro, not acute granule-cell mRNA abundance. The present RNA data therefore should not be expected to reproduce the conditioned-medium effect as a simple cerebellum-greater-than-dentate transcriptomic contrast (Fig. 5a,b).

Instead, the data point to stage- and context-dependent readiness. TGF-beta/BDNF signatures are stronger in dentate candidate-background contrasts than in cerebellar candidate-background contrasts, while SHH/PTCH/GLI is the clearest cerebellar pathway signal. In the stage-resolved mouse cerebellar dataset, TGFB2 rises from P0 to P8 and the combined TGF-beta/BDNF 2005 index peaks late, consistent with a postnatal cerebellar readiness window. BDNF mRNA, however, is sparse in cerebellar granule candidates and more evident in dentate datasets. This suggests that the BDNF component of the 2005 mechanism is supported more directly by the original antibody-neutralization, recombinant-factor, and signaling experiments than by cerebellar granule-cell transcript abundance.

The sender-receiver analysis extends this revision by nominating branch-specific niche interactions. Purkinje-to-granule precursor SHH and IGF1 interactions match established cerebellar developmental logic (Wechsler-Reya and Scott, 1999; Lewis et al., 2004). In the dentate SGZ, astrocyte, vascular, and immune-support interactions nominate SEMA6A-PLXNA2/4, APOE/LRP, C1QA/LRP1, and JAG1/NOTCH1 as candidates for stage-dependent maturation or permissive niche control. These remain hypotheses until tested by spatial, protein, or functional assays (Fig. 6).

### A transcriptomic assembly state, not a morphology gene list

The transcriptomic configuration score shows that morphology-relevant assembly is visible in RNA data as a stage- and identity-coupled state that biases neurite, synaptic, and excitability programs. Final morphology then depends on local inputs, niche signals, activity, spatial constraints, and developmental timing (Fig. 4; Fig. 5c,d).

The provisional epigenomic extension fits this view as a hypothesis-generating layer. If dentate and cerebellar granule cells use overlapping construction genes in different regional contexts, then promoters, enhancers, methylation state, and transcription-factor motif activity near fate and construction genes should help explain how the same toolkit is deployed differently (Fig. S5).

The hierarchical evidence model summarizes this mixed support across transcriptomic configuration, stage/niche context, provisional regulatory compatibility, morphology/activity calibration, and resource-constraint agreement. Mapping these terms onto the three hypotheses argues against a single hidden shared fate and favors the combined H2+H3 explanation: distinct lineages reuse related assembly programs, while stage, niche and circuit constraints help explain why compact granule-like designs are favored (Fig. 7; Fig. S5).

### Circuit constraints explain why compact granule-like designs can emerge

The sparse-coding model provides a mathematical reason why similar developmental assembly states may repeatedly resolve into compact granule-like designs. Dense expansion can maximize raw separation in a simplified model, but it is costly in wiring and activity. When morphology and resource costs are included, intermediate or sparse granule-like architectures become more favorable, whereas excessive sparsity loses useful information. This offers a quantitative language for a biological tradeoff: granule-like morphology may be a constrained solution for input expansion and pattern separation rather than an accidental resemblance between two brain regions (Fig. 5c,d; Leutgeb et al., 2007; Bird et al., 2024).

Public perturbation resources point to tractable future tests. SHH/Ptch-axis perturbations are most directly relevant to cerebellar granule precursor fate and proliferation, RBFOX-family perturbations to postmitotic hippocampal/dentate maturation and plasticity, and NFIA to regulatory control of construction programs. BDNF/TrkB and TGF-beta/SMAD resources can test pathway responsiveness, while HMGN2 remains a chromatin-competence candidate pending direct neural perturbation data (Fig. S6).

### Limitations

1. The study is computational and integrative, and should be read as a developmental model with prioritized experimental tests. Rank-meta modeling, hierarchical scoring, and sparse-coding simulation prioritize mechanisms, but direct developmental perturbation remains the highest-priority validation step.
2. Cross-dataset harmonization reduces but does not remove species, platform, stage, and annotation differences. Named-comparator tests depend on datasets with explicit pyramidal and Purkinje labels, and GSE325391 and GSE268609 remain selected-feature bridge trajectories rather than strict full-transcriptome trajectories in this implementation.
3. The epigenomic analysis is provisional. It defines regulatory targets and cerebellar multiome feasibility scores, but matched dentate peak-count scoring, methylome fitting, and direct enhancer perturbation remain future needs.
4. Ligand-receptor predictions are expression-based. They do not establish spatial contact, protein secretion, receptor activation, or conditioned-medium bioactivity.
5. Morphology and physiology validation are incomplete. NeuroMorpho reconstructions are not matched to transcriptomes, and DANDI evidence supports dentate spatial coding but not direct cerebellar granule-cell physiology.
6. Public perturbation signatures are heterogeneous and preliminary, combining bulk, cultured-neuron, neural-progenitor, and array resources. They prioritize but do not replace matched granule-cell perturb-seq or cell-type-specific perturbation experiments.

### Experimental tests

1. A stage-matched perturbation experiment should test whether TGF-beta/BDNF/SMAD/MAPK manipulation shifts dentate and cerebellar granule-lineage proliferation, neurite maturation, and synaptic/excitability modules in different developmental windows.
2. Perturbing Tier 1/Tier 2 assembly candidates should alter granule-cell neurite, synaptic, or excitability maturation more than regional fate identity.
3. Secretome/proteomic assays of cerebellar conditioned medium should identify factors beyond TGF-beta2 and BDNF that cooperate with maturation-stage state.
4. Sender-receiver candidates such as cerebellar SHH/IGF1 and dentate SGZ SEMA6A, APOE/LRP, C1QA/LRP1, and JAG1/NOTCH1 interactions should show spatial and protein-level support in matched developmental windows.
5. Matched perturb-seq or cell-type-specific perturbation should sharpen the preliminary public-resource result: SHH/Ptch-axis perturbations should alter cerebellar GNP fate, maturation, and synaptic/excitability modules; RBFOX3 perturbation should preferentially alter postmitotic neurite, synaptic, and plasticity modules; NFIA perturbation should affect regulatory construction targets more strongly than broad regional identity; and HMGN2 should be tested directly before being treated as causal.
6. Morphology-linked datasets should separate primary input-sampling geometry from dendritic-field complexity.
7. Additional DANDI or equivalent physiology datasets should show that sparse spatial or task-selective firing is more informative than long-timescale silence alone.

## Materials and Methods

### Dataset discovery, curation and primary-core design

Public datasets were selected to test whether dentate and cerebellar granule-cell lineages preserve regional identity while converging on later neuronal assembly programs. Candidate resources were collected from GEO, Allen-related resources, literature-linked human dentate/hippocampal studies, cerebellar developmental studies, candidate scATAC/multiome or methylome resources, NeuroMorpho.Org and DANDI. The strict primary core contained ten datasets: four mouse dentate datasets, three cerebellar datasets and three human dentate/hippocampal datasets. Mouse dentate resources included GSE104323/GSE95752, GSE292261 and GSE214309. Cerebellar resources included GSE122357, GSE165657 and GSE312658. Human dentate/hippocampal resources included GSE186538, GSE325391 and GSE268609. Supporting datasets were retained for label construction, marker validation, selected-feature bridge analyses, regulatory-compatibility extension or perturbation interpretation, and were kept analytically distinct from the strict primary core. The primary-core frame and supporting dataset tiers are provided in Tables S1 and S2.

### Human bridge-object construction and label harmonization

Human dentate/hippocampal bridge datasets were processed as selected-gene or reduced sparse objects when full transcriptome matrices were not locally available or were not needed for a specific validation layer. Source annotations, marker programs and reduced human-core label projections were used to harmonize candidate granule-cell, immature-neuron, neurogenic, pyramidal/comparator and background labels. Human scaffold datasets GSE185277 and GSE185553 were used for marker tuning and label construction from the human hippocampal immature-neuron lifespan resource, but they were not counted in the strict ten-dataset primary core.

### Candidate-cell annotation and pseudobulk construction

Candidate granule-cell, precursor, comparator and local background labels were assigned using source annotations, marker support and branch-specific label rules. Candidate and background groups were converted to pseudobulk summaries to stabilize cross-dataset comparisons and reduce sensitivity to differences in cell yield. Candidate-gene, expanded-gene, genome-wide symbol and full-MGI ortholog pseudobulk layers were generated separately. Per-cell module scores were used for selected validation checks, but primary cross-dataset tests used within-dataset or within-sample rank summaries.

### Ortholog-aware rank-meta modeling and candidate tiering

The primary transcriptomic screen used an ortholog-aware rank-meta design. For each dataset or contrast, candidate and local-background pseudobulk expression was converted to within-layer relative ranks, and candidate-background rank deltas were computed within the same dataset layer. Human and mouse genes were harmonized through one-to-one MGI orthologs where cross-species comparison was required. This design emphasizes reproducible within-dataset ordering rather than absolute expression magnitude, which is sensitive to species, chemistry, processing pipeline and sequencing depth.

Two evidence routes were tracked. The selected-feature route used curated candidate-gene, expanded-gene, module-gene and human bridge-feature matrices; this route maximized coverage for predefined granule-cell and mechanism features but was not treated as genome-wide discovery. The full-MGI/full-matrix route used broader pseudobulk matrices remapped to one-to-one MGI orthologs where full matrices were available. For each gene, screen and branch, nominal branch support required at least two datasets, median dataset rank delta \(\widetilde{\Delta}_{\mathrm{rank}} > 0\), positive dataset fraction \(f_{\Delta > 0} \ge 0.75\) and best one-sided dataset-level \(p_{\min} \le 0.25\) from the t, sign or Wilcoxon tests. Shared support was assigned without pooling raw expression across datasets: selected-feature shared support required nominal dentate and cerebellar support in the selected-feature route, full-matrix shared support required nominal dentate and cerebellar support in the full-MGI/full-matrix route, and both-screen shared support required both route-level calls.

Table S3 reports the formal-rank fields used for candidate ordering. Four possible screen-branch tests were tracked for each gene: selected-feature dentate, selected-feature cerebellar, full-matrix dentate and full-matrix cerebellar. The field formal_n_available_branches is the number of these screen-branches with sufficient data, and formal_n_nominal_branches is the number satisfying the nominal support rule above. FDR10 support required replication support plus the best within-screen/branch BH-adjusted one-sided q value \(\le 0.10\). The field formal_fisher_q_bh_all_available_branches is the BH-adjusted Fisher combined p value computed from the available one-sided branch-level \(p_{\min}\) values; it was used as a summary support statistic, not as a sole inclusion rule. Median-delta fields in Table S3 are median within-dataset candidate-minus-background rank deltas for the indicated screen and branch, so positive values indicate higher within-dataset rank in candidate granule populations than in local background populations. The formal-rank priority score was a deterministic ranking index:

$$
P_{\mathrm{formal}}
= 3N_{\mathrm{FDR10}}
+ 2N_{\mathrm{nominal}}
+ N_{\mathrm{replication}}
+ \sum_{b}\max(0,\widetilde{\Delta}_{\mathrm{rank},b})
+ \min_{b}(D_b).
$$

Here, \(b\) indexes available screen-branch tests, \(N_{\mathrm{FDR10}}\), \(N_{\mathrm{nominal}}\) and \(N_{\mathrm{replication}}\) are the numbers of branches meeting FDR10, nominal and replication support, respectively, \(\widetilde{\Delta}_{\mathrm{rank},b}\) is the median rank delta for branch \(b\), and \(D_b\) is median candidate-cell detection for branch \(b\). For mechanism-prioritized genes, a second unit-level check was fit within each gene, screen and branch:

$$
\Delta_{\mathrm{rank},u}
= \beta_0 + b_{d[u]} + \epsilon_u,\qquad
b_d \sim \mathcal{N}(0,\sigma_d^2),\qquad
\epsilon_u \sim \mathcal{N}(0,\sigma^2).
$$

Here, \(\Delta_{\mathrm{rank},u}\) is the candidate-minus-background rank delta for dataset/sample unit \(u\), \(d[u]\) identifies its dataset, \(b_{d[u]}\) is the dataset-specific random intercept, \(\epsilon_u\) is residual error, and \(\beta_0\) estimates the average directional rank shift for that gene in that screen-branch. The random intercept allows all units from the same dataset to share a baseline shift. The one-sided test was \(H_0:\beta_0 \le 0\) versus \(H_1:\beta_0 > 0\). Models were fit only when at least three units from at least two datasets were available. A linear mixed model with a dataset random intercept was used when possible; if the mixed model failed, an intercept-only ordinary least-squares model with dataset-cluster-robust standard errors was used as a fallback. formal_model_n_nominal_branches counts screen-branches with branch replication support and \(p_{\mathrm{greater}} \le 0.25\) in this unit-level model; the corresponding FDR10 count used branch replication support and Benjamini-Hochberg-adjusted \(q \le 0.10\).

Mechanism prioritization was used as an evidence-organization layer. Selected-feature and full-matrix mechanism triage classified genes into developmental regulatory control, neurite/cytoskeleton morphogenesis, axon guidance/adhesion, synaptic wiring, calcium signaling, potassium-channel/excitability and glutamate/GABA receptor or vesicle-release categories. Within each route, mechanism candidates required positive dentate and cerebellar candidate-background rank deltas, membership in a core mechanism class and BH-adjusted branch-level rank support at \(q < 0.10\) in both branches. Genes supported in both routes yielded 24 cross-screen consensus candidates. A dataset-robust candidate required at least two dataset/sample units in each available screen-branch, \(f_{\Delta > 0} \ge 0.75\) and median rank delta \(\widetilde{\Delta}_{\mathrm{rank}} > 0\). Six genes passed all four screen-branch tests: GABRA2, GPM6A, KCNK1, NFIA, NFIB and RFX3.

Final candidate tiers used five ordered evidence features: formal shared support in at least one route; same-symbol one-to-one MGI orthology within the mechanism-prioritized packet; both-screen formal support; central mechanism-priority support; and dataset-robust consensus support. Tier 1 comprised genes satisfying all five evidence features, including the four-of-four dataset-robust consensus criterion. Tier 2 comprised high-confidence wiring, synaptic, excitability, calcium, glutamate/GABA and guidance support genes with both-screen formal support and central mechanism-priority support but without the dataset-robust four-of-four designation. Tier 3 retained same-symbol mechanism genes with both-screen formal support but weaker central-priority status. Tier 4 retained same-symbol mechanism genes supported in only one route. Tier 5 retained exploratory non-identical-symbol one-to-one MGI ortholog hits from the broader formal shared-hit table. Candidate-gene outputs are provided in Tables S3-S8.

### Named-comparator specificity and granule-enrichment screens

Named-comparator analyses tested whether shared candidates were broad neuronal maturation genes rather than granule-lineage-enriched programs. Dentate granule-lineage groups were compared with GSE104323 pyramidal-comparator labels CA3-Pyr and Immature-Pyr. Cerebellar granule precursor/granule labels were compared with Purkinje cells in GSE122357 across P0, P8a and P8b. Mechanism-axis modules and the six Tier 1 seed genes were evaluated using within-sample ranks of group-level median log1p expression. For each dataset, sample, and module or gene, source groups were ordered by median log1p score and converted to percentile ranks using average ranks for ties; higher ranks therefore indicate higher relative expression within that local sample. Comparator and granule-enrichment outputs are provided in Tables S9, S10, S45 and S46.

### Niche, circuit and transcriptomic-configuration analyses

Curated gene modules were defined for regional fate, neurogenic niche/progenitor state, downstream neurite/morphology, downstream synaptic/excitability and targeted circuit or pathway axes. Candidate-versus-background support was evaluated as within-dataset rank deltas, and named-comparator support was evaluated separately to distinguish shared assembly machinery from broad neuronal identity.

The transcriptomic-configuration model tested whether granule-cell convergence is better described as an identity-coupled configuration than as a single shared fate program. Construction balance was defined as downstream neurite/synaptic module strength minus niche/progenitor module strength:

$$
B_C = \bar{R}_{NS} - R_{NP}.
$$

Regional fate polarity was defined as branch-matched fate strength minus opposed branch fate strength:

$$
P_F = R_{BM} - R_{OP}.
$$

The combined configuration score was the sum of construction balance and regional fate polarity:

$$
S_Q = B_C + P_F.
$$

Here, \(R_{NS}\) is the neurite/synaptic rank, \(R_{NP}\) is the niche/progenitor rank, \(R_{BM}\) is branch-matched fate rank, \(R_{OP}\) is opposed-branch fate rank, and \(S_Q\) is the configuration score.

Configuration contrasts were validated across primary-core full-MGI and selected-feature layers and then decomposed into construction and fate-polarity components. Associated outputs are provided in Tables S14-S18.

### Regulatory-compatibility extension and integrative model specification

The regulatory-compatibility extension used epigenomic resources as a supporting chromatin layer. Candidate resources were triaged for matched or comparable scATAC, multiome or methylome support, including human hippocampus/dentate multiome data from GSE268609, adult primate cerebellar cortex multiome data from GSE322785, mouse ventral hippocampus methylome/RNA/ATAC data from GSE245367 and murine hippocampal sci-ATAC data from GSE118987. GSE268609 was represented in the strict primary core through RNA/selected-gene projection. For GSE322785, regulatory target sets were built from candidate tiers, niche/circuit modules and developmental-origin modules, then mapped to gene-expression features and ATAC peaks overlapping the gene body or lying within 2 kb, 10 kb or 100 kb of the target gene in three human cerebellar H5 feature-barcode matrices. Provisional marker groups were assigned from selected marker panels and scored for each selected feature as:

$$
X_{g,f}=\log(1+\bar{c}_{g,f})+d_{g,f}.
$$

Here, \(\bar{c}_{g,f}\) is the mean count for feature \(f\) in marker group \(g\), and \(d_{g,f}\) is the fraction of barcodes in that group with nonzero counts for the feature. Target-set scores were summarized by averaging feature scores across features assigned to the same model term, target set and peak category. Granule-candidate versus comparator contrasts were then computed as:

$$
\Delta_{k}=\bar{X}_{\mathrm{granule},k}-\bar{X}_{\mathrm{comparator},k},
$$

where \(k\) denotes a feature type, model term, target set and peak-category combination. A stricter sensitivity layer retained marker calls only when they were supported by donor-specific selected-gene clusters, requiring at least 20 marker-call barcodes in the cluster, at least 2% cluster representation and at least twofold enrichment over the donor-level marker-call frequency; ambiguous and low-information calls were excluded. Robust-positive contrasts required both broad provisional and cluster-supported pooled granule-minus-comparator mean deltas to be at least 0.10, and strong robust-positive contrasts required both deltas to be at least 0.25. Regulatory target sets, marker-group scores, cluster-supported sensitivity scores and robust contrast summaries are provided in Tables S54-S88.

To state the causal mixed-model that could be fitted in a fully matched future dataset, we specified a matched-data equation in which granule design depends on branch-matched fate polarity, construction balance, their coupling term, stage/pseudotime, niche signal, regulatory compatibility, morphology sparse sampling, activity sparsity and circuit resource constraint, with dataset, species and assay random effects:

$$
Y_i = \beta_0 + \beta_F F_i + \beta_C C_i + \beta_I F_i C_i + \eta_i .
$$

$$
\eta_i =
\beta_T T_i + \beta_N N_i + \beta_E E_i + \beta_M M_i
+ \beta_A A_i + \beta_R R_i + \xi_i .
$$

$$
\xi_i = u_{d[i]} + u_{s[i]} + u_{a[i]} + \varepsilon_i .
$$

Here, \(Y_i\) is the granule-design score, \(F_i\) is branch-matched fate polarity, \(C_i\) is construction balance, \(F_iC_i\) is the fate-construction coupling term, \(T_i\) is stage or pseudotime, \(N_i\) is niche signal, \(E_i\) is regulatory compatibility, \(M_i\) is morphology sparse sampling, \(A_i\) is activity sparsity, \(R_i\) is circuit-resource constraint, \(\eta_i\) collects additional observed covariates and random effects, and \(u_{d[i]}\), \(u_{s[i]}\) and \(u_{a[i]}\) are dataset, species and assay random effects. Because the present datasets are not matched across all modalities, the implemented model was a weighted hierarchical evidence synthesis. An evidence unit was one scored observation from a sample-level contrast, branch summary, epigenomic sensitivity contrast, morphology or activity calibration, or simulation-calibration result; unit_id labels these scored observations, not genes. Each evidence unit was normalized to a bounded score in \([-1, 1]\), assigned a quality weight according to measurement level, and summarized first by hierarchy level and term:

$$
S_{\ell,t} =
\frac{\sum_{j \in \mathcal{U}_{\ell,t}} w_j x_j}
{\sum_{j \in \mathcal{U}_{\ell,t}} w_j},
$$

where \(\mathcal{U}_{\ell,t}\) is the set of evidence units assigned to hierarchy level \(\ell\) and term \(t\), \(x_j\) is the normalized evidence score for unit \(j\), and \(w_j\) is its quality weight. Component scores were then computed as term-balanced weighted averages so that large tables did not automatically dominate smaller calibration layers. Term definitions, evidence units, layer summaries, branch summaries and component scores are provided in Tables S89-S93.

For the hypothesis-comparison matrix, the same term-level evidence scores were mapped onto three prespecified conceptual hypotheses. This support-index step used a study-defined evidence-alignment summary based on multi-criteria decision analysis and weighted additive scoring (Thokala et al., 2016; Marsh et al., 2016). Each hypothesis was represented as a vector of expected directions across bounded evidence terms, and the observed evidence as the vector \(S_t\). Evidence contributes positively when its sign and magnitude agree with the hypothesis prediction, negatively when it points in the opposite direction, and not at all when that term was not used to discriminate that hypothesis. Because the evidence came from heterogeneous public datasets that could not be pooled into one matched meta-analysis, the score was used as a structured synthesis framework (Campbell et al., 2020). H1 represented a hidden shared granule-cell fate identity, H2 represented identity-coupled transcriptomic assembly convergence, and H3 represented stage/niche/circuit-constraint convergence. Each hypothesis \(h\) was assigned a prediction coefficient \(a_{h,t}\) for each evidence term \(t\), where positive values indicate expected support, negative values indicate mismatch, and zero indicates no discriminating role. The signed alignment score was:

$$
A_h =
\frac{\sum_t w_t a_{h,t} S_t}
{\sum_t w_t |a_{h,t}|},
$$

where \(S_t\) is the observed term-level evidence score and \(w_t\) is the term weight. The denominator is the maximum weighted relevance available to hypothesis \(h\), which normalizes \(A_h\) to the interval \([-1, 1]\) when \(S_t\) is bounded by \([-1, 1]\). For visualization, \(A_h\) was converted to a 0-100 support index:

$$
\mathrm{SupportIndex}_h = 50(1 + A_h).
$$

Thus, 50 is neutral, values above 50 indicate evidence alignment and values below 50 indicate opposition or mismatch. The H2+H3 row was included as a synthesis of the two complementary convergence layers. The hypothesis coefficient matrix and final support scores are provided in Tables S102 and S103.

Public perturbation resources were triaged for NFIA, BDNF/TrkB, TGF-beta/SMAD, SHH, RBFOX3 and HMGN2. Candidate datasets were retained if they included a direct perturbation, pathway agonist/inhibitor response, knockout/knockdown or closely related developmental model in neural, hippocampal, dentate, cerebellar, granule progenitor or niche-relevant cells. Tractable resources were converted into curated module-shift signatures. For each contrast and module, the signed module-shift score was:

$$
M_{c,m} =
\mathrm{median}_{g \in G_m \cap A_c}(e_{c,g})
\sqrt{\frac{|G_m \cap A_c|}{|G_m|}},
$$

where \(e_{c,g}\) is the gene-level perturbation effect for gene \(g\) in contrast \(c\), \(G_m\) is the curated module gene set, and \(A_c\) is the set of module genes available in that public dataset. Perturbation-resource triage and module-shift outputs are provided in Tables S94-S100.

### Diffusion, pseudotime and developmental-stage analyses

For local full-transcriptome matrices, counts were normalized as log1p counts per 10,000, highly variable genes were selected, TruncatedSVD embeddings were computed for sparse expression matrices, and k-nearest-neighbor diffusion or pseudotime structure was inferred from curated or marker-inferred roots. GSE325391 and GSE268609 were analyzed as selected-feature bridge trajectories because the locally constructed objects were designed for label and marker-panel validation rather than uniform whole-transcriptome trajectory inference.

Stage-aware analyses separated immature and mature/maturing granule-cell states where source annotations supported the distinction. Regional-origin versus shared-toolkit timing was tested by ordering GSE104323 dentate states from radial glia-like/progenitor groups through neuroblast, immature granule, juvenile granule and adult granule states and comparing them with GSE122357 P0, P8a and P8b cerebellar granule precursor/granule labels. The developmental-origin divergence audit rescored GSE104323 and GSE122357 with curated marker modules for deep neural progenitor competence, anterior/telencephalic patterning, medial pallium/dentate lineage, hindbrain/rhombic-lip patterning, shared postmitotic granule maturation and downstream neurite/synapse construction. Module scores were computed as group-level median log1p expression and converted to within-sample module ranks. Branch-origin polarity and shared convergence rank were defined as:

$$
P_O = R_{BO} - R_{OO},
\qquad
C_S = \frac{R_{PM} + R_{NS}}{2}.
$$

Here, \(R_{BO}\) and \(R_{OO}\) are branch-matched and opposed-origin ranks, \(R_{PM}\) is postmitotic-maturation rank and \(R_{NS}\) is neurite/synapse rank. Supporting developmental-stage and origin-divergence outputs are provided in Tables S50-S53.

### Pathway-readiness, stage-window and conditioned-medium analyses

Pathway-readiness analyses used curated modules for TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP/SMAD, Reelin, Semaphorin, SHH/PTCH/GLI, WNT/beta-catenin, FGF/MAPK, Notch/HES, differentiation/stop, neurogenic/permissive state, the historical TGF-beta/BDNF 2005 index and stop-minus-permissive balance. The historical index was derived from the conditioned-medium mechanism reported by Lu et al. (2005). Candidate-versus-background support was evaluated as branch-level rank deltas.

To test whether pathway behavior was stage-windowed, pathway signatures were modeled across normalized stage or pseudotime order. Normalized stage was scaled from the earliest available state within a dataset/axis to the latest or most mature/activity-linked state. For each signature, the fitted weighted least-squares model was:

$$
Y_i =
\beta_0 + \beta_1 t_i + \beta_2 t_i^2 + \beta_3 z_i + \eta_i .
$$

$$
\eta_i = \beta_4 t_i z_i + \beta_5 t_i^2 z_i + \varepsilon_i .
$$

Here, \(Y_i\) is the signature score, \(t_i\) is normalized stage, \(z_i\) is an indicator variable equal to 1 for cerebellar observations and 0 for dentate observations, \(\eta_i\) contains the cerebellum-by-stage interaction terms and residual error, and dentate is the reference branch. Observation weights were proportional to:

$$
w_i \propto
\sqrt{n_{\mathrm{cell},i}}
\times
\max(n_{\mathrm{feat},i}, 1),
$$

where \(n_{\mathrm{feat},i}\) is the number of pathway genes or features present for observation \(i\). Weights were median-normalized during fitting. Coefficients were fit by weighted least squares, and uncertainty was summarized with HC3 robust standard errors; term-wise \(p\) values were Benjamini-Hochberg adjusted within each signature. The quadratic term was used to capture stage-windowed pathway readiness without fitting higher-order curves to sparse stage data. Fitted branch peaks were reported only for concave quadratic trends:

$$
\hat{s}_{D}
= -\frac{\beta_1}{2\beta_2},
\qquad
\hat{s}_{C}
= -\frac{\beta_1+\beta_4}{2(\beta_2+\beta_5)} .
$$

Here, \(D\) and \(C\) denote dentate and cerebellar branches, respectively; \(\beta_1\) and \(\beta_2\) are the linear and quadratic stage coefficients for the dentate reference curve, whereas \(\beta_4\) and \(\beta_5\) are the cerebellar interaction terms. A peak was interpreted only when the branch-specific quadratic coefficient was negative, producing a downward-opening curve. Peak positions were clipped to the observed \([0,1]\) normalized stage interval. Secreted or extracellular candidates that could complement TGF-beta2 and BDNF were prioritized for future proteomic or functional assays. Pathway-readiness and stage-window outputs are provided in Tables S19-S28.

### Focused sender-receiver ligand-receptor analysis

Focused ligand-receptor prediction distinguished general pathway readiness from directional niche hypotheses. Curated ligand-receptor pairs were scored in GSE122357 cerebellum and GSE104323 dentate SGZ using annotated sender and granule-lineage receiver classes. Cerebellar sender classes included Purkinje, astroglial/Bergmann-proxy, microglial and endothelial groups; dentate sender classes included astrocyte, PVM/microglia-macrophage proxy and endothelial groups. Pyramidal neurons were used as named neuronal comparators in the transcriptomic specificity analyses, but they were not included as dentate sender classes in this focused SGZ niche screen. The ligand-receptor expression score was:

$$
Q_{srp}
=
\sqrt{\max\!\left(L^{r}_{sp}R^{r}_{rp},0\right)}
\sqrt{\max\!\left(L^{d}_{sp}R^{d}_{rp},0\right)},
$$

where \(L^{r}_{sp}\) and \(L^{d}_{sp}\) are sender-class relative expression and detection fraction for the ligand in pair \(p\), and \(R^{r}_{rp}\) and \(R^{d}_{rp}\) are the corresponding receiver-class receptor values. High expression support required \(Q_{srp} \ge 0.25\) and both detection fractions \(\ge 0.05\); moderate support required \(Q_{srp} \ge 0.10\) and both detection fractions \(\ge 0.02\). Ligand-receptor summaries and predictions are provided in Tables S38-S40.

### NeuroMorpho morphology validation

External morphology validation used NeuroMorpho records for dentate and cerebellar granule-cell reconstructions with dendritic morphometry, building on public digital neuronal reconstruction resources and morphometric mining methods (Ascoli et al., 2007; Halavi et al., 2012; Polavaram et al., 2014). Records were filtered for granule-cell identity, dentate or cerebellar anatomical assignment, non-missing dendritic stem/branch/length measurements and acceptable reconstruction integrity. Morphological comparisons focused on input-sampling architecture rather than assuming identical dendritic arbor geometry between regions. For hierarchical evidence synthesis, morphology comparisons were converted to bounded MorphologySparseSampling scores. Let \(r_m\) be the dentate-to-cerebellar median ratio for metric \(m\), and let \(\delta_m\) be the Cliff's delta for the dentate-versus-cerebellar comparison:

$$
S_{\mathrm{stems}} = |\delta_{\mathrm{stems}}|,
\qquad
S_{\mathrm{branch}} = 1-\min\!\left(|\log_2 r_{\mathrm{branch}}|,1\right),
$$

$$
S_{\mathrm{compact}} = 0.75|\delta_{\mathrm{compact}}|,
\qquad
S_{\mathrm{length}} = 0.50|\delta_{\mathrm{length}}|.
$$

Similar branch-count scale was scored as shared nontrivial branch complexity, whereas stem, compactness or length differences were scored as branch-specific implementations of a compact granule-like input-sampling principle. NeuroMorpho calibration targets were also converted to scores:

$$
S_{\mathrm{input}} =
1-\min\!\left(\frac{\max(v-1,0)}{15},1\right),
\qquad
S_{\mathrm{branch\ complexity}} = 0.75,
$$

where \(v\) is the median primary-stem target. Morphology summaries are reported in Tables S29 and S30, and morphology evidence units are included in Table S90.

### DANDI spatial-activity validation

Physiology validation used local DANDI 000003 NWB files as a dentate granule-cell activity layer. DANDI 000003 is the public archive version of the Senzai and Buzsaki hippocampal granule-cell and mossy-cell physiology dataset (Senzai and Buzsaki, 2017). NWB files were parsed for units, cell-type labels, behavior/position variables and awake-moving samples. Labeled granule units were analyzed for firing rate, active spatial-bin fraction, spatial information, spatial sparsity and population-vector structure. For hierarchical evidence synthesis, DANDI metrics were converted to bounded ActivitySparsity scores using \(\operatorname{clip}_{[-1,1]}(\cdot)\):

$$
S_{\mathrm{info}} =
\operatorname{clip}_{[-1,1]}\!\left(\frac{I_{\mathrm{GC}}}{0.80}\right),
\qquad
S_{\mathrm{active}} =
\operatorname{clip}_{[-1,1]}\!\left(\frac{1-f_{\mathrm{active,GC}}}{0.50}\right),
$$

$$
S_{\mathrm{info\ delta}} =
\operatorname{clip}_{[-1,1]}\!\left(\frac{I_{\mathrm{GC}}-I_{\mathrm{Pyr}}}{0.40}\right),
\qquad
S_{\mathrm{active\ delta}} =
\operatorname{clip}_{[-1,1]}\!\left(\frac{f_{\mathrm{active,Pyr}}-f_{\mathrm{active,GC}}}{0.30}\right),
$$

$$
S_{\mathrm{PV}} =
\operatorname{clip}_{[-1,1]}\!\left(\frac{\Delta_{\mathrm{far-near}}}{0.30}\right).
$$

Here, \(I\) denotes spatial information in bits per spike, \(f_{\mathrm{active}}\) denotes active spatial-bin fraction, Pyr denotes pyramidal-cell comparator units, and \(\Delta_{\mathrm{far-near}}\) is the median far-minus-near population-vector distance for labeled granule units. In the sparse-coding calibration, DANDI active spatial-bin fraction was used as a light activity proxy and upper-bound spatial-occupancy target. Activity summaries are reported in Tables S31-S33, and activity evidence units are included in Table S90.

### Sparse expansion-coding simulation and empirical calibration

The sparse expansion-coding analysis tested whether granule-like sparse expansion is computationally plausible under resource and activity constraints. The model was motivated by dentate pattern-separation physiology and by recent arguments that pattern-separation metrics must penalize information loss under extreme sparsity (Leutgeb et al., 2007; Bird et al., 2024). Correlated binary input patterns were generated from 18 prototype patterns across \(D_{\mathrm{in}}=64\) input dimensions, with input active fraction \(0.12\) and bit-flip noise rate \(0.08\), yielding 144 related input patterns per replicate. Each simulated architecture was defined by expansion ratio \(R\), input degree \(d\), and target output active fraction \(f\). The number of output units, nonzero projection weight and active output count were:

$$
M = \max\!\left(4,\operatorname{round}(D_{\mathrm{in}}R)\right),
\qquad
w_{\mathrm{nonzero}} = \frac{1}{\sqrt{d}},
\qquad
k = \min\!\left(M,\max\!\left(1,\operatorname{round}(fM)\right)\right).
$$

Each output unit sampled \(d\) input dimensions without replacement. For each input pattern, binary output activity was assigned by a top-\(k\) rule. The simulation grid crossed expansion ratios 0.5, 1, 2, 4, 8 and 16; input degrees 2, 4, 8, 16 and 32; and output active fractions 0.01, 0.03, 0.05, 0.10 and 0.20, with three replicates per grid point.

Simulation performance was quantified from pairwise input and output distances. A small numerical constant, \(\epsilon=10^{-9}\), was used to avoid division by zero. Separation gain was defined as median output Jaccard distance divided by median input Jaccard distance:

$$
G_{\mathrm{sep}}
=
\frac{\operatorname{median}\!\left(D^{J}_{\mathrm{out}}\right)}
{\operatorname{median}\!\left(D^{J}_{\mathrm{in}}\right)+\epsilon}.
$$

Near-pair separation gain, \(G_{\mathrm{near}}\), was computed the same way for input pairs in the lowest quartile of input Jaccard distance. Distance-structure preservation was measured as the Spearman correlation between pairwise input and output distances, \(\rho_{\mathrm{distance}}\). Collapse rate, \(C_{\mathrm{collapse}}\), was the fraction of non-identical input pairs with zero output distance. Mean output entropy, \(H_{\mathrm{out}}\), was the average binary entropy across output units. The activity penalty was:

$$
P_{\mathrm{activity}}
=
\exp\!\left[
-\frac{\left|\log\!\left((f_{\mathrm{obs}}+\epsilon)/0.05\right)\right|}{2}
\right].
$$

Information retention, useful pattern-separation score, relative wiring-activity cost and resource-adjusted useful score were:

$$
I_{\mathrm{retention}}
=
\max\!\left(0,\rho_{\mathrm{distance}}\right)
\left(1-C_{\mathrm{collapse}}\right),
$$

$$
U_{\mathrm{useful}}
=
G_{\mathrm{near}} I_{\mathrm{retention}} P_{\mathrm{activity}},
$$

$$
C_{\mathrm{resource}}
=
\frac{M d}{D_{\mathrm{in}}}
+ f_{\mathrm{obs}} M,
$$

$$
U_{\mathrm{resource\ constrained}}
=
\frac{U_{\mathrm{useful}}}{C_{\mathrm{resource}}+\epsilon}.
$$

Raw calibration and resource-constrained calibration were reported separately because they answer different questions: maximum separation/information performance versus biologically constrained sparse input-expansion plausibility. Simulation and calibration outputs are provided in Tables S34-S37.

## Acknowledgements

Not applicable.

## Competing interests

The author declares no competing or financial interests.

## Author contributions

J.L. conceived the study question, curated and interpreted the biological framework, supervised the computational analyses, reviewed the outputs, prepared the figures and wrote the manuscript.

## Funding

No specific funding was declared.

## Data availability

No new experimental datasets were generated for this study. All primary and supporting transcriptomic, regulatory-compatibility, perturbation, morphology, and physiology resources analyzed here are publicly available from GEO, NeuroMorpho.Org, DANDI, and related public repositories. The strict primary-core GEO accessions are GSE104323, GSE95752, GSE292261, GSE214309, GSE122357, GSE165657, GSE312658, GSE186538, GSE325391, and GSE268609. Additional supporting or validation resources include GSE185277, GSE185553, GSE322785, GSE245367, GSE118987, GSE84786, GSE71916, GSE242199, GSE81962, NeuroMorpho.Org granule-cell reconstruction records, and DANDI archive 000003. Processed summary tables, candidate-gene tables, model outputs, and result-level provenance are provided in the Supplementary Table packet. The accompanying machine-readable Supplementary Table Index separates reader-facing primary tables from supplementary control, audit, archive, and provenance tables.

## Code availability

Analysis code, run-order documentation, curated summary outputs, manuscript-facing figures, and supplementary table packets are available at GitHub (https://github.com/Ghostneuron/granule-cell-convergence) and archived on Zenodo (https://doi.org/10.5281/zenodo.21018501). The archive includes the final analysis scripts, environment information, and mappings between scripts, manuscript figures, supplementary tables, and result-level provenance in the Supplementary Table packet. Large public raw datasets are not redistributed with the code archive; download sources and accessions are provided in the Data Availability statement and supplementary tables.

## References

Artegiani, B., Lyubimova, A., Muraro, M., van Es, J. H., van Oudenaarden, A. and Clevers, H. (2017). A Single-Cell RNA Sequencing Study Reveals Cellular and Molecular Dynamics of the Hippocampal Neurogenic Niche. Cell Rep 21(11), 3271-3284. doi:10.1016/j.celrep.2017.11.050.

Ascoli, G. A., Donohue, D. E. and Halavi, M. (2007). NeuroMorpho.Org: a central resource for neuronal morphologies. J Neurosci 27(35), 9247-9251. doi:10.1523/JNEUROSCI.2055-07.2007.

Bird, A. D., Cuntz, H. and Jedlicka, P. (2024). Robust and consistent measures of pattern separation based on information theory and demonstrated in the dentate gyrus. PLoS Comput Biol 20(2), e1010706. doi:10.1371/journal.pcbi.1010706.

Boldrini, M., Fulmore, C. A., Tartt, A. N., Simeon, L. R., Pavlova, I., Poposka, V., Rosoklija, G. B., Stankov, A., Arango, V., Dwork, A. J., et al. (2018). Human Hippocampal Neurogenesis Persists throughout Aging. Cell Stem Cell 22(4), 589-599.e5. doi:10.1016/j.stem.2018.03.015.

Chen, P. B., Chen, R., LaPierre, N., Chen, Z., Mefford, J., Marcus, E., Heffel, M. G., Soto, D. C., Ernst, J., Luo, C., et al. (2024). Complementation testing identifies genes mediating effects at quantitative trait loci underlying fear-related behavior. Cell Genom 4(5), 100545. doi:10.1016/j.xgen.2024.100545.

Chen, X., Zhong, X., Yue, W., Wang, B., Woo, B., Goodarzi, H., Luo, Z., Lu, Q. R., Flamant, F., Reiter, J. F., et al. (2026). A transcription regulator atlas identifies TOX3 as an Atoh1 coactivator in cerebellar development and tumorigenesis. Proc Natl Acad Sci U S A 123(12), e2527163123. doi:10.1073/pnas.2527163123.

Disouky, A., Sanborn, M. A., Sabitha, K. R., Mostafa, M. M., Ayala, I. A., Bennett, D. A., Lu, Y., Zhou, Y., Keene, C. D., Weintraub, S., et al. (2026). Human hippocampal neurogenesis in adulthood, ageing and Alzheimer's disease. Nature 652(8112), 1264-1273. doi:10.1038/s41586-026-10169-4.

Eriksson, P. S., Perfilieva, E., Bjork-Eriksson, T., Alborn, A. M., Nordborg, C., Peterson, D. A. and Gage, F. H. (1998). Neurogenesis in the adult human hippocampus. Nat Med 4(11), 1313-1317. doi:10.1038/3305.

Flint, J., Heffel, M. G., Chen, Z., Mefford, J., Marcus, E., Chen, P. B., Ernst, J. and Luo, C. (2023). Single-cell methylation analysis of brain tissue prioritizes mutations that alter transcription. Cell Genom 3(12), 100454. doi:10.1016/j.xgen.2023.100454.

Franjic, D., Skarica, M., Ma, S., Arellano, J. I., Tebbenkamp, A. T. N., Choi, J., Xu, C., Li, Q., Morozov, Y. M., Andrijevic, D., et al. (2022). Transcriptomic taxonomy and neurogenic trajectories of adult human, macaque, and pig hippocampal and entorhinal cells. Neuron 110(3), 452-469.e14. doi:10.1016/j.neuron.2021.10.036.

Halavi, M., Hamilton, K. A., Parekh, R. and Ascoli, G. A. (2012). Digital reconstructions of neuronal morphology: three decades of research trends. Front Neurosci 6, 49. doi:10.3389/fnins.2012.00049.

Hochgerner, H., Zeisel, A., Lonnerberg, P. and Linnarsson, S. (2018). Conserved properties of dentate gyrus neurogenesis across postnatal development revealed by single-cell RNA sequencing. Nat Neurosci 21(2), 290-299. doi:10.1038/s41593-017-0056-2.

Jacobs, B., Johnson, N. L., Wahl, D., Schall, M., Maseko, B. C., Lewandowski, A., Raghanti, M. A., Wicinski, B., Butti, C., Hopkins, W. D., et al. (2014). Comparative neuronal morphology of the cerebellar cortex in afrotherians, carnivores, cetartiodactyls, and primates. Front Neuroanat 8, 24. doi:10.3389/fnana.2014.00024.

Kim, S. K., Cherskov, A., Sindhwani, A., Choi, S. H., Kim, H., Li, M. L., Zhang, M., Mato-Blanco, X., Liu, Y., Micali, N., et al. (2026). Human-specific features of the cerebellum and ZP2-regulated synapse development. Cell 189(6), 1802-1819.e28. doi:10.1016/j.cell.2026.02.014.

Kuhn, H. G., Dickinson-Anson, H. and Gage, F. H. (1996). Neurogenesis in the dentate gyrus of the adult rat: age-related decrease of neuronal progenitor proliferation. J Neurosci 16(6), 2027-2033. doi:10.1523/JNEUROSCI.16-06-02027.1996.

Leutgeb, J. K., Leutgeb, S., Moser, M. B. and Moser, E. I. (2007). Pattern separation in the dentate gyrus and CA3 of the hippocampus. Science 315(5814), 961-966. doi:10.1126/science.1135801.

Lewis, P. M., Gritli-Linde, A., Smeyne, R., Kottmann, A. and McMahon, A. P. (2004). Sonic hedgehog signaling is required for expansion of granule neuron precursors and patterning of the mouse cerebellum. Dev Biol 270(2), 393-410. doi:10.1016/j.ydbio.2004.03.007.

Lorente-Echeverria, B., Daaboul, D., Vandensteen, J., Marcassa, G., Naert, W., Vandenbempt, J., Leysen, E., Reverendo, M., Vlaeminck, I., Vervloessem, L., et al. (2025). A dynamic gene regulatory code drives synaptic development of hippocampal granule cells. Sci Adv 11(43), eadx5140. doi:10.1126/sciadv.adx5140.

Lu, J., Wu, Y., Sousa, N. and Almeida, O. F. (2005). SMAD pathway mediation of BDNF and TGF beta 2 regulation of proliferation and differentiation of hippocampal granule neurons. Development 132(14), 3231-3242. doi:10.1242/dev.01893.

Marco, A., Meharena, H. S., Dileep, V., Raju, R. M., Davila-Velderrain, J., Zhang, A. L., Adaikkan, C., Young, J. Z., Gao, F., Kellis, M., et al. (2020). Mapping the epigenomic and transcriptomic interplay during memory formation and recall in the hippocampal engram ensemble. Nat Neurosci 23(12), 1606-1617. doi:10.1038/s41593-020-00717-0.

Parylak, S. L., Qiu, F., Linker, S. B., Gallina, I. S., Lim, C. K., Preciado, D., McDonald, A. H., Zhou, X. and Gage, F. H. (2023). Neuronal activity-related transcription is blunted in immature compared to mature dentate granule cells. Hippocampus 33(4), 412-423. doi:10.1002/hipo.23515.

Peng, J., Sheng, A. L., Xiao, Q., Shen, L., Ju, X. C., Zhang, M., He, S. T., Wu, C. and Luo, Z. G. (2019). Single-cell transcriptomes reveal molecular specializations of neuronal cell types in the developing cerebellum. J Mol Cell Biol 11(8), 636-648. doi:10.1093/jmcb/mjy089.

Polavaram, S., Gillette, T. A., Parekh, R. and Ascoli, G. A. (2014). Statistical analysis and data mining of digital reconstructions of dendritic morphologies. Front Neuroanat 8, 138. doi:10.3389/fnana.2014.00138.

Ramnauth, A. D., Tippani, M., Divecha, H. R., Papariello, A. R., Miller, R. A., Nelson, E. D., Thompson, J. R., Pattie, E. A., Kleinman, J. E., Maynard, K. R., et al. (2025). Spatiotemporal analysis of gene expression in the human dentate gyrus reveals age-associated changes in cellular maturation and neuroinflammation. Cell Rep 44(2), 115300. doi:10.1016/j.celrep.2025.115300.

Senzai, Y. and Buzsaki, G. (2017). Physiological Properties and Behavioral Correlates of Hippocampal Granule Cells and Mossy Cells. Neuron 93(3), 691-704.e5. doi:10.1016/j.neuron.2016.12.011.

Sinnamon, J. R., Torkenczy, K. A., Linhoff, M. W., Vitak, S. A., Mulqueen, R. M., Pliner, H. A., Trapnell, C., Steemers, F. J., Mandel, G. and Adey, A. C. (2019). The accessible chromatin landscape of the murine hippocampus at single-cell resolution. Genome Res 29(5), 857-869. doi:10.1101/gr.243725.118.

Sorrells, S. F., Paredes, M. F., Cebrian-Silla, A., Sandoval, K., Qi, D., Kelley, K. W., James, D., Mayer, S., Chang, J., Auguste, K. I., et al. (2018). Human hippocampal neurogenesis drops sharply in children to undetectable levels in adults. Nature 555(7696), 377-381. doi:10.1038/nature25975.

Spalding, K. L., Bergmann, O., Alkass, K., Bernard, S., Salehpour, M., Huttner, H. B., Bostrom, E., Westerlund, I., Vial, C., Buchholz, B. A., et al. (2013). Dynamics of hippocampal neurogenesis in adult humans. Cell 153(6), 1219-1227. doi:10.1016/j.cell.2013.05.002.

Subash, P., Gray, A., Boswell, M., Cohen, S. L., Garner, R., Salehi, S., Fisher, C., Hobel, S., Ghosh, S., Halchenko, Y., et al. (2023). A comparison of neuroelectrophysiology databases. Sci Data 10(1), 719. doi:10.1038/s41597-023-02614-0.

Wechsler-Reya, R. J. and Scott, M. P. (1999). Control of neuronal precursor proliferation in the cerebellum by Sonic Hedgehog. Neuron 22(1), 103-114. doi:10.1016/s0896-6273(00)80682-0.

Wizeman, J. W., Guo, Q., Wilion, E. M. and Li, J. Y. (2019). Specification of diverse cell types during early neurogenesis of the mouse cerebellum. Elife 8, e42388. doi:10.7554/eLife.42388.

Yao, Z., van Velthoven, C. T. J., Kunst, M., Zhang, M., McMillen, D., Lee, C., Jung, W., Goldy, J., Abdelhak, A., Aitken, M., et al. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature 624(7991), 317-332. doi:10.1038/s41586-023-06812-z.

Yao, J., Dai, S., Zhu, R., Tan, J., Zhao, Q., Yin, Y., Sun, J., Du, X., Ge, L., Xu, J., et al. (2024). Deciphering molecular heterogeneity and dynamics of human hippocampal neural stem cells at different ages and injury states. Elife 12, RP89507. doi:10.7554/eLife.89507.

Zhong, S., Wang, M., Huang, L., Chen, Y., Ge, Y., Zhang, J., Shi, Y., Dong, H., Zhou, X., Wang, B., et al. (2023). Single-cell epigenomics and spatiotemporal transcriptomics reveal human cerebellar development. Nat Commun 14(1), 7613. doi:10.1038/s41467-023-43568-6.

Zhou, Y., Su, Y., Li, S., Kennedy, B. C., Zhang, D. Y., Bond, A. M., Sun, Y., Jacob, F., Lu, L., Hu, P., et al. (2022). Molecular landscapes of human hippocampal immature neurons across lifespan. Nature 607(7919), 527-533. doi:10.1038/s41586-022-04912-w.

## Figure Legends

### Figure 1. Biological question and primary-core study design.

**a,** Schematic comparison of cerebellar and dentate granule-cell morphology. Cerebellar granule cells are shown as compact neurons with short dendrites ending in claw-like input structures and a T-shaped parallel-fiber axon, whereas dentate granule cells are shown with a fan-shaped dendritic tree, dendritic spines, and a mossy-fiber axon. The shared feature is a compact excitatory input-expansion design despite distinct regional anatomy and circuit position. **b,** Working hypothesis: cerebellar granule cells arise through a rhombic-lip/SHH-associated lineage program, whereas dentate granule cells arise through a WNT/PROX1-associated lineage program. The model asks whether these distinct upstream fate programs converge on downstream assembly machinery for neurites, synapses, and excitability to generate related compact neuron designs. **c,** Strict primary-core dataset frame used for the primary analyses. The core includes mouse dentate datasets, cerebellar datasets, and human dentate/hippocampal datasets; scaffold and supporting resources were tracked separately from the strict primary core.

### Figure 2. Ortholog-aware candidate discovery and mechanism-axis organization.

**a,** Analysis workflow. Ten primary datasets were converted to within-sample/gene pseudobulk ranks, mapped through an MGI one-to-one ortholog frame, tested for branch-level and shared candidate effects, and organized into candidate tiers. **b,** Heatmap of prioritized Tier 1 and Tier 2 genes across selected-feature and full-MGI branch-analysis layers. The selected dentate layer includes mouse dentate plus selected human dentate/hippocampal bridge datasets, the full-MGI dentate layer is mouse dentate in the current implementation, and both cerebellar layers include mouse cerebellar datasets plus the human cerebellar validation dataset. Values represent median candidate-background rank deltas; positive values indicate higher relative expression in candidate granule populations. **c,** Mechanism-axis organization of prioritized genes. Candidate genes are grouped into developmental regulatory control, synaptic/excitatory maturation, neurite/cytoskeleton morphogenesis, and axon guidance/adhesion axes, with bars indicating Tier 1, Tier 2, and Tier 3/4 support. The conservative Tier 1 seed set includes GABRA2, GPM6A, KCNK1, NFIA, NFIB, and RFX3; Tier 2 genes add synaptic, excitability, calcium-signaling, and guidance candidates.

### Figure 3. Shared granule-cell programs are constrained by named-comparator specificity.

**a,** Named-comparator mechanism-axis audit. Median within-sample ranks are shown for dentate granule cells, pyramidal comparators, cerebellar granule cells, and Purkinje comparators across the four candidate-gene mechanism axes from Fig. 2c. **b,** Strict-core convergence and named-comparator specificity for the broader five-layer niche/circuit module set. The left plot reports shared cross-branch convergence, summarized from within-branch dentate and cerebellar candidate-background support, with the weaker branch limiting the shared score. Positive values indicate modules supported in both granule branches; near-zero or negative values indicate branch-specific or mixed support. The right plot reports local named-comparator deltas: dentate granule minus pyramidal comparator or cerebellar granule minus Purkinje comparator. Downstream neurite/morphology and synaptic/excitability modules show stronger cross-branch convergence than upstream fate/niche modules.

### Figure 4. Identity-coupled transcriptomic assembly configuration.

**a,** Definition of the transcriptomic configuration score. Construction balance is defined as downstream neurite/synaptic module strength minus niche/progenitor module strength. Regional fate polarity is defined as branch-matched fate strength minus opposed fate strength. The combined configuration score is the sum of construction balance and regional fate polarity. **b,** Local named-comparator configuration test. The heatmap shows module ranks for dentate granule, pyramidal comparator, cerebellar granule, and Purkinje comparator classes; the accompanying bar plot shows configuration-component deltas for dentate granule versus pyramidal and cerebellar granule versus Purkinje contrasts. **c,** Primary-core validation of the configuration score. Candidate-background configuration deltas are shown across full-MGI ortholog contrasts (F01-F08) and selected-feature contrasts (S01-S55); positive bars indicate higher configuration score in candidate granule populations than local background classes. The inset summarizes layer-level median effects for construction balance, regional fate polarity, and combined configuration. **d,** Driver decomposition. The stacked bars classify primary-core and local named-comparator contrasts by whether positive configuration is supported by both construction and fate components, primarily fate polarity, primarily construction balance, weak mixed support, or no positive configuration. The center panel shows construction and fate-polarity component deltas by contrast, and the right panel summarizes base module deltas. Together, these analyses support an identity-coupled assembly-configuration model rather than a single morphology-only gene signature.

### Figure 5. Stage-windowed niche signaling and computational constraints refine the convergence model.

**a,** Pathway-readiness audit. The three stacked plots show median candidate-background rank deltas for dentate and cerebellar branches across curated pathway modules, composite biological signatures, and ligand-receptor readiness contrasts. Positive bars indicate higher relative scores in candidate granule populations than local background classes. Dentate candidates show broad maturation and permissive-readiness signals, whereas cerebellar candidates show their clearest positive signal in SHH/PTCH/GLI. **b,** Fitted stage/pseudotime window analysis. Each small plot shows a signature score against normalized stage or pseudotime order, with points representing observed dataset/group scores and fitted branch trends summarizing dentate and cerebellar behavior using a quadratic stage-window model. The four signatures test the historical TGF-beta/BDNF 2005 mechanism, differentiation/stop signaling, neurogenic/permissive signaling, and stop-minus-permissive balance. **c,** Empirical calibration of the sparse expansion-coding model against NeuroMorpho morphology and DANDI dentate activity summaries. The upper scatter plot shows simulation grid points by model input degree and observed output active fraction; color encodes the empirical calibration score. The lower bar plot ranks named architecture families by median empirical calibration score. Raw empirical calibration favors dense expansion, whereas resource- and morphology-constrained calibration shifts the preferred family toward intermediate or granule-like sparse expansion. **d,** Revised final model. Distinct cerebellar and dentate lineages pass through stage/pseudotime-dependent maturation windows and identity-coupled assembly configuration. Branch-specific secreted cues and resource-constrained sparse expansion jointly shape convergence toward compact granule-cell morphology.

### Figure 6. Focused niche sender-receiver ligand-receptor prediction.

**a,** Cell coverage for the focused sender-receiver ligand-receptor screen. Cerebellar analyses used GSE122357 with granule-lineage receiver cells and Purkinje, astroglial/Bergmann-proxy, microglial, and endothelial sender classes. Dentate SGZ analyses used GSE104323 with granule-lineage receiver states and astrocyte, PVM/microglia-macrophage proxy, and endothelial sender classes. **b,** Top supported sender-pair summaries ranked by median supported ligand-receptor (LR) expression score. High-scoring examples include cerebellar Purkinje-to-granule precursor IGF1-IGF1R, SHH-PTCH1/SMO, guidance-related interactions, and dentate SGZ astrocyte/PVM/endothelial interactions involving SEMA6A-PLXNA2/4, APOE-LRP1/LDLR, C1QA-LRP1, and JAG1-NOTCH1. **c,** Heatmap of median sender-to-granule-lineage LR expression scores by sender class and pathway. Scores require ligand expression in the sender and receptor expression in the receiver, providing a directional niche hypothesis layer. Bergmann glia are represented by an astroglial/Bergmann proxy in GSE122357, and SGZ microglia by the PVM/microglia-macrophage proxy in GSE104323.

### Figure 7. Evidence-weighted hypothesis comparison supports an integrated convergence model.

**a,** Relationship among the three hypotheses tested in this study. H1 proposes a hidden shared granule-cell fate identity; H2 proposes identity-coupled transcriptomic assembly convergence; H3 proposes stage, niche and sparse-expansion circuit constraints. **b,** Observed evidence-term scores from the hierarchical synthesis. Y, granule configuration; F, branch-matched regional fate polarity; C, construction balance; I, fate-construction coupling; T, stage/pseudotime; N, niche signal; E, regulatory compatibility; M, morphology sparse sampling; A, activity sparsity; R, circuit-resource constraint. **c,** Hypothesis prediction coefficients used in the support-index calculation. Positive coefficients indicate expected support, negative coefficients indicate mismatch, and zero indicates that the term is not used to discriminate that hypothesis. **d,** Formula-calculated support indices. H1 scores below neutral, whereas H2, H3 and the H2+H3 synthesis score above neutral.

### Graphical Abstract. Integrated working model for granule-cell convergence.

The graphical abstract summarizes the evidence pipeline and biological model. The pipeline integrates strict primary-core transcriptomic evidence, ortholog rank-meta candidate tiers, trajectory/stage-window analysis, sparse-expansion computation, and external morphology/activity validation. The biological model proposes that cerebellar and dentate granule-cell lineages retain distinct upstream identities, pass through stage-windowed maturation-readiness states, and converge on identity-coupled transcriptomic assembly configuration. A shared construction toolkit for neurites, synapses, excitability, and guidance is filtered by circuit-level constraints favoring compact sparse expansion, yielding convergent compact excitatory granule-cell morphology despite distinct developmental origins.
