# Comparator-relative molecular recurrence in distinct dentate and cerebellar granule cells

Jie Lu

Affiliation: Independent Researcher, San Francisco, California, USA

Correspondence: ghostneuron@gmail.com

ORCID: 0000-0001-6843-9720

**Short title:** Molecular recurrence in granule cells

## Abstract

Neuronal names can group cells by morphology even when lineage and circuit context differ. Dentate and cerebellar granule cells are small, densely packed excitatory neurons in distinct regional systems. This study tested whether their shared designation corresponds to recurrent molecular organization relative to local projection-neuron comparators. Ten transcriptomic datasets were analyzed with within-dataset ranks and one-to-one ortholog mapping. A mouse-first screen prioritized six Tier 1 genes (GPM6A, NFIB, NFIA, KCNK1, RFX3 and GABRA2) and nine Tier 2 genes. All 15 retained directional support in every testable leave-one-dataset-out analysis. Collapsing nested comparisons to one median per independent dataset yielded seven of seven positive configuration deltas (median 1.000, bootstrap 95% confidence interval 0.625-1.500; one-sided exact P=0.0078). A five-module comparison was directional but not conventionally significant (exact P=0.10). In an external adult mouse Allen common-expression matrix, direct cerebellar-versus-dentate similarity did not favor downstream assembly over upstream fate and niche genes (Spearman correlations -0.379 and -0.329). Branch-local comparisons reproduced the predeclared direction for five of six Tier 1 genes, with GABRA2 as the exception, and for ten of 15 Tier 1 plus Tier 2 genes. Both fractions exceeded expression- and detection-matched null distributions (empirical P=0.00010 for each set), although exact sign tests across six or 15 genes were not conventionally significant. Only KCND2 exceeded every tested comparator in both targets. Thus, selected molecular features recur relative to regional comparators even though the two populations do not form a broad common adult molecular class.

## Keywords

Allen Brain Cell Atlas; cerebellum; comparative neuroanatomy; dentate gyrus; granule cell; single-cell transcriptomics; systems neuroscience

## Key points

- Local granule-versus-projection-neuron contrasts recover a small recurrent gene set across heterogeneous datasets.
- Adult Allen data reproduce candidate direction but not broad dentate-cerebellar molecular similarity.
- The data do not support a common adult molecular class for the two populations.

## 1. Introduction

"Granule cell" began as a histological description rather than a molecular or lineage classification. In the neuroanatomical tradition represented by Ramon y Cajal, the term was applied to small neurons whose densely packed somata gave a grain-like appearance. Dentate and cerebellar granule cells are both excitatory neurons positioned in large input-processing layers, but they participate in different systems. Cerebellar granule cells relay mossy-fiber input through parallel fibers within the Purkinje-cell circuit, whereas dentate granule cells transform entorhinal input before projecting to CA3 (Marr, 1969; Leutgeb et al., 2007). Their shared name therefore poses an organizational question: does repeated compact-neuron architecture carry a recurring molecular configuration, or do the two regions realize superficially similar forms through largely regional programs?

Their developmental histories are distinct. Cerebellar granule cells arise from the hindbrain rhombic lip, expand in a transient external granule layer and mature under local signaling that includes Sonic hedgehog (SHH) from Purkinje cells (Wechsler-Reya and Scott, 1999; Lewis et al., 2004; Wizeman et al., 2019). Dentate granule cells arise in a telencephalic hippocampal lineage and continue to be generated postnatally and during adulthood in rodents (Altman and Das, 1965; Altman and Das, 1966; Kuhn et al., 1996; Hochgerner et al., 2018). Morphological resemblance alone does not establish a common cell type.

The possibility of shared developmental elements nevertheless has a substantial history. Altman's work placed postnatal microneuron production and secondary germinal matrices at the center of granule-cell development. Yang, Zhong and Heintz subsequently identified RU49/Zipro1 expression in developing cerebellar, dentate and olfactory granule-cell populations and proposed that distinct lineages might reuse part of a developmental program (Yang et al., 1996). This proposal does not require a recent common migratory progenitor. It instead raises a more restricted question: can neurons specified in different regional lineages repeatedly recruit a subset of molecular features during differentiation?

Earlier culture work showed that cerebellar conditioned medium could suppress proliferation and promote differentiation of hippocampal granule-lineage cells under defined conditions, with TGF-beta2, BDNF, SMAD and MAPK signaling implicated as candidate mediators (Lu et al., 2005). This observation provides historical evidence of cross-tissue responsiveness, not a molecular validation target for the present recurrence analysis.

The comparative design addressed three progressively stricter questions. First, do candidate genes recur across independent dentate and cerebellar transcriptomic datasets after species and platform effects are reduced by within-dataset ranking? Second, does the result remain after nested contrasts are collapsed to independent datasets and individual datasets are withheld? Third, does a preselected candidate set generalize in a common adult mouse expression matrix containing regional projection-neuron comparators and additional compact or granule-named populations? The third question is necessary because genes associated with neuronal differentiation, synaptic maturation or excitability may be shared by many neuronal classes.

The discovery analyses test recurrence across the developing and adult resources that were available. The Allen analysis tests a narrower, predeclared prediction: whether the discovery-layer direction generalizes to an adult mouse common platform. Failure in Allen would reject adult generalization, but neither a positive nor a negative adult result can establish or exclude a transient developmental program.

The analyses support comparator-relative molecular recurrence rather than a shared adult cell identity. A limited candidate set recurs when cerebellar granule cells are compared with Purkinje cells and dentate granule cells with CA1/CA3 pyramidal cells. The same adult Allen data do not support preferential direct similarity of broad downstream assembly modules, and almost none of the candidates are exclusive to the two target populations. This distinction provides a molecular framework for comparing similarly named neurons without erasing their regional organization.

## 2. Materials and Methods

### 2.1. Dataset selection and analytical roles

Public transcriptomic resources were selected for source annotation, availability of dentate or cerebellar granule-lineage cells, usable regional background or comparator populations, age coverage and expression-matrix accessibility. The primary frame contained four mouse dentate, two mouse cerebellar, three human dentate or hippocampal, and one human cerebellar dataset. Supporting human hippocampal resources were used for label construction but were not counted as additional primary datasets. The evidentiary frame was restricted to transcriptomic datasets and the external Allen common matrix, allowing the same comparator-relative question to be evaluated across resources. Complete dataset metadata and machine-readable analysis tables are available in the public code repository.

### 2.2. Annotation, pseudobulk summaries and within-dataset ranks

Candidate granule cells, precursors, regional neuronal comparators and local background groups were assigned from source annotations and branch-specific marker rules. Expression was summarized by dataset/sample group before cross-group comparison. For each dataset layer, group-level median log1p expression was ordered within the local sample and converted to percentile ranks, using average ranks for ties. Candidate-minus-background rank deltas were then calculated within the same dataset and feature layer. Positive deltas therefore indicate higher relative expression in the candidate population; they are not log fold changes and are not comparable to absolute expression differences from another platform.

The transformation was chosen because raw expression magnitudes differ with chemistry, sequencing depth, genome annotation and preprocessing. It also avoids pooling cells from different studies as if they were biological replicates. Dataset-level summaries, rather than cell counts, were used for inferential tests.

### 2.3. Ortholog-aware screens and candidate tiering

Cross-species comparisons used one-to-one Mouse Genome Informatics ortholog mappings (RRID:SCR_006460). The selected-feature screen contained curated candidate, module and bridge-object genes. The full-matrix screen used broader available expression matrices remapped to the same ortholog frame. For a gene, screen and branch, nominal directional support required at least two datasets, median dataset rank delta greater than zero, at least 75% positive dataset deltas and a best one-sided dataset-level P value no greater than 0.25. The liberal nominal threshold was used for discovery across small dataset counts and was never interpreted as confirmatory significance. Shared support required the rule to be met separately in dentate and cerebellar branches.

Candidate ordering used the deterministic priority score

$$
P_{\mathrm{formal}}
=3N_{\mathrm{FDR10}}
+2N_{\mathrm{nominal}}
+N_{\mathrm{replication}}
+\sum_b \max(0,\widetilde{\Delta}_{\mathrm{rank},b})
+\min_b(D_b),
$$

where $b$ indexes available screen-branch tests, $N_{\mathrm{FDR10}}$, $N_{\mathrm{nominal}}$ and $N_{\mathrm{replication}}$ count tests meeting their respective criteria, $\widetilde{\Delta}_{\mathrm{rank},b}$ is the branch median rank delta and $D_b$ is median candidate-cell detection. FDR10 denotes a Benjamini-Hochberg false-discovery-rate threshold of 10%. Fisher-combined P values across available branches and their Benjamini-Hochberg-adjusted values were retained as summary fields, not sole inclusion criteria.

Tier 1 required five evidence features: formal shared support, same-symbol one-to-one orthology, both-screen support, central mechanism-priority membership and dataset-robust support in all four selected/full-matrix by dentate/cerebellar combinations. The same-symbol requirement was a conservative mapping-quality rule rather than a biological criterion. Among 64 one-to-one ortholog rows with nonidentical human and mouse symbols in the broader audit, none had nominal support in all four screen-branch combinations or in both screens. Tier 2 required both-screen formal support and central wiring, synaptic, excitability, calcium or guidance relevance but not the four-of-four Tier 1 designation. Later exploratory tiers are retained in the repository data tables for completeness but were not used in the external candidate test.

### 2.4. Species stratification and leave-one-dataset-out analysis

Tier 1 and Tier 2 rank deltas were summarized separately by species, screen and branch. The mouse-only audit required positive direction in selected-feature dentate, selected-feature cerebellar, full-matrix dentate and full-matrix cerebellar strata. For leave-one-dataset-out analysis, each dataset was removed in turn from a gene-screen-branch combination. A retained combination was called directionally robust when at least two datasets remained, the median rank delta was positive and at least 75% of remaining dataset deltas were positive. Only testable leave-one-dataset-out combinations were summarized.

Discovery-layer matched-null sets were drawn from noncandidate genes with the same number of available branches, a similar broad functional group when possible and nearest candidate-cell detection. Ten thousand matched sets were sampled. Empirical one-sided P values used $(1+k)/(1+N)$, where $k$ is the number of null statistics at least as large as the observed statistic and $N=10{,}000$. Because matching and candidate selection used discovery-layer information, this analysis was interpreted as selection-bias sensitivity rather than validation.

### 2.5. Dataset-level configuration and module analyses

Five curated modules represented cerebellar fate/rhombic-lip/SHH, dentate fate/WNT/PROX1, shared neurogenic niche/progenitor state, downstream neurite/morphology and downstream synaptic/excitability. Complete module definitions are available with the repository data tables. For eligible candidate-background contrasts, construction balance was

$$
B_C=\bar{R}_{\mathrm{downstream}}-R_{\mathrm{niche}},
$$

regional fate polarity was

$$
P_F=R_{\mathrm{branch\ matched\ fate}}-R_{\mathrm{opposed\ fate}},
$$

and the combined configuration score was

$$
S_Q=B_C+P_F.
$$

Nested contrasts were collapsed to the median $\Delta S_Q$ for each dataset. The overall direction was tested with a one-sided exact binomial sign test after zero values were removed. The 95% confidence interval for the median was obtained from 10,000 bootstrap samples of the seven datasets.

For the module-level comparison, the two downstream module medians were compared with the three upstream or niche module medians using an exact one-sided Mann-Whitney test. Cliff's delta was reported as an effect-size summary. Given the two-versus-three module design, the analysis was treated as descriptive.

### 2.6. Allen common-matrix external analysis

The external analysis used the official Consensus-WMB-Macosko-10X log2 expression matrices for cerebellum, hippocampal formation, isocortex and olfactory bulb, together with Allen cluster membership and integrated taxonomy tables from the Allen Brain Cell Atlas (RRID:SCR_024440; Yao et al., 2023). Allen data were not used to choose the candidates, but the tested positive direction was inherited from discovery. Subclass mappings defined cerebellar granule, Purkinje, mature and immature dentate granule, CA1/ProS, CA3, cortical L4/5 IT, mature olfactory-bulb GABA-proxy and immature olfactory-bulb GABA-proxy populations. The olfactory-bulb populations were included as falsifying comparators for recurrence across another granule-named or compact inhibitory class. Population names describe transcriptomic subclasses, and the proxy labels do not assert morphological identity for every cell.

For each matrix, cells were grouped by population and library. Library-population groups with fewer than 50 cells were excluded. Mean Allen log2 expression and detection fraction were calculated for the predeclared candidate and module genes. All uncertainty calculations used libraries as replicate units; cell counts were reported only as coverage.

For module scores, each gene was standardized across retained library-population groups,

$$
z_{g,u}=\frac{x_{g,u}-\bar{x}_g}{s_g},
$$

and the score for module $m$ in unit $u$ was the mean of its available gene z scores,

$$
M_{m,u}=\frac{1}{|G_m|}\sum_{g\in G_m}z_{g,u}.
$$

The cerebellar local contrast was cerebellar granule minus Purkinje. The dentate local contrast was mature dentate granule minus the equal-weighted mean of CA1/ProS and CA3. Equal weighting prevented the larger CA1 library count from defining the hippocampal reference. Sensitivity fields also compared mature dentate granule cells separately with CA1/ProS and CA3 and recorded whether each candidate remained positive against both references. Library bootstrap intervals used 2,000 resamples within each population.

Direct population similarity was calculated as the Spearman correlation between mean gene-z-score vectors. Upstream and downstream gene sets were evaluated separately. Nested bootstrap intervals resampled libraries within each population and genes within the tested set 2,000 times. Direct similarity and branch-local deltas were kept separate because the former asks whether two target profiles resemble each other, whereas the latter asks whether both differ in the same direction from their regional comparators.

For the external matched-null analysis, each Tier 1 or Tier 2 candidate was paired with its ten nearest noncandidate curated genes in two-dimensional standardized space defined by overall mean expression and mean detection. Ten thousand null sets were generated by drawing one matched gene for each candidate. The primary null statistic was the fraction of genes with both local deltas positive; the median of the smaller local delta was analyzed as a second statistic. Exact one-sided binomial sign tests across the six or 15 candidate genes were reported alongside empirical matched-null P values.

For cross-region transfer, only complete within-library population blocks were retained. The mature analysis included 63 cerebellar granule-Purkinje library pairs and 17 dentate granule-CA1-CA3 library triplets; the immature-dentate sensitivity analysis included nine triplets. For each predeclared feature panel, an L2-regularized, class-weighted logistic model was fitted in one region after standardizing genes using training-region means and standard deviations. The fitted model was applied to the other region without refitting or threshold adjustment. Performance was summarized by ROC AUC, fixed-threshold balanced accuracy and the within-library target-minus-mean-comparator decision margin. Training-label nulls used 2,000 permutations in which one target label was reassigned within each library block. Bidirectional support required both transfer directions; the intersection-union P value was the larger directional permutation P value.

Matched configuration nulls used 2,500 background genes present in all four Allen matrices. For each tested panel, 2,000 same-sized panels were assembled without replacement from genes nearest to the tested features in overall mean expression and detection. The test statistic was the smaller AUC across the two transfer directions. A complementary contrast-concordance analysis calculated the mean within-library cerebellar-granule-minus-Purkinje and dentate-granule-minus-mean-CA1/CA3 difference for every gene. Cosine concordance, Spearman correlation, same-sign fraction and dual-positive fraction were compared with the same matched panels. Benjamini-Hochberg correction was applied across the six feature panels. The adult immature-versus-mature sensitivity analysis used paired Wilcoxon tests across the nine libraries represented in both dentate states.

### 2.7. Statistical conventions

All tests were one-sided only when a positive directional hypothesis was specified in advance. Exact tests were used for the small dataset or module counts. Bootstrap intervals are percentile intervals and describe uncertainty under resampling of the stated unit; they are not corrections for dataset selection. Coverage may be imperfect when fewer than 20 library units are available, so those intervals were treated as descriptive. No cell-level P value was used as evidence of replication. Unless otherwise specified, statistical calculations used Python (RRID:SCR_008394), pandas (RRID:SCR_018214), NumPy, SciPy (RRID:SCR_008058) and scikit-learn (RRID:SCR_002577). Complete analysis code, environment information, source-accession metadata and intermediate tables are provided in the public code repository.

## 3. Results

### 3.1. Comparative framework separates regional identity from recurrent local shifts

Dentate and cerebellar granule cells were treated as distinct regional cell classes from the outset (Fig. 1a,b). The primary comparison was relational: cerebellar granule populations were evaluated against cerebellar references that included Purkinje cells, whereas dentate granule populations were evaluated against hippocampal references that included CA1/CA3 pyramidal cells. Direct similarity between the two granule populations was reserved for the external common-matrix analysis. The study did not test migration of a recent common progenitor.

The discovery frame contained ten datasets: four mouse dentate resources (GSE104323, GSE95752, GSE292261 and GSE214309), two mouse cerebellar resources (GSE122357 and GSE312658), one human cerebellar resource (GSE165657), and three human dentate or hippocampal resources (GSE186538, GSE325391 and GSE268609) (Fig. 1c). The mouse dentate arm included postnatal lineage and maturation datasets (Hochgerner et al., 2018), a synaptic-development resource (Lorente-Echeverria et al., 2025) and an activity-state comparison of immature and mature granule cells (Parylak et al., 2023). The cerebellar arm included developing mouse and human resources (Peng et al., 2019; Zhong et al., 2023) and a mouse perturbation dataset used only after within-dataset normalization (Chen et al., 2026). The human hippocampal arm supplied taxonomy, adult dentate and aging-context bridges (Franjic et al., 2022; Ramnauth et al., 2025; Disouky et al., 2026). The Zhou et al. lifespan resources were used for label construction but were not counted as additional primary datasets (Zhou et al., 2022). Because only one human cerebellar dataset met the primary criteria, human evidence was not used to make a species-level dentate-versus-cerebellar claim.

The analytical sequence comprised mouse-first discovery, dataset-level sensitivity analysis, external testing in a common expression matrix and interpretation constrained by both positive and negative findings (Fig. 1d). This design separated three questions that are often conflated by cell-type nomenclature: recurrence relative to local comparators, direct cross-region molecular similarity and transfer of a weighted molecular configuration.

### 3.2. Cross-dataset analysis prioritizes a small reproducible candidate set

Within each dataset, candidate granule populations were compared with local background populations using rank differences rather than pooled expression values. Two routes were retained. The selected-feature route maximized coverage of curated candidate and module genes, including reduced human bridge objects. The full-matrix route used broader matrices mapped through one-to-one MGI orthologs where available. These routes answer related but different questions: the first tests predefined features across more resources, whereas the second reduces dependence on feature preselection.

Tiering was used to organize discovery evidence, not to claim validation. Tier 1 comprised genes with formal shared support, one-to-one same-symbol orthology, support in both selected-feature and full-matrix routes, mechanism-priority membership and dataset-robust support in all four screen-branch combinations. This produced GPM6A, NFIB, NFIA, KCNK1, RFX3 and GABRA2. Tier 2 added nine genes with both-screen formal support and central wiring, synaptic, excitability, calcium or guidance relevance: ROBO2, PPP3CA, CACNA2D1, KCNJ6, GABRB3, GRIN2B, KCNJ3, KCND2 and STXBP5L.

All 30 mouse candidate-by-screen combinations, representing 15 genes across the selected-feature and full-matrix routes, had positive median rank deltas in both dentate and cerebellar branches (Fig. 2a). A mouse-only support audit likewise found positive direction for every candidate in all four selected/full-matrix by dentate/cerebellar combinations. The 15-gene set was therefore not dependent on human-only evidence. The value 0.5 for many cells reflects the coarse dataset-level rank scale and should not be interpreted as a common expression effect size. Human bridge results were less balanced and were retained as conservation context rather than as an independent cross-species replication.

The candidate direction did not depend on a single discovery dataset. Each Tier 1 or Tier 2 gene remained directionally supported in every leave-one-dataset-out combination that retained at least two datasets in the relevant branch and screen. Minimum leave-one-dataset-out median deltas ranged from 0.375 for ROBO2 to 0.5 for the other 14 candidates (Fig. 2b). As a selection-bias sensitivity analysis, not a validation test, the selected candidates were compared with genes matched for available branches, broad functional group and detection. For the mean minimum branch median delta, the observed value was 0.5 for Tier 1 and Tier 1 plus Tier 2, compared with null medians of 0.083 and 0.067, respectively (empirical P=0.00010 for each set; Fig. 2c). The small P values partly reflect the directional selection rule and quantify internal stability rather than external support.

The selected genes span regulatory and morphogenesis-associated factors, membrane or cytoskeletal features and synaptic or excitability genes (Fig. 2d). Their functional diversity is consistent with recurrence during differentiation or maturation, but the classification itself does not establish a coordinated program.

### 3.3. Dataset-level recurrence is robust, whereas broad module convergence is not

Multiple candidate-background contrasts could be derived from a single dataset. To preserve the independent dataset as the inferential unit, configuration deltas were collapsed to one median per eligible dataset. Seven datasets contained the module coverage and contrast structure required for this analysis. Every dataset median was positive (seven of seven), with an overall median of 1.000 and a dataset-bootstrap 95% confidence interval of 0.625-1.500 (one-sided exact sign test P=0.0078; Fig. 3a). With only seven datasets, the sign test depends on unanimity: one reversal would increase the exact one-sided P value to 0.0625, and the effect magnitude remains imprecisely estimated. The direction was positive in each regional and species stratum, but the stratum-specific exact tests were underpowered: three cerebellar datasets, four dentate datasets, three human datasets and four mouse datasets (Fig. 3b).

The configuration score combined branch-matched regional fate polarity with downstream-construction-over-progenitor balance. It therefore measures whether a candidate population simultaneously retains its regional identity and shifts toward selected neuronal-construction features relative to its local background. The positive dataset-level result supports reproducible organization of the predeclared features; it does not by itself show that cerebellar and dentate cells become directly similar.

A separate analysis compared five curated modules. Two downstream modules had a median convergence delta of 0.5, whereas three fate or niche modules had a median of -0.25. The resulting difference was 0.75 with Cliff's delta 1.0, but the exact one-sided Mann-Whitney test across two versus three modules was P=0.10 (Fig. 3c). This ordering is descriptive and does not establish downstream-module convergence at the module-level unit of inference.

### 3.4. Adult Allen data separate local recurrence from direct molecular similarity

External testing used the Allen Consensus-WMB-Macosko-10X adult mouse matrices and integrated taxonomy (Yao et al., 2023). Target and comparator populations were identified by Allen subclass labels. After requiring at least 50 cells per library-population group, the analysis included 121 cerebellar-granule libraries, 29 mature-dentate libraries, 12 immature-dentate libraries, 63 Purkinje libraries, 49 CA1/ProS libraries, 20 CA3 libraries, 98 cortical L4/5 IT libraries, 24 mature olfactory-bulb GABA-proxy libraries and 28 immature olfactory-bulb GABA-proxy libraries (Fig. 4a). The olfactory-bulb GABA subclasses were included as falsifying comparators to ask whether the signal simply extended to another granule-named or compact inhibitory population; their taxonomy does not establish classical granule-cell morphology for every included cell.

The common-matrix module profiles preserved strong regional distinctions (Fig. 4b). Cerebellar fate/rhombic-lip/SHH scores were higher in cerebellar granule cells than Purkinje cells (mean z-score delta 0.757), whereas dentate fate/WNT/PROX1 scores were higher in mature dentate granule cells than the equal-weighted CA1/CA3 reference (delta 0.495). Synaptic/excitability scores were also positive in both local contrasts (0.637 and 0.138). In contrast, the neurite/morphology module was nearly unchanged in the cerebellar comparison (0.039; bootstrap interval included zero) and lower in mature dentate granule cells than the CA1/CA3 reference (-0.383). The shared neurogenic-niche module was lower in cerebellar granule cells than Purkinje cells (-0.425) and near zero in the dentate comparison (0.015).

The gene-level branch-local test was more favorable than the broad module test. Five of six Tier 1 genes were higher in both cerebellar-granule-minus-Purkinje and dentate-granule-minus-mean-CA1/CA3 comparisons: GPM6A, NFIB, NFIA, KCNK1 and RFX3. GABRA2 was not positive in either local contrast. The same dual-positive criterion was met by ten of 15 Tier 1 plus Tier 2 genes (Fig. 4c). Exact one-sided sign tests across six and 15 preselected genes gave P=0.109 and P=0.151, respectively. Both fractions exceeded 10,000 expression- and detection-matched sets of noncandidate curated genes (empirical P=0.00010 for each set). This external comparison tests generalization of a direction selected in discovery, not a blind validation of an unselected gene set.

The equal-weighted CA1/CA3 reference did not create the overall result, although one Tier 1 gene depended on that composite. Four of six Tier 1 genes and nine of 15 Tier 1 plus Tier 2 genes were positive against Purkinje and against CA1 and CA3 considered separately. GPM6A was positive against the composite and CA3 but not against CA1.

Candidate recurrence was relative rather than exclusive. When both targets had to exceed every tested comparator, only KCND2 met the criterion. Quantitative expression profiles across all nine Allen populations show that NFIA, NFIB, RFX3 and the other prioritized genes were also strongly expressed in one or more comparator populations. The data therefore do not define a granule-cell barcode.

Finally, direct cerebellar-versus-dentate similarity provided a deliberately different test. Across the downstream assembly genes, the Spearman correlation between adult cerebellar and mature dentate profiles was -0.379; across upstream fate and niche genes, it was -0.329. The downstream-minus-upstream difference was -0.050 (Fig. 4d). Thus, the external common-matrix data do not reproduce preferential broad downstream convergence. The positive local-contrast result and negative direct-similarity result are compatible because they ask different questions: selected features can move in the same direction relative to regional comparators even when the complete adult target profiles remain dissimilar.

A stricter cross-region transfer test asked whether a multigene pattern learned from dentate granule versus CA1/CA3 libraries could identify cerebellar granule versus Purkinje libraries, and vice versa. Tier 1 ranked the held-out populations in the expected direction (area under the receiver-operating-characteristic curve, AUC, 1.000 and 0.901), but the minimum bidirectional AUC was not exceptional under library-blocked training-label permutations (intersection-union P=0.296) or expression- and detection-matched genomic panels (P=0.360). The matched-panel minimum AUC had a median of 0.844 and an interquartile range of 0.526-0.939, placing the observed value of 0.901 within a common null range. The fixed training-derived decision threshold also failed to transfer. Direct contrast analysis clarified the distinction: five of six Tier 1 genes were positive in both local contrasts, but their relative contrast magnitudes were weakly correlated across regions (Spearman rho=0.143; matched-panel P=0.557). The broad downstream panel showed reversed cross-region ranking in both directions (AUC 0.000 and 0.000) and anti-aligned regional contrasts (rho=-0.179). Same-direction candidate recurrence therefore does not establish a shared weighted configuration of the broad neuronal-construction toolkit.

## 4. Discussion

### 4.1. Comparator-relative recurrence, not shared molecular identity

The principal discovery is that molecular recurrence is detectable in a comparator-relative frame but not as broad direct adult similarity. Across heterogeneous discovery datasets, the Tier 1 and Tier 2 genes retained their direction when individual datasets were removed. In the external Allen common matrix, most candidates were higher in cerebellar granule cells than Purkinje cells and in dentate granule cells than CA1/CA3 pyramidal cells. The Allen data source was separate from candidate selection, but the predicted direction was inherited from discovery. The result is therefore an external generalization check rather than a blind validation.

The external analysis also sets the boundary of that conclusion. Cerebellar and mature dentate granule cells did not become directly similar across broad adult downstream modules, and only KCND2 was higher in both targets than every tested comparator. The data support comparator-relative recurrence, not a shared adult molecular identity. This distinction is biologically meaningful. A differentiating neuron can reuse a regulator, ion channel or membrane-associated factor relative to another cell class in its region without converging globally with a neuron from another brain territory.

The cross-region transfer analysis sharpens this boundary further. Tier 1 genes were enriched for positive local contrasts, but their relative effect sizes did not form a null-exceeding transferable configuration. Broad neurite and synaptic features were anti-aligned between the mature regional contrasts. The present evidence therefore concerns recurrence of selected components and their direction, not a distinctive molecular recipe that explains granule-cell morphology.

### 4.2. Implications for neuronal classification and regional organization

The analysis distinguishes two forms of molecular resemblance. Direct similarity asks whether two target populations occupy comparable expression states. Comparator-relative recurrence asks whether each target differs from a biologically relevant reference population in the same direction. The latter can recur while absolute profiles remain region-specific. For neurons named from histology, this distinction avoids inferring common identity from morphology while still allowing selected molecular relations to be compared across regions.

Regional comparators are therefore part of the biological definition of the contrast. Purkinje and pyramidal neurons differ in lineage, morphology and circuit role, so the two local contrasts are not interchangeable measures of an intrinsic granule-cell property. Their value is narrower: they ask whether selected features distinguish input-layer granule cells from major projection neurons within each regional system. The Allen analyses against CA1 and CA3 separately, cortical L4/5 intratelencephalic neurons and olfactory-bulb proxy populations show that this relation is neither comparator-free nor exclusive to the two targets.

### 4.3. Circuit architecture as a testable explanation for morphological recurrence

The failure to identify a transferable downstream molecular configuration does not make the morphological resemblance necessarily accidental. Cerebellar and dentate granule cells participate in different circuits, but both occur in large input-transformation layers composed of numerous compact neurons with restricted input sampling. Classic cerebellar theories assigned the mossy-fiber-granule-cell arrangement a role in pattern discrimination or expansion recoding (Marr, 1969; Albus, 1971), whereas dentate recordings and information-theoretic analyses support pattern separation in the hippocampal circuit (Leutgeb et al., 2007; Bird et al., 2024). More general analysis of expansion-layer networks shows that sparse synaptic connectivity can preserve useful high-dimensional representations while limiting connectivity costs (Litwin-Kumar et al., 2017). Analogous demands could narrow the range of viable architectures toward numerous compact, sparsely sampling neurons without specifying one shared transcriptional state.

This circuit-constraint explanation is testable but untested here. It predicts that matched developmental manipulations that reduce or redistribute afferent convergence should alter compact input-sampling morphology and sparse coding in the same direction in both lineages. If such manipulations change activity without coordinated morphological effects, or produce opposite structural responses in the two lineages, the proposed constraint-driven convergence would be weakened. The present transcriptomic data are compatible with this hypothesis but do not support it over other explanations.

### 4.4. Developmental interpretations remain hypotheses

The result can be compared with the restricted RU49/Zipro1 idea without treating that historical marker as supporting evidence. Yang and colleagues identified this factor, now annotated as ZSCAN21/Zscan21, in several developing granule-cell populations and proposed that anatomically distinct lineages might use elements of a common developmental mechanism (Yang et al., 1996). ZSCAN21 did not enter Tier 1 or Tier 2 because the selected-feature screen did not include it. In the broader full-matrix analysis, ZSCAN21 had a positive median rank delta of 0.5 in both branches: all three cerebellar datasets and three of four dentate datasets were positive, with GSE104323 as the exception. The combined evidence was nominal (Fisher's P=0.024) but did not survive Benjamini-Hochberg adjustment (q=0.254). ZSCAN21 was not externally tested and is retained only as historical context; this observation should not be interpreted as corroboration.

The externally tested GPM6A, NFIA, NFIB, RFX3, KCNK1 and Tier 2 results remain compatible with partial reuse of differentiation-associated features. They neither establish a shared specification mechanism nor reconstruct embryonic ancestry.

Altman's secondary-germinal-matrix framework supplies another useful interpretation. Both lineages include extensive precursor amplification and delayed postmitotic maturation, but they execute these processes in different regional structures and at different times. A shared dependence on proliferative amplification followed by compact-neuron differentiation could create partially similar molecular demands without erasing rhombic-lip versus telencephalic identity. The current data do not directly test whether transit amplification is the cause of candidate recurrence; lineage-resolved perturbations are required.

### 4.5. Relation to the 2005 conditioned-medium experiment

The 2005 culture study is historical motivation rather than a validation arm of the present analysis (Lu et al., 2005). It tested whether secreted cerebellar factors could change proliferation and differentiation of isolated hippocampal granule-lineage cells under a defined condition. The current study compares expression states across laboratories, ages, species and platforms, and therefore cannot establish which proteins were present at active concentrations in conditioned medium.

This distinction also explains why a simple expectation of higher cerebellar TGF-beta or BDNF transcript levels is not required by the 2005 result. A secreted-factor effect depends on producing cells, processing and release, extracellular stability, receptor competence and developmental stage. The decisive extension would therefore measure the conditioned-medium proteome and perturb candidate ligands and receptors in stage-matched dentate and cerebellar cultures. The transcriptomic candidate set can help prioritize readouts, but causal wording would be premature.

### 4.6. What the candidate genes may represent

The Tier 1 set spans several parts of neuronal differentiation. GPM6A has membrane and neurite-growth associations; NFIA and NFIB are developmental transcriptional regulators; RFX3 participates in transcriptional and ciliary programs; KCNK1 contributes to membrane excitability; and GABRA2 encodes a GABA-A receptor subunit. Their recurrence could therefore reflect a combination of neurite elaboration, state transition and excitability tuning. The Allen results caution against treating any of these genes as uniquely granule-cell-specific. GABRA2 did not reproduce in the local comparator test, while NFIA, NFIB and RFX3 were also prominent in other neuronal populations.

Tier 2 genes extend the set toward calcium handling, potassium conductance, glutamate or GABA signaling, guidance and vesicle-related functions. Their dual-positive enrichment is more informative than any single marker, but the cross-region transfer test did not identify a distinctive shared weighting of the set. The failure of the broad neurite/morphology module in the Allen analysis likewise argues against promoting this collection to a universal morphology program. Here, "widely used neuronal toolkits" means that most candidates are detectably or strongly expressed in non-target populations in the same common matrix. The matched-null analyses calibrate recurrence against expression- and detection-matched genes, but do not establish a uniquely defined universal toolkit.

### 4.7. Limitations

First, the discovery datasets differ in age, species, chemistry, annotation depth and available cell types. Within-dataset ranking limits, but does not eliminate, these differences. The human arm is asymmetric because there is only one primary human cerebellar resource and several human hippocampal resources. However, every Tier 1 and Tier 2 candidate retained positive mouse support in all four screen-branch combinations, so the reported set is not dependent on human-only direction.

Second, candidate selection, leave-one-dataset-out analysis and the discovery-layer matched-null analysis use related outputs. They assess stability and sensitivity, not validation. The Allen matrix is an external data source, but its directional hypothesis was selected in discovery; it tests adult generalization rather than an unselected candidate list.

Third, the Allen analysis is adult mouse data. It is well suited to cell-type specificity and common-platform comparison but cannot determine whether recurrence appears transiently during embryonic or postnatal differentiation. A sensitivity analysis using the adult immature-dentate subclass showed that broad downstream transfer differed sharply from the mature-dentate result, but only nine complete library blocks were available and the comparison represents adult cell state rather than developmental age. A developmental window therefore remains a testable explanation, not a rescue of the negative adult result.

Fourth, curated modules and the configuration score are investigator-defined summaries. Leave-one-gene-out stability and dataset-level inference reduce sensitivity to individual components, but these scores are not standardized biological scales. The five-module comparison is descriptive and underpowered.

Fifth, percentile bootstrap intervals are descriptive. Coverage can be imperfect for small library counts, especially the 12-library immature-dentate group, and those intervals should not be read as precise population bounds.

Finally, transcript abundance does not establish protein secretion, chromatin causality, morphology or circuit function. The study nominates experiments; it does not replace them.

### 4.8. Experimental tests of the organizational model

1. Perform matched developmental sampling of dentate and cerebellar granule lineages with single-cell RNA-seq and ATAC-seq, using the same chemistry and explicitly aligned proliferative, early postmitotic and mature states.

2. Perturb NFIA/NFIB/RFX3 and selected Tier 2 excitability genes singly and in combinations in lineage-specific cultures or conditional mouse models, then quantify proliferation exit, dendrite geometry, axon formation, synaptic physiology and candidate-gene recurrence.

3. Repeat the conditioned-medium experiment with quantitative secretome proteomics, ligand immunodepletion and receptor blockade across several recipient ages. This would directly test whether TGF-beta2, BDNF or additional secreted factors account for the cross-lineage effect.

4. Include Purkinje, CA1/CA3, olfactory-bulb interneuron and cortical compact-neuron controls in the same experiment. A true granule-associated differentiation feature should separate both targets from biologically relevant alternatives without requiring post hoc comparator selection.

5. Manipulate afferent convergence or input sparsity during matched developmental windows in each lineage, then measure dendritic sampling, firing sparsity and pattern separation. Concordant structural responses would support a circuit-constraint mechanism; absent or opposite responses would argue against it.

## 5. Conclusion

The analyses do not support classifying dentate and cerebellar granule cells as a shared adult molecular class. The two populations nevertheless show reproducible comparator-relative recurrence in a small gene set when each is evaluated against projection-neuron references from its own regional system. The absence of preferential broad direct similarity and null-exceeding cross-region transfer argues against a universal molecular configuration behind the shared histological name. The data instead identify a restricted recurring relation between regional identity and selected neuronal features. Matched developmental and circuit perturbations are required to determine whether that relation contributes to compact granule-cell organization.

## Acknowledgements

The author thanks the investigators who generated and shared the public datasets analyzed in this study.

## Use of artificial intelligence-assisted tools

OpenAI Codex (OpenAI) was used under the author's supervision to assist with analytical planning, code drafting, figure assembly, literature organization and language revision. The author reviewed the generated code and analytical outputs, verified numerical results and citations against source materials, revised the resulting text and takes full responsibility for the work. The tool was not used to generate or alter primary data.

## Funding

No specific funding was received for this work.

## Author contributions

J.L. conceived the study, curated the data, designed and performed the analyses, interpreted the results, prepared the figures and wrote the manuscript.

## Competing interests

The author declares no competing interests.

## Ethics statement

This study reanalyzed publicly available de-identified human and animal datasets and performed no new experiments involving human participants or animals.

## Data availability

All source data are publicly available from GEO under the accessions listed in Fig. 1c and in the repository dataset frame, and from the Allen Brain Cell Atlas (RRID:SCR_024440). The analysis uses the Allen Consensus-WMB-Macosko-10X expression matrices and associated taxonomy metadata. Machine-readable dataset metadata, candidate rankings, sensitivity analyses and Allen-derived summaries are available in the `dgd_reanalysis/Project/manuscript/source_tables` directory of the public code repository.

## Code availability

Analysis code and derived result tables required to reproduce the reported analyses are publicly available at https://github.com/Ghostneuron/granule-cell-convergence

## References

Altman, J. and Das, G. D. (1965). Autoradiographic and histological evidence of postnatal hippocampal neurogenesis in rats. J. Comp. Neurol. 124, 319-335.

Altman, J. and Das, G. D. (1966). Autoradiographic and histological studies of postnatal neurogenesis. I. A longitudinal investigation of the kinetics, migration and transformation of cells incorporating tritiated thymidine in neonate rats, with special reference to postnatal neurogenesis in some brain regions. J. Comp. Neurol. 126, 337-389.

Albus, J. S. (1971). A theory of cerebellar function. Math. Biosci. 10, 25-61.

Bird, A. D., Cuntz, H. and Jedlicka, P. (2024). Robust and consistent measures of pattern separation based on information theory and demonstrated in the dentate gyrus. PLoS Comput. Biol. 20, e1010706.

Chen, X., Zhong, X., Yue, W. et al. (2026). A transcription regulator atlas identifies TOX3 as an Atoh1 coactivator in cerebellar development and tumorigenesis. Proc. Natl. Acad. Sci. USA 123, e2527163123.

Disouky, A., Sanborn, M. A., Sabitha, K. R. et al. (2026). Human hippocampal neurogenesis in adulthood, ageing and Alzheimer's disease. Nature 652, 1264-1273.

Franjic, D., Skarica, M., Ma, S. et al. (2022). Transcriptomic taxonomy and neurogenic trajectories of adult human, macaque and pig hippocampal and entorhinal cells. Neuron 110, 452-469.e14.

Hochgerner, H., Zeisel, A., Lonnerberg, P. and Linnarsson, S. (2018). Conserved properties of dentate gyrus neurogenesis across postnatal development revealed by single-cell RNA sequencing. Nat. Neurosci. 21, 290-299.

Kuhn, H. G., Dickinson-Anson, H. and Gage, F. H. (1996). Neurogenesis in the dentate gyrus of the adult rat: age-related decrease of neuronal progenitor proliferation. J. Neurosci. 16, 2027-2033.

Leutgeb, J. K., Leutgeb, S., Moser, M. B. and Moser, E. I. (2007). Pattern separation in the dentate gyrus and CA3 of the hippocampus. Science 315, 961-966.

Lewis, P. M., Gritli-Linde, A., Smeyne, R., Kottmann, A. and McMahon, A. P. (2004). Sonic hedgehog signaling is required for expansion of granule neuron precursors and patterning of the mouse cerebellum. Dev. Biol. 270, 393-410.

Litwin-Kumar, A., Harris, K. D., Axel, R., Sompolinsky, H. and Abbott, L. F. (2017). Optimal degrees of synaptic connectivity. Neuron 93, 1153-1164.e7.

Lorente-Echeverria, B., Daaboul, D., Vandensteen, J. et al. (2025). A dynamic gene regulatory code drives synaptic development of hippocampal granule cells. Sci. Adv. 11, eadx5140.

Lu, J., Wu, Y., Sousa, N. and Almeida, O. F. X. (2005). SMAD pathway mediation of BDNF and TGF beta 2 regulation of proliferation and differentiation of hippocampal granule neurons. Development 132, 3231-3242.

Marr, D. (1969). A theory of cerebellar cortex. J. Physiol. 202, 437-470.

Parylak, S. L., Qiu, F., Linker, S. B. et al. (2023). Neuronal activity-related transcription is blunted in immature compared to mature dentate granule cells. Hippocampus 33, 412-423.

Peng, J., Sheng, A. L., Xiao, Q. et al. (2019). Single-cell transcriptomes reveal molecular specializations of neuronal cell types in the developing cerebellum. J. Mol. Cell Biol. 11, 636-648.

Ramnauth, A. D., Tippani, M., Divecha, H. R. et al. (2025). Spatiotemporal analysis of gene expression in the human dentate gyrus reveals age-associated changes in cellular maturation and neuroinflammation. Cell Rep. 44, 115300.

Ramon y Cajal, S. (1995). Histology of the Nervous System of Man and Vertebrates. New York: Oxford University Press. Original work published 1909-1911.

Wechsler-Reya, R. J. and Scott, M. P. (1999). Control of neuronal precursor proliferation in the cerebellum by Sonic Hedgehog. Neuron 22, 103-114.

Wizeman, J. W., Guo, Q., Wilion, E. M. and Li, J. Y. (2019). Specification of diverse cell types during early neurogenesis of the mouse cerebellum. eLife 8, e42388.

Yang, X. W., Zhong, R. and Heintz, N. (1996). Granule cell specification in the developing mouse brain as defined by expression of the zinc finger transcription factor RU49. Development 122, 555-566.

Yao, Z., van Velthoven, C. T. J., Kunst, M. et al. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature 624, 317-332.

Zhong, S., Wang, M., Huang, L. et al. (2023). Single-cell epigenomics and spatiotemporal transcriptomics reveal human cerebellar development. Nat. Commun. 14, 7613.

Zhou, Y., Su, Y., Li, S. et al. (2022). Molecular landscapes of human hippocampal immature neurons across lifespan. Nature 607, 527-533.

## Figure legends

### Figure 1. Comparative question, dataset frame and analytical sequence

**a,** Selected landmarks in the granule-cell question, from histological naming and Altman's postnatal microneuron studies to RU49/Zipro1, the 2005 conditioned-medium experiment and public single-cell data. **b,** Dentate and cerebellar granule cells arise from different regional lineages. The tested question is whether selected molecular features are reused relative to regional comparators; a recent shared migratory progenitor is not assumed. **c,** Composition of the ten-dataset discovery frame: four mouse dentate, two mouse cerebellar, three human dentate or hippocampal, and one human cerebellar dataset. **d,** Analytical sequence. Broad adult module convergence and developmental causality remained unproven in the external test.

### Figure 2. Candidate discovery and robustness analyses

**a,** Median dentate/cerebellar rank deltas for Tier 1 and Tier 2 candidates, separated by mouse selected-feature, mouse full-matrix and human selected-feature bridge screens. Values are the smaller of the dentate and cerebellar branch medians, so a positive value requires the same direction in both branches. **b,** Minimum median delta across all testable leave-one-dataset-out analyses. Green, Tier 1; purple, Tier 2. **c,** Observed mean minimum branch median delta compared with the median and 95% interval of selection-conditioned matched-null sets. The null analysis is used as a sensitivity test. **d,** Functional categories represented by Tier 1 and Tier 2 genes.

### Figure 3. Dataset-level recurrence and module-level analysis

**a,** One median configuration delta per eligible independent dataset. Blue, dentate; orange, cerebellum. Seven of seven medians were positive (one-sided exact sign test P=0.0078). **b,** Median dataset delta and bootstrap 95% confidence interval across all datasets and regional or species strata. **c,** Median convergence delta for five curated modules. Downstream modules were directionally greater than upstream or niche modules, but the exact two-versus-three module comparison was not conventionally significant (P=0.10).

### Figure 4. Adult Allen tests of local recurrence and direct similarity

**a,** Library coverage in the external Allen matrix for nine target or comparator populations after requiring at least 50 cells per library-population group. Cell counts are annotated as coverage and were not used as replicate counts. **b,** Mean library-level z scores for five predeclared modules. **c,** Branch-local expression contrasts for Tier 1 (green) and Tier 2 (purple) candidates. The x axis is cerebellar granule minus Purkinje mean Allen log2 expression. The y axis is mature dentate granule minus the equal-weighted CA1/ProS and CA3 mean. Genes in the upper-right quadrant are positive in both local contrasts. **d,** Downstream-minus-upstream Spearman similarity for the target pair and expanded comparator pairs. Orange marks the direct cerebellar-granule versus mature-dentate comparison. Its negative value indicates that downstream assembly genes were not more similar than upstream fate and niche genes in the adult common matrix.
