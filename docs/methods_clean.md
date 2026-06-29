## Methods

### Dataset discovery, curation and primary-core design

Public datasets were selected to test whether dentate and cerebellar granule-cell lineages preserve regional identity while converging on later neuronal assembly programs. Candidate resources were collected from GEO, Allen-related resources, literature-linked human dentate/hippocampal studies, cerebellar developmental studies, candidate scATAC/multiome or methylome resources, NeuroMorpho.Org and DANDI. The strict primary core contained ten datasets: four mouse dentate datasets, three cerebellar datasets and three human dentate/hippocampal datasets. Mouse dentate resources included GSE104323/GSE95752, GSE292261 and GSE214309. Cerebellar resources included GSE122357, GSE165657 and GSE312658. Human dentate/hippocampal resources included GSE186538, GSE325391 and GSE268609. Supporting datasets were retained for label construction, marker validation, selected-feature bridge analyses, regulatory-compatibility extension or perturbation interpretation, and were kept analytically distinct from the strict primary core. The primary-core frame and supporting dataset tiers are provided in Tables S1 and S2.

### Human bridge-object construction and label harmonization

Human dentate/hippocampal bridge datasets were processed as selected-gene or reduced sparse objects when full transcriptome matrices were not locally available or were not needed for a specific validation layer. Source annotations, marker programs and reduced human-core label projections were used to harmonize candidate granule-cell, immature-neuron, neurogenic, pyramidal/comparator and background labels. Human scaffold datasets GSE185277 and GSE185553 were used for marker tuning and label construction from the human hippocampal immature-neuron lifespan resource, but they were not counted in the strict ten-dataset primary core.

### Candidate-cell annotation and pseudobulk construction

Candidate granule-cell, precursor, comparator and local background labels were assigned using source annotations, marker support and branch-specific label rules. Candidate and background groups were converted to pseudobulk summaries to stabilize cross-dataset comparisons and reduce sensitivity to differences in cell yield. Candidate-gene, expanded-gene, genome-wide symbol and full-MGI ortholog pseudobulk layers were generated separately. Per-cell module scores were used for selected validation checks, but primary cross-dataset tests used within-dataset or within-sample rank summaries.

### Ortholog-aware rank-meta modeling and candidate tiering

The primary transcriptomic screen used an ortholog-aware rank-meta design. For each dataset or contrast, candidate and local-background pseudobulk expression was converted to within-layer relative ranks, and candidate-background rank deltas were computed within the same dataset layer. Human and mouse genes were harmonized through one-to-one MGI orthologs where cross-species comparison was required. This design emphasizes reproducible within-dataset ordering rather than absolute expression magnitude, which is sensitive to species, chemistry, processing pipeline and sequencing depth.

Two evidence routes were tracked. The selected-feature route used curated candidate-gene, expanded-gene, module-gene and human bridge-feature matrices; this route maximized coverage for predefined granule-cell and mechanism features but was not treated as genome-wide discovery. The full-MGI/full-matrix route used broader pseudobulk matrices remapped to one-to-one MGI orthologs where full matrices were available. For each gene, screen and branch, nominal branch support required at least two datasets, median dataset rank delta \(\widetilde{\Delta}_{\mathrm{rank}} > 0\), positive dataset fraction \(f_{\Delta > 0} \ge 0.75\) and best one-sided dataset-level \(p_{\min} \le 0.25\) from the t, sign or Wilcoxon tests. Shared support was assigned without pooling raw expression across datasets: selected-feature shared support required nominal dentate and cerebellar support in the selected-feature route, full-matrix shared support required nominal dentate and cerebellar support in the full-MGI/full-matrix route, and both-screen shared support required both route-level calls.

Table S3 reports the formal-rank fields used for candidate ordering. Four possible screen-branch tests were tracked for each gene: selected-feature dentate, selected-feature cerebellar, full-matrix dentate and full-matrix cerebellar. The field `formal_n_available_branches` is the number of these screen-branches with sufficient data, and `formal_n_nominal_branches` is the number satisfying the nominal support rule above. FDR10 support required replication support plus the best within-screen/branch BH-adjusted one-sided q value \(\le 0.10\). The field `formal_fisher_q_bh_all_available_branches` is the BH-adjusted Fisher combined p value computed from the available one-sided branch-level \(p_{\min}\) values; it was used as a summary support statistic, not as a sole inclusion rule. Median-delta fields in Table S3 are median within-dataset candidate-minus-background rank deltas for the indicated screen and branch, so positive values indicate higher within-dataset rank in candidate granule populations than in local background populations. The formal-rank priority score was a deterministic ranking index:

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

Here, \(\Delta_{\mathrm{rank},u}\) is the candidate-minus-background rank delta for dataset/sample unit \(u\), \(d[u]\) identifies its dataset, \(b_{d[u]}\) is the dataset-specific random intercept, \(\epsilon_u\) is residual error, and \(\beta_0\) estimates the average directional rank shift for that gene in that screen-branch. The random intercept allows all units from the same dataset to share a baseline shift. The one-sided test was \(H_0:\beta_0 \le 0\) versus \(H_1:\beta_0 > 0\). Models were fit only when at least three units from at least two datasets were available. A linear mixed model with a dataset random intercept was used when possible; if the mixed model failed, an intercept-only ordinary least-squares model with dataset-cluster-robust standard errors was used as a fallback. `formal_model_n_nominal_branches` counts screen-branches with branch replication support and \(p_{\mathrm{greater}} \le 0.25\) in this unit-level model; the corresponding FDR10 count used branch replication support and Benjamini-Hochberg-adjusted \(q \le 0.10\).

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

Here, \(Y_i\) is the granule-design score, \(F_i\) is branch-matched fate polarity, \(C_i\) is construction balance, \(F_iC_i\) is the fate-construction coupling term, \(T_i\) is stage or pseudotime, \(N_i\) is niche signal, \(E_i\) is regulatory compatibility, \(M_i\) is morphology sparse sampling, \(A_i\) is activity sparsity, \(R_i\) is circuit-resource constraint, \(\eta_i\) collects additional observed covariates and random effects, and \(u_{d[i]}\), \(u_{s[i]}\) and \(u_{a[i]}\) are dataset, species and assay random effects. Because the present datasets are not matched across all modalities, the implemented model was a weighted hierarchical evidence synthesis. An evidence unit was one scored observation from a sample-level contrast, branch summary, epigenomic sensitivity contrast, morphology or activity calibration, or simulation-calibration result; `unit_id` labels these scored observations, not genes. Each evidence unit was normalized to a bounded score in \([-1, 1]\), assigned a quality weight according to measurement level, and summarized first by hierarchy level and term:

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

External morphology validation used NeuroMorpho records for dentate and cerebellar granule-cell reconstructions with dendritic morphometry, building on public digital neuronal reconstruction resources and morphometric mining methods (Ascoli et al., 2007; Halavi et al., 2012; Polavaram et al., 2014). Records were filtered for granule-cell identity, dentate or cerebellar anatomical assignment, non-missing dendritic stem/branch/length measurements and acceptable reconstruction integrity. Morphological comparisons focused on input-sampling architecture rather than assuming identical dendritic arbor geometry between regions. For hierarchical evidence synthesis, morphology comparisons were converted to bounded `MorphologySparseSampling` scores. Let \(r_m\) be the dentate-to-cerebellar median ratio for metric \(m\), and let \(\delta_m\) be the Cliff's delta for the dentate-versus-cerebellar comparison:

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

Physiology validation used local DANDI 000003 NWB files as a dentate granule-cell activity layer. DANDI 000003 is the public archive version of the Senzai and Buzsaki hippocampal granule-cell and mossy-cell physiology dataset (Senzai and Buzsaki, 2017). NWB files were parsed for units, cell-type labels, behavior/position variables and awake-moving samples. Labeled granule units were analyzed for firing rate, active spatial-bin fraction, spatial information, spatial sparsity and population-vector structure. For hierarchical evidence synthesis, DANDI metrics were converted to bounded `ActivitySparsity` scores using \(\operatorname{clip}_{[-1,1]}(\cdot)\):

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
