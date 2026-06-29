# External Dataset Inventory for Granule Cell Convergence Project

Searched and cached on 2026-06-21.

## Purpose

This inventory complements the literature and local datasets already in the workspace for the question:

> Why do cerebellar granule cells and dentate gyrus/hippocampal granule cells share similar morphology despite arising in different brain regions?

The useful in silico project should not simply compare "granule cell" names. It should test whether these cells converge on shared morphogenetic programs: transcription factors, axon/dendrite/cytoskeleton modules, guidance/synapse genes, spatial neighborhood programs, and in the dentate gyrus, protein-level evidence.

## Search Sources

### NCBI GEO DataSets

Exact user-specified searches:

- Cerebellum, granule cell, single-cell sequencing: https://www.ncbi.nlm.nih.gov/gds?term=((cerebellum)%20AND%20granule%20cell)%20AND%20single%20cell%20sequencing
- Hippocampus, granule cell, single-cell sequencing: https://www.ncbi.nlm.nih.gov/gds/?term=((hippocampus)+AND+granule+cell)+AND+single+cell+sequencing
- GSE242688 proteomics/spatial transcriptomics: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242688

Expanded GEO searches were also run for:

- `dentate gyrus granule cell single cell`
- `cerebellar granule cell single cell`
- `cerebellum granule neuron single cell`
- `dentate granule single nucleus`
- rat-specific cerebellar/dentate granule single-cell searches

Local NCBI search cache:

- `Project/dataset_search_cache/ncbi_gds_granule_cell_merged_summary.tsv`
- `Project/dataset_search_cache/ncbi_gds_cerebellum_granule_singlecell_summary.tsv`
- `Project/dataset_search_cache/ncbi_gds_hippocampus_granule_singlecell_summary.tsv`
- `Project/dataset_search_cache/ncbi_gds_dentate_gyrus_granule_singlecell_summary.tsv`
- `Project/dataset_search_cache/ncbi_gds_rat_granule_loose_summary.tsv`

The merged table contains 123 unique GEO records plus a header. The strict rat organism searches returned no direct rat granule-cell single-cell GEO series; loose rat text searches mostly recovered rat platform records, not useful biological series. Rat should be treated as optional later support, not a core dataset axis.

### Allen Institute

Sources:

- Allen Brain Cell Atlas portal: https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas
- ABC Atlas WMB taxonomy documentation: https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-taxonomy.html
- Allen WMB taxonomy metadata bucket path referenced by docs: `metadata/WMB-taxonomy/20231215`

Local Allen files:

- `External_Data/Allen_Institute/WMB_taxonomy_cluster_annotation_CCN202307220.xlsx`
- `External_Data/Allen_Institute/WMB-taxonomy_access_page.html`
- `External_Data/Allen_Institute/WMB_granule_relevant_subclass_annotation.tsv`
- `External_Data/Allen_Institute/WMB_granule_relevant_supertype_annotation.tsv`
- `External_Data/Allen_Institute/WMB_granule_relevant_cluster_annotation.tsv`

The Allen taxonomy gives exact extraction keys:

| Cell set | Allen label | Useful markers | Extraction IDs |
|---|---|---|---|
| Dentate mature granule | `037 DG Glut` | `Prox1, Itpka, C1ql3`; TFs `Prox1, Glis3, Egr3` | supertype `0136-0139`; clusters `0502-0510` |
| Dentate immature/DG-PIR excitatory | `038 DG-PIR Ex IMN` | `Mex3a, Neurod1`; TFs `Prox1, Tbr1, Sox4, Pou3f3` | supertype `0140-0142`; clusters `0511-0517` |
| Cerebellar granule | `314 CB Granule Glut` | `Gabra6, Ror1`; TFs `Pax6, Neurod2, Etv1` | supertype `1154-1155`; clusters `5197-5201` |

These IDs are important because later we can pull only WMB cells/spots belonging to `DG Glut`, `DG-PIR Ex IMN`, and `CB Granule Glut`, instead of downloading the full atlas.

## Local Datasets Already Present Before This Search

### Cerebellum

| Dataset | Species | Local status | Project role |
|---|---:|---|---|
| GSE165657, molecular/spatial design of human cerebellar development | Human | `Single cell sequencing data/Cerebellum/GSE165657/` | Primary human developmental cerebellar anchor. |
| GSE122357, developing cerebellum single-cell transcriptomes | Mouse | `Single cell sequencing data/Cerebellum/GSE122357/` | Primary mouse cerebellar developmental anchor. |
| GSE246785, lifelong restructuring of cerebellar 3D genome | Mouse | `Single cell sequencing data/Cerebellum/GSE246785/` | Regulatory/3D genome support for granule-cell maturation. |

### Dentate Gyrus / Hippocampus

| Dataset | Species | Local status | Project role |
|---|---:|---|---|
| GSE104323, dentate gyrus neurogenesis across postnatal development | Mouse | `Single cell sequencing data/Dentate gyrus/GSE104323/` | Primary DG developmental trajectory anchor. |
| GSE95752, C1 single-cell dentate gyrus transcriptomes | Mouse | `Single cell sequencing data/Dentate gyrus/GSE95752/` | Single-cell DG support set. |
| GSE95315, 10x dentate gyrus transcriptomes | Mouse | GEO hit; related to the GSE104323/GSE95752 family | Useful if we need the original 10x GEO representation. |

## New Data Downloaded

### Dentate Gyrus / Hippocampus

| Dataset | Species | Local path | What was downloaded | Why it matters |
|---|---:|---|---|---|
| GSE292261, dynamic regulatory code in hippocampal granule cells, SS2 | Mouse | `External_Data/GEO/GSE292261/` | Filtered/unfiltered counts, sample metadata, filtered loom, GEO SOFT metadata | Very high priority: directly about hippocampal granule-cell synaptic development and gene regulation. |
| GSE214905, commissural and non-commissural dentate granule cells | Mouse | `External_Data/GEO/GSE214905/` | QC counts, no-QC counts, normalized counts, series matrix, listing | Tests whether DG granule-cell projection differences reuse morphology/guidance modules. |
| GSE214309, activity-related transcription in immature vs mature DGCs | Mouse | `External_Data/GEO/GSE214309/` | Counts, TPM, series matrix, listing | Excellent maturity/activity contrast for DGC morphology and plasticity modules. |
| GSE242688, MALDI/LC-MS/MS spatial proteomics plus spatial transcriptomics in hippocampus | Mouse | `External_Data/Proteomics/GSE242688/` | GEO RAW tar with Visium matrix files, series matrix | Protein/spatial support for DG granule-cell subpopulations. GEO summary reports two proteomically distinct DG granule-cell subpopulations. |

### Cerebellum

| Dataset | Species | Local path | What was downloaded | Why it matters |
|---|---:|---|---|---|
| GSE312658, TOX3/Atoh1 cerebellar development/tumorigenesis atlas | Mouse | `External_Data/GEO/GSE312658/` | Ctrl and cKO 10x barcodes/features/matrices, GEO SOFT metadata, filelist/listing | High priority regulatory perturbation dataset for cerebellar granule lineage. |
| GSE150153, human iPSC cerebellar organoids | Human | `External_Data/GEO/GSE150153/` | Processed RDS files, barcodes, genes, matrices, metadata for both samples, GEO SOFT metadata, filelist/listing | Useful human cerebellar development/organoid comparator, secondary to GSE165657. |

### Allen Institute

| Dataset | Species | Local path | What was downloaded | Why it matters |
|---|---:|---|---|---|
| Allen WMB taxonomy CCN202307220 | Mouse | `External_Data/Allen_Institute/` | Full taxonomy workbook plus DG/CB granule subset TSVs | Gives authoritative cell-set IDs, markers, and clusters for DG and CB granule cells across whole mouse brain. |

## Manifest-Only / Not Downloaded Yet

These were identified as useful but intentionally left as metadata or file listings because of size:

| Dataset | Species | Local manifest path | Size issue | Recommendation |
|---|---:|---|---|---|
| GSE292260, hippocampal granule-cell snMO multiome | Mouse | `External_Data/GEO/GSE292260/` | One RDS file is about 2.5 GB | Download when ready for chromatin/GEX multiome integration. High priority, but not needed for first transcriptomic pass. |
| GSE312741, Tox3 chromatin occupancy in cerebellum | Mouse | `External_Data/GEO/GSE312741/` | RAW tar about 274 MB; bigWigs about 360-369 MB each | Download if regulatory occupancy becomes a core aim. |
| GSE310490, PTEN variant cerebellar neuronal differentiation | Human | `External_Data/GEO/GSE310490/` | RAW tar about 20 GB; individual RDS files include multi-GB files | Useful human disease/organoid angle, but too large for broad exploratory download. |
| Full Allen WMB expression atlas | Mouse | `External_Data/Allen_Institute/` contains taxonomy only | Full expression atlas is much larger than the taxonomy | Next step should pull only `DG Glut`, `DG-PIR Ex IMN`, and `CB Granule Glut` subsets using Allen access tools. |

## Additional Candidate Resources From User-Suggested Papers

These are useful candidate resources but have not been downloaded or included in the current analysis matrix:

| Resource group | Component dataset/resource | Species | Status | Recommendation |
|---|---|---:|---|---|
| Human DG/hippocampal core construction pack | `GSE185277`, `GSE185553`, `GSE186538`, `GSE325391`, `GSE268609` | Human | `GSE185277` and `GSE185553` raw archives downloaded and sparse objects built; `GSE186538` human DG GC sparse subset built; remaining datasets not downloaded | Build/curate this branch first to correct the previous lack of a human dentate/hippocampal comparator for the human cerebellum core dataset. |
| Nature Neuroscience 2025 mouse DG aging/neurogenic-lineage resource | `GSE233363` | Mouse | Source listing saved; data not downloaded | Add after the human dentate branch as mouse DG aging, neurogenic-lineage maturation, and spatial niche-inflammation validation. |
| Cell Reports 2025 human DG spatial lifespan resource | `spatial_DG_lifespan`, Zenodo `10.5281/zenodo.10126687` | Human | Article PDF present in `Literature`; data not downloaded | Use as human spatial/lifespan context; the Zenodo data archive is large, so acquire only when spatial validation is needed. |
| Nature 2022 human hippocampal immature-neuron disease context | `GSE198323` | Human | Source listing saved; data not downloaded | Keep as AD/disease context after healthy/reference human imGC construction. |

## Highest-Value Project Data Matrix

Use this as the practical analysis stack:

| Layer | Dentate/Hippocampus | Cerebellum | Purpose |
|---|---|---|---|
| Developmental scRNA | GSE104323, GSE95752/GSE95315 | GSE122357, GSE165657, GSE150153 | Compare maturation trajectories and shared morphogenetic modules. |
| Human dentate/hippocampal core extension | GSE185277, GSE185553, GSE186538, GSE325391, GSE268609 | GSE165657 | Correct the current human-side asymmetry before cross-species conclusions. |
| Modern whole-brain taxonomy | Allen `DG Glut`, `DG-PIR Ex IMN` | Allen `CB Granule Glut` | Independent adult mouse validation and marker extraction. |
| Regulatory perturbation/multiome | GSE292261 now; GSE292260 later | GSE246785, GSE312658; GSE312741 later | Test whether shared morphology maps onto shared regulatory logic. |
| Projection/activity states | GSE214905, GSE214309 | Not direct; use cerebellar developmental/maturation contrasts | Ask whether DG granule-cell morphology modules vary with projection/activity. |
| Proteomics/spatial support | GSE242688 | None equivalent yet | Protein/transcript concordance in DG; useful as an orthogonal validation layer. |

## Recommendation

The project is reasonable and now more data-supported than before. The strongest framing is not "same morphology means same cell type," but "distant excitatory microgranular neurons may converge on a partially shared morphogenetic program under different regional patterning constraints."

Recommended first paper-shaped analysis:

1. Define a granule-cell morphology module.
   - Start with curated genes for dendrite/axon growth, cytoskeleton, synaptogenesis, axon guidance, ion channels, and immediate early/plasticity.
   - Add unbiased modules from DEG/WGCNA or NMF across DG and CB datasets.

2. Build matched developmental trajectories.
   - Cerebellum: GSE122357/GSE165657/GSE312658.
   - Dentate: GSE104323/GSE292261/GSE214309.
   - Compare pseudotime modules, not only terminal markers.

3. Test convergence vs regional identity.
   - Shared morphogenesis: `Neurod`, cytoskeletal, guidance, synaptic-vesicle, mitochondrial/redox, microtubule genes.
   - Region-specific identity: CB `Atoh1/Pax6/Zic/Barhl/Gabra6`; DG `Prox1/Tbr/Neurod/Mex3a/Glis3`.
   - The most interesting result would be partial convergence: common structural modules embedded in different upstream GRNs.

4. Validate externally.
   - Use Allen WMB cell-set IDs to validate adult mouse marker/module expression in `DG Glut` vs `CB Granule Glut`.
   - Use GSE242688 to ask whether DG transcript modules have protein-level support.
   - Add morphology reconstructions later from NeuroMorpho.Org or published reconstructions to connect gene modules to actual shape metrics.

5. Keep scope publication-friendly.
   - A strong review/meta-analysis is feasible now.
   - A computational paper is feasible if the analysis adds a real quantitative convergence score, module conservation analysis, and external validation with Allen WMB plus GSE242688.

## Immediate Next Computational Steps

1. Build a sample/dataset manifest with species, age/stage, region, platform, and local path.
2. Convert all matrix formats into a common AnnData/Seurat-compatible format.
3. Annotate dentate and cerebellar granule cells using local markers and Allen labels.
4. Compute cross-dataset ortholog-mapped gene modules for mouse-human comparison.
5. Rank candidate convergence genes by:
   - shared DG/CB expression,
   - developmental upregulation,
   - morphology/guidance/synapse ontology,
   - regulatory evidence from GSE246785/GSE312658/GSE292260,
   - protein support in GSE242688 where available.
