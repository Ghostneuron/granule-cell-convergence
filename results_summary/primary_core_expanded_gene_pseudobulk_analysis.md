# Expanded Primary-Core Gene Pseudobulk Screen

Date built: 2026-06-22

## Scope

This analysis expands the 67-gene candidate pseudobulk screen to the 2,169-gene selected human-core universe. It is not a final whole-transcriptome DE model, but it is the first broad discovery layer that keeps all 10 primary datasets in the same target-gene frame.

Minimum broad-class size for rank statistics: 20 cells/nuclei.

## Coverage

- Selected-gene universe: 2,169 genes.
- Pseudobulk expression rows: 537,070.
- Primary datasets represented: 10/10.
- Genes tested in statistics: 2,169.
- Shared-positive rank genes: 1,279.
- Shared-positive genes passing BH<0.10 in both branches: 579.
- Shared-positive genes passing BH<0.20 in both branches: 787.
- Original 67-gene packet genes recovered among shared-positive hits: 42.

- `GSE104323` / `10X_all_cells`: 1937/2169 selected genes, 24185/24185 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464549_P0`: 1919/2169 selected genes, 7468/7468 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464550_P8a`: 1925/2169 selected genes, 6088/6088 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464551_P8b`: 1926/2169 selected genes, 8168/8168 labeled observations (`full_raw_matrix`).
- `GSE165657` / `Cerebellum_aggr`: 2031/2169 selected genes, 73550/73550 labeled observations (`full_raw_matrix`).
- `GSE186538` / `selected_bridge_all_samples`: 2169/2169 selected genes, 32004/32004 labeled observations (`selected_norm_bridge`).
- `GSE214309` / `snRNA_counts`: 1800/2169 selected genes, 401/401 labeled observations (`full_raw_matrix`).
- `GSE268609` / `selected_bridge_all_samples`: 2169/2169 selected genes, 366175/366175 labeled observations (`selected_gene_bridge`).
- `GSE292261` / `SS2_filtered_counts`: 1838/2169 selected genes, 509/509 labeled observations (`full_raw_matrix`).
- `GSE312658` / `Ctrl`: 1949/2169 selected genes, 4857/4857 labeled observations (`full_raw_matrix`).
- `GSE312658` / `cKO`: 1949/2169 selected genes, 6946/6946 labeled observations (`full_raw_matrix`).
- `GSE325391` / `selected_bridge_all_samples`: 2169/2169 selected genes, 59075/59075 labeled observations (`selected_gene_bridge`).
- `GSE95752` / `C1_all_cells`: 1889/2169 selected genes, 2303/2303 labeled observations (`full_raw_matrix`).

## Original Candidate Recovery

- `Map2` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Syt1` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Stmn2` (shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500.
- `Gap43` (shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500.
- `Ndufa4` (supporting_metabolic_validation): dentate delta 0.500, cerebellar delta 0.500.
- `Robo2` (shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500.
- `Cfl1` (shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500.
- `Ndufb8` (supporting_metabolic_validation): dentate delta 0.500, cerebellar delta 0.500.
- `Stmn3` (shared_structural_executor): dentate delta 0.500, cerebellar delta 0.500.
- `Neurod2` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Tubb3` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Prdx2` (supporting_metabolic_validation): dentate delta 0.500, cerebellar delta 0.500.
- `Ror1` (regional_cerebellar_identity): dentate delta 0.500, cerebellar delta 0.500.
- `Neurod1` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Slc17a7` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Tbr1` (shared_granule_neuronal_state): dentate delta 0.500, cerebellar delta 0.500.
- `Gabra6` (regional_cerebellar_identity): dentate delta 0.500, cerebellar delta 0.500.
- `St18` (shared_granule_neuronal_state): dentate delta 0.333, cerebellar delta 0.500.
- `Etv1` (regional_cerebellar_identity): dentate delta 0.333, cerebellar delta 0.500.
- `Zic1` (regional_cerebellar_identity): dentate delta 0.333, cerebellar delta 0.500.

## New Shared-Positive Candidates

- `CELF2` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GPM6A` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `TCF4` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `MAP1B` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `TTC3` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIA` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP3CA` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CALM1` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIB` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `RTN3` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `RTN1` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NRXN1` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `APP` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `BASP1` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ANK3` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CALM2` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `LUC7L3` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `DDX5` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `FRMD4A` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CAMTA1` (high_information_gene): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.

## Interpretation

- This is a discovery screen over the selected-gene bridge universe, not a final genome-wide model.
- Genes that are shared-positive here should be divided into broad neuronal state, structural/morphogenesis executors, metabolic/supporting genes, and regional identity leakage before entering a manuscript mechanism figure.
- The most important use is prioritization: it identifies new candidate genes to inspect in the future whole-transcriptome ortholog-aware DE model and tests whether the 67-gene packet is recovered in a broader feature space.
- Branch-specific hits are useful too: they help separate shared morphology executors from region-specific wiring or maturation programs.

## Outputs

- Expression table: `Project/results/primary_core_expanded_gene_pseudobulk_expression.tsv.gz`
- Coverage table: `Project/results/primary_core_expanded_gene_pseudobulk_coverage.tsv`
- Statistics table: `Project/results/primary_core_expanded_gene_pseudobulk_statistics.tsv`
- Shared hits: `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.tsv`
- Branch-specific hits: `Project/results/primary_core_expanded_gene_pseudobulk_branch_specific.tsv`
- Shared-hit plot: `Project/results/primary_core_expanded_gene_pseudobulk_shared_hits.png`
