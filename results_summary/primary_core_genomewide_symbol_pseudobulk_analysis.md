# Genome-Wide Same-Symbol Full-Matrix Pseudobulk Screen

Date built: 2026-06-22

## Scope

This analysis streams local full matrices for primary-core datasets and aggregates a full human-symbol target universe from the GSE186538 human DG taxonomy gene list. Mouse genes are mapped by same-root upper-case symbol, so this is an ortholog-ready screen rather than a curated ortholog model.

Contrast statistics require at least two eligible broad classes within a dataset/sample/gene. Single-class human DG anchor samples are retained for expression/detection but excluded from candidate-versus-background rank tests.

## Coverage

- Full-symbol target universe: 33,939 genes.
- Pseudobulk expression rows: 628,339.
- Datasets with full-matrix expression represented: 8/10 primary datasets.
- Datasets contributing to rank contrasts: 7 (GSE104323, GSE122357, GSE165657, GSE214309, GSE292261, GSE312658, GSE95752).
- Genes tested in contrast statistics: 21,253.
- Shared-positive rank genes: 6,440.
- Shared-positive genes passing BH<0.10 in both branches: 283.
- Branch-specific genes: 5,286.

- `GSE104323` / `10X_all_cells`: 15596/33939 target symbols, 24185/24185 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464549_P0`: 13327/33939 target symbols, 7468/7468 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464550_P8a`: 13551/33939 target symbols, 6088/6088 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464551_P8b`: 13539/33939 target symbols, 8168/8168 labeled observations (`full_raw_matrix`).
- `GSE165657` / `Cerebellum_aggr`: 20308/33939 target symbols, 73550/73550 labeled observations (`full_raw_matrix`).
- `GSE186538` / `selected_bridge_all_samples`: 33939/33939 target symbols, 32067/32067 labeled observations (`full_sparse_subset`).
- `GSE214309` / `snRNA_counts`: 14120/33939 target symbols, 401/401 labeled observations (`full_raw_matrix`).
- `GSE292261` / `SS2_filtered_counts`: 10886/33939 target symbols, 509/509 labeled observations (`full_raw_matrix`).
- `GSE312658` / `Ctrl`: 16239/33939 target symbols, 4857/4857 labeled observations (`full_raw_matrix`).
- `GSE312658` / `cKO`: 16239/33939 target symbols, 6946/6946 labeled observations (`full_raw_matrix`).
- `GSE95752` / `C1_all_cells`: 12962/33939 target symbols, 2303/2303 labeled observations (`full_raw_matrix`).

## Top Shared Hits

- `CELF2`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GPM6A`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIB`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIA`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP3CA`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `FTH1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ANK3`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CAMTA1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `YWHAH`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NBEA`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GRIA2`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP2R2C`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NRN1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CELF1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ATP2B1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CD47`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `BEX2`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CELF4`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `MAPK1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `SYT1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `SDCBP`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GATAD2B`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPM1H`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `SPOP`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `STXBP1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `MRFAP1`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NEUROD2`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CDC42SE2`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP1CB`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ASH1L`: dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.

## Interpretation

- This is the first full-matrix discovery pass, but it still uses same-symbol mapping rather than a curated mouse-human ortholog table.
- The result should be used to prioritize candidates and to design the final mixed-effect/ortholog DE model, not as the final gene-level claim.
- The current full-matrix contrast layer excludes `GSE325391` and `GSE268609` from genome-wide tests because they are represented locally by selected-gene bridge objects or very large source objects requiring dedicated export.
- `GSE186538` contributes full human DG expression/detection but not rank contrast because the extracted local object is a DG GC subset without local non-DG background.

## Outputs

- Expression table: `Project/results/primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz`
- Coverage table: `Project/results/primary_core_genomewide_symbol_pseudobulk_coverage.tsv`
- Statistics table: `Project/results/primary_core_genomewide_symbol_pseudobulk_statistics.tsv`
- Shared hits: `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.tsv`
- Branch-specific hits: `Project/results/primary_core_genomewide_symbol_pseudobulk_branch_specific.tsv`
- Shared-hit plot: `Project/results/primary_core_genomewide_symbol_pseudobulk_shared_hits.png`
