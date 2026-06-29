# MGI One-to-One Ortholog Full-Matrix Pseudobulk Screen

Date built: 2026-06-22

## Scope

This screen expands the previous full-matrix same-symbol analysis by resolving mouse matrix rows through the MGI mouse-symbol side of one-to-one human-mouse homology classes. The canonical output gene is the human symbol.

It remains a rank-based pseudobulk screen, not the final mixed-effect DE model. It is designed to make the current gene-level evidence ortholog-aware before formal modeling.

## Coverage

- MGI one-to-one target genes: 17,611.
- Non-identical human/mouse symbol targets: 1,366.
- Pseudobulk expression rows: 508,917.
- Datasets with expression represented: 8/10 primary datasets.
- Datasets contributing to rank contrasts: 7 (GSE104323, GSE122357, GSE165657, GSE214309, GSE292261, GSE312658, GSE95752).
- Genes tested in contrast statistics: 16,704.
- Shared-positive ortholog genes: 6,413.
- Shared-positive same-symbol genes: 6,167.
- Shared-positive non-identical-symbol genes: 246.
- Shared-positive genes passing BH<0.10 in both branches: 286.
- Branch-specific genes: 5,191.

- `GSE104323` / `10X_all_cells`: 15624/17611 MGI targets, 905/1366 non-identical targets, 24185/24185 labeled observations (`mouse`).
- `GSE122357` / `GSM3464549_P0`: 13320/17611 MGI targets, 608/1366 non-identical targets, 7468/7468 labeled observations (`mouse`).
- `GSE122357` / `GSM3464550_P8a`: 13556/17611 MGI targets, 628/1366 non-identical targets, 6088/6088 labeled observations (`mouse`).
- `GSE122357` / `GSM3464551_P8b`: 13541/17611 MGI targets, 624/1366 non-identical targets, 8168/8168 labeled observations (`mouse`).
- `GSE165657` / `Cerebellum_aggr`: 15558/17611 MGI targets, 909/1366 non-identical targets, 73550/73550 labeled observations (`human`).
- `GSE186538` / `selected_bridge_all_samples`: 16407/17611 MGI targets, 990/1366 non-identical targets, 32067/32067 labeled observations (`human`).
- `GSE214309` / `snRNA_counts`: 14059/17611 MGI targets, 721/1366 non-identical targets, 401/401 labeled observations (`mouse`).
- `GSE292261` / `SS2_filtered_counts`: 10925/17611 MGI targets, 473/1366 non-identical targets, 509/509 labeled observations (`mouse`).
- `GSE312658` / `Ctrl`: 16243/17611 MGI targets, 1001/1366 non-identical targets, 4857/4857 labeled observations (`mouse`).
- `GSE312658` / `cKO`: 16243/17611 MGI targets, 1001/1366 non-identical targets, 6946/6946 labeled observations (`mouse`).
- `GSE95752` / `C1_all_cells`: 12964/17611 MGI targets, 560/1366 non-identical targets, 2303/2303 labeled observations (`mouse`).

## Top Shared Hits

- `CELF2` (same_symbol; mouse `Celf2`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GPM6A` (same_symbol; mouse `Gpm6a`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIB` (same_symbol; mouse `Nfib`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NFIA` (same_symbol; mouse `Nfia`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP3CA` (same_symbol; mouse `Ppp3ca`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `FTH1` (same_symbol; mouse `Fth1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ANK3` (same_symbol; mouse `Ank3`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CAMTA1` (same_symbol; mouse `Camta1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `YWHAH` (same_symbol; mouse `Ywhah`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NBEA` (same_symbol; mouse `Nbea`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GRIA2` (same_symbol; mouse `Gria2`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP2R2C` (same_symbol; mouse `Ppp2r2c`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NRN1` (same_symbol; mouse `Nrn1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CELF1` (same_symbol; mouse `Celf1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ATP2B1` (same_symbol; mouse `Atp2b1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CD47` (same_symbol; mouse `Cd47`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CELF4` (same_symbol; mouse `Celf4`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `MAPK1` (same_symbol; mouse `Mapk1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `SYT1` (same_symbol; mouse `Syt1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `SDCBP` (same_symbol; mouse `Sdcbp`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `GATAD2B` (same_symbol; mouse `Gatad2b`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPM1H` (same_symbol; mouse `Ppm1h`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `STXBP1` (same_symbol; mouse `Stxbp1`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `TMEM178A` (nonidentical_symbol; mouse `Tmem178`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `NEUROD2` (same_symbol; mouse `Neurod2`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `CDC42SE2` (same_symbol; mouse `Cdc42se2`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `PPP1CB` (same_symbol; mouse `Ppp1cb`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ASH1L` (same_symbol; mouse `Ash1l`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `RUNDC3A` (same_symbol; mouse `Rundc3a`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.
- `ST18` (same_symbol; mouse `St18`): dentate delta 0.500, cerebellar delta 0.500, BH<0.10 both branches=True.

## Non-Identical Symbol Shared Hits

- `TMEM178A` / mouse `Tmem178`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF148` / mouse `Zfp148`: dentate delta 0.500, cerebellar delta 0.500.
- `C1orf21` / mouse `1700025G04Rik`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF536` / mouse `Zfp536`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF226` / mouse `Zfp61`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF263` / mouse `Zfp263`: dentate delta 0.500, cerebellar delta 0.500.
- `C5orf22` / mouse `6030458C11Rik`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF667` / mouse `Zfp667`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF12` / mouse `Zfp12`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF189` / mouse `Zfp189`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF706` / mouse `Zfp706`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF207` / mouse `Zfp207`: dentate delta 0.500, cerebellar delta 0.500.
- `PRP4K` / mouse `Prpf4b`: dentate delta 0.500, cerebellar delta 0.500.
- `MIR124-1HG` / mouse `Mir124a-1hg`: dentate delta 0.500, cerebellar delta 0.500.
- `C11orf58` / mouse `1110004F10Rik`: dentate delta 0.500, cerebellar delta 0.500.
- `RAB7A` / mouse `Rab7`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF292` / mouse `Zfp292`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF644` / mouse `Zfp644`: dentate delta 0.500, cerebellar delta 0.500.
- `TUBB` / mouse `Tubb5`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF326` / mouse `Zfp326`: dentate delta 0.500, cerebellar delta 0.500.
- `TOMM70` / mouse `Tomm70a`: dentate delta 0.500, cerebellar delta 0.500.
- `GPRASP3` / mouse `Bhlhb9`: dentate delta 0.500, cerebellar delta 0.500.
- `TMEM167A` / mouse `Tmem167`: dentate delta 0.500, cerebellar delta 0.500.
- `CENPC` / mouse `Cenpc1`: dentate delta 0.500, cerebellar delta 0.500.
- `ZNF280D` / mouse `Zfp280d`: dentate delta 0.500, cerebellar delta 0.500.
- `CAMLG` / mouse `Caml`: dentate delta 0.500, cerebellar delta 0.500.
- `C9orf85` / mouse `1110059E24Rik`: dentate delta 0.500, cerebellar delta 0.500.
- `NT5C3A` / mouse `Nt5c3`: dentate delta 0.500, cerebellar delta 0.500.
- `C1orf52` / mouse `2410004B18Rik`: dentate delta 0.500, cerebellar delta 0.500.
- `C16orf87` / mouse `4921524J17Rik`: dentate delta 0.500, cerebellar delta 0.500.

## Interpretation

- This pass corrects the main limitation of the earlier same-symbol full-matrix screen: mouse genes with non-identical human ortholog symbols can now contribute.
- `GSE214309` contributes symbol-resolved rows in this pass, although the source file starts with Ensembl-style rows. A dedicated Ensembl-to-symbol bridge would still be useful to rescue any residual Ensembl-only rows.
- The next step is to re-run the dataset-aware meta-model using this MGI ortholog full-matrix expression layer instead of the same-symbol layer.

## Outputs

- Expression table: `Project/results/primary_core_mgi_ortholog_full_matrix_expression.tsv.gz`
- Coverage table: `Project/results/primary_core_mgi_ortholog_full_matrix_coverage.tsv`
- Statistics table: `Project/results/primary_core_mgi_ortholog_full_matrix_statistics.tsv`
- Shared hits: `Project/results/primary_core_mgi_ortholog_full_matrix_shared_hits.tsv`
- Branch-specific hits: `Project/results/primary_core_mgi_ortholog_full_matrix_branch_specific.tsv`
- Non-identical-symbol hits: `Project/results/primary_core_mgi_ortholog_full_matrix_nonidentical_symbol_hits.tsv`
- Shared-hit plot: `Project/results/primary_core_mgi_ortholog_full_matrix_shared_hits.png`
