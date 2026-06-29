# Primary-Core Candidate-Gene Pseudobulk Analysis

Date built: 2026-06-22

## Scope

This layer aggregates the 67 candidate genes from the human bridge packet by refined candidate-cell class across the strict 10-dataset primary core.

It is DE-adjacent but not yet a genome-wide mixed-effect differential-expression model. It tests whether the proposed structural-executor genes are elevated in dentate and cerebellar candidate granule populations relative to local backgrounds.

Minimum broad-class size used for rank statistics: 20 cells/nuclei.

## Coverage

- Primary datasets represented by at least one candidate gene: 10/10.
- Pseudobulk expression rows: 16,719.
- Candidate genes tested: 67.
- Structural-executor genes with positive dentate and cerebellar candidate deltas: 14.
- Structural-executor genes passing the stricter exploratory BH<0.2 rule in both branches: 12.

- `GSE104323` / `10X_all_cells`: 66/67 candidate genes, 24185/24185 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464549_P0`: 64/67 candidate genes, 7468/7468 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464550_P8a`: 64/67 candidate genes, 6088/6088 labeled observations (`full_raw_matrix`).
- `GSE122357` / `GSM3464551_P8b`: 64/67 candidate genes, 8168/8168 labeled observations (`full_raw_matrix`).
- `GSE165657` / `Cerebellum_aggr`: 66/67 candidate genes, 73550/73550 labeled observations (`full_raw_matrix`).
- `GSE186538` / `selected_bridge_all_samples`: 67/67 candidate genes, 32004/32004 labeled observations (`selected_norm_bridge`).
- `GSE214309` / `snRNA_counts`: 64/67 candidate genes, 401/401 labeled observations (`full_raw_matrix`).
- `GSE268609` / `selected_bridge_all_samples`: 67/67 candidate genes, 366175/366175 labeled observations (`selected_gene_bridge`).
- `GSE292261` / `SS2_filtered_counts`: 56/67 candidate genes, 509/509 labeled observations (`full_raw_matrix`).
- `GSE312658` / `Ctrl`: 66/67 candidate genes, 4857/4857 labeled observations (`full_raw_matrix`).
- `GSE312658` / `cKO`: 66/67 candidate genes, 6946/6946 labeled observations (`full_raw_matrix`).
- `GSE325391` / `selected_bridge_all_samples`: 67/67 candidate genes, 59075/59075 labeled observations (`selected_gene_bridge`).
- `GSE95752` / `C1_all_cells`: 61/67 candidate genes, 2303/2303 labeled observations (`full_raw_matrix`).

## Top Shared-Executor Signals

- `Cfl1` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.500, shared=True.
- `Gap43` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.500, shared=True.
- `Robo2` (axon_guidance_synapse): dentate delta 0.500, cerebellar delta 0.500, shared=True.
- `Stmn2` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.500, shared=True.
- `Stmn3` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.500, shared=True.
- `Cdk5r1` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.333, shared=True.
- `Dpysl2` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.333, shared=True.
- `Ephb2` (axon_guidance_synapse): dentate delta 0.500, cerebellar delta 0.333, shared=True.
- `L1cam` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.333, shared=True.
- `Mapt` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.167, shared=True.
- `Dpysl3` (morphogenesis_cytoskeleton): dentate delta 0.333, cerebellar delta 0.167, shared=True.
- `Elavl4` (morphogenesis_cytoskeleton): dentate delta 0.167, cerebellar delta 0.500, shared=True.
- `Ncam1` (morphogenesis_cytoskeleton): dentate delta 0.500, cerebellar delta 0.167, shared=True.
- `Cntn6` (axon_guidance_synapse): dentate delta 0.417, cerebellar delta 0.167, shared=True.

## Interpretation

- The useful positive result is not that dentate and cerebellar granule cells are transcriptionally identical. The module analysis already argues against that. The stronger claim is that both lineages repeatedly use an elevated structural-executor gene axis on top of distinct regional identity programs.
- This candidate-gene pass supports prioritizing shared structural genes whose rank is positive in both dentate-candidate and cerebellar-candidate pseudobulks.
- `GSE214309` contributes gene-symbol-resolved candidate-gene pseudobulk evidence in this focused pass. Its broader whole-transcriptome cleanup should still be revisited before genome-wide DE.
- The next stricter analysis should expand from these 67 genes to genome-wide ortholog-aware pseudobulk DE, with donor/sample/stage modeled explicitly.

## Outputs

- Expression table: `Project/results/primary_core_candidate_gene_pseudobulk_expression.tsv`
- Coverage table: `Project/results/primary_core_candidate_gene_pseudobulk_coverage.tsv`
- Gene statistics: `Project/results/primary_core_candidate_gene_pseudobulk_statistics.tsv`
- Ranked hits: `Project/results/primary_core_candidate_gene_pseudobulk_hits.tsv`
- Effect plot: `Project/results/primary_core_candidate_gene_pseudobulk_effects.png`
