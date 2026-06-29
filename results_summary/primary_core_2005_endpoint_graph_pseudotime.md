# 2005 Endpoint Graph Pseudotime

## Purpose

This is the first cell-level trajectory layer for the 2005 paper endpoints. It uses paper-endpoint genes, PCA, kNN graph construction, and shortest-path distance from biologically defined root cells. It is not yet a full-transcriptome Scanpy/Monocle trajectory.

## Roots

- `GSE104323`: RGL_young/RGL dentate lineage cells.
- `GSE292261`: P5 candidate dentate granule cells.
- `GSE122357`: P0 candidate cerebellar granule cells.

## Main Findings

- `GSE104323` is the cleanest trajectory: proliferation versus graph pseudotime is rho 0.149, p 0.000; neuronal maturation is rho 0.637, p 0.000.
- `GSE292261` shows postnatal dentate timing but not a simple age-only line: proliferation is rho -0.117, p 0.102; TGF/SMAD is rho -0.153, p 0.033.
- `GSE122357` shows a cerebellar candidate trajectory from P0 roots, but P8a/P8b differ: TGF/SMAD is rho 0.347, p 0.000; secreted stop-factor axis is rho 0.121, p 0.000.

## Interpretation

The graph layer supports the user's concern: the 2005 readouts are stage/trajectory dependent. The cleanest evidence is the adult dentate lineage, where neuronal maturation rises along graph pseudotime and proliferation behaves as a trajectory-windowed state rather than a simple late-versus-early marker. The postnatal dentate and cerebellar datasets show trajectory windows rather than simple monotonic age effects, so future manuscript claims should avoid treating age as pseudotime.

The next stronger analysis should use full-transcriptome highly variable genes and a dedicated trajectory package or an equivalent diffusion/PAGA workflow, then overlay the 2005 endpoint modules, TGF-beta/BDNF, and conditioned-medium secretome candidates.

## Outputs

- Cell scores: `Project/results/primary_core_2005_endpoint_graph_pseudotime_cell_scores.tsv.gz` (19,192 cells).
- Group summary: `Project/results/primary_core_2005_endpoint_graph_pseudotime_group_summary.tsv` (16 rows).
- Module correlations: `Project/results/primary_core_2005_endpoint_graph_pseudotime_module_correlations.tsv` (27 rows).
- Plot: `Project/results/primary_core_2005_endpoint_graph_pseudotime.png`.
