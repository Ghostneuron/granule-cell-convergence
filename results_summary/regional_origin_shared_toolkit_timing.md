# Regional-Origin Versus Shared-Toolkit Timing Analysis

## Purpose

This analysis tests whether dentate and cerebellar granule cells look like a single recent progenitor lineage or instead retain different regional fate polarity while recruiting shared granule/maturation genes during postmitotic development.

## Main Findings

- Dentate progenitor states do not yet show strong dentate-over-cerebellar fate polarity (median -0.177), but dentate postmitotic/mature granule states do (median 0.479; shift 0.656).
- Cerebellar precursor/P0 states show cerebellar-over-dentate fate polarity (median 0.375), whereas labeled maturing granule cells are lower or near tied (median -0.062).
- Dentate `NEUROD1`/`RBFOX3` postmitotic-special rank increases from progenitor to postmitotic/mature states by 0.375; construction rank increases by 0.302.
- Cerebellar `NEUROD1`/`RBFOX3` postmitotic-special rank changes from P0/precursor to maturing states by 0.062; construction rank increases by 0.156.
- The four-gene shared-special average changes by 0.057 in dentate and -0.156 in cerebellum because `NFIA`/`HMGN2` behave more like early regulatory/chromatin-state genes while `NEUROD1`/`RBFOX3` track postmitotic maturation more directly.

## Interpretation

- The timing pattern argues against a single recent shared dentate/cerebellar progenitor in the sampled states. Branch-specific fate polarity appears in different lineage windows rather than as a shared root state.
- `NFIA`, `NEUROD1`, `RBFOX3`, and `HMGN2` are best interpreted as a reused stage-composite toolkit: `NFIA`/`HMGN2` mark regulatory/chromatin competence, while `NEUROD1`/`RBFOX3` better mark postmitotic neuronal maturation.
- The result supports the manuscript model: distinct regional origins, followed by partial convergence through shared neurogenic/maturation and downstream construction layers.
- The cerebellar timing axis remains limited by three postnatal `GSE122357` samples, so the direction of cerebellar changes should be treated as supportive rather than definitive.

## Outputs

- Module timing units: `Project/results/regional_origin_shared_toolkit_timing_module_units.tsv`
- Gene timing units: `Project/results/regional_origin_shared_toolkit_timing_gene_units.tsv`
- State summary: `Project/results/regional_origin_shared_toolkit_timing_state_summary.tsv`
- Timing metrics: `Project/results/regional_origin_shared_toolkit_timing_metrics.tsv`
- Plot: `Project/results/regional_origin_shared_toolkit_timing.png`
