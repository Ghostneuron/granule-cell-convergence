# Seed Gene Named-Comparator Specificity

This focused analysis scores the six conservative Tier 1 seed genes in the two datasets that retain named local comparator labels: `GSE104323` for dentate granule-lineage versus pyramidal comparators and `GSE122357` for cerebellar granule-lineage versus Purkinje comparators.

Ranks are within-sample ranks of group-level median log1p expression for each gene. A positive delta means the granule-lineage groups are higher than the named comparator in that branch.

## Summary

- `GABRA2`: mixed_or_tied; dentate-vs-pyramidal delta -0.208, cerebellar-vs-Purkinje delta 0.000.
- `GPM6A`: named_comparator_enriched_in_both_branches; dentate-vs-pyramidal delta -0.135, cerebellar-vs-Purkinje delta -0.125.
- `KCNK1`: mixed_or_tied; dentate-vs-pyramidal delta 0.000, cerebellar-vs-Purkinje delta -0.375.
- `NFIA`: granule_enriched_vs_named_comparators; dentate-vs-pyramidal delta 0.219, cerebellar-vs-Purkinje delta 0.719.
- `NFIB`: mixed_or_tied; dentate-vs-pyramidal delta 0.000, cerebellar-vs-Purkinje delta 0.562.
- `RFX3`: mixed_or_tied; dentate-vs-pyramidal delta 0.000, cerebellar-vs-Purkinje delta 0.000.

## Interpretation

- The seed set is not a Purkinje- or pyramidal-cell enriched signature as a group.
- Only `NFIA` is granule-enriched against both named comparator branches in this focused rank test.
- `GPM6A` is relatively higher in the named comparators in both branches, `GABRA2` is higher in pyramidal comparators but tied with Purkinje cells, `KCNK1` is tied in dentate but higher in Purkinje cells, and `NFIB`/`RFX3` are mixed or tied.
- This supports treating the seed set as a shared neuronal assembly/configuration program rather than a uniquely granule-cell-exclusive marker panel.

## Outputs

- Unit table: `Project/results/primary_core_seed_gene_named_comparator_units.tsv`
- Summary table: `Project/results/primary_core_seed_gene_named_comparator_summary.tsv`

Total gene/group units: 288.
