# Granule-cell comparator-relative recurrence

This repository contains analysis code and derived tables for the manuscript
**"Limited molecular recurrence in distinct dentate and cerebellar granule
cells relative to regional comparators"**.

The current, manuscript-matched analysis is in [`dgd_reanalysis/`](dgd_reanalysis/).
The older project directories are retained as an exploratory archive and are
not the evidentiary basis of the current manuscript.

## Main result

Dentate and cerebellar granule cells do not form a single molecular identity.
Instead, a limited set of genes recurs in the same direction when each granule
population is compared with its own regional neuronal comparator. This
comparator-relative signal was positive in seven independent datasets, while
an independent adult mouse Allen common-matrix test did not support broad
direct adult module convergence. The result therefore identifies constrained
molecular reuse, not a shared recent lineage or a demonstrated causal
developmental mechanism.

## Current packet

- `dgd_reanalysis/Project/scripts/`: robustness, Allen common-matrix, figure,
  and supplementary-table scripts.
- `dgd_reanalysis/Project/results/`: compact inputs and derived result tables
  used by the current figures and supplementary tables.
- `dgd_reanalysis/Project/manuscript/source_tables/`: Tables S1-S12 as TSV
  files.
- `dgd_reanalysis/requirements.txt`: Python dependencies.

Raw public single-cell matrices and large intermediate expression matrices are
not redistributed. Dataset accessions are listed in Table S1, and Allen input
requirements are documented in the packet README.

## Citation

Please cite the associated manuscript and this repository. Author metadata are
provided in `CITATION.cff`.

## License status

No software license has yet been assigned. Contact the author before reusing or
redistributing the code.
