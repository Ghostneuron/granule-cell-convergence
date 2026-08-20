# Granule-cell comparator-relative recurrence

This repository contains analysis code and derived tables for the manuscript
**"Comparator-relative molecular recurrence in distinct dentate and
cerebellar granule cells"**.

The current, manuscript-matched analysis is in [`dgd_reanalysis/`](dgd_reanalysis/).
The historical directory name is retained to preserve stable paths. Other
top-level project directories are an exploratory archive and are not the
evidentiary basis of the current manuscript.

## Main result

Dentate and cerebellar granule cells do not form a single molecular identity.
Instead, a limited set of genes recurs in the same direction when each granule
population is compared with its own regional neuronal comparator. This
comparator-relative signal was positive in seven independent datasets, while
an independent adult mouse Allen common-matrix test did not support broad
direct adult module convergence or a null-exceeding transferable multigene
configuration. The result identifies limited molecular recurrence relative to
regional comparators, not a shared recent lineage, a common adult molecular
class or a demonstrated causal developmental mechanism.

## Current packet

- `dgd_reanalysis/Project/scripts/`: robustness, Allen common-matrix,
  cross-region transfer, figure and repository-table scripts.
- `dgd_reanalysis/Project/results/`: compact inputs and derived result tables
  used by the current figures and machine-readable tables.
- `dgd_reanalysis/Project/manuscript/JCN_manuscript.md`: Markdown copy of the
  current Journal of Comparative Neurology submission.
- `dgd_reanalysis/Project/manuscript/main_figures/`: the four current main
  figures.
- `dgd_reanalysis/Project/manuscript/source_tables/`: Tables S1-S13 as TSV
  files. These are repository data-table identifiers, not a journal
  supplementary upload.
- `dgd_reanalysis/requirements.txt`: Python dependencies.

Raw public single-cell matrices and large intermediate expression matrices are
not redistributed. Dataset accessions are listed in Table S1, and Allen input
requirements are documented in the packet README. The older supplementary
figures and 101-table exploratory collection remain archived but are not part
of the current JCN evidence packet.

## Citation

Please cite the associated manuscript and this repository. Author metadata are
provided in `CITATION.cff`.

## License status

No software license has yet been assigned. Contact the author before reusing or
redistributing the code.
