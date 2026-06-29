#!/usr/bin/env python3
"""Prepare GSE186538 human DG granule-cell candidate metadata."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "External_Data/GEO/GSE186538/GSE186538_Human_cell_meta.txt.gz"
RESULTS = ROOT / "Project/results"

OUT_CELLS = RESULTS / "gse186538_human_dg_gc_candidate_cells.tsv"
OUT_SUMMARY = RESULTS / "gse186538_human_dg_gc_candidate_summary.tsv"
OUT_MD = RESULTS / "gse186538_human_dg_gc_candidate_summary.md"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(META, sep="\t")
    meta.insert(0, "cell_index_1based", range(1, len(meta) + 1))

    is_dg = meta["region"].eq("DG")
    is_dg_gc = meta["cluster"].astype(str).str.startswith("DG GC")
    candidates = meta.loc[is_dg & is_dg_gc].copy()
    candidates.to_csv(OUT_CELLS, sep="\t", index=False)

    summary = (
        candidates.groupby(["samplename", "cluster"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["cluster", "samplename"])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    cluster_counts = candidates["cluster"].value_counts().to_dict()
    sample_counts = candidates["samplename"].value_counts().to_dict()

    lines = [
        "# GSE186538 human DG granule-cell candidates",
        "",
        "Date prepared: 2026-06-21",
        "",
        "## Selection Rule",
        "",
        "Cells were selected from the human metadata using `region == DG` and `cluster` names beginning with `DG GC`.",
        "",
        "## Counts",
        "",
        f"- Candidate DG granule cells: {len(candidates)}",
        f"- Source human metadata cells: {len(meta)}",
        "",
        "## Cluster Counts",
        "",
    ]
    for cluster, n_cells in cluster_counts.items():
        lines.append(f"- `{cluster}`: {n_cells}")
    lines.extend(["", "## Sample Counts", ""])
    for sample, n_cells in sample_counts.items():
        lines.append(f"- `{sample}`: {n_cells}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Candidate cell table: `{OUT_CELLS.relative_to(ROOT)}`",
            f"- Candidate summary table: `{OUT_SUMMARY.relative_to(ROOT)}`",
            "",
            "The count matrix is already local as Matrix Market format and is oriented genes-by-cells. The next compute-heavy step is to stream-extract these 32,067 candidate columns from the full matrix.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote {OUT_CELLS}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_MD}")
    print(f"candidate_cells={len(candidates)}")


if __name__ == "__main__":
    main()
