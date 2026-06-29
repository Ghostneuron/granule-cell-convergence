#!/usr/bin/env python3
"""Gene-level named-comparator specificity for the six Tier 1 seed genes."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

import build_primary_core_granule_specificity_named_comparators as base


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

SEED_GENES = ["GPM6A", "NFIB", "NFIA", "KCNK1", "RFX3", "GABRA2"]
SEED_MOUSE = {
    "GPM6A": "Gpm6a",
    "NFIB": "Nfib",
    "NFIA": "Nfia",
    "KCNK1": "Kcnk1",
    "RFX3": "Rfx3",
    "GABRA2": "Gabra2",
}

OUT_UNITS = RESULTS / "primary_core_seed_gene_named_comparator_units.tsv"
OUT_SUMMARY = RESULTS / "primary_core_seed_gene_named_comparator_summary.tsv"
OUT_MD = RESULTS / "primary_core_seed_gene_named_comparator_specificity.md"


def role_order(role: str) -> int:
    return {
        "dentate_granule": 0,
        "pyramidal_comparator": 1,
        "cerebellar_granule": 2,
        "purkinje_comparator": 3,
        "other_local_cell_type": 4,
    }.get(role, 9)


def compute_gene_units(
    *,
    dataset: str,
    sample: str,
    expression: pd.DataFrame,
    cell_groups: pd.Series,
) -> list[dict[str, object]]:
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    index_by_canon = {base.canon(gene): gene for gene in expression.index}
    common_cells = [cell for cell in expression.columns if cell in set(cell_groups.index)]
    if not common_cells:
        return []
    expression = expression[common_cells]
    cell_groups = cell_groups.loc[common_cells]

    records: list[dict[str, object]] = []
    for human_gene in SEED_GENES:
        mouse_gene = SEED_MOUSE[human_gene]
        present_gene = index_by_canon.get(base.canon(mouse_gene))
        if present_gene is None:
            continue
        values = np.log1p(expression.loc[present_gene].to_numpy(dtype=float))
        score_df = pd.DataFrame(
            {
                "cell_id": common_cells,
                "source_group": cell_groups.to_numpy(),
                "log1p_expression": values,
                "detected": values > 0,
            }
        )
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("log1p_expression", "size"),
                median_log1p_expression=("log1p_expression", "median"),
                mean_log1p_expression=("log1p_expression", "mean"),
                detection_fraction=("detected", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_gene_rank"] = group_summary["median_log1p_expression"].rank(
            pct=True, method="average"
        )
        for _, row in group_summary.iterrows():
            group = str(row["source_group"])
            records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "source_group": group,
                    "specificity_role": base.group_role(dataset, group),
                    "human_gene": human_gene,
                    "mouse_gene": mouse_gene,
                    "source_gene_symbol": present_gene,
                    "n_cells": int(row["n_cells"]),
                    "median_log1p_expression": float(row["median_log1p_expression"]),
                    "mean_log1p_expression": float(row["mean_log1p_expression"]),
                    "detection_fraction": float(row["detection_fraction"]),
                    "within_sample_gene_rank": float(row["within_sample_gene_rank"]),
                }
            )
    return records


def gse104323_gene_units() -> list[dict[str, object]]:
    meta = pd.read_csv(base.GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    expression = base.read_selected_rows(base.GSE104323_EXPR, "\t", set(SEED_MOUSE.values()))
    return compute_gene_units(dataset="GSE104323", sample="10X_all_cells", expression=expression, cell_groups=cell_groups)


def gse122357_gene_units() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in base.GSE122357_FILES.items():
            expression = base.read_selected_rows_from_tar(base.GSE122357_TAR, member, set(SEED_MOUSE.values()), tmpdir)
            labels = base.load_gse122357_label_map(prefix)
            records.extend(compute_gene_units(dataset="GSE122357", sample=sample, expression=expression, cell_groups=labels))
    return records


def role_median(sub: pd.DataFrame, role: str) -> float:
    values = sub.loc[sub["specificity_role"].eq(role), "within_sample_gene_rank"].dropna()
    return float(values.median()) if not values.empty else np.nan


def build_summary(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for human_gene, sub in units.groupby("human_gene", sort=False):
        dentate = role_median(sub, "dentate_granule")
        pyramidal = role_median(sub, "pyramidal_comparator")
        cerebellar = role_median(sub, "cerebellar_granule")
        purkinje = role_median(sub, "purkinje_comparator")
        dentate_delta = dentate - pyramidal if pd.notna(dentate) and pd.notna(pyramidal) else np.nan
        cerebellar_delta = cerebellar - purkinje if pd.notna(cerebellar) and pd.notna(purkinje) else np.nan
        if pd.notna(dentate_delta) and pd.notna(cerebellar_delta) and dentate_delta > 0 and cerebellar_delta > 0:
            call = "granule_enriched_vs_named_comparators"
        elif pd.notna(dentate_delta) and pd.notna(cerebellar_delta) and dentate_delta < 0 and cerebellar_delta < 0:
            call = "named_comparator_enriched_in_both_branches"
        elif pd.notna(dentate_delta) and dentate_delta < 0 and pd.notna(cerebellar_delta) and cerebellar_delta > 0:
            call = "pyramidal_enriched_but_cerebellar_granule_enriched"
        elif pd.notna(dentate_delta) and dentate_delta > 0 and pd.notna(cerebellar_delta) and cerebellar_delta < 0:
            call = "dentate_granule_enriched_but_purkinje_enriched"
        else:
            call = "mixed_or_tied"
        rows.append(
            {
                "human_gene": human_gene,
                "mouse_gene": SEED_MOUSE[human_gene],
                "dentate_granule_median_rank": dentate,
                "pyramidal_comparator_median_rank": pyramidal,
                "cerebellar_granule_median_rank": cerebellar,
                "purkinje_comparator_median_rank": purkinje,
                "dentate_vs_pyramidal_rank_delta": dentate_delta,
                "cerebellar_vs_purkinje_rank_delta": cerebellar_delta,
                "specificity_call": call,
            }
        )
    return pd.DataFrame(rows).sort_values("human_gene")


def write_report(summary: pd.DataFrame, units: pd.DataFrame) -> None:
    lines = [
        "# Seed Gene Named-Comparator Specificity",
        "",
        "This focused analysis scores the six conservative Tier 1 seed genes in the two datasets that retain named local comparator labels: `GSE104323` for dentate granule-lineage versus pyramidal comparators and `GSE122357` for cerebellar granule-lineage versus Purkinje comparators.",
        "",
        "Ranks are within-sample ranks of group-level median log1p expression for each gene. A positive delta means the granule-lineage groups are higher than the named comparator in that branch.",
        "",
        "## Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- `{row['human_gene']}`: {row['specificity_call']}; "
            f"dentate-vs-pyramidal delta {row['dentate_vs_pyramidal_rank_delta']:.3f}, "
            f"cerebellar-vs-Purkinje delta {row['cerebellar_vs_purkinje_rank_delta']:.3f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The seed set is not a Purkinje- or pyramidal-cell enriched signature as a group.",
            "- Only `NFIA` is granule-enriched against both named comparator branches in this focused rank test.",
            "- `GPM6A` is relatively higher in the named comparators in both branches, `GABRA2` is higher in pyramidal comparators but tied with Purkinje cells, `KCNK1` is tied in dentate but higher in Purkinje cells, and `NFIB`/`RFX3` are mixed or tied.",
            "- This supports treating the seed set as a shared neuronal assembly/configuration program rather than a uniquely granule-cell-exclusive marker panel.",
            "",
            "## Outputs",
            "",
            f"- Unit table: `{OUT_UNITS.relative_to(ROOT)}`",
            f"- Summary table: `{OUT_SUMMARY.relative_to(ROOT)}`",
            "",
            f"Total gene/group units: {len(units):,}.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    records = []
    records.extend(gse104323_gene_units())
    records.extend(gse122357_gene_units())
    units = pd.DataFrame(records)
    units["role_order"] = units["specificity_role"].map(role_order)
    units = units.sort_values(["human_gene", "dataset", "sample", "role_order", "source_group"]).drop(columns=["role_order"])
    summary = build_summary(units)
    units.to_csv(OUT_UNITS, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    write_report(summary, units)
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
