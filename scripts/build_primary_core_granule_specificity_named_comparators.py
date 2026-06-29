#!/usr/bin/env python3
"""Test whether mechanism axes are enriched in granule cells versus named comparators.

This is a first-pass named-comparator specificity analysis using local datasets
that retain explicit source labels:

- GSE104323: dentate granule-lineage groups versus CA3/immature pyramidal cells.
- GSE122357: cerebellar granule/precursor groups versus Purkinje cells.

The analysis uses the Tier 1-4 mechanism-axis genes from the manuscript packet.
It is pathway/module-level, not raw differential expression.
"""

from __future__ import annotations

import os
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

AXIS_GENES = RESULTS / "primary_core_mechanism_axis_gene_table.tsv"

GSE104323_EXPR = EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz"
GSE104323_META = EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz"

GSE122357_TAR = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
GSE122357_LABELS = EXTERNAL / "GEO/GSE122357/GSE122357_cell_number.xlsx"

OUT_UNITS = RESULTS / "primary_core_granule_specificity_named_comparator_units.tsv"
OUT_AXIS = RESULTS / "primary_core_granule_specificity_named_comparator_axis_summary.tsv"
OUT_COVERAGE = RESULTS / "primary_core_granule_specificity_named_comparator_gene_coverage.tsv"
OUT_PLOT = RESULTS / "primary_core_granule_specificity_named_comparators.png"
OUT_MD = RESULTS / "primary_core_granule_specificity_named_comparators.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


AXIS_ORDER = [
    "developmental_regulatory_control",
    "neurite_cytoskeleton_morphogenesis",
    "axon_guidance_adhesion",
    "synaptic_excitability_maturation",
]

AXIS_LABELS = {
    "developmental_regulatory_control": "Developmental regulatory",
    "neurite_cytoskeleton_morphogenesis": "Neurite/cytoskeleton",
    "axon_guidance_adhesion": "Axon guidance/adhesion",
    "synaptic_excitability_maturation": "Synaptic/excitability",
}

GSE104323_GRANULE_GROUPS = {"GC-adult", "GC-juv", "Immature-GC", "Neuroblast"}
GSE104323_PYRAMIDAL_GROUPS = {"CA3-Pyr", "Immature-Pyr"}

GSE122357_GRANULE_GROUPS = {"Granule cells", "Granule precursor"}
GSE122357_PURKINJE_GROUPS = {"Purkinje cells"}

GSE122357_FILES = {
    "P0": ("GSM3464549_P0.csv.gz", "p0"),
    "P8a": ("GSM3464550_P8a.csv.gz", "p8a"),
    "P8b": ("GSM3464551_P8b.csv.gz", "p8b"),
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().upper()


def group_role(dataset: str, group: str) -> str:
    if dataset == "GSE104323":
        if group in GSE104323_GRANULE_GROUPS:
            return "dentate_granule"
        if group in GSE104323_PYRAMIDAL_GROUPS:
            return "pyramidal_comparator"
        return "other_local_cell_type"
    if dataset == "GSE122357":
        if group in GSE122357_GRANULE_GROUPS:
            return "cerebellar_granule"
        if group in GSE122357_PURKINJE_GROUPS:
            return "purkinje_comparator"
        return "other_local_cell_type"
    return "other_local_cell_type"


def load_axis_gene_sets() -> tuple[dict[str, list[str]], pd.DataFrame]:
    df = pd.read_csv(AXIS_GENES, sep="\t")
    df = df.loc[df["tier_rank"].le(4)].copy()
    df["mouse_symbol_clean"] = df["mouse_symbol"].fillna(df["gene"]).astype(str)
    axis_sets: dict[str, list[str]] = {}
    for axis in AXIS_ORDER:
        genes = sorted(set(df.loc[df["mechanism_axis"].eq(axis), "mouse_symbol_clean"].dropna().astype(str)))
        axis_sets[axis] = genes
    return axis_sets, df


def read_selected_rows(path: Path, sep: str, wanted_genes: set[str], chunksize: int = 500) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    wanted_canon = {canon(gene) for gene in wanted_genes}
    for chunk in pd.read_csv(path, sep=sep, chunksize=chunksize, low_memory=False):
        gene_col = chunk.columns[0]
        chunk[gene_col] = chunk[gene_col].astype(str)
        sub = chunk.loc[chunk[gene_col].map(canon).isin(wanted_canon)].copy()
        if not sub.empty:
            sub = sub.drop_duplicates(gene_col).set_index(gene_col)
            pieces.append(sub)
    if not pieces:
        return pd.DataFrame()
    out = pd.concat(pieces, axis=0)
    out = out.loc[~out.index.duplicated(keep="first")]
    return out.apply(pd.to_numeric, errors="coerce").fillna(0)


def read_selected_rows_from_tar(
    tar_path: Path, member_name: str, wanted_genes: set[str], tmpdir: Path
) -> pd.DataFrame:
    with tarfile.open(tar_path) as tar:
        member = tar.getmember(member_name)
        source = tar.extractfile(member)
        if source is None:
            raise FileNotFoundError(f"Could not read {member_name} from {tar_path}")
        path = tmpdir / member_name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(source.read())
    return read_selected_rows(path, ",", wanted_genes)


def compute_group_axis_units(
    *,
    dataset: str,
    sample: str,
    expression: pd.DataFrame,
    cell_groups: pd.Series,
    axis_sets: dict[str, list[str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    records: list[dict[str, object]] = []
    coverage: list[dict[str, object]] = []
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    index_by_canon = {canon(gene): gene for gene in expression.index}

    common_cells = [cell for cell in expression.columns if cell in set(cell_groups.index)]
    if not common_cells:
        return records, coverage
    expression = expression[common_cells]
    cell_groups = cell_groups.loc[common_cells]

    for axis, genes in axis_sets.items():
        present = [index_by_canon[canon(gene)] for gene in genes if canon(gene) in index_by_canon]
        coverage.append(
            {
                "dataset": dataset,
                "sample": sample,
                "mechanism_axis": axis,
                "mechanism_axis_label": AXIS_LABELS.get(axis, axis),
                "n_axis_genes": len(genes),
                "n_present_genes": len(present),
                "present_genes": ",".join(present),
                "missing_genes": ",".join(gene for gene in genes if canon(gene) not in index_by_canon),
            }
        )
        if not present:
            continue
        cell_scores = np.log1p(expression.loc[present].to_numpy(dtype=float)).mean(axis=0)
        score_df = pd.DataFrame(
            {
                "cell_id": common_cells,
                "source_group": cell_groups.to_numpy(),
                "axis_score": cell_scores,
            }
        )
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("axis_score", "size"),
                median_axis_score=("axis_score", "median"),
                mean_axis_score=("axis_score", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_axis_rank"] = group_summary["median_axis_score"].rank(pct=True, method="average")
        for _, row in group_summary.iterrows():
            group = str(row["source_group"])
            records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "source_group": group,
                    "specificity_role": group_role(dataset, group),
                    "mechanism_axis": axis,
                    "mechanism_axis_label": AXIS_LABELS.get(axis, axis),
                    "n_present_axis_genes": len(present),
                    "n_cells": int(row["n_cells"]),
                    "median_axis_score": float(row["median_axis_score"]),
                    "mean_axis_score": float(row["mean_axis_score"]),
                    "within_sample_axis_rank": float(row["within_sample_axis_rank"]),
                }
            )
    return records, coverage


def gse104323_units(axis_sets: dict[str, list[str]], wanted_genes: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    meta = pd.read_csv(GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    expression = read_selected_rows(GSE104323_EXPR, "\t", wanted_genes)
    return compute_group_axis_units(
        dataset="GSE104323",
        sample="10X_all_cells",
        expression=expression,
        cell_groups=cell_groups,
        axis_sets=axis_sets,
    )


def load_gse122357_label_map(sample_prefix: str) -> pd.Series:
    labels = pd.read_excel(GSE122357_LABELS, sheet_name="Sheet1", dtype=str)
    records: dict[str, str] = {}
    for col in labels.columns:
        for value in labels[col].dropna().astype(str):
            value = value.strip()
            prefix = f"{sample_prefix}_"
            if not value.lower().startswith(prefix):
                continue
            barcode = value[len(prefix) :]
            if barcode:
                records[barcode] = col
    return pd.Series(records, dtype=str)


def gse122357_units(axis_sets: dict[str, list[str]], wanted_genes: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    all_records: list[dict[str, object]] = []
    all_coverage: list[dict[str, object]] = []
    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in GSE122357_FILES.items():
            expression = read_selected_rows_from_tar(GSE122357_TAR, member, wanted_genes, tmpdir)
            labels = load_gse122357_label_map(prefix)
            records, coverage = compute_group_axis_units(
                dataset="GSE122357",
                sample=sample,
                expression=expression,
                cell_groups=labels,
                axis_sets=axis_sets,
            )
            all_records.extend(records)
            all_coverage.extend(coverage)
    return all_records, all_coverage


def build_axis_summary(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for axis, sub in units.groupby("mechanism_axis", sort=False):
        def med(role: str) -> float:
            values = sub.loc[sub["specificity_role"].eq(role), "within_sample_axis_rank"].dropna()
            return float(values.median()) if not values.empty else np.nan

        dentate = med("dentate_granule")
        pyramidal = med("pyramidal_comparator")
        cereb = med("cerebellar_granule")
        purkinje = med("purkinje_comparator")
        granule_min = np.nanmin([dentate, cereb]) if pd.notna(dentate) or pd.notna(cereb) else np.nan
        comparator_max = np.nanmax([pyramidal, purkinje]) if pd.notna(pyramidal) or pd.notna(purkinje) else np.nan
        dentate_vs_pyramidal = dentate - pyramidal if pd.notna(dentate) and pd.notna(pyramidal) else np.nan
        cerebellar_vs_purkinje = cereb - purkinje if pd.notna(cereb) and pd.notna(purkinje) else np.nan
        index = granule_min - comparator_max if pd.notna(granule_min) and pd.notna(comparator_max) else np.nan
        if pd.notna(dentate_vs_pyramidal) and pd.notna(cerebellar_vs_purkinje) and dentate_vs_pyramidal > 0 and cerebellar_vs_purkinje > 0:
            call = "granule_enriched_vs_named_comparators"
        elif pd.notna(dentate_vs_pyramidal) and dentate_vs_pyramidal > 0:
            call = "dentate_granule_enriched_but_not_cerebellar_vs_purkinje"
        elif pd.notna(cerebellar_vs_purkinje) and cerebellar_vs_purkinje > 0:
            call = "cerebellar_granule_enriched_but_not_dentate_vs_pyramidal"
        else:
            call = "not_granule_specific_vs_named_comparators"

        rows.append(
            {
                "mechanism_axis": axis,
                "mechanism_axis_label": AXIS_LABELS.get(axis, axis),
                "dentate_granule_median_rank": dentate,
                "pyramidal_comparator_median_rank": pyramidal,
                "cerebellar_granule_median_rank": cereb,
                "purkinje_comparator_median_rank": purkinje,
                "dentate_vs_pyramidal_rank_delta": dentate_vs_pyramidal,
                "cerebellar_vs_purkinje_rank_delta": cerebellar_vs_purkinje,
                "granule_min_rank": granule_min,
                "comparator_max_rank": comparator_max,
                "granule_specificity_index": index,
                "specificity_call": call,
            }
        )
    out = pd.DataFrame(rows)
    out["axis_order"] = out["mechanism_axis"].map({axis: i for i, axis in enumerate(AXIS_ORDER)}).fillna(99)
    return out.sort_values("axis_order").drop(columns=["axis_order"])


def plot_units(units: pd.DataFrame) -> None:
    plot_roles = ["dentate_granule", "pyramidal_comparator", "cerebellar_granule", "purkinje_comparator"]
    pivot = (
        units.loc[units["specificity_role"].isin(plot_roles)]
        .groupby(["mechanism_axis", "specificity_role"], sort=False)["within_sample_axis_rank"]
        .median()
        .unstack("specificity_role")
        .reindex(index=AXIS_ORDER, columns=plot_roles)
    )
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(plot_roles)))
    ax.set_xticklabels([role.replace("_", " ") for role in plot_roles], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([AXIS_LABELS.get(axis, axis) for axis in pivot.index])
    ax.set_title("Named-comparator module specificity")
    for i, axis in enumerate(pivot.index):
        for j, role in enumerate(plot_roles):
            value = pivot.loc[axis, role]
            text = "NA" if pd.isna(value) else f"{value:.2f}"
            color = "#f6f6f6" if pd.notna(value) and value < 0.35 else "#202020"
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Median within-sample rank")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(units: pd.DataFrame, axis_summary: pd.DataFrame, coverage: pd.DataFrame) -> None:
    n_gse104323 = units.loc[units["dataset"].eq("GSE104323"), "source_group"].nunique()
    n_gse122357 = units.loc[units["dataset"].eq("GSE122357"), "source_group"].nunique()
    calls = axis_summary["specificity_call"].value_counts().to_dict()
    lines = [
        "# Granule Specificity Named-Comparator Analysis",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This first-pass analysis asks whether the mechanism-axis similarity between dentate and cerebellar granule cells is higher than in named non-granule neuronal comparators. It uses only local datasets with explicit source labels for the relevant named groups.",
        "",
        "## Local Named Comparators Used",
        "",
        "- `GSE104323`: dentate granule-lineage groups (`GC-adult`, `GC-juv`, `Immature-GC`, `Neuroblast`) versus pyramidal comparators (`CA3-Pyr`, `Immature-Pyr`).",
        "- `GSE122357`: cerebellar `Granule cells` and `Granule precursor` versus `Purkinje cells` across P0, P8a, and P8b.",
        "",
        "The analysis scores Tier 1-4 mechanism-axis genes, then ranks source groups within each sample and axis. A positive specificity result requires both dentate granule > pyramidal comparator and cerebellar granule > Purkinje comparator.",
        "",
        "## Scope",
        "",
        f"- GSE104323 source groups scored: {n_gse104323}.",
        f"- GSE122357 source groups scored: {n_gse122357}.",
        f"- Mechanism axes tested: {axis_summary['mechanism_axis'].nunique()}.",
        f"- Specificity calls: {', '.join(f'{key}: {value}' for key, value in calls.items())}.",
        "",
        "## Axis Results",
        "",
    ]
    for _, row in axis_summary.iterrows():
        lines.append(
            f"- {row['mechanism_axis_label']}: {row['specificity_call']}; "
            f"dentate-vs-pyramidal delta {row['dentate_vs_pyramidal_rank_delta']:.3f}, "
            f"cerebellar-vs-Purkinje delta {row['cerebellar_vs_purkinje_rank_delta']:.3f}, "
            f"granule specificity index {row['granule_specificity_index']:.3f}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This analysis directly addresses the user's specificity question for the named local comparators that are currently available.",
            "- A pathway is treated as granule-enriched only if both granule branches exceed their named local comparator.",
            "- If an axis fails this test, it may still be biologically important, but it should be phrased as a broader neuronal maturation or morphology pathway rather than a unique granule-cell pathway.",
            "- This is still a module-level analysis. A stronger next step would add more explicit pyramidal/Purkinje datasets or Allen expression matrices and test gene-level specificity with raw-count models.",
            "",
            "## Outputs",
            "",
            f"- Unit table: `{rel(OUT_UNITS)}`",
            f"- Axis summary: `{rel(OUT_AXIS)}`",
            f"- Gene coverage: `{rel(OUT_COVERAGE)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    axis_sets, axis_gene_df = load_axis_gene_sets()
    wanted_genes = set(axis_gene_df["mouse_symbol_clean"].dropna().astype(str))

    records: list[dict[str, object]] = []
    coverage_records: list[dict[str, object]] = []
    rec, cov = gse104323_units(axis_sets, wanted_genes)
    records.extend(rec)
    coverage_records.extend(cov)
    rec, cov = gse122357_units(axis_sets, wanted_genes)
    records.extend(rec)
    coverage_records.extend(cov)

    units = pd.DataFrame(records)
    coverage = pd.DataFrame(coverage_records)
    axis_summary = build_axis_summary(units)

    units.to_csv(OUT_UNITS, sep="\t", index=False)
    axis_summary.to_csv(OUT_AXIS, sep="\t", index=False)
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    plot_units(units)
    write_report(units, axis_summary, coverage)

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Units: {len(units):,}")
    print(axis_summary[["mechanism_axis", "specificity_call", "granule_specificity_index"]].to_string(index=False))


if __name__ == "__main__":
    main()
