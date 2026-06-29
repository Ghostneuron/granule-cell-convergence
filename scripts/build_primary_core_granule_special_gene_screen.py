#!/usr/bin/env python3
"""Genome-wide local screen for genes enriched in granule cells versus named comparators.

The screen uses the two local full-matrix datasets with explicit comparator
labels:

- GSE104323: dentate granule-lineage groups versus pyramidal comparators.
- GSE122357: cerebellar granule-lineage groups versus Purkinje cells.

Two contrasts are reported:

- granule_lineage: includes neuroblast/precursor states.
- postmitotic_granule: focuses on granule-cell states and excludes the most
  precursor-like groups where possible.
"""

from __future__ import annotations

import gzip
import os
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

import build_primary_core_granule_specificity_named_comparators as base


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
SUPP_FIGURES = ROOT / "Project/manuscript/Supplementary figures"

ORTHOLOG_MAP = RESULTS / "primary_core_mgi_ortholog_meta_model_map.tsv"

OUT_UNITS = RESULTS / "primary_core_granule_special_gene_named_comparator_units.tsv.gz"
OUT_SUMMARY = RESULTS / "primary_core_granule_special_gene_named_comparator_summary.tsv"
OUT_TOP = RESULTS / "primary_core_granule_special_gene_named_comparator_top_candidates.tsv"
OUT_MD = RESULTS / "primary_core_granule_special_gene_named_comparator_screen.md"
OUT_PLOT = RESULTS / "primary_core_granule_special_gene_named_comparator_top_candidates.png"
OUT_SUPP_PLOT = SUPP_FIGURES / "Fig.S2_primary_core_granule_special_gene_named_comparator_top_candidates.png"

MIN_GRANULE_DETECTION = 0.05
MIN_POSITIVE_CEREBELLAR_SAMPLES = 2
MIN_STRONG_DELTA = 0.05
CHUNKSIZE = 250

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONTRASTS = {
    "granule_lineage": {
        "label": "Granule lineage",
        "GSE104323_granule": {"Neuroblast", "Immature-GC", "GC-juv", "GC-adult"},
        "GSE104323_comparator": {"CA3-Pyr", "Immature-Pyr"},
        "GSE122357_granule": {"Granule precursor", "Granule cells"},
        "GSE122357_comparator": {"Purkinje cells"},
    },
    "postmitotic_granule": {
        "label": "Postmitotic/granule-cell state",
        "GSE104323_granule": {"Immature-GC", "GC-juv", "GC-adult"},
        "GSE104323_comparator": {"CA3-Pyr", "Immature-Pyr"},
        "GSE122357_granule": {"Granule cells"},
        "GSE122357_comparator": {"Purkinje cells"},
    },
}


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().strip('"').upper()


def clean_symbol(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().strip('"')


def load_ortholog_map() -> dict[str, dict[str, object]]:
    if not ORTHOLOG_MAP.exists():
        return {}
    ortho = pd.read_csv(ORTHOLOG_MAP, sep="\t")
    ortho = ortho.loc[ortho["mgi_one_to_one_human_mouse"].astype(bool)].copy()
    out: dict[str, dict[str, object]] = {}
    for _, row in ortho.iterrows():
        key = canon(row.get("canonical_mouse_symbol", row.get("mouse_symbol", "")))
        if key:
            out[key] = row.to_dict()
    return out


def group_indices(columns: list[str], cell_groups: pd.Series, wanted_groups: set[str]) -> dict[str, np.ndarray]:
    index_by_cell = {cell: i for i, cell in enumerate(columns)}
    out: dict[str, list[int]] = {group: [] for group in wanted_groups}
    for cell, group in cell_groups.items():
        if group in wanted_groups and cell in index_by_cell:
            out[group].append(index_by_cell[cell])
    return {group: np.array(indices, dtype=int) for group, indices in out.items() if indices}


def role_stat(
    log_values: np.ndarray,
    detected: np.ndarray,
    indices_by_group: dict[str, np.ndarray],
    groups: set[str],
) -> tuple[np.ndarray, np.ndarray, int]:
    group_medians: list[np.ndarray] = []
    group_detection: list[np.ndarray] = []
    for group in groups:
        idx = indices_by_group.get(group)
        if idx is None or len(idx) == 0:
            continue
        group_medians.append(np.median(log_values[:, idx], axis=1))
        group_detection.append(np.mean(detected[:, idx], axis=1))
    if not group_medians:
        n = log_values.shape[0]
        return np.full(n, np.nan), np.full(n, np.nan), 0
    return np.median(np.vstack(group_medians), axis=0), np.median(np.vstack(group_detection), axis=0), len(group_medians)


def process_expression_file(
    *,
    path: Path,
    sep: str,
    dataset: str,
    sample: str,
    cell_groups: pd.Series,
    ortholog: dict[str, dict[str, object]],
    chunksize: int = CHUNKSIZE,
) -> list[dict[str, object]]:
    wanted_groups: set[str] = set()
    for contrast in CONTRASTS.values():
        wanted_groups.update(contrast[f"{dataset}_granule"])
        wanted_groups.update(contrast[f"{dataset}_comparator"])

    header = pd.read_csv(path, sep=sep, nrows=0).columns.tolist()
    gene_col = header[0]
    wanted_cells = set(cell_groups.loc[cell_groups.isin(wanted_groups)].index)
    common_columns = [column for column in header[1:] if column in wanted_cells]
    if not common_columns:
        return []
    usecols = [gene_col] + common_columns
    indices_by_group = group_indices(common_columns, cell_groups.loc[cell_groups.index.intersection(common_columns)], wanted_groups)
    if not indices_by_group:
        return []

    records: list[dict[str, object]] = []
    for chunk in pd.read_csv(path, sep=sep, usecols=usecols, chunksize=chunksize, low_memory=False):
        genes = chunk[gene_col].map(clean_symbol)
        chunk = chunk.drop(columns=[gene_col])
        values = chunk[common_columns].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(dtype=float)
        log_values = np.log1p(values)
        detected = values > 0

        canonical = genes.map(canon).to_numpy()
        gene_symbols = genes.to_numpy()
        for contrast_id, contrast in CONTRASTS.items():
            granule_median, granule_detection, n_granule_groups = role_stat(
                log_values,
                detected,
                indices_by_group,
                contrast[f"{dataset}_granule"],
            )
            comp_median, comp_detection, n_comp_groups = role_stat(
                log_values,
                detected,
                indices_by_group,
                contrast[f"{dataset}_comparator"],
            )
            delta = granule_median - comp_median
            for i, canonical_mouse_symbol in enumerate(canonical):
                if not canonical_mouse_symbol or pd.isna(delta[i]):
                    continue
                ortho = ortholog.get(canonical_mouse_symbol, {})
                records.append(
                    {
                        "dataset": dataset,
                        "sample": sample,
                        "contrast_id": contrast_id,
                        "contrast_label": contrast["label"],
                        "source_gene_symbol": gene_symbols[i],
                        "canonical_mouse_symbol": canonical_mouse_symbol,
                        "human_symbol": ortho.get("human_symbol", ""),
                        "mouse_symbol": ortho.get("mouse_symbol", gene_symbols[i]),
                        "mgi_one_to_one_human_mouse": bool(ortho),
                        "granule_group_median_log1p": float(granule_median[i]),
                        "comparator_group_median_log1p": float(comp_median[i]),
                        "granule_vs_comparator_log1p_delta": float(delta[i]),
                        "granule_detection_fraction": float(granule_detection[i]),
                        "comparator_detection_fraction": float(comp_detection[i]),
                        "n_granule_groups_present": int(n_granule_groups),
                        "n_comparator_groups_present": int(n_comp_groups),
                    }
                )
    return records


def process_gse104323(ortholog: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    meta = pd.read_csv(base.GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    return process_expression_file(
        path=base.GSE104323_EXPR,
        sep="\t",
        dataset="GSE104323",
        sample="10X_all_cells",
        cell_groups=cell_groups,
        ortholog=ortholog,
    )


def extract_tar_member(tar_path: Path, member_name: str, tmpdir: Path) -> Path:
    out = tmpdir / member_name
    out.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        source = tar.extractfile(member_name)
        if source is None:
            raise FileNotFoundError(member_name)
        out.write_bytes(source.read())
    return out


def process_gse122357(ortholog: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in base.GSE122357_FILES.items():
            path = extract_tar_member(base.GSE122357_TAR, member, tmpdir)
            labels = base.load_gse122357_label_map(prefix)
            records.extend(
                process_expression_file(
                    path=path,
                    sep=",",
                    dataset="GSE122357",
                    sample=sample,
                    cell_groups=labels,
                    ortholog=ortholog,
                )
            )
    return records


def collapse_duplicate_units(units: pd.DataFrame) -> pd.DataFrame:
    sort_cols = ["dataset", "sample", "contrast_id", "canonical_mouse_symbol"]
    numeric_cols = [
        "granule_group_median_log1p",
        "comparator_group_median_log1p",
        "granule_vs_comparator_log1p_delta",
        "granule_detection_fraction",
        "comparator_detection_fraction",
    ]
    first_cols = [
        "contrast_label",
        "source_gene_symbol",
        "human_symbol",
        "mouse_symbol",
        "mgi_one_to_one_human_mouse",
        "n_granule_groups_present",
        "n_comparator_groups_present",
    ]
    agg = {col: "max" for col in numeric_cols}
    agg.update({col: "first" for col in first_cols})
    return units.groupby(sort_cols, as_index=False, sort=False).agg(agg)


def build_summary(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (contrast_id, canonical_mouse_symbol), sub in units.groupby(["contrast_id", "canonical_mouse_symbol"], sort=False):
        dentate = sub.loc[sub["dataset"].eq("GSE104323")]
        cereb = sub.loc[sub["dataset"].eq("GSE122357")]
        if dentate.empty or cereb.empty:
            continue
        drow = dentate.iloc[0]
        cereb_delta = cereb["granule_vs_comparator_log1p_delta"].to_numpy(dtype=float)
        cereb_detection = cereb["granule_detection_fraction"].to_numpy(dtype=float)
        cereb_support_samples = cereb.loc[cereb["granule_vs_comparator_log1p_delta"].gt(0), "sample"].tolist()
        dentate_delta = float(drow["granule_vs_comparator_log1p_delta"])
        cereb_median_delta = float(np.median(cereb_delta))
        support_n = len(cereb_support_samples)
        shared_positive = (
            dentate_delta > 0
            and cereb_median_delta > 0
            and support_n >= MIN_POSITIVE_CEREBELLAR_SAMPLES
            and float(drow["granule_detection_fraction"]) >= MIN_GRANULE_DETECTION
            and float(np.median(cereb_detection)) >= MIN_GRANULE_DETECTION
        )
        strong_shared = (
            shared_positive
            and dentate_delta >= MIN_STRONG_DELTA
            and cereb_median_delta >= MIN_STRONG_DELTA
            and support_n == cereb["sample"].nunique()
        )
        specificity_score = min(dentate_delta, cereb_median_delta) + 0.05 * support_n
        rows.append(
            {
                "contrast_id": contrast_id,
                "contrast_label": drow["contrast_label"],
                "source_gene_symbol": drow["source_gene_symbol"],
                "canonical_mouse_symbol": canonical_mouse_symbol,
                "human_symbol": drow["human_symbol"],
                "mouse_symbol": drow["mouse_symbol"],
                "mgi_one_to_one_human_mouse": bool(drow["mgi_one_to_one_human_mouse"]),
                "dentate_granule_vs_pyramidal_log1p_delta": dentate_delta,
                "dentate_granule_detection_fraction": float(drow["granule_detection_fraction"]),
                "dentate_pyramidal_detection_fraction": float(drow["comparator_detection_fraction"]),
                "cerebellar_granule_vs_purkinje_median_log1p_delta": cereb_median_delta,
                "cerebellar_granule_vs_purkinje_min_log1p_delta": float(np.min(cereb_delta)),
                "cerebellar_positive_sample_count": support_n,
                "cerebellar_samples_positive": ",".join(cereb_support_samples),
                "cerebellar_granule_median_detection_fraction": float(np.median(cereb_detection)),
                "cerebellar_purkinje_median_detection_fraction": float(np.median(cereb["comparator_detection_fraction"])),
                "shared_positive_named_comparator": shared_positive,
                "strong_shared_named_comparator": strong_shared,
                "specificity_score": float(specificity_score),
                "combined_log1p_delta": float(dentate_delta + cereb_median_delta),
            }
        )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary
    return summary.sort_values(
        [
            "contrast_id",
            "strong_shared_named_comparator",
            "shared_positive_named_comparator",
            "specificity_score",
            "combined_log1p_delta",
        ],
        ascending=[True, False, False, False, False],
    )


def build_top(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return summary
    top = summary.loc[summary["shared_positive_named_comparator"]].copy()
    return top.sort_values(
        ["contrast_id", "strong_shared_named_comparator", "specificity_score", "combined_log1p_delta"],
        ascending=[True, False, False, False],
    )


def plot_top(top: pd.DataFrame) -> None:
    if top.empty:
        return
    panels = []
    for contrast_id in CONTRASTS:
        sub = top.loc[top["contrast_id"].eq(contrast_id)].head(20).copy()
        if not sub.empty:
            panels.append((contrast_id, sub))
    if not panels:
        return
    fig, axes = plt.subplots(len(panels), 1, figsize=(9.5, 4.2 * len(panels)), squeeze=False)
    for panel_i, (ax, (contrast_id, sub)) in enumerate(zip(axes[:, 0], panels)):
        labels = sub["human_symbol"].replace("", np.nan).fillna(sub["canonical_mouse_symbol"]).tolist()
        y = np.arange(len(sub))
        ax.barh(y - 0.18, sub["dentate_granule_vs_pyramidal_log1p_delta"], height=0.34, label="Dentate > pyramidal", color="#4c78a8")
        ax.barh(y + 0.18, sub["cerebellar_granule_vs_purkinje_median_log1p_delta"], height=0.34, label="Cerebellar > Purkinje", color="#f58518")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.set_xlabel("Log1p expression delta")
        ax.set_title(
            f"{chr(ord('a') + panel_i)}. {CONTRASTS[contrast_id]['label']}",
            loc="left",
            fontsize=11,
            fontweight="bold",
        )
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    SUPP_FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUPP_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, top: pd.DataFrame) -> None:
    lines = [
        "# Granule-Special Gene Named-Comparator Screen",
        "",
        "This screen searches for genes that are higher in dentate granule-lineage groups than pyramidal comparators in `GSE104323` and higher in cerebellar granule-lineage groups than Purkinje cells in `GSE122357`.",
        "",
        "A gene is called shared-positive when dentate delta is positive, median cerebellar delta is positive, at least two cerebellar samples are positive, and median granule detection is at least 5% in both branches. A strong shared call additionally requires deltas >= 0.05 and all three cerebellar samples positive.",
        "",
        "## Screen Counts",
        "",
    ]
    for contrast_id in CONTRASTS:
        sub = summary.loc[summary["contrast_id"].eq(contrast_id)]
        top_sub = top.loc[top["contrast_id"].eq(contrast_id)]
        strong = int(sub["strong_shared_named_comparator"].sum()) if not sub.empty else 0
        lines.append(
            f"- {CONTRASTS[contrast_id]['label']}: tested {len(sub):,} genes; "
            f"shared-positive {len(top_sub):,}; strong shared {strong:,}."
        )

    lines.extend(["", "## Top Shared-Positive Candidates", ""])
    for contrast_id in CONTRASTS:
        lines.append(f"### {CONTRASTS[contrast_id]['label']}")
        sub = top.loc[top["contrast_id"].eq(contrast_id)].head(25)
        if sub.empty:
            lines.append("")
            lines.append("No shared-positive candidates passed the current thresholds.")
            lines.append("")
            continue
        names = [
            str(row["human_symbol"]) if str(row["human_symbol"]).strip() else str(row["canonical_mouse_symbol"])
            for _, row in sub.iterrows()
        ]
        lines.append("")
        lines.append(", ".join(f"`{name}`" for name in names))
        lines.append("")

    nfia = summary.loc[summary["canonical_mouse_symbol"].eq("NFIA")]
    if not nfia.empty:
        lines.extend(["## NFIA Position", ""])
        for _, row in nfia.iterrows():
            rank = (
                top.loc[top["contrast_id"].eq(row["contrast_id"]), "canonical_mouse_symbol"]
                .reset_index(drop=True)
                .eq("NFIA")
            )
            rank_text = "not shared-positive"
            if rank.any():
                rank_text = f"rank {int(np.flatnonzero(rank.to_numpy())[0] + 1)} among shared-positive genes"
            lines.append(
                f"- {row['contrast_label']}: {rank_text}; dentate delta {row['dentate_granule_vs_pyramidal_log1p_delta']:.3f}, "
                f"cerebellar median delta {row['cerebellar_granule_vs_purkinje_median_log1p_delta']:.3f}."
            )
        lines.append("")

    lines.extend(
        [
            "## Interpretation",
            "",
            "- The screen identifies local named-comparator granule-enriched genes, not universal granule-cell-specific markers across the whole brain.",
            "- The granule-lineage contrast is broader and includes precursor/neuroblast biology.",
            "- The postmitotic/granule-cell-state contrast is more conservative for mature or differentiating granule-cell identity.",
            "- Candidates should be cross-checked against broader cell atlases before being described as granule-cell-specific.",
            "",
            "## Outputs",
            "",
            f"- Unit table: `{OUT_UNITS.relative_to(ROOT)}`",
            f"- Summary table: `{OUT_SUMMARY.relative_to(ROOT)}`",
            f"- Top candidates: `{OUT_TOP.relative_to(ROOT)}`",
            f"- Plot: `{OUT_PLOT.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    ortholog = load_ortholog_map()
    records: list[dict[str, object]] = []
    records.extend(process_gse104323(ortholog))
    records.extend(process_gse122357(ortholog))
    units = collapse_duplicate_units(pd.DataFrame(records))
    summary = build_summary(units)
    top = build_top(summary)
    units.to_csv(OUT_UNITS, sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    top.to_csv(OUT_TOP, sep="\t", index=False)
    plot_top(top)
    write_report(summary, top)
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    for contrast_id in CONTRASTS:
        sub = summary.loc[summary["contrast_id"].eq(contrast_id)]
        top_sub = top.loc[top["contrast_id"].eq(contrast_id)]
        print(
            f"{CONTRASTS[contrast_id]['label']}: tested {len(sub):,}, "
            f"shared-positive {len(top_sub):,}, strong {int(sub['strong_shared_named_comparator'].sum())}"
        )
        show = top_sub.head(15)
        if not show.empty:
            cols = [
                "human_symbol",
                "canonical_mouse_symbol",
                "dentate_granule_vs_pyramidal_log1p_delta",
                "cerebellar_granule_vs_purkinje_median_log1p_delta",
                "cerebellar_positive_sample_count",
                "specificity_score",
            ]
            print(show[cols].to_string(index=False))


if __name__ == "__main__":
    main()
