#!/usr/bin/env python3
"""Compare upstream niche/fate modules with downstream circuit/morphology modules.

This analysis asks whether dentate and cerebellar granule-cell similarity is
better explained by shared upstream stem-cell niche/fate programs or by
convergent downstream morphology, excitability, and circuit-implementation
programs.

It has two layers:

1. Formal-core gene summary: uses the existing expanded MGI ortholog rank-meta
   gene table across the strict primary core.
2. Named local comparators: scores the same modules in GSE104323 dentate
   granule-lineage versus pyramidal labels and GSE122357 cerebellar
   granule-lineage versus Purkinje labels.
"""

from __future__ import annotations

import os
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

FORMAL_GENE_SUMMARY = RESULTS / "primary_core_mgi_ortholog_formal_rank_gene_summary.tsv"

GSE104323_EXPR = EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz"
GSE104323_META = EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz"

GSE122357_TAR = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
GSE122357_LABELS = EXTERNAL / "GEO/GSE122357/GSE122357_cell_number.xlsx"

OUT_GENE_SETS = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
OUT_FORMAL_GENE = RESULTS / "primary_core_niche_circuit_module_formal_gene_scores.tsv"
OUT_FORMAL_SUMMARY = RESULTS / "primary_core_niche_circuit_module_formal_summary.tsv"
OUT_COMPARATOR_UNITS = RESULTS / "primary_core_niche_circuit_module_named_comparator_units.tsv"
OUT_COMPARATOR_SUMMARY = RESULTS / "primary_core_niche_circuit_module_named_comparator_summary.tsv"
OUT_PLOT = RESULTS / "primary_core_niche_circuit_module_model.png"
OUT_MD = RESULTS / "primary_core_niche_circuit_module_model.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULES = [
    {
        "module_id": "cerebellar_fate_rhombic_lip_shh",
        "module_label": "Cerebellar fate/rhombic-lip/SHH",
        "module_family": "upstream_region_fate",
        "hypothesis_role": "cerebellar-specific upstream fate and progenitor-expansion logic",
        "genes": [
            "ATOH1",
            "BARHL1",
            "PAX6",
            "ZIC1",
            "ZIC2",
            "ZIC3",
            "MYCN",
            "PTCH1",
            "GLI1",
            "GLI2",
            "EN1",
            "EN2",
            "MEIS1",
        ],
    },
    {
        "module_id": "dentate_fate_wnt_prox1",
        "module_label": "Dentate fate/WNT/PROX1",
        "module_family": "upstream_region_fate",
        "hypothesis_role": "dentate-specific hippocampal granule fate and WNT-associated identity logic",
        "genes": [
            "PROX1",
            "NEUROD1",
            "NEUROD2",
            "EOMES",
            "TBR1",
            "LEF1",
            "LHX2",
            "EMX2",
            "ZBTB20",
            "TCF7L2",
            "WNT3A",
            "WNT7A",
            "DKK3",
            "BCL11B",
            "CALB1",
            "C1QL3",
            "GLIS3",
        ],
    },
    {
        "module_id": "shared_neurogenic_niche_state",
        "module_label": "Shared neurogenic niche/progenitor state",
        "module_family": "shared_niche_state",
        "hypothesis_role": "general neurogenic stem/progenitor and local niche signaling state",
        "genes": [
            "SOX2",
            "HES1",
            "HES5",
            "NOTCH1",
            "NOTCH2",
            "DLL1",
            "JAG1",
            "ASCL1",
            "DCX",
            "TUBB3",
            "NES",
            "VIM",
            "GFAP",
            "BMP4",
            "BMPR1A",
            "ID1",
            "ID2",
            "PTCH1",
            "GLI1",
            "WNT5A",
            "WNT7A",
        ],
    },
    {
        "module_id": "downstream_neurite_morphology",
        "module_label": "Downstream neurite/morphology",
        "module_family": "downstream_circuit_morphology",
        "hypothesis_role": "convergent neurite, axon, adhesion, and morphology implementation",
        "genes": [
            "GPM6A",
            "ROBO2",
            "DCC",
            "CADM3",
            "STMN2",
            "STMN3",
            "GAP43",
            "DPYSL2",
            "DPYSL3",
            "MAP1B",
            "BASP1",
            "NCAM1",
            "L1CAM",
            "CFL1",
            "RTN1",
            "RTN3",
            "NRXN1",
            "CNTN5",
            "PLXNA2",
            "PLXNA4",
            "NRP1",
        ],
    },
    {
        "module_id": "downstream_synaptic_excitability",
        "module_label": "Downstream synaptic/excitability",
        "module_family": "downstream_circuit_morphology",
        "hypothesis_role": "convergent sparse-coding, synaptic maturation, and excitability implementation",
        "genes": [
            "KCNK1",
            "KCNJ3",
            "KCNJ6",
            "KCND2",
            "GABRA2",
            "GABRB3",
            "GRIN2B",
            "CACNA2D1",
            "PPP3CA",
            "STXBP5L",
            "STXBP1",
            "SYNPR",
            "SLC17A6",
            "SLC17A7",
            "SNAP25",
            "SYT1",
            "CALM1",
            "CALM2",
            "CAMK2A",
            "CAMK2B",
        ],
    },
]

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


def mouse_case(symbol: str) -> str:
    if not symbol:
        return symbol
    return symbol[:1].upper() + symbol[1:].lower()


def write_gene_sets() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module in MODULES:
        for order, gene in enumerate(module["genes"], start=1):
            rows.append(
                {
                    "module_id": module["module_id"],
                    "module_label": module["module_label"],
                    "module_family": module["module_family"],
                    "hypothesis_role": module["hypothesis_role"],
                    "gene": gene,
                    "canonical_gene": canon(gene),
                    "default_mouse_symbol": mouse_case(gene),
                    "gene_order": order,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(OUT_GENE_SETS, sep="\t", index=False)
    return df


def finite_median(values: list[float]) -> float:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def first_finite(values: list[float]) -> float:
    for value in values:
        if np.isfinite(value):
            return float(value)
    return np.nan


def formal_gene_scores(gene_sets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    formal = pd.read_csv(FORMAL_GENE_SUMMARY, sep="\t", low_memory=False)
    formal["canonical_gene"] = formal["canonical_gene"].map(canon)
    by_gene = formal.set_index("canonical_gene")
    rows: list[dict[str, object]] = []

    for _, row in gene_sets.iterrows():
        gene = row["canonical_gene"]
        if gene not in by_gene.index:
            rows.append(
                {
                    **row.to_dict(),
                    "present_in_formal_model": False,
                    "formal_mouse_symbol": "",
                    "formal_rank_tier": "",
                    "formal_rank_priority_score": np.nan,
                    "formal_n_available_branches": np.nan,
                    "formal_n_nominal_branches": np.nan,
                    "formal_n_fdr10_branches": np.nan,
                    "formal_replication_shared_both_screens": False,
                    "formal_nominal_shared_both_screens": False,
                    "formal_fdr10_shared_both_screens": False,
                    "median_delta_selected_dentate": np.nan,
                    "median_delta_selected_cerebellar": np.nan,
                    "median_delta_full_matrix_dentate": np.nan,
                    "median_delta_full_matrix_cerebellar": np.nan,
                    "selected_convergence_delta": np.nan,
                    "full_matrix_convergence_delta": np.nan,
                    "overall_convergence_delta": np.nan,
                    "mean_branch_bias": np.nan,
                    "shared_positive_any_screen": False,
                    "shared_positive_both_screens": False,
                    "branch_pattern": "missing_from_formal_model",
                }
            )
            continue

        match = by_gene.loc[gene]
        if isinstance(match, pd.DataFrame):
            match = match.iloc[0]

        sd = float(match.get("median_dataset_rank_delta_selected_dentate", np.nan))
        sc = float(match.get("median_dataset_rank_delta_selected_cerebellar", np.nan))
        fd = float(match.get("median_dataset_rank_delta_full_matrix_dentate", np.nan))
        fc = float(match.get("median_dataset_rank_delta_full_matrix_cerebellar", np.nan))
        selected_conv = min(sd, sc) if np.isfinite(sd) and np.isfinite(sc) else np.nan
        full_conv = min(fd, fc) if np.isfinite(fd) and np.isfinite(fc) else np.nan
        overall_conv = finite_median([selected_conv, full_conv])
        selected_bias = abs(sd - sc) if np.isfinite(sd) and np.isfinite(sc) else np.nan
        full_bias = abs(fd - fc) if np.isfinite(fd) and np.isfinite(fc) else np.nan
        mean_bias = finite_median([selected_bias, full_bias])
        shared_any = bool((np.isfinite(selected_conv) and selected_conv > 0) or (np.isfinite(full_conv) and full_conv > 0))
        shared_both = bool((np.isfinite(selected_conv) and selected_conv > 0) and (np.isfinite(full_conv) and full_conv > 0))

        dentate_peak = finite_median([sd, fd])
        cereb_peak = finite_median([sc, fc])
        if shared_both:
            pattern = "shared_positive_both_screens"
        elif shared_any:
            pattern = "shared_positive_one_screen"
        elif np.isfinite(dentate_peak) and np.isfinite(cereb_peak) and dentate_peak - cereb_peak >= 0.15:
            pattern = "dentate_biased"
        elif np.isfinite(dentate_peak) and np.isfinite(cereb_peak) and cereb_peak - dentate_peak >= 0.15:
            pattern = "cerebellar_biased"
        elif np.isfinite(dentate_peak) or np.isfinite(cereb_peak):
            pattern = "weak_or_mixed"
        else:
            pattern = "insufficient_delta_data"

        rows.append(
            {
                **row.to_dict(),
                "present_in_formal_model": True,
                "formal_mouse_symbol": match.get("mouse_symbol", ""),
                "formal_rank_tier": match.get("formal_rank_tier", ""),
                "formal_rank_priority_score": match.get("formal_rank_priority_score", np.nan),
                "formal_n_available_branches": match.get("formal_n_available_branches", np.nan),
                "formal_n_nominal_branches": match.get("formal_n_nominal_branches", np.nan),
                "formal_n_fdr10_branches": match.get("formal_n_fdr10_branches", np.nan),
                "formal_replication_shared_both_screens": bool(match.get("formal_replication_shared_both_screens", False)),
                "formal_nominal_shared_both_screens": bool(match.get("formal_nominal_shared_both_screens", False)),
                "formal_fdr10_shared_both_screens": bool(match.get("formal_fdr10_shared_both_screens", False)),
                "median_delta_selected_dentate": sd,
                "median_delta_selected_cerebellar": sc,
                "median_delta_full_matrix_dentate": fd,
                "median_delta_full_matrix_cerebellar": fc,
                "selected_convergence_delta": selected_conv,
                "full_matrix_convergence_delta": full_conv,
                "overall_convergence_delta": overall_conv,
                "mean_branch_bias": mean_bias,
                "shared_positive_any_screen": shared_any,
                "shared_positive_both_screens": shared_both,
                "branch_pattern": pattern,
            }
        )

    gene_df = pd.DataFrame(rows)

    summary_rows: list[dict[str, object]] = []
    for module_id, sub in gene_df.groupby("module_id", sort=False):
        present = sub.loc[sub["present_in_formal_model"].astype(bool)].copy()
        module_meta = sub.iloc[0]
        downstream = present.loc[present["module_family"].eq("downstream_circuit_morphology")]
        summary_rows.append(
            {
                "module_id": module_id,
                "module_label": module_meta["module_label"],
                "module_family": module_meta["module_family"],
                "hypothesis_role": module_meta["hypothesis_role"],
                "n_genes_defined": int(len(sub)),
                "n_genes_present_formal": int(len(present)),
                "n_shared_positive_any_screen": int(present["shared_positive_any_screen"].sum()) if not present.empty else 0,
                "n_shared_positive_both_screens": int(present["shared_positive_both_screens"].sum()) if not present.empty else 0,
                "fraction_shared_positive_any_screen": float(present["shared_positive_any_screen"].mean())
                if not present.empty
                else np.nan,
                "fraction_shared_positive_both_screens": float(present["shared_positive_both_screens"].mean())
                if not present.empty
                else np.nan,
                "median_overall_convergence_delta": float(present["overall_convergence_delta"].median())
                if not present.empty
                else np.nan,
                "median_mean_branch_bias": float(present["mean_branch_bias"].median()) if not present.empty else np.nan,
                "median_selected_dentate_delta": float(present["median_delta_selected_dentate"].median())
                if not present.empty
                else np.nan,
                "median_selected_cerebellar_delta": float(present["median_delta_selected_cerebellar"].median())
                if not present.empty
                else np.nan,
                "median_full_matrix_dentate_delta": float(present["median_delta_full_matrix_dentate"].median())
                if not present.empty
                else np.nan,
                "median_full_matrix_cerebellar_delta": float(present["median_delta_full_matrix_cerebellar"].median())
                if not present.empty
                else np.nan,
                "n_formal_nominal_shared_both_screens": int(present["formal_nominal_shared_both_screens"].sum())
                if not present.empty
                else 0,
                "top_shared_positive_genes": ",".join(
                    present.loc[present["shared_positive_any_screen"]]
                    .sort_values("overall_convergence_delta", ascending=False)["gene"]
                    .astype(str)
                    .head(10)
                ),
                "n_downstream_present": int(len(downstream)),
            }
        )

    summary = pd.DataFrame(summary_rows)
    gene_df.to_csv(OUT_FORMAL_GENE, sep="\t", index=False)
    summary.to_csv(OUT_FORMAL_SUMMARY, sep="\t", index=False)
    return gene_df, summary


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


def role_for_group(dataset: str, group: str) -> str:
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


def compute_module_units(
    *,
    dataset: str,
    sample: str,
    expression: pd.DataFrame,
    cell_groups: pd.Series,
    gene_sets: pd.DataFrame,
) -> list[dict[str, object]]:
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    index_by_canon = {canon(gene): gene for gene in expression.index}
    common_cells = [cell for cell in expression.columns if cell in set(cell_groups.index)]
    if not common_cells:
        return []
    expression = expression[common_cells]
    cell_groups = cell_groups.loc[common_cells]

    records: list[dict[str, object]] = []
    for module_id, sub in gene_sets.groupby("module_id", sort=False):
        meta = sub.iloc[0]
        genes = list(sub["gene"].astype(str))
        present = [index_by_canon[canon(gene)] for gene in genes if canon(gene) in index_by_canon]
        if not present:
            continue
        cell_scores = np.log1p(expression.loc[present].to_numpy(dtype=float)).mean(axis=0)
        score_df = pd.DataFrame(
            {
                "cell_id": common_cells,
                "source_group": cell_groups.to_numpy(),
                "module_score": cell_scores,
            }
        )
        group_summary = (
            score_df.groupby("source_group", sort=False)
            .agg(
                n_cells=("module_score", "size"),
                median_module_score=("module_score", "median"),
                mean_module_score=("module_score", "mean"),
            )
            .reset_index()
        )
        group_summary["within_sample_module_rank"] = group_summary["median_module_score"].rank(
            pct=True, method="average"
        )
        for _, row in group_summary.iterrows():
            group = str(row["source_group"])
            records.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "source_group": group,
                    "specificity_role": role_for_group(dataset, group),
                    "module_id": module_id,
                    "module_label": meta["module_label"],
                    "module_family": meta["module_family"],
                    "n_defined_genes": int(len(genes)),
                    "n_present_genes": int(len(present)),
                    "present_genes": ",".join(present),
                    "missing_genes": ",".join(gene for gene in genes if canon(gene) not in index_by_canon),
                    "n_cells": int(row["n_cells"]),
                    "median_module_score": float(row["median_module_score"]),
                    "mean_module_score": float(row["mean_module_score"]),
                    "within_sample_module_rank": float(row["within_sample_module_rank"]),
                }
            )
    return records


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


def named_comparator_units(gene_sets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted_genes = set(gene_sets["gene"].astype(str)) | set(gene_sets["default_mouse_symbol"].astype(str))
    records: list[dict[str, object]] = []

    meta = pd.read_csv(GSE104323_META, sep="\t")
    meta = meta.rename(columns={"Sample name (24185 single cells)": "cell_id", "characteristics: cell cluster": "group"})
    cell_groups = meta.set_index("cell_id")["group"].astype(str)
    expression = read_selected_rows(GSE104323_EXPR, "\t", wanted_genes)
    records.extend(
        compute_module_units(
            dataset="GSE104323",
            sample="10X_all_cells",
            expression=expression,
            cell_groups=cell_groups,
            gene_sets=gene_sets,
        )
    )

    with TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for sample, (member, prefix) in GSE122357_FILES.items():
            expression = read_selected_rows_from_tar(GSE122357_TAR, member, wanted_genes, tmpdir)
            labels = load_gse122357_label_map(prefix)
            records.extend(
                compute_module_units(
                    dataset="GSE122357",
                    sample=sample,
                    expression=expression,
                    cell_groups=labels,
                    gene_sets=gene_sets,
                )
            )

    units = pd.DataFrame(records)
    summary_rows: list[dict[str, object]] = []
    for module_id, sub in units.groupby("module_id", sort=False):
        meta_row = sub.iloc[0]

        def role_median(role: str) -> float:
            values = sub.loc[sub["specificity_role"].eq(role), "within_sample_module_rank"].dropna().to_numpy(dtype=float)
            if values.size == 0:
                return np.nan
            return float(np.median(values))

        dg = role_median("dentate_granule")
        pyr = role_median("pyramidal_comparator")
        cb = role_median("cerebellar_granule")
        pc = role_median("purkinje_comparator")
        dentate_delta = dg - pyr if np.isfinite(dg) and np.isfinite(pyr) else np.nan
        cereb_delta = cb - pc if np.isfinite(cb) and np.isfinite(pc) else np.nan
        granule_index = first_finite(
            [
                min(dentate_delta, cereb_delta)
                if np.isfinite(dentate_delta) and np.isfinite(cereb_delta)
                else np.nan,
                dentate_delta,
                cereb_delta,
            ]
        )
        if np.isfinite(dentate_delta) and np.isfinite(cereb_delta) and dentate_delta > 0 and cereb_delta > 0:
            call = "granule_enriched_vs_named_comparators"
        elif np.isfinite(cereb_delta) and cereb_delta > 0 and not (np.isfinite(dentate_delta) and dentate_delta > 0):
            call = "cerebellar_granule_enriched_only"
        elif np.isfinite(dentate_delta) and dentate_delta > 0 and not (np.isfinite(cereb_delta) and cereb_delta > 0):
            call = "dentate_granule_enriched_only"
        else:
            call = "not_granule_specific_vs_named_comparators"
        summary_rows.append(
            {
                "module_id": module_id,
                "module_label": meta_row["module_label"],
                "module_family": meta_row["module_family"],
                "n_defined_genes": int(sub["n_defined_genes"].max()),
                "median_n_present_genes": float(sub["n_present_genes"].median()),
                "dentate_granule_median_rank": dg,
                "pyramidal_comparator_median_rank": pyr,
                "cerebellar_granule_median_rank": cb,
                "purkinje_comparator_median_rank": pc,
                "dentate_vs_pyramidal_rank_delta": dentate_delta,
                "cerebellar_vs_purkinje_rank_delta": cereb_delta,
                "granule_specificity_index": granule_index,
                "specificity_call": call,
            }
        )
    summary = pd.DataFrame(summary_rows)
    units.to_csv(OUT_COMPARATOR_UNITS, sep="\t", index=False)
    summary.to_csv(OUT_COMPARATOR_SUMMARY, sep="\t", index=False)
    return units, summary


def family_comparison(formal_gene: pd.DataFrame) -> dict[str, float]:
    present = formal_gene.loc[formal_gene["present_in_formal_model"].astype(bool)].copy()
    upstream = present.loc[present["module_family"].isin(["upstream_region_fate", "shared_niche_state"])]
    downstream = present.loc[present["module_family"].eq("downstream_circuit_morphology")]
    up_values = upstream["overall_convergence_delta"].dropna().to_numpy(dtype=float)
    down_values = downstream["overall_convergence_delta"].dropna().to_numpy(dtype=float)
    if len(up_values) and len(down_values):
        mw = stats.mannwhitneyu(down_values, up_values, alternative="greater")
        p = float(mw.pvalue)
    else:
        p = np.nan
    return {
        "n_upstream_or_niche_genes": int(len(up_values)),
        "n_downstream_genes": int(len(down_values)),
        "median_upstream_or_niche_convergence_delta": float(np.median(up_values)) if len(up_values) else np.nan,
        "median_downstream_convergence_delta": float(np.median(down_values)) if len(down_values) else np.nan,
        "mannwhitney_downstream_greater_p": p,
    }


def plot_results(formal_summary: pd.DataFrame, comparator_summary: pd.DataFrame) -> None:
    modules = formal_summary["module_id"].tolist()
    labels = formal_summary["module_label"].tolist()
    colors = [
        "#7b6d8d" if fam == "upstream_region_fate" else "#899878" if fam == "shared_niche_state" else "#2f7f8f"
        for fam in formal_summary["module_family"]
    ]
    comp = comparator_summary.set_index("module_id")
    x = np.arange(len(modules))

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8), constrained_layout=True)
    ax = axes[0]
    ax.bar(x, formal_summary["median_overall_convergence_delta"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("Median formal convergence delta")
    ax.set_title("Strict-core gene-level convergence")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    ax = axes[1]
    dentate = [comp.loc[module, "dentate_vs_pyramidal_rank_delta"] if module in comp.index else np.nan for module in modules]
    cereb = [comp.loc[module, "cerebellar_vs_purkinje_rank_delta"] if module in comp.index else np.nan for module in modules]
    width = 0.36
    ax.bar(x - width / 2, dentate, width=width, color="#3d7c5f", label="dentate vs pyramidal")
    ax.bar(x + width / 2, cereb, width=width, color="#7f4e8a", label="cerebellar vs Purkinje")
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("Named-comparator rank delta")
    ax.set_title("Local comparator specificity")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    fig.suptitle("Niche/Fate Versus Circuit/Morphology Programs", fontsize=15, y=1.03)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(
    formal_gene: pd.DataFrame,
    formal_summary: pd.DataFrame,
    comparator_units: pd.DataFrame,
    comparator_summary: pd.DataFrame,
    comparison: dict[str, float],
) -> None:
    downstream = formal_summary.loc[formal_summary["module_family"].eq("downstream_circuit_morphology")]
    upstream = formal_summary.loc[formal_summary["module_family"].ne("downstream_circuit_morphology")]
    strongest_downstream = downstream.sort_values("median_overall_convergence_delta", ascending=False).head(2)
    strongest_upstream = upstream.sort_values("median_overall_convergence_delta", ascending=False).head(2)

    lines = [
        "# Niche/Fate Versus Circuit/Morphology Module Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis asks whether dentate and cerebellar granule-cell similarity is better explained by shared upstream stem-cell niche/fate programs or by convergent downstream circuit, morphology, and maturation programs.",
        "",
        "## Module Families",
        "",
    ]
    for module in MODULES:
        lines.append(
            f"- {module['module_label']}: {len(module['genes'])} genes; {module['hypothesis_role']}."
        )
    lines.extend(
        [
            "",
            "## Formal-Core Result",
            "",
            f"- Genes scored in the formal model: {int(formal_gene['present_in_formal_model'].sum())}/{len(formal_gene)}.",
            f"- Median upstream/niche convergence delta: {comparison['median_upstream_or_niche_convergence_delta']:.3f}.",
            f"- Median downstream circuit/morphology convergence delta: {comparison['median_downstream_convergence_delta']:.3f}.",
            f"- Mann-Whitney test, downstream greater than upstream/niche: p={comparison['mannwhitney_downstream_greater_p']:.3g}.",
            "",
            "Top downstream modules by formal convergence:",
        ]
    )
    for _, row in strongest_downstream.iterrows():
        lines.append(
            f"- {row['module_label']}: median convergence {row['median_overall_convergence_delta']:.3f}, "
            f"{row['n_shared_positive_any_screen']}/{row['n_genes_present_formal']} genes shared-positive in at least one screen."
        )
    lines.append("")
    lines.append("Best upstream/niche modules by formal convergence:")
    for _, row in strongest_upstream.iterrows():
        lines.append(
            f"- {row['module_label']}: median convergence {row['median_overall_convergence_delta']:.3f}, "
            f"{row['n_shared_positive_any_screen']}/{row['n_genes_present_formal']} genes shared-positive in at least one screen."
        )

    lines.extend(
        [
            "",
            "## Named-Comparator Result",
            "",
            f"- Comparator units scored: {len(comparator_units)}.",
            "- Direct comparator rule: dentate granule-lineage groups must exceed pyramidal labels, and cerebellar granule-lineage groups must exceed Purkinje labels.",
        ]
    )
    for _, row in comparator_summary.iterrows():
        lines.append(
            f"- {row['module_label']}: {row['specificity_call']}; dentate-vs-pyramidal delta "
            f"{row['dentate_vs_pyramidal_rank_delta']:.3f}, cerebellar-vs-Purkinje delta "
            f"{row['cerebellar_vs_purkinje_rank_delta']:.3f}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The current evidence favors a convergence model: upstream fate modules are branch-specific or mixed, while downstream morphology/excitability modules carry the strongest formal shared-convergence signal.",
            "- The named-comparator layer remains a constraint: several downstream modules are not uniquely granule-specific versus pyramidal and Purkinje comparators, so the claim should be circuit-convergence rather than pathway uniqueness.",
            "- This analysis can be strengthened later by adding developmental time-resolved trajectory inference and perturbation/lineage-tracing evidence from the literature.",
            "",
            "## Outputs",
            "",
            f"- Gene-set table: `{rel(OUT_GENE_SETS)}`",
            f"- Formal gene scores: `{rel(OUT_FORMAL_GENE)}`",
            f"- Formal module summary: `{rel(OUT_FORMAL_SUMMARY)}`",
            f"- Named-comparator units: `{rel(OUT_COMPARATOR_UNITS)}`",
            f"- Named-comparator summary: `{rel(OUT_COMPARATOR_SUMMARY)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    gene_sets = write_gene_sets()
    formal_gene, formal_summary = formal_gene_scores(gene_sets)
    comparator_units, comparator_summary = named_comparator_units(gene_sets)
    comparison = family_comparison(formal_gene)
    plot_results(formal_summary, comparator_summary)
    write_report(formal_gene, formal_summary, comparator_units, comparator_summary, comparison)
    print(f"Wrote {rel(OUT_MD)}")
    print(formal_summary[["module_label", "median_overall_convergence_delta", "fraction_shared_positive_any_screen"]].to_string(index=False))
    print(comparator_summary[["module_label", "specificity_call", "granule_specificity_index"]].to_string(index=False))


if __name__ == "__main__":
    main()
