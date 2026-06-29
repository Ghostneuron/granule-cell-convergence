#!/usr/bin/env python3
"""Aim 2 niche/pathway signaling audit for granule-cell convergence.

This analysis finishes the first computable version of Specific Aim 2:

    Test whether niche-signaling programs recapitulate the historical
    TGF-beta2/BDNF/SMAD-MAPK mechanism and whether cerebellar and dentate
    granule-cell contexts differ in differentiation/stop versus persistent
    neurogenic/permissive signaling.

The available primary-core matrices are broad-class pseudobulk summaries, not
spatially resolved ligand-sender/receiver assays. Therefore the output is a
conservative pathway-readiness and module-rank audit, not a direct proof of
cell-cell communication.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

FULL_EXPR = RESULTS / "primary_core_mgi_ortholog_full_matrix_expression.tsv.gz"
SELECTED_EXPR = RESULTS / "primary_core_expanded_gene_pseudobulk_expression.tsv.gz"

OUT_GENE_SETS = RESULTS / "primary_core_aim2_pathway_gene_sets.tsv"
OUT_UNITS = RESULTS / "primary_core_aim2_pathway_units.tsv.gz"
OUT_CONTRASTS = RESULTS / "primary_core_aim2_pathway_contrasts.tsv"
OUT_SUMMARY = RESULTS / "primary_core_aim2_pathway_summary.tsv"
OUT_COVERAGE = RESULTS / "primary_core_aim2_pathway_coverage.tsv"
OUT_SIGNATURE_UNITS = RESULTS / "primary_core_aim2_signature_units.tsv.gz"
OUT_SIGNATURE_CONTRASTS = RESULTS / "primary_core_aim2_signature_contrasts.tsv"
OUT_SIGNATURE_SUMMARY = RESULTS / "primary_core_aim2_signature_summary.tsv"
OUT_LR_PAIRS = RESULTS / "primary_core_aim2_ligand_receptor_pairs.tsv"
OUT_LR_UNITS = RESULTS / "primary_core_aim2_ligand_receptor_units.tsv.gz"
OUT_LR_CONTRASTS = RESULTS / "primary_core_aim2_ligand_receptor_contrasts.tsv"
OUT_LR_SUMMARY = RESULTS / "primary_core_aim2_ligand_receptor_summary.tsv"
OUT_PLOT = RESULTS / "primary_core_aim2_niche_pathway_model.png"
OUT_MD = RESULTS / "primary_core_aim2_niche_pathway_model.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPRESSION_LAYERS = [
    {
        "expression_layer": "full_mgi_ortholog_matrix",
        "path": FULL_EXPR,
    },
    {
        "expression_layer": "selected_feature_matrix",
        "path": SELECTED_EXPR,
    },
]


PATHWAY_MODULES = [
    {
        "pathway_id": "tgf_beta_smad",
        "pathway_label": "TGF-beta/SMAD",
        "pathway_family": "differentiation_stop",
        "hypothesis_role": "historical TGF-beta2-associated differentiation/stop signaling",
        "genes": {
            "TGFB2": "ligand",
            "TGFB1": "ligand",
            "TGFB3": "ligand",
            "TGFBR1": "receptor",
            "TGFBR2": "receptor",
            "TGFBR3": "receptor",
            "SMAD2": "effector",
            "SMAD3": "effector",
            "SMAD4": "effector",
            "SMAD7": "feedback",
            "TGFBI": "response",
            "SERPINE1": "response",
            "ID1": "response",
            "ID2": "response",
            "ID3": "response",
        },
    },
    {
        "pathway_id": "bdnf_trkb_mapk",
        "pathway_label": "BDNF/TrkB/MAPK",
        "pathway_family": "differentiation_stop",
        "hypothesis_role": "historical BDNF/TrkB plus MAPK maturation pathway",
        "genes": {
            "BDNF": "ligand",
            "NTRK2": "receptor",
            "SORT1": "coreceptor",
            "NGFR": "coreceptor",
            "SHC1": "adaptor",
            "GRB2": "adaptor",
            "SOS1": "adaptor",
            "RAF1": "effector",
            "MAP2K1": "effector",
            "MAPK1": "effector",
            "MAPK3": "effector",
            "MAPK14": "effector",
            "CREB1": "response",
            "ELK1": "response",
            "EGR1": "response",
            "FOS": "response",
            "JUN": "response",
            "DUSP6": "feedback",
        },
    },
    {
        "pathway_id": "bmp_smad",
        "pathway_label": "BMP/SMAD",
        "pathway_family": "differentiation_stop",
        "hypothesis_role": "BMP-associated differentiation and neurogenic brake signaling",
        "genes": {
            "BMP2": "ligand",
            "BMP4": "ligand",
            "BMP7": "ligand",
            "BMPR1A": "receptor",
            "BMPR1B": "receptor",
            "BMPR2": "receptor",
            "ACVR1": "receptor",
            "SMAD1": "effector",
            "SMAD5": "effector",
            "SMAD9": "effector",
            "ID1": "response",
            "ID2": "response",
            "ID3": "response",
            "ID4": "response",
            "NOG": "antagonist",
        },
    },
    {
        "pathway_id": "reelin_migration_stop",
        "pathway_label": "Reelin/migration-stop",
        "pathway_family": "differentiation_stop",
        "hypothesis_role": "lamination, migration-stop, and maturation positioning cues",
        "genes": {
            "RELN": "ligand",
            "VLDLR": "receptor",
            "LRP8": "receptor",
            "DAB1": "effector",
            "CRK": "adaptor",
            "CRKL": "adaptor",
            "RAP1A": "effector",
            "PAFAH1B1": "effector",
            "CDK5": "effector",
            "CDK5R1": "effector",
        },
    },
    {
        "pathway_id": "semaphorin_guidance",
        "pathway_label": "Semaphorin/plexin guidance",
        "pathway_family": "differentiation_stop",
        "hypothesis_role": "axon guidance, pruning, fasciculation, and local circuit architecture",
        "genes": {
            "SEMA3A": "ligand",
            "SEMA3C": "ligand",
            "SEMA3E": "ligand",
            "SEMA5A": "ligand",
            "SEMA6A": "ligand",
            "SEMA6D": "ligand",
            "SEMA7A": "ligand",
            "PLXNA1": "receptor",
            "PLXNA2": "receptor",
            "PLXNA3": "receptor",
            "PLXNA4": "receptor",
            "PLXNB1": "receptor",
            "PLXNC1": "receptor",
            "NRP1": "coreceptor",
            "NRP2": "coreceptor",
            "RHOA": "effector",
            "RAC1": "effector",
        },
    },
    {
        "pathway_id": "shh_granule_expansion",
        "pathway_label": "SHH/PTCH/GLI",
        "pathway_family": "progenitor_expansion_or_regional_fate",
        "hypothesis_role": "cerebellar granule precursor expansion and regional fate context",
        "genes": {
            "SHH": "ligand",
            "PTCH1": "receptor",
            "PTCH2": "receptor",
            "SMO": "effector",
            "GLI1": "response",
            "GLI2": "response",
            "GLI3": "response",
            "HHIP": "feedback",
            "MYCN": "response",
            "ATOH1": "regional_fate",
        },
    },
    {
        "pathway_id": "wnt_beta_catenin",
        "pathway_label": "WNT/beta-catenin",
        "pathway_family": "neurogenic_permissive",
        "hypothesis_role": "dentate neurogenic permissive and granule fate signaling",
        "genes": {
            "WNT3A": "ligand",
            "WNT5A": "ligand",
            "WNT7A": "ligand",
            "WNT7B": "ligand",
            "FZD1": "receptor",
            "FZD3": "receptor",
            "FZD5": "receptor",
            "FZD7": "receptor",
            "LRP5": "coreceptor",
            "LRP6": "coreceptor",
            "CTNNB1": "effector",
            "LEF1": "response",
            "TCF7L2": "response",
            "AXIN2": "feedback",
            "DKK1": "antagonist",
            "DKK3": "modulator",
        },
    },
    {
        "pathway_id": "fgf_mapk",
        "pathway_label": "FGF/MAPK",
        "pathway_family": "neurogenic_permissive",
        "hypothesis_role": "neurogenic niche maintenance and MAPK-linked growth signaling",
        "genes": {
            "FGF2": "ligand",
            "FGF8": "ligand",
            "FGF9": "ligand",
            "FGF10": "ligand",
            "FGF17": "ligand",
            "FGFR1": "receptor",
            "FGFR2": "receptor",
            "FGFR3": "receptor",
            "FRS2": "adaptor",
            "MAP2K1": "effector",
            "MAPK1": "effector",
            "MAPK3": "effector",
            "SPRY1": "feedback",
            "SPRY2": "feedback",
            "SPRY4": "feedback",
            "ETV4": "response",
            "ETV5": "response",
        },
    },
    {
        "pathway_id": "notch_hes",
        "pathway_label": "Notch/HES",
        "pathway_family": "neurogenic_permissive",
        "hypothesis_role": "persistent neurogenic progenitor-state and differentiation timing control",
        "genes": {
            "NOTCH1": "receptor",
            "NOTCH2": "receptor",
            "NOTCH3": "receptor",
            "DLL1": "ligand",
            "DLL3": "ligand",
            "DLL4": "ligand",
            "JAG1": "ligand",
            "JAG2": "ligand",
            "RBPJ": "effector",
            "HES1": "response",
            "HES5": "response",
            "HEY1": "response",
            "HEY2": "response",
            "LFNG": "modulator",
        },
    },
]


SIGNATURES = [
    {
        "signature_id": "differentiation_stop_index",
        "signature_label": "Differentiation/stop signaling",
        "pathways": [
            "tgf_beta_smad",
            "bdnf_trkb_mapk",
            "bmp_smad",
            "reelin_migration_stop",
            "semaphorin_guidance",
        ],
    },
    {
        "signature_id": "neurogenic_permissive_index",
        "signature_label": "Neurogenic/permissive signaling",
        "pathways": ["wnt_beta_catenin", "fgf_mapk", "notch_hes"],
    },
    {
        "signature_id": "tgf_bdnf_2005_index",
        "signature_label": "TGF-beta/BDNF 2005 mechanism",
        "pathways": ["tgf_beta_smad", "bdnf_trkb_mapk"],
    },
    {
        "signature_id": "stop_minus_permissive_index",
        "signature_label": "Stop minus permissive balance",
        "formula": "differentiation_stop_index - neurogenic_permissive_index",
    },
]


LR_PAIRS = [
    ("tgf_beta_smad", "TGFB2", "TGFBR1", "TGF-beta2 to TGFBR1"),
    ("tgf_beta_smad", "TGFB2", "TGFBR2", "TGF-beta2 to TGFBR2"),
    ("tgf_beta_smad", "TGFB1", "TGFBR2", "TGF-beta1 to TGFBR2"),
    ("bdnf_trkb_mapk", "BDNF", "NTRK2", "BDNF to TrkB"),
    ("shh_granule_expansion", "SHH", "PTCH1", "SHH to PTCH1"),
    ("shh_granule_expansion", "SHH", "SMO", "SHH to SMO readiness"),
    ("wnt_beta_catenin", "WNT3A", "FZD5", "WNT3A to FZD5"),
    ("wnt_beta_catenin", "WNT7A", "FZD3", "WNT7A to FZD3"),
    ("wnt_beta_catenin", "WNT7A", "LRP6", "WNT7A to LRP6 readiness"),
    ("fgf_mapk", "FGF2", "FGFR1", "FGF2 to FGFR1"),
    ("fgf_mapk", "FGF8", "FGFR1", "FGF8 to FGFR1"),
    ("fgf_mapk", "FGF9", "FGFR3", "FGF9 to FGFR3"),
    ("bmp_smad", "BMP4", "BMPR1A", "BMP4 to BMPR1A"),
    ("bmp_smad", "BMP7", "BMPR2", "BMP7 to BMPR2"),
    ("notch_hes", "DLL1", "NOTCH1", "DLL1 to NOTCH1"),
    ("notch_hes", "JAG1", "NOTCH1", "JAG1 to NOTCH1"),
    ("notch_hes", "JAG1", "NOTCH2", "JAG1 to NOTCH2"),
    ("reelin_migration_stop", "RELN", "VLDLR", "Reelin to VLDLR"),
    ("reelin_migration_stop", "RELN", "LRP8", "Reelin to LRP8"),
    ("semaphorin_guidance", "SEMA6A", "PLXNA2", "SEMA6A to PLXNA2"),
    ("semaphorin_guidance", "SEMA6A", "PLXNA4", "SEMA6A to PLXNA4"),
    ("semaphorin_guidance", "SEMA3A", "NRP1", "SEMA3A to NRP1"),
    ("semaphorin_guidance", "SEMA7A", "PLXNC1", "SEMA7A to PLXNC1"),
]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def canon(symbol: object) -> str:
    if pd.isna(symbol):
        return ""
    return str(symbol).strip().upper()


def branch_kind(core_branch: str) -> str:
    branch = str(core_branch).lower()
    if "cerebell" in branch:
        return "cerebellar"
    if "dentate" in branch or "hippocamp" in branch:
        return "dentate"
    return "unknown"


def target_class_for_branch(core_branch: str) -> str:
    kind = branch_kind(core_branch)
    if kind == "cerebellar":
        return "cerebellar_candidate"
    if kind == "dentate":
        return "dentate_candidate"
    return ""


def finite_median(values: pd.Series | list[float]) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.median(arr))


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def write_pathway_gene_sets() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module in PATHWAY_MODULES:
        for order, (gene, gene_role) in enumerate(module["genes"].items(), start=1):
            rows.append(
                {
                    "pathway_id": module["pathway_id"],
                    "pathway_label": module["pathway_label"],
                    "pathway_family": module["pathway_family"],
                    "hypothesis_role": module["hypothesis_role"],
                    "gene": gene,
                    "canonical_gene": canon(gene),
                    "gene_role": gene_role,
                    "gene_order": order,
                }
            )
    gene_sets = pd.DataFrame(rows)
    gene_sets.to_csv(OUT_GENE_SETS, sep="\t", index=False)
    return gene_sets


def write_lr_pair_table() -> pd.DataFrame:
    pathway_lookup = {
        module["pathway_id"]: (module["pathway_label"], module["pathway_family"]) for module in PATHWAY_MODULES
    }
    rows = []
    for pair_id, (pathway_id, ligand, receptor, label) in enumerate(LR_PAIRS, start=1):
        pathway_label, pathway_family = pathway_lookup[pathway_id]
        rows.append(
            {
                "pair_id": f"LR{pair_id:03d}",
                "pathway_id": pathway_id,
                "pathway_label": pathway_label,
                "pathway_family": pathway_family,
                "pair_label": label,
                "ligand_gene": ligand,
                "ligand_canonical_gene": canon(ligand),
                "receptor_gene": receptor,
                "receptor_canonical_gene": canon(receptor),
            }
        )
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT_LR_PAIRS, sep="\t", index=False)
    return pairs


def read_gene_expression(path: Path, genes: set[str], expression_layer: str) -> pd.DataFrame:
    use_cols = [
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "n_cells",
        "canonical_gene",
        "gene",
        "human_symbol",
        "mouse_symbol",
        "detection_fraction",
        "mean_log1p_expression",
        "eligible_class",
        "mean_log1p_rank_within_sample_gene",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", usecols=lambda col: col in use_cols, chunksize=100_000, low_memory=False):
        chunk["canonical_gene"] = chunk["canonical_gene"].map(canon)
        sub = chunk.loc[chunk["canonical_gene"].isin(genes)].copy()
        if not sub.empty:
            sub["expression_layer"] = expression_layer
            pieces.append(sub)
    if not pieces:
        return pd.DataFrame(columns=[*use_cols, "expression_layer"])
    out = pd.concat(pieces, ignore_index=True)
    out["eligible_class"] = bool_series(out["eligible_class"])
    out = out.loc[out["eligible_class"]].copy()
    return out


def load_expression_layers(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_expr: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    gene_sets = pd.read_csv(OUT_GENE_SETS, sep="\t")

    for layer in EXPRESSION_LAYERS:
        expr = read_gene_expression(layer["path"], genes, layer["expression_layer"])
        if not expr.empty:
            all_expr.append(expr)
        for pathway_id, sub in gene_sets.groupby("pathway_id", sort=False):
            present = sorted(set(expr["canonical_gene"]).intersection(set(sub["canonical_gene"])))
            coverage_rows.append(
                {
                    "expression_layer": layer["expression_layer"],
                    "pathway_id": pathway_id,
                    "pathway_label": sub["pathway_label"].iloc[0],
                    "pathway_family": sub["pathway_family"].iloc[0],
                    "n_defined_genes": int(sub["canonical_gene"].nunique()),
                    "n_present_genes": int(len(present)),
                    "present_genes": ",".join(present),
                    "missing_genes": ",".join(sorted(set(sub["canonical_gene"]) - set(present))),
                }
            )
    if not all_expr:
        raise RuntimeError("No Aim 2 pathway genes found in expression layers.")
    expr_all = pd.concat(all_expr, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(OUT_COVERAGE, sep="\t", index=False)
    return expr_all, coverage


def build_pathway_units(expr_all: pd.DataFrame, gene_sets: pd.DataFrame) -> pd.DataFrame:
    module_map = gene_sets[
        ["canonical_gene", "pathway_id", "pathway_label", "pathway_family", "hypothesis_role", "gene_role"]
    ].drop_duplicates()
    expr = expr_all.merge(module_map, on="canonical_gene", how="left")
    group_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "pathway_id",
        "pathway_label",
        "pathway_family",
        "hypothesis_role",
    ]
    units = (
        expr.groupby(group_cols, sort=False)
        .agg(
            n_cells=("n_cells", "max"),
            n_pathway_genes_present=("canonical_gene", "nunique"),
            median_pathway_gene_rank=("mean_log1p_rank_within_sample_gene", "median"),
            mean_pathway_gene_rank=("mean_log1p_rank_within_sample_gene", "mean"),
            median_detection_fraction=("detection_fraction", "median"),
            genes_present=("canonical_gene", lambda values: ",".join(sorted(set(values)))),
            gene_roles_present=("gene_role", lambda values: ",".join(sorted(set(values)))),
        )
        .reset_index()
    )
    units["branch_kind"] = units["core_branch"].map(branch_kind)
    units["target_class_for_branch"] = units["core_branch"].map(target_class_for_branch)
    units["is_branch_candidate"] = units["broad_class"].eq(units["target_class_for_branch"])
    units["is_background_class"] = ~units["is_branch_candidate"]
    units = units.sort_values(["expression_layer", "dataset", "sample", "broad_class", "pathway_id"])
    units.to_csv(OUT_UNITS, sep="\t", index=False, compression="gzip")
    return units


def build_pathway_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "pathway_id",
    ]
    for keys, sub in units.groupby(group_cols, sort=False):
        expression_layer, dataset, core_branch, sample, source_layer, expression_scope, pathway_id = keys
        target = target_class_for_branch(core_branch)
        if not target:
            continue
        candidate = sub.loc[sub["broad_class"].eq(target)]
        background = sub.loc[
            ~sub["broad_class"].eq(target)
            & ~sub["broad_class"].astype(str).str.contains("low_support", case=False, na=False)
        ]
        if candidate.empty or background.empty:
            continue
        cand = finite_median(candidate["median_pathway_gene_rank"])
        bg = finite_median(background["median_pathway_gene_rank"])
        meta = sub.iloc[0]
        rows.append(
            {
                "expression_layer": expression_layer,
                "dataset": dataset,
                "core_branch": core_branch,
                "branch_kind": branch_kind(core_branch),
                "sample": sample,
                "source_layer": source_layer,
                "expression_scope": expression_scope,
                "pathway_id": pathway_id,
                "pathway_label": meta["pathway_label"],
                "pathway_family": meta["pathway_family"],
                "hypothesis_role": meta["hypothesis_role"],
                "target_class": target,
                "background_classes": ",".join(sorted(set(background["broad_class"].astype(str)))),
                "candidate_median_pathway_rank": cand,
                "background_median_pathway_rank": bg,
                "delta_pathway_rank": cand - bg if np.isfinite(cand) and np.isfinite(bg) else np.nan,
                "candidate_median_n_genes_present": finite_median(candidate["n_pathway_genes_present"]),
                "background_median_n_genes_present": finite_median(background["n_pathway_genes_present"]),
                "n_background_classes": int(len(background)),
            }
        )
    contrasts = pd.DataFrame(rows).sort_values(["expression_layer", "branch_kind", "dataset", "sample", "pathway_id"])
    contrasts.to_csv(OUT_CONTRASTS, sep="\t", index=False)
    return contrasts


def positive_test(values: pd.Series) -> tuple[int, int, float, float]:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return 0, 0, np.nan, np.nan
    positive = int(np.sum(arr > 0))
    n = int(arr.size)
    sign_p = float(stats.binomtest(positive, n=n, p=0.5, alternative="greater").pvalue)
    wilcoxon_p = float(stats.wilcoxon(arr, alternative="greater").pvalue) if n > 1 and np.any(arr != 0) else np.nan
    return positive, n, sign_p, wilcoxon_p


def summarize_delta_table(
    table: pd.DataFrame,
    *,
    value_col: str,
    id_cols: list[str],
    out_path: Path,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    groupings: list[tuple[str, list[str]]] = [
        ("all", []),
        ("by_branch", ["branch_kind"]),
        ("by_expression_layer", ["expression_layer"]),
        ("by_id", id_cols),
        ("by_id_and_branch", [*id_cols, "branch_kind"]),
        ("by_id_layer_branch", [*id_cols, "expression_layer", "branch_kind"]),
    ]
    for level, cols in groupings:
        iterator = [((), table)] if not cols else table.groupby(cols, sort=False)
        for keys, sub in iterator:
            if not isinstance(keys, tuple):
                keys = (keys,)
            row: dict[str, object] = {"summary_level": level}
            for col, key in zip(cols, keys):
                row[col] = key
            positive, n, sign_p, wilcoxon_p = positive_test(sub[value_col])
            row["n_contrasts"] = n
            row["n_positive"] = positive
            row["fraction_positive"] = positive / n if n else np.nan
            row["median_delta"] = finite_median(sub[value_col])
            row["mean_delta"] = float(pd.to_numeric(sub[value_col], errors="coerce").mean()) if n else np.nan
            row["sign_test_p"] = sign_p
            row["wilcoxon_p"] = wilcoxon_p
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, sep="\t", index=False)
    return summary


def build_signature_units(pathway_units: pd.DataFrame) -> pd.DataFrame:
    idx_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
        "branch_kind",
        "target_class_for_branch",
        "is_branch_candidate",
        "is_background_class",
    ]
    pivot = pathway_units.pivot_table(
        index=idx_cols,
        columns="pathway_id",
        values="median_pathway_gene_rank",
        aggfunc="median",
    ).reset_index()
    pivot.columns.name = None
    rows: list[dict[str, object]] = []
    for _, row in pivot.iterrows():
        base = {col: row[col] for col in idx_cols}
        computed: dict[str, float] = {}
        for sig in SIGNATURES:
            sig_id = sig["signature_id"]
            if "pathways" in sig:
                values = [row[pathway] for pathway in sig["pathways"] if pathway in row.index and pd.notna(row[pathway])]
                score = float(np.mean(values)) if values else np.nan
            else:
                score = computed.get("differentiation_stop_index", np.nan) - computed.get(
                    "neurogenic_permissive_index", np.nan
                )
            computed[sig_id] = score
            rows.append(
                {
                    **base,
                    "signature_id": sig_id,
                    "signature_label": sig["signature_label"],
                    "signature_score": score,
                }
            )
    units = pd.DataFrame(rows)
    units.to_csv(OUT_SIGNATURE_UNITS, sep="\t", index=False, compression="gzip")
    return units


def build_signature_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "signature_id",
    ]
    for keys, sub in units.groupby(group_cols, sort=False):
        expression_layer, dataset, core_branch, sample, source_layer, expression_scope, signature_id = keys
        target = target_class_for_branch(core_branch)
        if not target:
            continue
        candidate = sub.loc[sub["broad_class"].eq(target)]
        background = sub.loc[
            ~sub["broad_class"].eq(target)
            & ~sub["broad_class"].astype(str).str.contains("low_support", case=False, na=False)
        ]
        if candidate.empty or background.empty:
            continue
        cand = finite_median(candidate["signature_score"])
        bg = finite_median(background["signature_score"])
        rows.append(
            {
                "expression_layer": expression_layer,
                "dataset": dataset,
                "core_branch": core_branch,
                "branch_kind": branch_kind(core_branch),
                "sample": sample,
                "source_layer": source_layer,
                "expression_scope": expression_scope,
                "signature_id": signature_id,
                "signature_label": sub["signature_label"].iloc[0],
                "target_class": target,
                "background_classes": ",".join(sorted(set(background["broad_class"].astype(str)))),
                "candidate_median_signature_score": cand,
                "background_median_signature_score": bg,
                "delta_signature_score": cand - bg if np.isfinite(cand) and np.isfinite(bg) else np.nan,
                "n_background_classes": int(len(background)),
            }
        )
    contrasts = pd.DataFrame(rows).sort_values(["expression_layer", "branch_kind", "dataset", "sample", "signature_id"])
    contrasts.to_csv(OUT_SIGNATURE_CONTRASTS, sep="\t", index=False)
    return contrasts


def build_lr_units(expr_all: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    pair_genes = set(pairs["ligand_canonical_gene"]) | set(pairs["receptor_canonical_gene"])
    expr = expr_all.loc[expr_all["canonical_gene"].isin(pair_genes)].copy()
    idx_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "broad_class",
    ]
    rank_pivot = expr.pivot_table(
        index=idx_cols,
        columns="canonical_gene",
        values="mean_log1p_rank_within_sample_gene",
        aggfunc="median",
    ).reset_index()
    rank_pivot.columns.name = None
    det_pivot = expr.pivot_table(
        index=idx_cols,
        columns="canonical_gene",
        values="detection_fraction",
        aggfunc="median",
    ).reset_index()
    det_pivot.columns.name = None

    rows: list[dict[str, object]] = []
    for _, base in rank_pivot.iterrows():
        match = det_pivot
        for col in idx_cols:
            match = match.loc[match[col].eq(base[col])]
        det_row = match.iloc[0] if not match.empty else None
        for _, pair in pairs.iterrows():
            ligand = pair["ligand_canonical_gene"]
            receptor = pair["receptor_canonical_gene"]
            ligand_rank = base[ligand] if ligand in rank_pivot.columns else np.nan
            receptor_rank = base[receptor] if receptor in rank_pivot.columns else np.nan
            ligand_det = det_row[ligand] if det_row is not None and ligand in det_pivot.columns else np.nan
            receptor_det = det_row[receptor] if det_row is not None and receptor in det_pivot.columns else np.nan
            if not (np.isfinite(ligand_rank) or np.isfinite(receptor_rank)):
                continue
            rows.append(
                {
                    **{col: base[col] for col in idx_cols},
                    "pair_id": pair["pair_id"],
                    "pathway_id": pair["pathway_id"],
                    "pathway_label": pair["pathway_label"],
                    "pathway_family": pair["pathway_family"],
                    "pair_label": pair["pair_label"],
                    "ligand_gene": pair["ligand_gene"],
                    "receptor_gene": pair["receptor_gene"],
                    "ligand_rank": ligand_rank,
                    "receptor_rank": receptor_rank,
                    "ligand_detection_fraction": ligand_det,
                    "receptor_detection_fraction": receptor_det,
                    "pair_mean_rank": np.nanmean([ligand_rank, receptor_rank]),
                    "pair_min_rank": np.nanmin([ligand_rank, receptor_rank])
                    if np.isfinite(ligand_rank) and np.isfinite(receptor_rank)
                    else np.nan,
                    "pair_min_detection": np.nanmin([ligand_det, receptor_det])
                    if np.isfinite(ligand_det) and np.isfinite(receptor_det)
                    else np.nan,
                }
            )
    units = pd.DataFrame(rows)
    units["branch_kind"] = units["core_branch"].map(branch_kind)
    units["target_class_for_branch"] = units["core_branch"].map(target_class_for_branch)
    units["is_branch_candidate"] = units["broad_class"].eq(units["target_class_for_branch"])
    units["is_background_class"] = ~units["is_branch_candidate"]
    units.to_csv(OUT_LR_UNITS, sep="\t", index=False, compression="gzip")
    return units


def build_lr_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [
        "expression_layer",
        "dataset",
        "core_branch",
        "sample",
        "source_layer",
        "expression_scope",
        "pair_id",
    ]
    for keys, sub in units.groupby(group_cols, sort=False):
        expression_layer, dataset, core_branch, sample, source_layer, expression_scope, pair_id = keys
        target = target_class_for_branch(core_branch)
        if not target:
            continue
        candidate = sub.loc[sub["broad_class"].eq(target)]
        background = sub.loc[
            ~sub["broad_class"].eq(target)
            & ~sub["broad_class"].astype(str).str.contains("low_support", case=False, na=False)
        ]
        if candidate.empty or background.empty:
            continue
        cand = finite_median(candidate["pair_min_rank"])
        bg = finite_median(background["pair_min_rank"])
        if not (np.isfinite(cand) and np.isfinite(bg)):
            continue
        meta = sub.iloc[0]
        rows.append(
            {
                "expression_layer": expression_layer,
                "dataset": dataset,
                "core_branch": core_branch,
                "branch_kind": branch_kind(core_branch),
                "sample": sample,
                "source_layer": source_layer,
                "expression_scope": expression_scope,
                "pair_id": pair_id,
                "pathway_id": meta["pathway_id"],
                "pathway_label": meta["pathway_label"],
                "pathway_family": meta["pathway_family"],
                "pair_label": meta["pair_label"],
                "target_class": target,
                "background_classes": ",".join(sorted(set(background["broad_class"].astype(str)))),
                "candidate_median_pair_min_rank": cand,
                "background_median_pair_min_rank": bg,
                "delta_pair_min_rank": cand - bg,
                "candidate_median_pair_mean_rank": finite_median(candidate["pair_mean_rank"]),
                "background_median_pair_mean_rank": finite_median(background["pair_mean_rank"]),
                "n_background_classes": int(len(background)),
            }
        )
    contrasts = pd.DataFrame(rows).sort_values(["expression_layer", "branch_kind", "dataset", "sample", "pair_id"])
    contrasts.to_csv(OUT_LR_CONTRASTS, sep="\t", index=False)
    return contrasts


def plot_results(pathway_summary: pd.DataFrame, signature_summary: pd.DataFrame, lr_summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12.8, 12.8), constrained_layout=True)

    branch_colors = {"dentate": "#2f7f8f", "cerebellar": "#7f4e8a"}

    ax = axes[0]
    pathway_plot = pathway_summary.loc[pathway_summary["summary_level"].eq("by_id_and_branch")].copy()
    pathway_order = [module["pathway_id"] for module in PATHWAY_MODULES]
    pathway_plot["pathway_id"] = pd.Categorical(pathway_plot["pathway_id"], pathway_order, ordered=True)
    pathway_plot = pathway_plot.sort_values(["pathway_id", "branch_kind"])
    x = np.arange(len(pathway_order))
    width = 0.38
    for offset, branch in [(-width / 2, "dentate"), (width / 2, "cerebellar")]:
        vals = []
        labels = []
        for pathway in pathway_order:
            sub = pathway_plot.loc[pathway_plot["pathway_id"].eq(pathway) & pathway_plot["branch_kind"].eq(branch)]
            vals.append(float(sub["median_delta"].iloc[0]) if not sub.empty else np.nan)
            labels.append(pathway)
        ax.bar(x + offset, vals, width=width, label=branch, color=branch_colors[branch])
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n") for p in pathway_order], fontsize=8)
    ax.set_ylabel("Median candidate-background\npathway-rank delta")
    ax.set_title("Aim 2 pathway module contrasts")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    ax = axes[1]
    sig_plot = signature_summary.loc[signature_summary["summary_level"].eq("by_id_and_branch")].copy()
    sig_order = [sig["signature_id"] for sig in SIGNATURES]
    x = np.arange(len(sig_order))
    for offset, branch in [(-width / 2, "dentate"), (width / 2, "cerebellar")]:
        vals = []
        for sig in sig_order:
            sub = sig_plot.loc[sig_plot["signature_id"].eq(sig) & sig_plot["branch_kind"].eq(branch)]
            vals.append(float(sub["median_delta"].iloc[0]) if not sub.empty else np.nan)
        ax.bar(x + offset, vals, width=width, label=branch, color=branch_colors[branch])
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", "\n") for s in sig_order], fontsize=8)
    ax.set_ylabel("Median candidate-background\nsignature delta")
    ax.set_title("Differentiation/stop versus neurogenic/permissive balance")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    ax = axes[2]
    lr_plot = lr_summary.loc[lr_summary["summary_level"].eq("by_id_and_branch")].copy()
    lr_pathways = [module["pathway_id"] for module in PATHWAY_MODULES]
    lr_plot = (
        lr_plot.groupby(["pathway_id", "branch_kind"], sort=False)
        .agg(median_delta=("median_delta", "median"))
        .reset_index()
    )
    x = np.arange(len(lr_pathways))
    for offset, branch in [(-width / 2, "dentate"), (width / 2, "cerebellar")]:
        vals = []
        for pathway in lr_pathways:
            sub = lr_plot.loc[lr_plot["pathway_id"].eq(pathway) & lr_plot["branch_kind"].eq(branch)]
            vals.append(float(sub["median_delta"].iloc[0]) if not sub.empty else np.nan)
        ax.bar(x + offset, vals, width=width, label=branch, color=branch_colors[branch])
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n") for p in lr_pathways], fontsize=8)
    ax.set_ylabel("Median candidate-background\nLR min-rank delta")
    ax.set_title("Ligand-receptor readiness contrasts")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    fig.suptitle("Aim 2 Niche/Pathway Signaling Audit", fontsize=16, y=1.01)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def row_for(summary: pd.DataFrame, level: str, **filters: object) -> pd.Series | None:
    sub = summary.loc[summary["summary_level"].eq(level)].copy()
    for col, value in filters.items():
        if col not in sub.columns:
            return None
        sub = sub.loc[sub[col].eq(value)]
    if sub.empty:
        return None
    return sub.iloc[0]


def format_p(value: object) -> str:
    try:
        val = float(value)
    except Exception:
        return "NA"
    if not np.isfinite(val):
        return "NA"
    return f"{val:.3g}"


def write_report(
    *,
    units: pd.DataFrame,
    contrasts: pd.DataFrame,
    summary: pd.DataFrame,
    coverage: pd.DataFrame,
    signature_units: pd.DataFrame,
    signature_contrasts: pd.DataFrame,
    signature_summary: pd.DataFrame,
    lr_units: pd.DataFrame,
    lr_contrasts: pd.DataFrame,
    lr_summary: pd.DataFrame,
) -> None:
    all_pathway = row_for(summary, "all")
    all_sig = row_for(signature_summary, "all")
    all_lr = row_for(lr_summary, "all")
    tgf_bdnf_all = row_for(signature_summary, "by_id", signature_id="tgf_bdnf_2005_index")
    sig_full = row_for(signature_summary, "by_expression_layer", expression_layer="full_mgi_ortholog_matrix")
    sig_selected = row_for(signature_summary, "by_expression_layer", expression_layer="selected_feature_matrix")
    stop_dentate = row_for(signature_summary, "by_id_and_branch", signature_id="differentiation_stop_index", branch_kind="dentate")
    stop_cereb = row_for(signature_summary, "by_id_and_branch", signature_id="differentiation_stop_index", branch_kind="cerebellar")
    perm_dentate = row_for(signature_summary, "by_id_and_branch", signature_id="neurogenic_permissive_index", branch_kind="dentate")
    perm_cereb = row_for(signature_summary, "by_id_and_branch", signature_id="neurogenic_permissive_index", branch_kind="cerebellar")
    balance_dentate = row_for(signature_summary, "by_id_and_branch", signature_id="stop_minus_permissive_index", branch_kind="dentate")
    balance_cereb = row_for(signature_summary, "by_id_and_branch", signature_id="stop_minus_permissive_index", branch_kind="cerebellar")
    tgf_bdnf_dentate = row_for(signature_summary, "by_id_and_branch", signature_id="tgf_bdnf_2005_index", branch_kind="dentate")
    tgf_bdnf_cereb = row_for(signature_summary, "by_id_and_branch", signature_id="tgf_bdnf_2005_index", branch_kind="cerebellar")
    shh_cereb = row_for(summary, "by_id_and_branch", pathway_id="shh_granule_expansion", branch_kind="cerebellar")
    tgf_path_cereb = row_for(summary, "by_id_and_branch", pathway_id="tgf_beta_smad", branch_kind="cerebellar")
    tgf_path_dentate = row_for(summary, "by_id_and_branch", pathway_id="tgf_beta_smad", branch_kind="dentate")

    lines = [
        "# Aim 2 Niche/Pathway Signaling Audit",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis addresses Specific Aim 2 by scoring targeted niche, ligand-receptor, and pathway-readiness programs in the primary-core pseudobulk expression layers.",
        "",
        "The analysis is intentionally conservative: the available primary-core matrices are broad-class expression summaries, not spatial sender-receiver assays. Therefore the results test pathway/module readiness and candidate-versus-background enrichment, not direct cell-cell communication.",
        "",
        "## Inputs",
        "",
        "- Full MGI one-to-one ortholog expression layer.",
        "- 2,169-gene selected-feature expression layer.",
        "- Curated pathway modules for TGF-beta/SMAD, BDNF/TrkB/MAPK, BMP/SMAD, Reelin, Semaphorin, SHH, WNT, FGF, and Notch.",
        "- Curated ligand-receptor pairs centered on the 2005 TGF-beta2/BDNF mechanism and related niche pathways.",
        "",
        "## Scale",
        "",
        f"- Pathway class units: {len(units):,} across {units['dataset'].nunique():,} datasets.",
        f"- Pathway candidate-background contrasts: {len(contrasts):,}.",
        f"- Signature class units: {len(signature_units):,}.",
        f"- Signature candidate-background contrasts: {len(signature_contrasts):,}.",
        f"- Ligand-receptor readiness units: {len(lr_units):,}.",
        f"- Ligand-receptor candidate-background contrasts: {len(lr_contrasts):,}.",
        "",
        "## Main Results",
        "",
    ]

    if all_pathway is not None:
        lines.append(
            f"- Pathway modules were candidate-enriched in {int(all_pathway['n_positive'])}/{int(all_pathway['n_contrasts'])} contrasts "
            f"(median delta {float(all_pathway['median_delta']):.3f}, Wilcoxon p={format_p(all_pathway['wilcoxon_p'])})."
        )
    if all_sig is not None:
        lines.append(
            f"- Composite signatures were candidate-enriched in {int(all_sig['n_positive'])}/{int(all_sig['n_contrasts'])} contrasts "
            f"(median delta {float(all_sig['median_delta']):.3f}, Wilcoxon p={format_p(all_sig['wilcoxon_p'])})."
        )
    if all_lr is not None:
        lines.append(
            f"- Ligand-receptor readiness pairs were candidate-enriched in {int(all_lr['n_positive'])}/{int(all_lr['n_contrasts'])} contrasts "
            f"(median delta {float(all_lr['median_delta']):.3f}, Wilcoxon p={format_p(all_lr['wilcoxon_p'])})."
        )
    if tgf_bdnf_all is not None:
        lines.append(
            f"- The TGF-beta/BDNF 2005 mechanism index was positive in {int(tgf_bdnf_all['n_positive'])}/{int(tgf_bdnf_all['n_contrasts'])} contrasts "
            f"(median delta {float(tgf_bdnf_all['median_delta']):.3f}, Wilcoxon p={format_p(tgf_bdnf_all['wilcoxon_p'])})."
        )
    if sig_full is not None and sig_selected is not None:
        lines.append(
            f"- Layer caveat: the composite-signature signal is driven mainly by the selected-feature matrix "
            f"({int(sig_selected['n_positive'])}/{int(sig_selected['n_contrasts'])} positive, median delta {float(sig_selected['median_delta']):.3f}); "
            f"the full MGI layer is not broadly positive "
            f"({int(sig_full['n_positive'])}/{int(sig_full['n_contrasts'])} positive, median delta {float(sig_full['median_delta']):.3f})."
        )

    lines.extend(["", "Branch-level signature results:"])
    for label, row in [
        ("Dentate differentiation/stop", stop_dentate),
        ("Cerebellar differentiation/stop", stop_cereb),
        ("Dentate neurogenic/permissive", perm_dentate),
        ("Cerebellar neurogenic/permissive", perm_cereb),
        ("Dentate stop-minus-permissive", balance_dentate),
        ("Cerebellar stop-minus-permissive", balance_cereb),
    ]:
        if row is not None:
            lines.append(
                f"- {label}: {int(row['n_positive'])}/{int(row['n_contrasts'])} positive, "
                f"median delta {float(row['median_delta']):.3f}, Wilcoxon p={format_p(row['wilcoxon_p'])}."
            )

    lines.extend(["", "## Aim 2 Hypothesis Test Outcome", ""])
    lines.append(
        "The starting prediction was that cerebellar granule-cell context might show stronger differentiation/stop signaling, while the dentate context might retain a more persistent neurogenic/permissive signature."
    )
    if stop_cereb is not None and tgf_bdnf_cereb is not None and shh_cereb is not None:
        lines.append(
            f"- Cerebellar candidates do not show broad differentiation/stop enrichment in this pseudobulk audit: "
            f"differentiation/stop {int(stop_cereb['n_positive'])}/{int(stop_cereb['n_contrasts'])} positive, "
            f"median delta {float(stop_cereb['median_delta']):.3f}; TGF-beta/BDNF index "
            f"{int(tgf_bdnf_cereb['n_positive'])}/{int(tgf_bdnf_cereb['n_contrasts'])} positive, "
            f"median delta {float(tgf_bdnf_cereb['median_delta']):.3f}."
        )
        lines.append(
            f"- The clearest cerebellar pathway signal is SHH/PTCH/GLI: "
            f"{int(shh_cereb['n_positive'])}/{int(shh_cereb['n_contrasts'])} positive, "
            f"median delta {float(shh_cereb['median_delta']):.3f}."
        )
    if stop_dentate is not None and perm_dentate is not None and tgf_bdnf_dentate is not None:
        lines.append(
            f"- Dentate candidates show broad pathway readiness: differentiation/stop "
            f"{int(stop_dentate['n_positive'])}/{int(stop_dentate['n_contrasts'])} positive, "
            f"neurogenic/permissive {int(perm_dentate['n_positive'])}/{int(perm_dentate['n_contrasts'])} positive, "
            f"and TGF-beta/BDNF index {int(tgf_bdnf_dentate['n_positive'])}/{int(tgf_bdnf_dentate['n_contrasts'])} positive."
        )
    if tgf_path_cereb is not None and tgf_path_dentate is not None:
        lines.append(
            f"- The TGF-beta/SMAD pathway itself is dentate-enriched rather than cerebellar-enriched in candidate-background contrasts: "
            f"dentate median delta {float(tgf_path_dentate['median_delta']):.3f}, cerebellar median delta {float(tgf_path_cereb['median_delta']):.3f}."
        )
    lines.append(
        "Conclusion: Aim 2 supports context-dependent niche/maturation signaling, but not the simple version of a cerebellar-biased TGF-beta/BDNF stop-signaling program."
    )

    lines.extend(["", "Pathway-level highlights:"])
    pathway_highlights = summary.loc[summary["summary_level"].eq("by_id_and_branch")].copy()
    pathway_highlights = pathway_highlights.sort_values("median_delta", ascending=False)
    for _, row in pathway_highlights.head(10).iterrows():
        lines.append(
            f"- {row['pathway_id']} / {row['branch_kind']}: median delta {float(row['median_delta']):.3f}; "
            f"{int(row['n_positive'])}/{int(row['n_contrasts'])} positive."
        )

    coverage_min = coverage.groupby("pathway_label")["n_present_genes"].min().sort_values()
    lines.extend(["", "Minimum gene coverage across expression layers:"])
    for pathway_label, n_present in coverage_min.items():
        lines.append(f"- {pathway_label}: {int(n_present)} genes present.")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Aim 2 is now computably tested at the pathway-readiness level.",
            "- The historical TGF-beta2/BDNF/SMAD-MAPK mechanism can be discussed as a prioritized niche/maturation axis, but in this primary-core audit it behaves as context-dependent and stronger in the dentate branch than in cerebellar candidate-background contrasts.",
            "- SHH/PTCH/GLI is the strongest cerebellar pathway signal, consistent with cerebellar granule precursor biology.",
            "- The key manuscript-safe interpretation should distinguish pathway readiness from true spatial niche signaling.",
            "- If stronger Aim 2 evidence is required, the next layer should use spatial datasets, ligand-sender/receiver cell labels, or raw object-level cell-type communication tools.",
            "",
            "## Outputs",
            "",
            f"- Pathway gene sets: `{rel(OUT_GENE_SETS)}`",
            f"- Pathway units: `{rel(OUT_UNITS)}`",
            f"- Pathway contrasts: `{rel(OUT_CONTRASTS)}`",
            f"- Pathway summary: `{rel(OUT_SUMMARY)}`",
            f"- Pathway coverage: `{rel(OUT_COVERAGE)}`",
            f"- Signature units: `{rel(OUT_SIGNATURE_UNITS)}`",
            f"- Signature contrasts: `{rel(OUT_SIGNATURE_CONTRASTS)}`",
            f"- Signature summary: `{rel(OUT_SIGNATURE_SUMMARY)}`",
            f"- Ligand-receptor pairs: `{rel(OUT_LR_PAIRS)}`",
            f"- Ligand-receptor units: `{rel(OUT_LR_UNITS)}`",
            f"- Ligand-receptor contrasts: `{rel(OUT_LR_CONTRASTS)}`",
            f"- Ligand-receptor summary: `{rel(OUT_LR_SUMMARY)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    gene_sets = write_pathway_gene_sets()
    lr_pairs = write_lr_pair_table()
    wanted_genes = set(gene_sets["canonical_gene"]) | set(lr_pairs["ligand_canonical_gene"]) | set(
        lr_pairs["receptor_canonical_gene"]
    )
    expr_all, coverage = load_expression_layers(wanted_genes)
    units = build_pathway_units(expr_all, gene_sets)
    contrasts = build_pathway_contrasts(units)
    summary = summarize_delta_table(
        contrasts, value_col="delta_pathway_rank", id_cols=["pathway_id"], out_path=OUT_SUMMARY
    )
    signature_units = build_signature_units(units)
    signature_contrasts = build_signature_contrasts(signature_units)
    signature_summary = summarize_delta_table(
        signature_contrasts,
        value_col="delta_signature_score",
        id_cols=["signature_id"],
        out_path=OUT_SIGNATURE_SUMMARY,
    )
    lr_units = build_lr_units(expr_all, lr_pairs)
    lr_contrasts = build_lr_contrasts(lr_units)
    lr_summary = summarize_delta_table(
        lr_contrasts,
        value_col="delta_pair_min_rank",
        id_cols=["pathway_id", "pair_id"],
        out_path=OUT_LR_SUMMARY,
    )
    plot_results(summary, signature_summary, lr_summary)
    write_report(
        units=units,
        contrasts=contrasts,
        summary=summary,
        coverage=coverage,
        signature_units=signature_units,
        signature_contrasts=signature_contrasts,
        signature_summary=signature_summary,
        lr_units=lr_units,
        lr_contrasts=lr_contrasts,
        lr_summary=lr_summary,
    )

    all_sig = row_for(signature_summary, "all")
    tgf_bdnf = row_for(signature_summary, "by_id", signature_id="tgf_bdnf_2005_index")
    print(f"Wrote {rel(OUT_MD)}")
    if all_sig is not None:
        print(
            "signature_all",
            int(all_sig["n_positive"]),
            "/",
            int(all_sig["n_contrasts"]),
            "median_delta",
            f"{float(all_sig['median_delta']):.3f}",
            "wilcoxon_p",
            format_p(all_sig["wilcoxon_p"]),
        )
    if tgf_bdnf is not None:
        print(
            "tgf_bdnf_2005",
            int(tgf_bdnf["n_positive"]),
            "/",
            int(tgf_bdnf["n_contrasts"]),
            "median_delta",
            f"{float(tgf_bdnf['median_delta']):.3f}",
            "wilcoxon_p",
            format_p(tgf_bdnf["wilcoxon_p"]),
        )


if __name__ == "__main__":
    main()
