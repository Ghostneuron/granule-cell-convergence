#!/usr/bin/env python3
"""MGI-filtered dataset-aware meta-model for primary-core pseudobulk screens.

This model is intentionally conservative. The local full-matrix extraction was
made in a same-symbol frame, so this script keeps only MGI one_to_one
human-mouse ortholog classes where the human and mouse symbols collapse to the
same canonical symbol. That avoids many-to-many/paralog ambiguity while making
the current same-symbol screen more ortholog-aware.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

MGI_HOMOLOGY = EXTERNAL / "Orthology/HOM_MouseHumanSequence.rpt"
SELECTED_EXPR = RESULTS / "primary_core_expanded_gene_pseudobulk_expression.tsv.gz"
GENOME_EXPR = RESULTS / "primary_core_genomewide_symbol_pseudobulk_expression.tsv.gz"
SELECTED_TRIAGE = RESULTS / "primary_core_expanded_gene_mechanism_triage.tsv"
GENOME_TRIAGE = RESULTS / "primary_core_genomewide_symbol_mechanism_triage.tsv"
CONSENSUS = RESULTS / "primary_core_cross_screen_mechanism_consensus.tsv"
CONSENSUS_VALIDATION = RESULTS / "primary_core_consensus_candidate_dataset_validation.tsv"

OUT_ORTHOLOG_MAP = RESULTS / "primary_core_mgi_ortholog_meta_model_map.tsv"
OUT_DELTAS = RESULTS / "primary_core_mgi_ortholog_meta_model_unit_deltas.tsv.gz"
OUT_BRANCH = RESULTS / "primary_core_mgi_ortholog_meta_model_branch_summary.tsv"
OUT_GENE = RESULTS / "primary_core_mgi_ortholog_meta_model_gene_summary.tsv"
OUT_HITS = RESULTS / "primary_core_mgi_ortholog_meta_model_shared_hits.tsv"
OUT_MECHANISM_HITS = RESULTS / "primary_core_mgi_ortholog_meta_model_mechanism_hits.tsv"
OUT_PLOT = RESULTS / "primary_core_mgi_ortholog_meta_model_top_hits.png"
OUT_MD = RESULTS / "primary_core_mgi_ortholog_meta_model.md"

MIN_CLASS_CELLS = 20
MIN_DATASETS_PER_BRANCH = 2
SUPPORT_POSITIVE_DATASET_FRACTION = 0.75
STRICT_DATASET_SIGN_P = 0.25

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DENTATE_BACKGROUNDS = {
    "non_dentate_background",
    "dentate_low_support",
    "other_or_ambiguous",
    "broad_neuronal_structural_warning",
}
CEREBELLAR_BACKGROUNDS = {"other_or_ambiguous", "broad_neuronal_structural_warning"}


def canon_gene(gene: object) -> str:
    if pd.isna(gene):
        return ""
    return str(gene).strip().strip('"').strip("'").upper()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    if not valid.any():
        return out
    idx = p_values.index[valid].to_numpy()
    order = np.argsort(p_values.loc[idx].to_numpy(dtype=float))
    p_sorted = p_values.loc[idx[order]].to_numpy(dtype=float)
    adjusted = np.minimum.accumulate((p_sorted * len(p_sorted) / np.arange(1, len(p_sorted) + 1))[::-1])[::-1]
    out.loc[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def sign_p_greater(n_positive: int, n_total: int) -> float:
    if n_total <= 0:
        return np.nan
    return float(stats.binomtest(n_positive, n_total, 0.5, alternative="greater").pvalue)


def wilcoxon_p_greater(values: pd.Series) -> float:
    vals = values.dropna().to_numpy(dtype=float)
    if len(vals) < 3 or np.allclose(vals, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(vals, alternative="greater", zero_method="wilcox").pvalue)
    except ValueError:
        return np.nan


def bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].map(lambda value: False if pd.isna(value) else str(value).lower() in {"true", "1", "yes"})


def load_mgi_ortholog_map() -> tuple[pd.DataFrame, dict[str, int]]:
    df = pd.read_csv(MGI_HOMOLOGY, sep="\t", dtype=str).fillna("")
    df["db_class_key"] = df["DB Class Key"].astype(str)
    df["taxon_id"] = df["NCBI Taxon ID"].astype(str)
    df["canonical_symbol"] = df["Symbol"].map(canon_gene)

    rows: list[dict[str, object]] = []
    classes_with_human_mouse = 0
    for key, sub in df.groupby("db_class_key", sort=False):
        human = sub.loc[sub["taxon_id"].eq("9606")].copy()
        mouse = sub.loc[sub["taxon_id"].eq("10090")].copy()
        if human.empty or mouse.empty:
            continue
        classes_with_human_mouse += 1
        human_symbols = sorted(set(human["Symbol"].astype(str)))
        mouse_symbols = sorted(set(mouse["Symbol"].astype(str)))
        human_canon = sorted(set(human["canonical_symbol"]))
        mouse_canon = sorted(set(mouse["canonical_symbol"]))
        one_to_one = len(human_symbols) == 1 and len(mouse_symbols) == 1
        same_symbol = one_to_one and human_canon[0] == mouse_canon[0]
        if one_to_one:
            human_row = human.iloc[0]
            mouse_row = mouse.iloc[0]
            rows.append(
                {
                    "db_class_key": key,
                    "canonical_gene": human_canon[0],
                    "human_symbol": human_symbols[0],
                    "mouse_symbol": mouse_symbols[0],
                    "canonical_mouse_symbol": mouse_canon[0],
                    "mgi_one_to_one_human_mouse": True,
                    "same_canonical_symbol": bool(same_symbol),
                    "strict_same_symbol_one_to_one": bool(same_symbol),
                    "human_entrez_id": human_row.get("EntrezGene ID", ""),
                    "mouse_entrez_id": mouse_row.get("EntrezGene ID", ""),
                    "human_hgnc_id": human_row.get("HGNC ID", ""),
                    "mouse_mgi_id": mouse_row.get("Mouse MGI ID", ""),
                    "human_location": human_row.get("Genetic Location", ""),
                    "mouse_location": mouse_row.get("Genetic Location", ""),
                }
            )

    ortho = pd.DataFrame(rows).drop_duplicates("canonical_gene")
    summary = {
        "mgi_rows": len(df),
        "mgi_classes_with_human_mouse": classes_with_human_mouse,
        "one_to_one_pairs": int(ortho["mgi_one_to_one_human_mouse"].sum()),
        "strict_same_symbol_one_to_one_pairs": int(ortho["strict_same_symbol_one_to_one"].sum()),
    }
    return ortho.sort_values("canonical_gene"), summary


def load_expression(path: Path, screen: str, strict_genes: set[str]) -> pd.DataFrame:
    usecols = [
        "dataset",
        "core_branch",
        "sample",
        "broad_class",
        "n_cells",
        "gene",
        "canonical_gene",
        "detection_fraction",
        "mean_log1p_expression",
        "eligible_class",
        "mean_log1p_rank_within_sample_gene",
    ]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    df = df.loc[df["canonical_gene"].isin(strict_genes)].copy()
    df["screen"] = screen
    df["eligible_class"] = df["eligible_class"].astype(str).str.lower().isin(["true", "1", "yes"])
    df["n_cells"] = pd.to_numeric(df["n_cells"], errors="coerce").fillna(0).astype(int)
    df["detection_fraction"] = pd.to_numeric(df["detection_fraction"], errors="coerce")
    df["mean_log1p_expression"] = pd.to_numeric(df["mean_log1p_expression"], errors="coerce")
    df["mean_log1p_rank_within_sample_gene"] = pd.to_numeric(
        df["mean_log1p_rank_within_sample_gene"], errors="coerce"
    )
    return df


def summarize_sample_gene(sub: pd.DataFrame, branch: str) -> dict[str, object] | None:
    if branch == "dentate":
        target_class = "dentate_candidate"
        backgrounds = DENTATE_BACKGROUNDS
    elif branch == "cerebellar":
        target_class = "cerebellar_candidate"
        backgrounds = CEREBELLAR_BACKGROUNDS
    else:
        raise ValueError(branch)

    eligible = sub.loc[sub["eligible_class"] & sub["n_cells"].ge(MIN_CLASS_CELLS)].copy()
    target = eligible.loc[eligible["broad_class"].eq(target_class)]
    background = eligible.loc[eligible["broad_class"].isin(backgrounds)]
    if target.empty or background.empty:
        return None
    candidate_rank = target["mean_log1p_rank_within_sample_gene"].median()
    background_rank = background["mean_log1p_rank_within_sample_gene"].median()
    candidate_expr = target["mean_log1p_expression"].median()
    background_expr = background["mean_log1p_expression"].median()
    candidate_detection = target["detection_fraction"].median()
    background_detection = background["detection_fraction"].median()
    return {
        "candidate_rank": candidate_rank,
        "background_rank": background_rank,
        "rank_delta": candidate_rank - background_rank,
        "candidate_mean_log1p_expression": candidate_expr,
        "background_mean_log1p_expression": background_expr,
        "mean_log1p_expression_delta": candidate_expr - background_expr,
        "candidate_detection": candidate_detection,
        "background_detection": background_detection,
        "detection_delta": candidate_detection - background_detection,
        "candidate_n_cells": int(target["n_cells"].sum()),
        "background_n_cells": int(background["n_cells"].sum()),
        "n_background_classes": int(background["broad_class"].nunique()),
    }


def build_unit_deltas(expr: pd.DataFrame, ortho: pd.DataFrame) -> pd.DataFrame:
    eligible = expr.loc[expr["eligible_class"] & expr["n_cells"].ge(MIN_CLASS_CELLS)].copy()
    eligible["branch_tested"] = np.select(
        [
            eligible["core_branch"].isin(["mouse_dentate", "human_dentate_hippocampus"]),
            eligible["core_branch"].eq("cerebellum"),
        ],
        ["dentate", "cerebellar"],
        default="",
    )
    eligible = eligible.loc[eligible["branch_tested"].ne("")].copy()

    dentate = eligible["branch_tested"].eq("dentate")
    cerebellar = eligible["branch_tested"].eq("cerebellar")
    candidate = (dentate & eligible["broad_class"].eq("dentate_candidate")) | (
        cerebellar & eligible["broad_class"].eq("cerebellar_candidate")
    )
    background = (dentate & eligible["broad_class"].isin(DENTATE_BACKGROUNDS)) | (
        cerebellar & eligible["broad_class"].isin(CEREBELLAR_BACKGROUNDS)
    )
    eligible["delta_role"] = np.select([candidate, background], ["candidate", "background"], default="")
    eligible = eligible.loc[eligible["delta_role"].ne("")].copy()

    group_cols = ["screen", "dataset", "sample", "core_branch", "branch_tested", "canonical_gene"]
    grouped = (
        eligible.groupby(group_cols + ["delta_role"], sort=False)
        .agg(
            source_gene=("gene", "first"),
            rank=("mean_log1p_rank_within_sample_gene", "median"),
            mean_log1p_expression=("mean_log1p_expression", "median"),
            detection=("detection_fraction", "median"),
            n_cells=("n_cells", "sum"),
            n_classes=("broad_class", "nunique"),
        )
        .reset_index()
    )

    target = grouped.loc[grouped["delta_role"].eq("candidate")].drop(columns=["delta_role", "n_classes"])
    target = target.rename(
        columns={
            "source_gene": "source_gene",
            "rank": "candidate_rank",
            "mean_log1p_expression": "candidate_mean_log1p_expression",
            "detection": "candidate_detection",
            "n_cells": "candidate_n_cells",
        }
    )
    bg = grouped.loc[grouped["delta_role"].eq("background")].drop(columns=["delta_role", "source_gene"])
    bg = bg.rename(
        columns={
            "rank": "background_rank",
            "mean_log1p_expression": "background_mean_log1p_expression",
            "detection": "background_detection",
            "n_cells": "background_n_cells",
            "n_classes": "n_background_classes",
        }
    )
    merged = target.merge(bg, on=group_cols, how="inner")
    merged["rank_delta"] = merged["candidate_rank"] - merged["background_rank"]
    merged["mean_log1p_expression_delta"] = (
        merged["candidate_mean_log1p_expression"] - merged["background_mean_log1p_expression"]
    )
    merged["detection_delta"] = merged["candidate_detection"] - merged["background_detection"]
    merged["positive_delta"] = merged["rank_delta"].gt(0)

    meta = ortho[["canonical_gene", "human_symbol", "mouse_symbol", "db_class_key"]].drop_duplicates("canonical_gene")
    merged = merged.merge(meta, on="canonical_gene", how="left")
    merged["gene"] = merged["human_symbol"].fillna(merged["source_gene"])
    cols = [
        "screen",
        "dataset",
        "sample",
        "core_branch",
        "branch_tested",
        "gene",
        "canonical_gene",
        "human_symbol",
        "mouse_symbol",
        "db_class_key",
        "candidate_rank",
        "background_rank",
        "rank_delta",
        "candidate_mean_log1p_expression",
        "background_mean_log1p_expression",
        "mean_log1p_expression_delta",
        "candidate_detection",
        "background_detection",
        "detection_delta",
        "candidate_n_cells",
        "background_n_cells",
        "n_background_classes",
        "positive_delta",
    ]
    return merged[cols].copy()


def summarize_branches(deltas: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["canonical_gene", "screen", "branch_tested"]
    meta_cols = ["gene", "human_symbol", "mouse_symbol", "db_class_key"]

    unit = (
        deltas.groupby(group_cols, sort=False)
        .agg(
            **{col: (col, "first") for col in meta_cols},
            n_units=("rank_delta", "size"),
            n_positive_units=("positive_delta", "sum"),
            median_unit_delta=("rank_delta", "median"),
            min_unit_delta=("rank_delta", "min"),
            median_candidate_detection=("candidate_detection", "median"),
            median_detection_delta=("detection_delta", "median"),
        )
        .reset_index()
    )
    unit["n_positive_units"] = unit["n_positive_units"].astype(int)
    unit["positive_unit_fraction"] = unit["n_positive_units"] / unit["n_units"]
    unit["unit_sign_p_greater"] = [
        sign_p_greater(pos, total) for pos, total in zip(unit["n_positive_units"], unit["n_units"], strict=False)
    ]

    dataset_delta = (
        deltas.groupby(group_cols + ["dataset"], sort=False)
        .agg(
            dataset_rank_delta=("rank_delta", "median"),
            dataset_detection_delta=("detection_delta", "median"),
            dataset_candidate_detection=("candidate_detection", "median"),
        )
        .reset_index()
    )
    dataset_delta["dataset_positive"] = dataset_delta["dataset_rank_delta"].gt(0)
    dataset = (
        dataset_delta.groupby(group_cols, sort=False)
        .agg(
            n_datasets=("dataset", "nunique"),
            n_positive_datasets=("dataset_positive", "sum"),
            median_dataset_delta=("dataset_rank_delta", "median"),
            min_dataset_delta=("dataset_rank_delta", "min"),
            median_dataset_detection_delta=("dataset_detection_delta", "median"),
            median_dataset_candidate_detection=("dataset_candidate_detection", "median"),
        )
        .reset_index()
    )
    dataset["n_positive_datasets"] = dataset["n_positive_datasets"].astype(int)
    dataset["positive_dataset_fraction"] = dataset["n_positive_datasets"] / dataset["n_datasets"]
    dataset["dataset_sign_p_greater"] = [
        sign_p_greater(pos, total)
        for pos, total in zip(dataset["n_positive_datasets"], dataset["n_datasets"], strict=False)
    ]
    dataset["dataset_wilcoxon_p_greater"] = (
        dataset_delta.groupby(group_cols, sort=False)["dataset_rank_delta"].apply(wilcoxon_p_greater).to_numpy()
    )
    positive_lists = dataset_delta.loc[dataset_delta["dataset_positive"]].groupby(group_cols, sort=False)["dataset"].apply(
        lambda values: ",".join(sorted(map(str, values)))
    )
    nonpositive_lists = dataset_delta.loc[~dataset_delta["dataset_positive"]].groupby(group_cols, sort=False)[
        "dataset"
    ].apply(lambda values: ",".join(sorted(map(str, values))))
    list_index = pd.MultiIndex.from_frame(dataset[group_cols])
    dataset["datasets_positive"] = positive_lists.reindex(list_index).fillna("").to_numpy()
    dataset["datasets_nonpositive"] = nonpositive_lists.reindex(list_index).fillna("").to_numpy()

    summary = unit.merge(dataset, on=group_cols, how="left")
    summary["branch_meta_support"] = (
        summary["n_datasets"].ge(MIN_DATASETS_PER_BRANCH)
        & summary["positive_dataset_fraction"].ge(SUPPORT_POSITIVE_DATASET_FRACTION)
        & summary["median_dataset_delta"].gt(0)
    )
    summary["branch_meta_strict"] = summary["branch_meta_support"] & summary["dataset_sign_p_greater"].le(
        STRICT_DATASET_SIGN_P
    )
    summary["dataset_sign_p_adj_bh"] = np.nan
    summary["unit_sign_p_adj_bh"] = np.nan
    for (_, _), idx in summary.groupby(["screen", "branch_tested"]).groups.items():
        summary.loc[idx, "dataset_sign_p_adj_bh"] = bh_adjust(summary.loc[idx, "dataset_sign_p_greater"])
        summary.loc[idx, "unit_sign_p_adj_bh"] = bh_adjust(summary.loc[idx, "unit_sign_p_greater"])
    return summary.sort_values(
        ["branch_meta_strict", "branch_meta_support", "median_dataset_delta", "positive_dataset_fraction"],
        ascending=[False, False, False, False],
    )


def load_annotations() -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    if SELECTED_TRIAGE.exists():
        selected = pd.read_csv(
            SELECTED_TRIAGE,
            sep="\t",
            usecols=[
                "canonical_gene",
                "mechanism_class",
                "manuscript_use",
                "mechanism_priority_score",
                "is_original_candidate_gene",
            ],
        ).drop_duplicates("canonical_gene")
        selected = selected.rename(
            columns={
                "mechanism_class": "selected_mechanism_class",
                "manuscript_use": "selected_manuscript_use",
                "mechanism_priority_score": "selected_mechanism_priority_score",
                "is_original_candidate_gene": "selected_is_original_candidate_gene",
            }
        )
        pieces.append(selected)
    if GENOME_TRIAGE.exists():
        genome = pd.read_csv(
            GENOME_TRIAGE,
            sep="\t",
            usecols=[
                "canonical_gene",
                "mechanism_class",
                "manuscript_use",
                "mechanism_priority_score",
                "is_original_candidate_gene",
            ],
        ).drop_duplicates("canonical_gene")
        genome = genome.rename(
            columns={
                "mechanism_class": "genome_mechanism_class",
                "manuscript_use": "genome_manuscript_use",
                "mechanism_priority_score": "genome_mechanism_priority_score",
                "is_original_candidate_gene": "genome_is_original_candidate_gene",
            }
        )
        pieces.append(genome)
    if CONSENSUS.exists():
        consensus = pd.read_csv(
            CONSENSUS,
            sep="\t",
            usecols=["canonical_gene", "consensus_tier", "combined_priority_score"],
        ).drop_duplicates("canonical_gene")
        pieces.append(consensus)
    if CONSENSUS_VALIDATION.exists():
        validation = pd.read_csv(
            CONSENSUS_VALIDATION,
            sep="\t",
            usecols=["canonical_gene", "robust_all_available_branches", "n_robust_screen_branches"],
        ).drop_duplicates("canonical_gene")
        validation = validation.rename(
            columns={
                "robust_all_available_branches": "consensus_candidate_dataset_robust_all_available",
                "n_robust_screen_branches": "consensus_candidate_n_robust_screen_branches",
            }
        )
        pieces.append(validation)

    if not pieces:
        return pd.DataFrame(columns=["canonical_gene"])
    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece, on="canonical_gene", how="outer")
    return merged


def build_gene_summary(branch_summary: pd.DataFrame, ortho: pd.DataFrame) -> pd.DataFrame:
    value_cols = [
        "n_units",
        "n_datasets",
        "positive_dataset_fraction",
        "median_dataset_delta",
        "min_dataset_delta",
        "median_dataset_candidate_detection",
        "dataset_sign_p_greater",
        "dataset_sign_p_adj_bh",
        "branch_meta_support",
        "branch_meta_strict",
    ]
    wide = branch_summary.pivot_table(
        index=["canonical_gene", "gene", "human_symbol", "mouse_symbol", "db_class_key"],
        columns=["screen", "branch_tested"],
        values=value_cols,
        aggfunc="first",
    ).reset_index()
    wide.columns = [
        "_".join(str(part) for part in col if str(part) != "") if isinstance(col, tuple) else col for col in wide.columns
    ]

    bool_cols = [col for col in wide.columns if col.startswith("branch_meta_")]
    for col in bool_cols:
        wide[col] = bool_col(wide, col)

    support_cols = [col for col in wide.columns if col.startswith("branch_meta_support_")]
    strict_cols = [col for col in wide.columns if col.startswith("branch_meta_strict_")]
    available_cols = [col for col in wide.columns if col.startswith("n_datasets_")]
    median_delta_cols = [col for col in wide.columns if col.startswith("median_dataset_delta_")]
    detection_cols = [col for col in wide.columns if col.startswith("median_dataset_candidate_detection_")]

    wide["n_supported_screen_branches"] = wide[support_cols].sum(axis=1) if support_cols else 0
    wide["n_strict_screen_branches"] = wide[strict_cols].sum(axis=1) if strict_cols else 0
    wide["n_available_screen_branches"] = (
        wide[available_cols].notna().sum(axis=1) if available_cols else 0
    )
    wide["sum_positive_median_dataset_delta"] = (
        wide[median_delta_cols].clip(lower=0).sum(axis=1) if median_delta_cols else 0
    )
    wide["minimum_available_candidate_detection"] = (
        wide[detection_cols].min(axis=1) if detection_cols else np.nan
    )

    def col(name: str) -> pd.Series:
        return wide[name] if name in wide.columns else pd.Series(False, index=wide.index)

    selected_shared = col("branch_meta_support_selected_dentate") & col("branch_meta_support_selected_cerebellar")
    genome_shared = col("branch_meta_support_full_matrix_dentate") & col("branch_meta_support_full_matrix_cerebellar")
    selected_strict = col("branch_meta_strict_selected_dentate") & col("branch_meta_strict_selected_cerebellar")
    genome_strict = col("branch_meta_strict_full_matrix_dentate") & col("branch_meta_strict_full_matrix_cerebellar")
    wide["shared_meta_support_selected"] = selected_shared
    wide["shared_meta_support_full_matrix"] = genome_shared
    wide["shared_meta_strict_selected"] = selected_strict
    wide["shared_meta_strict_full_matrix"] = genome_strict
    wide["shared_meta_support_both_screens"] = selected_shared & genome_shared
    wide["shared_meta_strict_both_screens"] = selected_strict & genome_strict

    conditions = [
        wide["shared_meta_strict_both_screens"],
        wide["shared_meta_support_both_screens"],
        genome_strict,
        selected_strict,
        genome_shared,
        selected_shared,
    ]
    tiers = [
        "strict_shared_both_screens",
        "supported_shared_both_screens",
        "strict_shared_full_matrix_only",
        "strict_shared_selected_only",
        "supported_shared_full_matrix_only",
        "supported_shared_selected_only",
    ]
    wide["ortholog_meta_tier"] = np.select(conditions, tiers, default="not_shared_or_incomplete")
    wide["ortholog_meta_priority_score"] = (
        wide["n_strict_screen_branches"] * 2.0
        + wide["n_supported_screen_branches"]
        + wide["sum_positive_median_dataset_delta"].fillna(0)
        + wide["minimum_available_candidate_detection"].fillna(0)
    )

    annotations = load_annotations()
    if not annotations.empty:
        wide = wide.merge(annotations, on="canonical_gene", how="left")

    ortho_cols = [
        "canonical_gene",
        "mgi_one_to_one_human_mouse",
        "same_canonical_symbol",
        "strict_same_symbol_one_to_one",
        "human_entrez_id",
        "mouse_entrez_id",
        "human_hgnc_id",
        "mouse_mgi_id",
    ]
    wide = wide.merge(ortho[ortho_cols], on="canonical_gene", how="left")
    return wide.sort_values(
        [
            "shared_meta_strict_both_screens",
            "shared_meta_support_both_screens",
            "n_strict_screen_branches",
            "n_supported_screen_branches",
            "ortholog_meta_priority_score",
        ],
        ascending=[False, False, False, False, False],
    )


def mechanism_hit_mask(df: pd.DataFrame) -> pd.Series:
    return (
        df.get("selected_manuscript_use", pd.Series("", index=df.index)).eq("mechanism_figure_candidate")
        | df.get("genome_manuscript_use", pd.Series("", index=df.index)).eq("mechanism_figure_candidate")
        | df.get("consensus_tier", pd.Series("", index=df.index)).eq("consensus_figure_candidate")
        | bool_col(df, "consensus_candidate_dataset_robust_all_available")
    )


def add_mechanism_priority(gene_summary: pd.DataFrame) -> pd.DataFrame:
    out = gene_summary.loc[gene_summary["ortholog_meta_tier"].ne("not_shared_or_incomplete")].copy()
    if out.empty:
        out["mechanism_hit_tier"] = []
        return out
    robust = bool_col(out, "consensus_candidate_dataset_robust_all_available")
    consensus_figure = out.get("consensus_tier", pd.Series("", index=out.index)).eq("consensus_figure_candidate")
    selected_figure = out.get("selected_manuscript_use", pd.Series("", index=out.index)).eq(
        "mechanism_figure_candidate"
    )
    genome_figure = out.get("genome_manuscript_use", pd.Series("", index=out.index)).eq("mechanism_figure_candidate")
    out["mechanism_hit_tier"] = np.select(
        [
            robust,
            consensus_figure,
            selected_figure & genome_figure,
            selected_figure,
            genome_figure,
        ],
        [
            "dataset_robust_consensus_figure",
            "consensus_figure_candidate",
            "dual_screen_mechanism_figure",
            "selected_screen_mechanism_figure",
            "full_matrix_mechanism_figure",
        ],
        default="not_mechanism_prioritized",
    )
    priority = {
        "dataset_robust_consensus_figure": 0,
        "consensus_figure_candidate": 1,
        "dual_screen_mechanism_figure": 2,
        "selected_screen_mechanism_figure": 3,
        "full_matrix_mechanism_figure": 4,
        "not_mechanism_prioritized": 9,
    }
    out["mechanism_hit_rank"] = out["mechanism_hit_tier"].map(priority).fillna(9).astype(int)
    return out


def mechanism_hits_from_gene_summary(gene_summary: pd.DataFrame) -> pd.DataFrame:
    prioritized = add_mechanism_priority(gene_summary)
    hits = prioritized.loc[mechanism_hit_mask(prioritized)].copy()
    return hits.sort_values(
        [
            "mechanism_hit_rank",
            "shared_meta_strict_both_screens",
            "shared_meta_support_both_screens",
            "n_strict_screen_branches",
            "n_supported_screen_branches",
            "ortholog_meta_priority_score",
        ],
        ascending=[True, False, False, False, False, False],
    )


def plot_top_hits(gene_summary: pd.DataFrame) -> None:
    plot_df = add_mechanism_priority(gene_summary)
    if plot_df.empty:
        return
    mechanism = mechanism_hits_from_gene_summary(gene_summary)
    if not mechanism.empty:
        plot_df = mechanism
    plot_df = plot_df.sort_values(
        [
            "mechanism_hit_rank",
            "shared_meta_strict_both_screens",
            "shared_meta_support_both_screens",
            "n_strict_screen_branches",
            "n_supported_screen_branches",
            "ortholog_meta_priority_score",
        ],
        ascending=[True, False, False, False, False, False],
    ).head(32)

    cols = [
        ("selected dentate", "median_dataset_delta_selected_dentate"),
        ("selected cerebellar", "median_dataset_delta_selected_cerebellar"),
        ("full dentate", "median_dataset_delta_full_matrix_dentate"),
        ("full cerebellar", "median_dataset_delta_full_matrix_cerebellar"),
    ]
    data = []
    labels = []
    for _, row in plot_df.iterrows():
        labels.append(row["gene"])
        data.append([row.get(col, np.nan) for _, col in cols])
    arr = np.array(data, dtype=float)

    fig_h = max(6.0, 0.32 * len(labels) + 1.8)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    im = ax.imshow(arr, aspect="auto", cmap="PiYG", vmin=-0.5, vmax=0.5)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([label for label, _ in cols], rotation=30, ha="right")
    ax.set_title("MGI one_to_one same-symbol ortholog meta-model")
    ax.set_xlabel("Screen and branch")
    ax.set_ylabel("Gene")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            value = arr[i, j]
            text = "NA" if np.isnan(value) else f"{value:.2f}"
            color = "#f6f6f6" if not np.isnan(value) and abs(value) > 0.28 else "#202020"
            ax.text(j, i, text, ha="center", va="center", fontsize=7, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Median dataset rank delta")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(
    *,
    ortho: pd.DataFrame,
    mgi_summary: dict[str, int],
    selected_expr: pd.DataFrame,
    genome_expr: pd.DataFrame,
    deltas: pd.DataFrame,
    branch_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    hits: pd.DataFrame,
    mechanism_hits: pd.DataFrame,
) -> None:
    strict_both = hits.loc[hits["ortholog_meta_tier"].eq("strict_shared_both_screens")]
    supported_both = hits.loc[hits["ortholog_meta_tier"].eq("supported_shared_both_screens")]
    full_only = hits.loc[hits["ortholog_meta_tier"].str.contains("full_matrix_only", na=False)]
    selected_only = hits.loc[hits["ortholog_meta_tier"].str.contains("selected_only", na=False)]
    top = hits.head(20)
    mechanism_top = mechanism_hits.head(24)
    robust_consensus = gene_summary.loc[bool_col(gene_summary, "consensus_candidate_dataset_robust_all_available")]

    lines = [
        "# Primary-Core MGI Ortholog Meta-Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis adds a conservative ortholog-aware layer to the selected-feature and full-matrix pseudobulk screens. It uses the official MGI human-mouse homology report, keeps one_to_one human-mouse classes, and further restricts the strict model to classes where the human and mouse symbols have the same canonical symbol.",
        "",
        "Because the current local matrix extraction was done in a same-symbol frame, non-identical one_to_one orthologs are intentionally deferred until a mouse-symbol-aware extraction is built.",
        "",
        "## Ortholog Scope",
        "",
        f"- MGI report rows: {mgi_summary['mgi_rows']:,}.",
        f"- MGI human-mouse homology classes: {mgi_summary['mgi_classes_with_human_mouse']:,}.",
        f"- One_to_one human-mouse pairs: {mgi_summary['one_to_one_pairs']:,}.",
        f"- Strict same-symbol one_to_one pairs: {mgi_summary['strict_same_symbol_one_to_one_pairs']:,}.",
        f"- Strict pairs represented in selected-feature expression rows: {selected_expr['canonical_gene'].nunique():,}.",
        f"- Strict pairs represented in full-matrix expression rows: {genome_expr['canonical_gene'].nunique():,}.",
        "",
        "## Meta-Model Summary",
        "",
        f"- Unit delta rows: {len(deltas):,}.",
        f"- Branch summary rows: {len(branch_summary):,}.",
        f"- Gene summary rows: {len(gene_summary):,}.",
        f"- Shared strict both-screen hits: {len(strict_both):,}.",
        f"- Shared supported both-screen hits: {len(supported_both):,}.",
        f"- Shared full-matrix-only hits: {len(full_only):,}.",
        f"- Shared selected-only hits: {len(selected_only):,}.",
        f"- Mechanism-prioritized shared hits: {len(mechanism_hits):,}.",
        "",
        "A branch is supported when at least two datasets contribute, at least 75% of datasets have positive candidate-versus-background deltas, and the median dataset delta is positive. A branch is strict when it also passes the dataset-level sign-test threshold p<=0.25, which is intentionally permissive because some branches have only two or three independent datasets.",
        "",
        "## Top Shared Ortholog Hits",
        "",
    ]
    for _, row in top.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['ortholog_meta_tier']}): "
            f"{int(row['n_supported_screen_branches'])}/{int(row['n_available_screen_branches'])} supported screen/branches, "
            f"score {row['ortholog_meta_priority_score']:.2f}."
        )

    if not mechanism_top.empty:
        lines.extend(["", "## Mechanism-Prioritized Hits", ""])
        for _, row in mechanism_top.iterrows():
            mechanism_class = row.get("genome_mechanism_class") or row.get("selected_mechanism_class") or "unclassified"
            lines.append(
                f"- `{row['gene']}` ({row['mechanism_hit_tier']}; {mechanism_class}): "
                f"{int(row['n_supported_screen_branches'])}/{int(row['n_available_screen_branches'])} supported screen/branches."
            )

    if not robust_consensus.empty:
        lines.extend(["", "## Consensus Candidate Check", ""])
        genes = ", ".join(f"`{g}`" for g in robust_consensus["gene"].head(20))
        lines.append(
            "The dataset-robust six-gene consensus shortlist remains inside this strict MGI one_to_one same-symbol frame: "
            f"{genes}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This is stronger than the same-symbol screen alone because it removes many-to-many and non-MGI-supported symbol matches.",
            "- It is still not the final ortholog DE model because non-identical one_to_one orthologs are absent from the current same-symbol matrix extraction.",
            "- The strongest manuscript-ready tier should come from genes that are supported in both selected-feature and full-matrix screens, and especially from genes that also passed the 24-candidate dataset-aware validation.",
            "- Rat is not represented in the MGI report used here; a rat extension should be added only when rat primary datasets enter the core.",
            "",
            "## Outputs",
            "",
            f"- Ortholog map: `{rel(OUT_ORTHOLOG_MAP)}`",
            f"- Unit deltas: `{rel(OUT_DELTAS)}`",
            f"- Branch summary: `{rel(OUT_BRANCH)}`",
            f"- Gene summary: `{rel(OUT_GENE)}`",
            f"- Shared hits: `{rel(OUT_HITS)}`",
            f"- Mechanism-prioritized hits: `{rel(OUT_MECHANISM_HITS)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    ortho, mgi_summary = load_mgi_ortholog_map()
    ortho.to_csv(OUT_ORTHOLOG_MAP, sep="\t", index=False)
    strict_genes = set(ortho.loc[ortho["strict_same_symbol_one_to_one"], "canonical_gene"])

    print(f"Strict same-symbol one_to_one pairs: {len(strict_genes):,}", flush=True)
    selected_expr = load_expression(SELECTED_EXPR, "selected", strict_genes)
    print(f"Selected expression genes: {selected_expr['canonical_gene'].nunique():,}", flush=True)
    genome_expr = load_expression(GENOME_EXPR, "full_matrix", strict_genes)
    print(f"Full-matrix expression genes: {genome_expr['canonical_gene'].nunique():,}", flush=True)
    expr = pd.concat([selected_expr, genome_expr], ignore_index=True, sort=False)
    deltas = build_unit_deltas(expr, ortho)
    print(f"Built unit deltas: {len(deltas):,}", flush=True)
    branch_summary = summarize_branches(deltas)
    print(f"Built branch summaries: {len(branch_summary):,}", flush=True)
    gene_summary = build_gene_summary(branch_summary, ortho)
    hits = gene_summary.loc[gene_summary["ortholog_meta_tier"].ne("not_shared_or_incomplete")].copy()
    hits = hits.sort_values(
        [
            "shared_meta_strict_both_screens",
            "shared_meta_support_both_screens",
            "n_strict_screen_branches",
            "n_supported_screen_branches",
            "ortholog_meta_priority_score",
        ],
        ascending=[False, False, False, False, False],
    )
    mechanism_hits = mechanism_hits_from_gene_summary(gene_summary)

    deltas.to_csv(OUT_DELTAS, sep="\t", index=False, compression="gzip")
    branch_summary.to_csv(OUT_BRANCH, sep="\t", index=False)
    gene_summary.to_csv(OUT_GENE, sep="\t", index=False)
    hits.to_csv(OUT_HITS, sep="\t", index=False)
    mechanism_hits.to_csv(OUT_MECHANISM_HITS, sep="\t", index=False)
    plot_top_hits(gene_summary)
    write_report(
        ortho=ortho,
        mgi_summary=mgi_summary,
        selected_expr=selected_expr,
        genome_expr=genome_expr,
        deltas=deltas,
        branch_summary=branch_summary,
        gene_summary=gene_summary,
        hits=hits,
        mechanism_hits=mechanism_hits,
    )

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Strict same-symbol one_to_one pairs: {len(strict_genes):,}")
    print(f"Unit deltas: {len(deltas):,}")
    print(f"Shared hits: {len(hits):,}")


if __name__ == "__main__":
    main()
