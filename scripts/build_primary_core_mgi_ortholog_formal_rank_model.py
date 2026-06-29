#!/usr/bin/env python3
"""Formal replication tests for the expanded MGI ortholog rank meta-model.

This is a dataset-level rank/pseudobulk validation layer, not raw-count DE.
It tests whether each gene has a positive granule-cell-versus-background
rank delta across independent datasets, then runs mixed/intercept models for
the smaller mechanism-prioritized candidate set.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import build_primary_core_mgi_ortholog_meta_model as meta


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

UNIT_DELTAS = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_unit_deltas.tsv.gz"
EXPANDED_GENE_SUMMARY = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_gene_summary.tsv"
EXPANDED_MECHANISM = RESULTS / "primary_core_mgi_ortholog_expanded_meta_model_mechanism_hits.tsv"

OUT_DATASET_DELTAS = RESULTS / "primary_core_mgi_ortholog_formal_rank_dataset_deltas.tsv.gz"
OUT_BRANCH_TESTS = RESULTS / "primary_core_mgi_ortholog_formal_rank_branch_tests.tsv"
OUT_GENE_SUMMARY = RESULTS / "primary_core_mgi_ortholog_formal_rank_gene_summary.tsv"
OUT_SHARED_HITS = RESULTS / "primary_core_mgi_ortholog_formal_rank_shared_hits.tsv"
OUT_MECHANISM_LONG = RESULTS / "primary_core_mgi_ortholog_formal_rank_mechanism_model_long.tsv"
OUT_MECHANISM = RESULTS / "primary_core_mgi_ortholog_formal_rank_mechanism_hits.tsv"
OUT_PLOT = RESULTS / "primary_core_mgi_ortholog_formal_rank_mechanism_hits.png"
OUT_MD = RESULTS / "primary_core_mgi_ortholog_formal_rank_model.md"

MIN_DATASETS = 2
REPLICATION_FRACTION = 0.75
NOMINAL_P = 0.25
FDR_Q = 0.10

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def rel(path: Path) -> str:
    return meta.rel(path)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    return meta.bool_col(df, col)


def one_sample_t_greater(values: pd.Series) -> tuple[float, float, float]:
    vals = numeric(values).dropna().to_numpy(dtype=float)
    if len(vals) < 2:
        return np.nan, np.nan, np.nan
    mean = float(np.mean(vals))
    sd = float(np.std(vals, ddof=1))
    if not np.isfinite(sd) or sd <= 1e-12:
        return np.nan, np.nan, np.nan
    t_stat = mean / (sd / np.sqrt(len(vals)))
    p_value = float(stats.t.sf(t_stat, len(vals) - 1))
    effect_size = mean / sd
    return float(t_stat), p_value, float(effect_size)


def sign_p_greater(values: pd.Series) -> float:
    vals = numeric(values).dropna()
    if vals.empty:
        return np.nan
    return meta.sign_p_greater(int(vals.gt(0).sum()), int(len(vals)))


def wilcoxon_p_greater(values: pd.Series) -> float:
    return meta.wilcoxon_p_greater(numeric(values))


def min_valid(values: list[float]) -> float:
    vals = [float(v) for v in values if pd.notna(v)]
    return float(min(vals)) if vals else np.nan


def fisher_p(values: list[float]) -> float:
    vals = [float(v) for v in values if pd.notna(v) and 0 <= float(v) <= 1]
    if len(vals) < 2:
        return np.nan
    clipped = np.clip(vals, 1e-300, 1.0)
    return float(stats.combine_pvalues(clipped, method="fisher").pvalue)


def load_unit_deltas() -> pd.DataFrame:
    usecols = [
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
        "rank_delta",
        "mean_log1p_expression_delta",
        "candidate_detection",
        "detection_delta",
        "candidate_n_cells",
        "background_n_cells",
        "positive_delta",
    ]
    df = pd.read_csv(UNIT_DELTAS, sep="\t", usecols=usecols, low_memory=False)
    for col in [
        "rank_delta",
        "mean_log1p_expression_delta",
        "candidate_detection",
        "detection_delta",
        "candidate_n_cells",
        "background_n_cells",
    ]:
        df[col] = numeric(df[col])
    df["positive_delta"] = df["rank_delta"].gt(0)
    return df


def build_dataset_deltas(unit_deltas: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "canonical_gene",
        "gene",
        "human_symbol",
        "mouse_symbol",
        "db_class_key",
        "screen",
        "branch_tested",
        "dataset",
    ]
    out = (
        unit_deltas.groupby(group_cols, sort=False)
        .agg(
            n_units=("rank_delta", "size"),
            n_positive_units=("positive_delta", "sum"),
            dataset_rank_delta=("rank_delta", "median"),
            dataset_mean_log1p_expression_delta=("mean_log1p_expression_delta", "median"),
            dataset_candidate_detection=("candidate_detection", "median"),
            dataset_detection_delta=("detection_delta", "median"),
            dataset_candidate_n_cells=("candidate_n_cells", "sum"),
            dataset_background_n_cells=("background_n_cells", "sum"),
            samples=("sample", lambda values: ",".join(sorted(set(map(str, values))))),
            core_branches=("core_branch", lambda values: ",".join(sorted(set(map(str, values))))),
        )
        .reset_index()
    )
    out["n_positive_units"] = out["n_positive_units"].astype(int)
    out["dataset_positive"] = out["dataset_rank_delta"].gt(0)
    return out


def build_branch_tests(dataset_deltas: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "canonical_gene",
        "gene",
        "human_symbol",
        "mouse_symbol",
        "db_class_key",
        "screen",
        "branch_tested",
    ]
    records: list[dict[str, object]] = []
    for keys, sub in dataset_deltas.groupby(group_cols, sort=False):
        values = numeric(sub["dataset_rank_delta"]).dropna()
        expr_values = numeric(sub["dataset_mean_log1p_expression_delta"]).dropna()
        detection_values = numeric(sub["dataset_detection_delta"]).dropna()
        n_datasets = int(values.size)
        n_positive = int(values.gt(0).sum())
        t_stat, t_p, standardized_effect = one_sample_t_greater(values)
        sign_p = sign_p_greater(values)
        wilcoxon_p = wilcoxon_p_greater(values)
        best_p = min_valid([t_p, sign_p, wilcoxon_p])
        record = dict(zip(group_cols, keys, strict=False))
        record.update(
            {
                "n_datasets": n_datasets,
                "n_positive_datasets": n_positive,
                "positive_dataset_fraction": float(n_positive / n_datasets) if n_datasets else np.nan,
                "n_units": int(sub["n_units"].sum()),
                "median_dataset_rank_delta": float(values.median()) if n_datasets else np.nan,
                "mean_dataset_rank_delta": float(values.mean()) if n_datasets else np.nan,
                "min_dataset_rank_delta": float(values.min()) if n_datasets else np.nan,
                "dataset_rank_delta_iqr": float(values.quantile(0.75) - values.quantile(0.25))
                if n_datasets
                else np.nan,
                "median_dataset_mean_log1p_expression_delta": float(expr_values.median())
                if not expr_values.empty
                else np.nan,
                "median_dataset_detection_delta": float(detection_values.median())
                if not detection_values.empty
                else np.nan,
                "median_dataset_candidate_detection": float(numeric(sub["dataset_candidate_detection"]).median()),
                "dataset_t_stat_greater": t_stat,
                "dataset_t_p_greater": t_p,
                "dataset_standardized_effect": standardized_effect,
                "dataset_sign_p_greater": sign_p,
                "dataset_wilcoxon_p_greater": wilcoxon_p,
                "dataset_best_nominal_p_greater": best_p,
                "datasets_positive": ",".join(sorted(map(str, sub.loc[sub["dataset_positive"], "dataset"]))),
                "datasets_nonpositive": ",".join(sorted(map(str, sub.loc[~sub["dataset_positive"], "dataset"]))),
            }
        )
        records.append(record)

    out = pd.DataFrame(records)
    out["branch_replication_support"] = (
        out["n_datasets"].ge(MIN_DATASETS)
        & out["median_dataset_rank_delta"].gt(0)
        & out["positive_dataset_fraction"].ge(REPLICATION_FRACTION)
    )
    out["branch_nominal_support"] = out["branch_replication_support"] & out[
        "dataset_best_nominal_p_greater"
    ].le(NOMINAL_P)

    for (_, _), idx in out.groupby(["screen", "branch_tested"]).groups.items():
        for p_col in [
            "dataset_t_p_greater",
            "dataset_sign_p_greater",
            "dataset_wilcoxon_p_greater",
            "dataset_best_nominal_p_greater",
        ]:
            out.loc[idx, p_col.replace("_p_", "_q_bh_")] = meta.bh_adjust(out.loc[idx, p_col])

    q_cols = [
        "dataset_t_q_bh_greater",
        "dataset_sign_q_bh_greater",
        "dataset_wilcoxon_q_bh_greater",
        "dataset_best_nominal_q_bh_greater",
    ]
    out["dataset_best_q_bh_greater"] = out[q_cols].min(axis=1, skipna=True)
    out["branch_fdr10_support"] = out["branch_replication_support"] & out["dataset_best_q_bh_greater"].le(FDR_Q)

    return out.sort_values(
        [
            "branch_fdr10_support",
            "branch_nominal_support",
            "branch_replication_support",
            "median_dataset_rank_delta",
            "positive_dataset_fraction",
            "dataset_best_nominal_p_greater",
        ],
        ascending=[False, False, False, False, False, True],
    )


def flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        "_".join(str(part) for part in col if str(part) != "") if isinstance(col, tuple) else col for col in df.columns
    ]
    return df


def wide_col(df: pd.DataFrame, col: str, default: object = False) -> pd.Series:
    if col in df.columns:
        return df[col]
    return pd.Series(default, index=df.index)


def add_shared_fisher_p(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()

    def branch_p(row: pd.Series, screen: str, branch: str) -> float:
        return row.get(f"dataset_best_nominal_p_greater_{screen}_{branch}", np.nan)

    selected_p: list[float] = []
    full_p: list[float] = []
    all_p: list[float] = []
    for _, row in out.iterrows():
        selected_vals = [branch_p(row, "selected", "dentate"), branch_p(row, "selected", "cerebellar")]
        full_vals = [branch_p(row, "full_matrix", "dentate"), branch_p(row, "full_matrix", "cerebellar")]
        selected_p.append(fisher_p(selected_vals))
        full_p.append(fisher_p(full_vals))
        all_p.append(fisher_p(selected_vals + full_vals))
    out["formal_fisher_p_selected_shared"] = selected_p
    out["formal_fisher_p_full_matrix_shared"] = full_p
    out["formal_fisher_p_all_available_branches"] = all_p
    out["formal_fisher_q_bh_selected_shared"] = meta.bh_adjust(out["formal_fisher_p_selected_shared"])
    out["formal_fisher_q_bh_full_matrix_shared"] = meta.bh_adjust(out["formal_fisher_p_full_matrix_shared"])
    out["formal_fisher_q_bh_all_available_branches"] = meta.bh_adjust(out["formal_fisher_p_all_available_branches"])
    return out


def build_gene_summary(branch_tests: pd.DataFrame) -> pd.DataFrame:
    value_cols = [
        "n_datasets",
        "n_units",
        "positive_dataset_fraction",
        "median_dataset_rank_delta",
        "mean_dataset_rank_delta",
        "min_dataset_rank_delta",
        "median_dataset_candidate_detection",
        "median_dataset_mean_log1p_expression_delta",
        "median_dataset_detection_delta",
        "dataset_t_p_greater",
        "dataset_t_q_bh_greater",
        "dataset_sign_p_greater",
        "dataset_sign_q_bh_greater",
        "dataset_wilcoxon_p_greater",
        "dataset_wilcoxon_q_bh_greater",
        "dataset_best_nominal_p_greater",
        "dataset_best_nominal_q_bh_greater",
        "dataset_best_q_bh_greater",
        "branch_replication_support",
        "branch_nominal_support",
        "branch_fdr10_support",
    ]
    wide = branch_tests.pivot_table(
        index=["canonical_gene", "gene", "human_symbol", "mouse_symbol", "db_class_key"],
        columns=["screen", "branch_tested"],
        values=value_cols,
        aggfunc="first",
    ).reset_index()
    wide = flatten_columns(wide)

    for col in [c for c in wide.columns if c.startswith("branch_")]:
        wide[col] = bool_col(wide, col)

    replication_cols = [c for c in wide.columns if c.startswith("branch_replication_support_")]
    nominal_cols = [c for c in wide.columns if c.startswith("branch_nominal_support_")]
    fdr_cols = [c for c in wide.columns if c.startswith("branch_fdr10_support_")]
    available_cols = [c for c in wide.columns if c.startswith("n_datasets_")]
    median_cols = [c for c in wide.columns if c.startswith("median_dataset_rank_delta_")]
    detection_cols = [c for c in wide.columns if c.startswith("median_dataset_candidate_detection_")]

    wide["formal_n_available_branches"] = wide[available_cols].notna().sum(axis=1) if available_cols else 0
    wide["formal_n_replication_branches"] = wide[replication_cols].sum(axis=1) if replication_cols else 0
    wide["formal_n_nominal_branches"] = wide[nominal_cols].sum(axis=1) if nominal_cols else 0
    wide["formal_n_fdr10_branches"] = wide[fdr_cols].sum(axis=1) if fdr_cols else 0
    wide["formal_sum_positive_median_rank_delta"] = (
        wide[median_cols].clip(lower=0).sum(axis=1) if median_cols else 0
    )
    wide["formal_min_candidate_detection"] = wide[detection_cols].min(axis=1) if detection_cols else np.nan

    for support_type in ["replication", "nominal", "fdr10"]:
        prefix = f"branch_{support_type}_support"
        selected = bool_col(wide, f"{prefix}_selected_dentate") & bool_col(
            wide, f"{prefix}_selected_cerebellar"
        )
        full = bool_col(wide, f"{prefix}_full_matrix_dentate") & bool_col(
            wide, f"{prefix}_full_matrix_cerebellar"
        )
        wide[f"formal_{support_type}_shared_selected"] = selected
        wide[f"formal_{support_type}_shared_full_matrix"] = full
        wide[f"formal_{support_type}_shared_any_screen"] = selected | full
        wide[f"formal_{support_type}_shared_both_screens"] = selected & full

    wide = add_shared_fisher_p(wide)

    conditions = [
        wide["formal_fdr10_shared_both_screens"],
        wide["formal_nominal_shared_both_screens"],
        wide["formal_replication_shared_both_screens"],
        wide["formal_fdr10_shared_full_matrix"],
        wide["formal_fdr10_shared_selected"],
        wide["formal_nominal_shared_full_matrix"],
        wide["formal_nominal_shared_selected"],
        wide["formal_replication_shared_full_matrix"],
        wide["formal_replication_shared_selected"],
    ]
    tiers = [
        "formal_fdr10_shared_both_screens",
        "formal_nominal_shared_both_screens",
        "formal_replication_shared_both_screens",
        "formal_fdr10_shared_full_matrix",
        "formal_fdr10_shared_selected",
        "formal_nominal_shared_full_matrix",
        "formal_nominal_shared_selected",
        "formal_replication_shared_full_matrix",
        "formal_replication_shared_selected",
    ]
    wide["formal_rank_tier"] = np.select(conditions, tiers, default="not_formally_shared")
    wide["formal_rank_priority_score"] = (
        wide["formal_n_fdr10_branches"] * 3.0
        + wide["formal_n_nominal_branches"] * 2.0
        + wide["formal_n_replication_branches"]
        + wide["formal_sum_positive_median_rank_delta"].fillna(0)
        + wide["formal_min_candidate_detection"].fillna(0)
    )

    expanded = pd.read_csv(EXPANDED_GENE_SUMMARY, sep="\t", low_memory=False)
    keep_cols = [
        "canonical_gene",
        "ortholog_meta_tier",
        "ortholog_symbol_class",
        "n_supported_screen_branches",
        "n_strict_screen_branches",
        "n_available_screen_branches",
        "consensus_candidate_dataset_robust_all_available",
        "consensus_candidate_n_robust_screen_branches",
        "consensus_tier",
        "combined_priority_score",
        "selected_mechanism_class",
        "selected_manuscript_use",
        "genome_mechanism_class",
        "genome_manuscript_use",
    ]
    keep_cols = [col for col in keep_cols if col in expanded.columns]
    wide = wide.merge(expanded[keep_cols].drop_duplicates("canonical_gene"), on="canonical_gene", how="left")

    return wide.sort_values(
        [
            "formal_fdr10_shared_both_screens",
            "formal_nominal_shared_both_screens",
            "formal_replication_shared_both_screens",
            "formal_n_fdr10_branches",
            "formal_n_nominal_branches",
            "formal_n_replication_branches",
            "formal_rank_priority_score",
        ],
        ascending=[False, False, False, False, False, False, False],
    )


def fit_intercept_model(sub: pd.DataFrame) -> dict[str, object]:
    clean = sub.loc[sub["rank_delta"].notna(), ["rank_delta", "dataset"]].copy()
    clean["dataset"] = clean["dataset"].astype(str)
    n_units = int(len(clean))
    n_datasets = int(clean["dataset"].nunique())
    if n_units < 3 or n_datasets < 2:
        return {
            "formal_model_method": "not_fit_too_few_units",
            "formal_model_estimate": np.nan,
            "formal_model_se": np.nan,
            "formal_model_z_greater": np.nan,
            "formal_model_p_greater": np.nan,
            "formal_model_converged": False,
            "formal_model_message": "too_few_units_or_datasets",
        }

    try:
        from statsmodels.regression.mixed_linear_model import MixedLM

        y = clean["rank_delta"].astype(float).to_numpy()
        exog = pd.DataFrame({"intercept": np.ones(n_units)})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = MixedLM(y, exog, groups=clean["dataset"]).fit(
                reml=False, method="lbfgs", maxiter=200, disp=False
            )
        estimate = float(result.params["intercept"])
        se = float(result.bse["intercept"])
        z_value = estimate / se if np.isfinite(se) and se > 0 else np.nan
        p_value = float(stats.norm.sf(z_value)) if pd.notna(z_value) else np.nan
        return {
            "formal_model_method": "mixedlm_dataset_random_intercept",
            "formal_model_estimate": estimate,
            "formal_model_se": se,
            "formal_model_z_greater": z_value,
            "formal_model_p_greater": p_value,
            "formal_model_converged": bool(getattr(result, "converged", False)),
            "formal_model_message": "",
        }
    except Exception as exc:
        try:
            import statsmodels.api as sm

            y = clean["rank_delta"].astype(float).to_numpy()
            x = np.ones((n_units, 1))
            fit = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": clean["dataset"]})
            estimate = float(fit.params[0])
            se = float(fit.bse[0])
            z_value = estimate / se if np.isfinite(se) and se > 0 else np.nan
            p_value = float(stats.norm.sf(z_value)) if pd.notna(z_value) else np.nan
            return {
                "formal_model_method": "cluster_robust_intercept_fallback",
                "formal_model_estimate": estimate,
                "formal_model_se": se,
                "formal_model_z_greater": z_value,
                "formal_model_p_greater": p_value,
                "formal_model_converged": True,
                "formal_model_message": str(exc)[:180],
            }
        except Exception as fallback_exc:
            return {
                "formal_model_method": "model_failed",
                "formal_model_estimate": np.nan,
                "formal_model_se": np.nan,
                "formal_model_z_greater": np.nan,
                "formal_model_p_greater": np.nan,
                "formal_model_converged": False,
                "formal_model_message": f"{exc}; fallback: {fallback_exc}"[:240],
            }


def build_mechanism_models(unit_deltas: pd.DataFrame, branch_tests: pd.DataFrame) -> pd.DataFrame:
    mechanism = pd.read_csv(EXPANDED_MECHANISM, sep="\t", low_memory=False)
    genes = set(mechanism["canonical_gene"].astype(str))
    sub = unit_deltas.loc[unit_deltas["canonical_gene"].astype(str).isin(genes)].copy()
    group_cols = [
        "canonical_gene",
        "gene",
        "human_symbol",
        "mouse_symbol",
        "db_class_key",
        "screen",
        "branch_tested",
    ]
    records: list[dict[str, object]] = []
    for keys, group in sub.groupby(group_cols, sort=False):
        record = dict(zip(group_cols, keys, strict=False))
        record.update(
            {
                "formal_model_n_units": int(len(group)),
                "formal_model_n_datasets": int(group["dataset"].nunique()),
            }
        )
        record.update(fit_intercept_model(group))
        records.append(record)
    out = pd.DataFrame(records)

    branch_cols = [
        "canonical_gene",
        "screen",
        "branch_tested",
        "dataset_best_nominal_p_greater",
        "dataset_best_q_bh_greater",
        "branch_replication_support",
        "branch_nominal_support",
        "branch_fdr10_support",
    ]
    out = out.merge(branch_tests[branch_cols], on=["canonical_gene", "screen", "branch_tested"], how="left")
    for (_, _), idx in out.groupby(["screen", "branch_tested"]).groups.items():
        out.loc[idx, "formal_model_q_bh_greater"] = meta.bh_adjust(out.loc[idx, "formal_model_p_greater"])
    out["formal_model_nominal_support"] = out["branch_replication_support"] & out["formal_model_p_greater"].le(
        NOMINAL_P
    )
    out["formal_model_fdr10_support"] = out["branch_replication_support"] & out["formal_model_q_bh_greater"].le(FDR_Q)
    return out.sort_values(
        [
            "formal_model_fdr10_support",
            "formal_model_nominal_support",
            "branch_nominal_support",
            "formal_model_p_greater",
        ],
        ascending=[False, False, False, True],
    )


def build_mechanism_summary(gene_summary: pd.DataFrame, mechanism_long: pd.DataFrame) -> pd.DataFrame:
    mechanism = pd.read_csv(EXPANDED_MECHANISM, sep="\t", low_memory=False)
    summary = mechanism.merge(gene_summary, on="canonical_gene", how="left", suffixes=("", "_formal"))

    value_cols = [
        "formal_model_estimate",
        "formal_model_se",
        "formal_model_p_greater",
        "formal_model_q_bh_greater",
        "formal_model_nominal_support",
        "formal_model_fdr10_support",
        "formal_model_method",
        "formal_model_n_units",
        "formal_model_n_datasets",
    ]
    wide = mechanism_long.pivot_table(
        index="canonical_gene",
        columns=["screen", "branch_tested"],
        values=value_cols,
        aggfunc="first",
    ).reset_index()
    wide = flatten_columns(wide)
    for col in [c for c in wide.columns if c.startswith("formal_model_") and c.endswith("_support")]:
        wide[col] = bool_col(wide, col)
    summary = summary.merge(wide, on="canonical_gene", how="left")

    model_nominal_cols = [c for c in summary.columns if c.startswith("formal_model_nominal_support_")]
    model_fdr_cols = [c for c in summary.columns if c.startswith("formal_model_fdr10_support_")]
    summary["formal_model_n_nominal_branches"] = (
        summary[model_nominal_cols].sum(axis=1) if model_nominal_cols else 0
    )
    summary["formal_model_n_fdr10_branches"] = summary[model_fdr_cols].sum(axis=1) if model_fdr_cols else 0

    return summary.sort_values(
        [
            "mechanism_hit_rank",
            "formal_fdr10_shared_both_screens",
            "formal_nominal_shared_both_screens",
            "formal_replication_shared_both_screens",
            "formal_model_n_fdr10_branches",
            "formal_model_n_nominal_branches",
            "formal_rank_priority_score",
        ],
        ascending=[True, False, False, False, False, False, False],
    )


def plot_mechanism_hits(mechanism_summary: pd.DataFrame) -> None:
    if mechanism_summary.empty:
        return
    plot_df = mechanism_summary.head(36).copy()
    cols = [
        ("selected dentate", "median_dataset_rank_delta_selected_dentate", "branch_nominal_support_selected_dentate"),
        (
            "selected cerebellar",
            "median_dataset_rank_delta_selected_cerebellar",
            "branch_nominal_support_selected_cerebellar",
        ),
        ("full dentate", "median_dataset_rank_delta_full_matrix_dentate", "branch_nominal_support_full_matrix_dentate"),
        (
            "full cerebellar",
            "median_dataset_rank_delta_full_matrix_cerebellar",
            "branch_nominal_support_full_matrix_cerebellar",
        ),
    ]
    labels = plot_df["gene"].astype(str).tolist()
    arr = np.array([[row.get(value_col, np.nan) for _, value_col, _ in cols] for _, row in plot_df.iterrows()])

    fig_h = max(6.0, 0.30 * len(labels) + 1.8)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    im = ax.imshow(arr, aspect="auto", cmap="PiYG", vmin=-0.5, vmax=0.5)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([label for label, _, _ in cols], rotation=30, ha="right")
    ax.set_title("Formal MGI ortholog rank-meta validation")
    ax.set_xlabel("Screen and branch")
    ax.set_ylabel("Mechanism-prioritized gene")
    for i, (_, row) in enumerate(plot_df.iterrows()):
        for j, (_, value_col, support_col) in enumerate(cols):
            value = row.get(value_col, np.nan)
            support = bool(row.get(support_col, False))
            suffix = "*" if support else ""
            text = "NA" if pd.isna(value) else f"{float(value):.2f}{suffix}"
            color = "#f6f6f6" if not pd.isna(value) and abs(float(value)) > 0.28 else "#202020"
            ax.text(j, i, text, ha="center", va="center", fontsize=7, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Median dataset rank delta")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(
    *,
    unit_deltas: pd.DataFrame,
    dataset_deltas: pd.DataFrame,
    branch_tests: pd.DataFrame,
    gene_summary: pd.DataFrame,
    shared_hits: pd.DataFrame,
    mechanism_long: pd.DataFrame,
    mechanism_summary: pd.DataFrame,
) -> None:
    replication_shared = shared_hits.loc[bool_col(shared_hits, "formal_replication_shared_any_screen")]
    nominal_shared = shared_hits.loc[bool_col(shared_hits, "formal_nominal_shared_any_screen")]
    fdr_shared = shared_hits.loc[bool_col(shared_hits, "formal_fdr10_shared_any_screen")]
    both_screen = shared_hits.loc[
        bool_col(shared_hits, "formal_replication_shared_both_screens")
        | bool_col(shared_hits, "formal_nominal_shared_both_screens")
        | bool_col(shared_hits, "formal_fdr10_shared_both_screens")
    ]
    robust = mechanism_summary.loc[
        bool_col(mechanism_summary, "consensus_candidate_dataset_robust_all_available")
    ].copy()

    model_methods = mechanism_long["formal_model_method"].value_counts(dropna=False).to_dict()
    method_text = ", ".join(f"{key}: {value}" for key, value in model_methods.items())

    lines = [
        "# Formal MGI Ortholog Rank-Meta Validation",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis adds a stricter statistical validation layer to the expanded MGI ortholog meta-model. It uses the existing pseudobulk rank deltas and tests whether granule-cell candidate classes are consistently above branch-specific background classes across independent datasets.",
        "",
        "This remains a rank-meta pseudobulk model, not raw-count DESeq2/edgeR differential expression. Its strength is cross-dataset replication; its limitation is small independent dataset count in several branches.",
        "",
        "## Input Scope",
        "",
        f"- Unit delta rows: {len(unit_deltas):,}.",
        f"- Dataset-level delta rows: {len(dataset_deltas):,}.",
        f"- Branch tests: {len(branch_tests):,}.",
        f"- Gene summaries: {len(gene_summary):,}.",
        f"- Mechanism-prioritized genes modeled: {mechanism_summary['canonical_gene'].nunique():,}.",
        "",
        "## Test Definitions",
        "",
        f"- Replication support: at least {MIN_DATASETS} datasets, median dataset rank delta > 0, and positive dataset fraction >= {REPLICATION_FRACTION:.2f}.",
        f"- Nominal support: replication support plus best one-sided dataset-level p <= {NOMINAL_P:.2f} from t, sign, or Wilcoxon tests.",
        f"- FDR10 support: replication support plus best within-screen/branch BH q <= {FDR_Q:.2f}.",
        "- Shared support requires both dentate and cerebellar branches to pass in the selected screen, full-matrix screen, or both.",
        "",
        "## Summary",
        "",
        f"- Branches with replication support: {int(branch_tests['branch_replication_support'].sum()):,}.",
        f"- Branches with nominal support: {int(branch_tests['branch_nominal_support'].sum()):,}.",
        f"- Branches with FDR10 support: {int(branch_tests['branch_fdr10_support'].sum()):,}.",
        f"- Formally shared hits: {len(shared_hits):,}.",
        f"- Replication-shared hits: {len(replication_shared):,}.",
        f"- Nominal-shared hits: {len(nominal_shared):,}.",
        f"- FDR10-shared hits: {len(fdr_shared):,}.",
        f"- Both-screen shared hits: {len(both_screen):,}.",
        f"- Mechanism model methods: {method_text}.",
        "",
        "## Mechanism-Prioritized Result",
        "",
    ]
    for _, row in mechanism_summary.head(30).iterrows():
        mechanism_class = row.get("genome_mechanism_class") or row.get("selected_mechanism_class") or "unclassified"
        lines.append(
            f"- `{row['gene']}` ({row.get('mechanism_hit_tier', 'mechanism_candidate')}; {mechanism_class}): "
            f"{row.get('formal_rank_tier', 'not_formally_shared')}, "
            f"{int(row.get('formal_n_nominal_branches', 0))}/"
            f"{int(row.get('formal_n_available_branches', 0))} nominal branches, "
            f"{int(row.get('formal_model_n_nominal_branches', 0))} model-supported branches."
        )

    if not robust.empty:
        genes = ", ".join(f"`{gene}`" for gene in robust["gene"].astype(str))
        lines.extend(
            [
                "",
                "## Consensus Six-Gene Check",
                "",
                f"- The dataset-robust consensus genes retained after formal rank-meta validation are: {genes}.",
            ]
        )

    lines.extend(["", "## Top Formal Shared Hits", ""])
    for _, row in shared_hits.head(25).iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['formal_rank_tier']}; {row.get('ortholog_symbol_class', 'ortholog')}): "
            f"{int(row['formal_n_nominal_branches'])}/{int(row['formal_n_available_branches'])} nominal branches, "
            f"score {row['formal_rank_priority_score']:.2f}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This formal layer supports the project strategy: the strongest genes are not just morphology-plausible, they replicate as dentate-plus-cerebellar granule-cell enriched across independent datasets.",
            "- The six-gene consensus (`GPM6A`, `NFIB`, `NFIA`, `KCNK1`, `RFX3`, `GABRA2`) remains the safest manuscript seed set.",
            "- Additional synaptic and wiring genes are useful as pathway/context support, while non-identical-symbol ortholog hits should remain secondary until raw-count or external validation is added.",
            "",
            "## Outputs",
            "",
            f"- Dataset deltas: `{rel(OUT_DATASET_DELTAS)}`",
            f"- Branch tests: `{rel(OUT_BRANCH_TESTS)}`",
            f"- Gene summary: `{rel(OUT_GENE_SUMMARY)}`",
            f"- Formal shared hits: `{rel(OUT_SHARED_HITS)}`",
            f"- Mechanism model long table: `{rel(OUT_MECHANISM_LONG)}`",
            f"- Mechanism summary: `{rel(OUT_MECHANISM)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    unit_deltas = load_unit_deltas()
    dataset_deltas = build_dataset_deltas(unit_deltas)
    branch_tests = build_branch_tests(dataset_deltas)
    gene_summary = build_gene_summary(branch_tests)
    shared_hits = gene_summary.loc[gene_summary["formal_rank_tier"].ne("not_formally_shared")].copy()
    shared_hits = shared_hits.sort_values(
        [
            "formal_fdr10_shared_both_screens",
            "formal_nominal_shared_both_screens",
            "formal_replication_shared_both_screens",
            "formal_n_fdr10_branches",
            "formal_n_nominal_branches",
            "formal_n_replication_branches",
            "formal_rank_priority_score",
        ],
        ascending=[False, False, False, False, False, False, False],
    )
    mechanism_long = build_mechanism_models(unit_deltas, branch_tests)
    mechanism_summary = build_mechanism_summary(gene_summary, mechanism_long)

    dataset_deltas.to_csv(OUT_DATASET_DELTAS, sep="\t", index=False, compression="gzip")
    branch_tests.to_csv(OUT_BRANCH_TESTS, sep="\t", index=False)
    gene_summary.to_csv(OUT_GENE_SUMMARY, sep="\t", index=False)
    shared_hits.to_csv(OUT_SHARED_HITS, sep="\t", index=False)
    mechanism_long.to_csv(OUT_MECHANISM_LONG, sep="\t", index=False)
    mechanism_summary.to_csv(OUT_MECHANISM, sep="\t", index=False)
    plot_mechanism_hits(mechanism_summary)
    write_report(
        unit_deltas=unit_deltas,
        dataset_deltas=dataset_deltas,
        branch_tests=branch_tests,
        gene_summary=gene_summary,
        shared_hits=shared_hits,
        mechanism_long=mechanism_long,
        mechanism_summary=mechanism_summary,
    )

    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Dataset deltas: {len(dataset_deltas):,}")
    print(f"Branch tests: {len(branch_tests):,}")
    print(f"Formal shared hits: {len(shared_hits):,}")
    print(f"Mechanism hits: {len(mechanism_summary):,}")


if __name__ == "__main__":
    main()
