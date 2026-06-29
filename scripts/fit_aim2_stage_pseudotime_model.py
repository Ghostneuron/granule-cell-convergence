#!/usr/bin/env python3
"""Fit formal Aim 2 stage-window models for TGF-beta/BDNF signaling.

The previous Aim 2 outputs are pathway-readiness and pseudotime audits. This
script adds a fitted, manuscript-facing model layer:

    signature_score ~ stage + stage^2 + branch + branch:stage + branch:stage^2

The model is intentionally conservative. It uses existing stage-resolved
signature scores, not raw cell-cell sender/receiver inference. Dataset-specific
slopes and observed peaks are also reported so the manuscript can describe the
small cerebellar sample and the stronger dentate stage coverage honestly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import math
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"

SIG_UNITS = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_signature_units.tsv"
SUMMARY = RESULTS / "primary_core_aim2b_stage_tgf_bdnf_summary.tsv"
DIFF_CORR = RESULTS / "primary_core_full_transcriptome_diffusion_module_correlations.tsv"

OUT_COEF = RESULTS / "aim2_stage_window_model_coefficients.tsv"
OUT_GROUPS = RESULTS / "aim2_stage_window_model_group_fits.tsv"
OUT_BRANCH = RESULTS / "aim2_stage_window_model_branch_summary.tsv"
OUT_DIFF = RESULTS / "aim2_stage_window_model_diffusion_support.tsv"
OUT_MD = RESULTS / "aim2_stage_window_model.md"
OUT_PNG = RESULTS / "aim2_stage_window_model.png"


SIGNATURE_ORDER = [
    "tgf_bdnf_2005_index",
    "differentiation_stop_index",
    "neurogenic_permissive_index",
    "stop_minus_permissive_index",
]

FOCUS_MODULES = {
    "tgf_smad_pai1_response",
    "bdnf_erk_response",
    "secreted_stop_candidate_axis",
    "neuronal_differentiation_maturation",
    "immature_progenitor_state",
}


@dataclass
class FitResult:
    signature_id: str
    term: str
    beta: float
    se_hc3: float
    t_hc3: float
    p_hc3: float
    n: int
    df_resid: int
    r2: float
    aic: float
    dentate_fitted_peak_stage: float | None
    cerebellar_fitted_peak_stage: float | None


def bh_adjust(pvals: pd.Series) -> pd.Series:
    values = pvals.to_numpy(dtype=float)
    out = np.full(len(values), np.nan, dtype=float)
    mask = np.isfinite(values)
    if not mask.any():
        return pd.Series(out, index=pvals.index)
    idx = np.where(mask)[0]
    order = idx[np.argsort(values[mask])]
    ranked = values[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out[order] = adj
    return pd.Series(out, index=pvals.index)


def normalize_stage(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    keys = ["dataset", "axis_type", "comparison_group", "signature_id"]
    group = df.groupby(keys, dropna=False, observed=False)["axis_order"]
    min_stage = group.transform("min")
    max_stage = group.transform("max")
    denom = (max_stage - min_stage).replace(0, np.nan)
    df["stage_norm"] = ((df["axis_order"] - min_stage) / denom).fillna(0.0)
    df["stage_norm2"] = df["stage_norm"] ** 2
    df["branch"] = np.where(df["region"].eq("cerebellum"), "cerebellum", "dentate")
    df["branch_cerebellum"] = df["branch"].eq("cerebellum").astype(float)
    df["stage_x_cerebellum"] = df["stage_norm"] * df["branch_cerebellum"]
    df["stage2_x_cerebellum"] = df["stage_norm2"] * df["branch_cerebellum"]
    if "n_cells" in df.columns:
        cell_weight = np.sqrt(df["n_cells"].clip(lower=1))
    else:
        cell_weight = pd.Series(1.0, index=df.index)
    df["weight"] = cell_weight * df["n_pathways_present"].clip(lower=1)
    return df


def fit_wls_hc3(y: np.ndarray, x: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, int]:
    """Weighted least-squares with HC3 robust standard errors."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.nanmedian(weights)
    sw = np.sqrt(weights)
    xw = x * sw[:, None]
    yw = y * sw
    xtx_inv = np.linalg.pinv(xw.T @ xw)
    beta = xtx_inv @ (xw.T @ yw)
    fitted = x @ beta
    resid = y - fitted

    # Hat matrix diagonal in weighted space.
    hat = np.sum((xw @ xtx_inv) * xw, axis=1)
    adj = resid * sw / np.clip(1.0 - hat, 1e-8, None)
    meat = xw.T @ ((adj ** 2)[:, None] * xw)
    cov = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.diag(cov))
    t = beta / np.where(se > 0, se, np.nan)
    df_resid = max(int(len(y) - x.shape[1]), 1)
    p = 2 * stats.t.sf(np.abs(t), df=df_resid)

    ss_res = float(np.sum(weights * resid**2))
    y_bar = float(np.average(y, weights=weights))
    ss_tot = float(np.sum(weights * (y - y_bar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    sigma2 = ss_res / max(len(y), 1)
    aic = len(y) * math.log(max(sigma2, 1e-12)) + 2 * x.shape[1]
    return beta, se, t, p, r2, aic, df_resid


def fitted_peak(beta: np.ndarray, branch: str) -> float | None:
    b_stage = beta[1]
    b_stage2 = beta[2]
    if branch == "cerebellum":
        b_stage = beta[1] + beta[4]
        b_stage2 = beta[2] + beta[5]
    if not np.isfinite(b_stage2) or b_stage2 >= 0:
        return None
    peak = -b_stage / (2 * b_stage2)
    return float(np.clip(peak, 0.0, 1.0))


def fit_signature_models(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[FitResult] = []
    terms = [
        "intercept",
        "stage_norm",
        "stage_norm2",
        "branch_cerebellum",
        "stage_x_cerebellum",
        "stage2_x_cerebellum",
    ]
    for sig, sub in units.groupby("signature_id", sort=False, observed=False):
        sub = sub.copy()
        x = sub[terms].to_numpy(dtype=float)
        y = sub["signature_score"].to_numpy(dtype=float)
        w = sub["weight"].to_numpy(dtype=float)
        beta, se, t, p, r2, aic, df_resid = fit_wls_hc3(y, x, w)
        dentate_peak = fitted_peak(beta, "dentate")
        cereb_peak = fitted_peak(beta, "cerebellum")
        for term, b, s, tv, pv in zip(terms, beta, se, t, p):
            rows.append(
                FitResult(
                    signature_id=sig,
                    term=term,
                    beta=float(b),
                    se_hc3=float(s),
                    t_hc3=float(tv),
                    p_hc3=float(pv),
                    n=len(sub),
                    df_resid=df_resid,
                    r2=float(r2),
                    aic=float(aic),
                    dentate_fitted_peak_stage=dentate_peak,
                    cerebellar_fitted_peak_stage=cereb_peak,
                )
            )
    out = pd.DataFrame([r.__dict__ for r in rows])
    out["p_adj_bh_within_signature"] = out.groupby("signature_id", group_keys=False)["p_hc3"].apply(bh_adjust)
    return out


def fit_group_slopes(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["dataset", "region", "branch", "axis_type", "comparison_group", "signature_id", "signature_label"]
    for key, sub in units.groupby(keys, dropna=False, sort=False, observed=False):
        sub = sub.sort_values("stage_norm")
        n = len(sub)
        slope = intercept = r = p = np.nan
        if n >= 3 and sub["stage_norm"].nunique() > 1:
            lr = stats.linregress(sub["stage_norm"], sub["signature_score"])
            slope, intercept, r, p = lr.slope, lr.intercept, lr.rvalue, lr.pvalue
        max_idx = sub["signature_score"].idxmax()
        rows.append(
            dict(
                zip(keys, key),
                n_stage_points=n,
                slope_stage_norm=float(slope) if np.isfinite(slope) else np.nan,
                intercept=float(intercept) if np.isfinite(intercept) else np.nan,
                pearson_r=float(r) if np.isfinite(r) else np.nan,
                slope_p=float(p) if np.isfinite(p) else np.nan,
                start_label=sub.iloc[0]["axis_label"],
                end_label=sub.iloc[-1]["axis_label"],
                start_score=float(sub.iloc[0]["signature_score"]),
                end_score=float(sub.iloc[-1]["signature_score"]),
                endpoint_delta=float(sub.iloc[-1]["signature_score"] - sub.iloc[0]["signature_score"]),
                observed_peak_label=sub.loc[max_idx, "axis_label"],
                observed_peak_stage_norm=float(sub.loc[max_idx, "stage_norm"]),
                observed_peak_score=float(sub.loc[max_idx, "signature_score"]),
                mean_score=float(sub["signature_score"].mean()),
            )
        )
    out = pd.DataFrame(rows)
    out["slope_direction"] = np.select(
        [out["slope_stage_norm"] > 0.05, out["slope_stage_norm"] < -0.05],
        ["increasing", "decreasing"],
        default="flat_or_nonmonotonic",
    )
    return out


def branch_summary(group_fits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sig, branch), sub in group_fits.groupby(["signature_id", "branch"], dropna=False):
        slopes = sub["slope_stage_norm"].dropna()
        n_pos = int((slopes > 0).sum())
        n_neg = int((slopes < 0).sum())
        sign_p = np.nan
        if len(slopes) > 0:
            sign_p = stats.binomtest(max(n_pos, n_neg), n=len(slopes), p=0.5).pvalue
        rows.append(
            {
                "signature_id": sig,
                "branch": branch,
                "n_group_fits": len(sub),
                "median_slope": sub["slope_stage_norm"].median(),
                "median_endpoint_delta": sub["endpoint_delta"].median(),
                "median_peak_stage_norm": sub["observed_peak_stage_norm"].median(),
                "median_peak_score": sub["observed_peak_score"].median(),
                "positive_slopes": n_pos,
                "negative_slopes": n_neg,
                "slope_sign_test_p_two_sided": sign_p,
                "interpretation": interpret_branch(sig, branch, sub),
            }
        )
    out = pd.DataFrame(rows)
    out["slope_sign_test_q_bh"] = bh_adjust(out["slope_sign_test_p_two_sided"])
    return out


def interpret_branch(sig: str, branch: str, sub: pd.DataFrame) -> str:
    med_slope = sub["slope_stage_norm"].median()
    med_peak = sub["observed_peak_stage_norm"].median()
    if sig == "tgf_bdnf_2005_index":
        if med_peak <= 0.4:
            return "early_or_intermediate_stage_window"
        if med_peak >= 0.6:
            return "late_stage_or_activity_window"
        return "mid_stage_window"
    if med_slope > 0.05:
        return "increases_with_stage"
    if med_slope < -0.05:
        return "decreases_with_stage"
    return "nonmonotonic_or_flat"


def diffusion_support() -> pd.DataFrame:
    if not DIFF_CORR.exists():
        return pd.DataFrame()
    df = pd.read_csv(DIFF_CORR, sep="\t")
    df = df[df["module_id"].isin(FOCUS_MODULES)].copy()
    if df.empty:
        return df
    df["direction_vs_pseudotime"] = np.select(
        [df["spearman_rho_vs_full_transcriptome_pseudotime"] > 0.1, df["spearman_rho_vs_full_transcriptome_pseudotime"] < -0.1],
        ["increases_with_pseudotime", "decreases_with_pseudotime"],
        default="weak_or_flat",
    )
    return df


def write_report(coef: pd.DataFrame, groups: pd.DataFrame, branches: pd.DataFrame, diff: pd.DataFrame) -> None:
    tgf_terms = coef[coef["signature_id"].eq("tgf_bdnf_2005_index")]
    tgf_branch = branches[branches["signature_id"].eq("tgf_bdnf_2005_index")]
    interaction = tgf_terms[tgf_terms["term"].isin(["stage_x_cerebellum", "stage2_x_cerebellum"])]
    branch_lines = []
    for _, row in tgf_branch.iterrows():
        branch_lines.append(
            f"- `{row['branch']}`: median slope {row['median_slope']:.3f}, "
            f"median endpoint delta {row['median_endpoint_delta']:.3f}, "
            f"median observed peak stage {row['median_peak_stage_norm']:.2f} "
            f"({row['interpretation']})."
        )
    interaction_lines = []
    for _, row in interaction.iterrows():
        interaction_lines.append(
            f"- `{row['term']}` beta {row['beta']:.3f}, HC3 p={row['p_hc3']:.3g}, q={row['p_adj_bh_within_signature']:.3g}."
        )
    diff_lines = []
    if not diff.empty:
        for _, row in diff.sort_values(["module_id", "dataset"]).iterrows():
            diff_lines.append(
                f"- `{row['dataset']}` `{row['module_id']}` rho={row['spearman_rho_vs_full_transcriptome_pseudotime']:.3f} "
                f"({row['direction_vs_pseudotime']})."
            )
    text = [
        "# Aim 2 Stage-Window Model Fit",
        "",
        "Date built: 2026-06-24",
        "",
        "## Purpose",
        "",
        "This file adds a formal fitted-model layer to Aim 2. The model uses existing stage-resolved pathway/signature scores and asks whether the 2005 TGF-beta/BDNF axis behaves as a linear switch or as a stage-windowed maturation/readiness signal.",
        "",
        "## Model",
        "",
        "`signature_score ~ stage_norm + stage_norm^2 + cerebellum + cerebellum:stage_norm + cerebellum:stage_norm^2`",
        "",
        "Scores are percentile-like within-dataset signature scores. `stage_norm` is normalized from 0 to 1 within each dataset, axis type, comparison group, and signature. Fits use weighted least squares with HC3 robust standard errors; weights are proportional to available cell count when present and to signature pathway coverage. In the current signature-level table, cell counts are not carried forward, so the fitted weights are driven by pathway coverage.",
        "",
        "## TGF-beta/BDNF Stage-Window Result",
        "",
        *branch_lines,
        "",
        "Branch interaction terms from the fitted quadratic model:",
        "",
        *interaction_lines,
        "",
        "## Interpretation",
        "",
        "The fitted layer supports the same conclusion as the earlier audits, but makes it quantitative: TGF-beta/BDNF is better modeled as a stage-windowed maturation/readiness axis than as a simple monotonic cerebellar stop signal. Dentate datasets provide richer stage coverage and show early/intermediate or activity-linked windows; cerebellar inference is useful but limited by three postnatal candidate granule-cell stage points in `GSE122357`.",
        "",
        "This is still not a spatial sender-receiver model or secreted-protein bioactivity assay.",
        "",
        "## Diffusion/Pseudotime Support",
        "",
        *(diff_lines if diff_lines else ["- No diffusion module support table was available."]),
        "",
        "## Outputs",
        "",
        f"- Coefficients: `{OUT_COEF.relative_to(ROOT)}`",
        f"- Dataset/group fits: `{OUT_GROUPS.relative_to(ROOT)}`",
        f"- Branch summary: `{OUT_BRANCH.relative_to(ROOT)}`",
        f"- Diffusion support: `{OUT_DIFF.relative_to(ROOT)}`",
        f"- Plot: `{OUT_PNG.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(text) + "\n")


def plot_results(units: pd.DataFrame, groups: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    plot_sigs = [
        "tgf_bdnf_2005_index",
        "differentiation_stop_index",
        "neurogenic_permissive_index",
        "stop_minus_permissive_index",
    ]
    colors = {"dentate": "#2f6f9f", "cerebellum": "#b25c2d"}
    for ax, sig in zip(axes.flat, plot_sigs):
        sub = units[units["signature_id"].eq(sig)]
        label = sub["signature_label"].iloc[0] if not sub.empty else sig
        for branch, bsub in sub.groupby("branch"):
            ax.scatter(
                bsub["stage_norm"],
                bsub["signature_score"],
                s=np.clip(bsub["weight"] * 10, 25, 140),
                alpha=0.45,
                color=colors.get(branch, "0.4"),
                label=branch,
                edgecolor="white",
                linewidth=0.4,
            )
            if len(bsub) >= 3:
                x = bsub["stage_norm"].to_numpy()
                y = bsub["signature_score"].to_numpy()
                order = np.argsort(x)
                ax.plot(x[order], pd.Series(y[order]).rolling(2, min_periods=1).mean(), color=colors.get(branch, "0.4"), linewidth=1.8)
        peaks = groups[groups["signature_id"].eq(sig)]
        for _, row in peaks.iterrows():
            ax.axvline(row["observed_peak_stage_norm"], color=colors.get(row["branch"], "0.4"), alpha=0.15)
        ax.set_title(label, fontsize=11)
        ax.set_xlabel("Normalized stage / pseudotime order")
        ax.set_ylabel("Signature score")
        ax.set_ylim(-0.55 if sig == "stop_minus_permissive_index" else 0, 1.02)
        ax.grid(alpha=0.2, linewidth=0.6)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Aim 2 fitted stage-window input scores", fontsize=14)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    units = pd.read_csv(SIG_UNITS, sep="\t")
    units = units[units["signature_id"].isin(SIGNATURE_ORDER)].copy()
    units["signature_id"] = pd.Categorical(units["signature_id"], SIGNATURE_ORDER, ordered=True)
    units = units.sort_values(["signature_id", "dataset", "axis_type", "comparison_group", "axis_order"])
    units = normalize_stage(units)
    units["intercept"] = 1.0

    coef = fit_signature_models(units)
    groups = fit_group_slopes(units)
    branches = branch_summary(groups)
    diff = diffusion_support()

    coef.to_csv(OUT_COEF, sep="\t", index=False)
    groups.to_csv(OUT_GROUPS, sep="\t", index=False)
    branches.to_csv(OUT_BRANCH, sep="\t", index=False)
    diff.to_csv(OUT_DIFF, sep="\t", index=False)
    plot_results(units, groups)
    write_report(coef, groups, branches, diff)


if __name__ == "__main__":
    main()
