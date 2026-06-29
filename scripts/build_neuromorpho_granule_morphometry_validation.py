#!/usr/bin/env python3
"""Build a first NeuroMorpho granule-cell morphometry validation layer.

The goal is not to exhaustively analyze every reconstruction. The public
dentate set is much larger than the cerebellar granule set, so this script
uses all available cerebellar granule hits under the current strict query and
a reproducible species-stratified dentate sample.
"""

from __future__ import annotations

import json
import math
import os
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CACHE = ROOT / "Project" / "dataset_search_cache" / "phys_morph" / "neuromorpho_validation"

OUT_META = RESULTS / "neuromorpho_granule_morphometry_metadata.tsv"
OUT_MORPH = RESULTS / "neuromorpho_granule_morphometry_sample.tsv"
OUT_SUMMARY = RESULTS / "neuromorpho_granule_morphometry_summary.tsv"
OUT_COMPARISON = RESULTS / "neuromorpho_granule_morphometry_comparison.tsv"
OUT_PRIORS = RESULTS / "neuromorpho_aim3_input_degree_priors.tsv"
OUT_PLOT = RESULTS / "neuromorpho_granule_morphometry_validation.png"
OUT_MD = RESULTS / "neuromorpho_granule_morphometry_validation.md"

BASE = "https://neuromorpho.org/api"
PAGE_SIZE = 50
MAX_WORKERS = 6
SLEEP_BETWEEN_PAGE_FETCHES = 0.05

GROUPS = [
    {
        "sampling_group": "DG_mouse_sample",
        "analysis_region": "dentate_gyrus",
        "query": "cell_type:granule",
        "filters": ["brain_region:dentate gyrus", "species:mouse"],
        "target_n": 250,
        "sample_mode": "even_pages",
    },
    {
        "sampling_group": "DG_rat_sample",
        "analysis_region": "dentate_gyrus",
        "query": "cell_type:granule",
        "filters": ["brain_region:dentate gyrus", "species:rat"],
        "target_n": 250,
        "sample_mode": "even_pages",
    },
    {
        "sampling_group": "DG_human_all",
        "analysis_region": "dentate_gyrus",
        "query": "cell_type:granule",
        "filters": ["brain_region:dentate gyrus", "species:human"],
        "target_n": None,
        "sample_mode": "all",
    },
    {
        "sampling_group": "CB_all_species_all",
        "analysis_region": "cerebellum",
        "query": "brain_region:cerebellum",
        "filters": ["cell_type:granule"],
        "target_n": None,
        "sample_mode": "all",
    },
]

MORPH_FIELDS = [
    "n_stems",
    "n_bifs",
    "n_branch",
    "length",
    "surface",
    "volume",
    "width",
    "height",
    "depth",
    "diameter",
    "eucDistance",
    "pathDistance",
    "branch_Order",
    "contraction",
    "fragmentation",
    "partition_asymmetry",
    "fractal_Dim",
    "soma_Surface",
]

SUMMARY_METRICS = [
    "n_stems",
    "n_bifs",
    "n_branch",
    "length",
    "surface",
    "volume",
    "pathDistance",
    "branch_Order",
    "fractal_Dim",
]


def slug(text: str) -> str:
    return (
        text.replace(" ", "_")
        .replace(":", "-")
        .replace("/", "_")
        .replace(";", "_")
        .replace(",", "_")
    )


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def api_get_json(url: str, cache_path: Path, retries: int = 3) -> dict[str, Any]:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text())
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read().decode("utf-8")
            cache_path.write_text(payload)
            return json.loads(payload)
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def select_url(query: str, filters: list[str], page: int, size: int) -> str:
    params: list[tuple[str, str | int]] = [("q", query)]
    params.extend(("fq", f) for f in filters)
    params.extend([("page", page), ("size", size)])
    return f"{BASE}/neuron/select?{urllib.parse.urlencode(params)}"


def get_total(group: dict[str, Any]) -> tuple[int, int]:
    cache_path = CACHE / "pages" / group["sampling_group"] / "page_0_size_1.json"
    data = api_get_json(select_url(group["query"], group["filters"], 0, 1), cache_path)
    if data.get("status") == 404:
        return 0, 0
    page = data.get("page") or {}
    total = int(page.get("totalElements", 0))
    total_pages = int(math.ceil(total / PAGE_SIZE)) if total else 0
    return total, total_pages


def page_numbers_for_group(total: int, target_n: int | None, sample_mode: str) -> list[int]:
    if total == 0:
        return []
    total_pages = int(math.ceil(total / PAGE_SIZE))
    if sample_mode == "all" or target_n is None or target_n >= total:
        return list(range(total_pages))
    n_pages = min(total_pages, int(math.ceil(target_n / PAGE_SIZE)))
    if n_pages <= 1:
        return [0]
    return sorted({int(round(v)) for v in np.linspace(0, total_pages - 1, n_pages)})


def resources_from_page(data: dict[str, Any]) -> list[dict[str, Any]]:
    embedded = data.get("_embedded") or {}
    return embedded.get("neuronResources") or []


def collect_metadata() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    group_stats: list[dict[str, Any]] = []
    for group in GROUPS:
        total, total_pages = get_total(group)
        pages = page_numbers_for_group(total, group["target_n"], group["sample_mode"])
        seen: set[int] = set()
        fetched = 0
        for page in pages:
            cache_path = CACHE / "pages" / group["sampling_group"] / f"page_{page}_size_{PAGE_SIZE}.json"
            data = api_get_json(select_url(group["query"], group["filters"], page, PAGE_SIZE), cache_path)
            time.sleep(SLEEP_BETWEEN_PAGE_FETCHES)
            for item in resources_from_page(data):
                neuron_id = int(item["neuron_id"])
                if neuron_id in seen:
                    continue
                seen.add(neuron_id)
                fetched += 1
                rows.append(
                    {
                        "analysis_region": group["analysis_region"],
                        "sampling_group": group["sampling_group"],
                        "neuron_id": neuron_id,
                        "neuron_name": item.get("neuron_name", ""),
                        "archive": item.get("archive", ""),
                        "species": item.get("species", ""),
                        "strain": item.get("strain", ""),
                        "brain_region": ";".join(item.get("brain_region") or []),
                        "cell_type": ";".join(item.get("cell_type") or []),
                        "age_classification": item.get("age_classification", ""),
                        "age_scale": item.get("age_scale", ""),
                        "min_age": item.get("min_age", ""),
                        "max_age": item.get("max_age", ""),
                        "gender": item.get("gender", ""),
                        "experiment_condition": ";".join(item.get("experiment_condition") or []),
                        "protocol": item.get("protocol", ""),
                        "domain": item.get("domain", ""),
                        "attributes": item.get("attributes", ""),
                        "physical_integrity": item.get("physical_Integrity", ""),
                        "reconstruction_software": item.get("reconstruction_software", ""),
                        "slicing_thickness": item.get("slicing_thickness", ""),
                        "reference_pmid": ";".join(item.get("reference_pmid") or []),
                        "reference_doi": ";".join(item.get("reference_doi") or []),
                        "measurements_url": ((item.get("_links") or {}).get("measurements") or {}).get("href", ""),
                    }
                )
        group_stats.append(
            {
                "sampling_group": group["sampling_group"],
                "query_total": total,
                "query_total_pages": total_pages,
                "selected_pages": ",".join(map(str, pages)),
                "metadata_rows": fetched,
            }
        )

    meta = pd.DataFrame(rows).drop_duplicates(subset=["neuron_id"]).reset_index(drop=True)
    stats = pd.DataFrame(group_stats)
    stats.to_csv(RESULTS / "neuromorpho_granule_morphometry_sampling_plan.tsv", sep="\t", index=False)
    return meta


def fetch_one_morphometry(neuron_id: int) -> dict[str, Any]:
    cache_path = CACHE / "morphometry" / f"{neuron_id}.json"
    url = f"{BASE}/morphometry/id/{neuron_id}"
    data = api_get_json(url, cache_path)
    if data.get("status") == 404:
        data = {"neuron_id": neuron_id}
    data.setdefault("neuron_id", neuron_id)
    return data


def collect_morphometry(neuron_ids: list[int]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        future_map = {pool.submit(fetch_one_morphometry, int(nid)): int(nid) for nid in neuron_ids}
        for future in as_completed(future_map):
            rows.append(future.result())
    return pd.DataFrame(rows)


def numericize(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    out = df.copy()
    for field in fields:
        if field in out.columns:
            out[field] = pd.to_numeric(out[field], errors="coerce")
    return out


def add_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    condition = out["experiment_condition"].fillna("").str.lower()
    integrity = out["physical_integrity"].fillna("").str.lower()
    domain = out["domain"].fillna("").str.lower()
    out["is_control_like"] = condition.str.contains("control|not reported|healthy|normal", regex=True)
    out["dendrites_present"] = domain.str.contains("dendrite")
    out["moderate_or_complete_dendrites"] = integrity.str.contains("complete|moderate", regex=True)
    out["has_core_morphometry"] = out[["n_stems", "n_branch", "length"]].notna().all(axis=1)
    out["analysis_inclusion"] = (
        out["dendrites_present"] & out["moderate_or_complete_dendrites"] & out["has_core_morphometry"]
    )
    out["compact_branch_index"] = out["n_branch"] / out["length"].replace(0, np.nan)
    return out


def summarize_metric(values: pd.Series) -> dict[str, float]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "q25": np.nan,
            "q75": np.nan,
            "min": np.nan,
            "max": np.nan,
        }
    return {
        "n": int(clean.shape[0]),
        "mean": float(clean.mean()),
        "median": float(clean.median()),
        "q25": float(clean.quantile(0.25)),
        "q75": float(clean.quantile(0.75)),
        "min": float(clean.min()),
        "max": float(clean.max()),
    }


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    included = df[df["analysis_inclusion"]].copy()
    rows: list[dict[str, Any]] = []
    for keys, sub in included.groupby(["analysis_region", "species"], dropna=False):
        analysis_region, species = keys
        for metric in SUMMARY_METRICS + ["compact_branch_index"]:
            stats = summarize_metric(sub[metric])
            rows.append(
                {
                    "analysis_region": analysis_region,
                    "species": species,
                    "metric": metric,
                    **stats,
                }
            )
    for analysis_region, sub in included.groupby("analysis_region", dropna=False):
        for metric in SUMMARY_METRICS + ["compact_branch_index"]:
            stats = summarize_metric(sub[metric])
            rows.append(
                {
                    "analysis_region": analysis_region,
                    "species": "ALL_INCLUDED",
                    "metric": metric,
                    **stats,
                }
            )
    return pd.DataFrame(rows)


def cliffs_delta(x: pd.Series, y: pd.Series) -> float:
    x_vals = pd.to_numeric(x, errors="coerce").dropna().to_numpy()
    y_vals = pd.to_numeric(y, errors="coerce").dropna().to_numpy()
    if len(x_vals) == 0 or len(y_vals) == 0:
        return float("nan")
    gt = 0
    lt = 0
    for value in x_vals:
        gt += int(np.sum(value > y_vals))
        lt += int(np.sum(value < y_vals))
    return float((gt - lt) / (len(x_vals) * len(y_vals)))


def mannwhitney_p_approx(x: pd.Series, y: pd.Series) -> float:
    """Two-sided Mann-Whitney p-value using a tie-corrected normal approximation."""
    x_vals = pd.to_numeric(x, errors="coerce").dropna().to_numpy()
    y_vals = pd.to_numeric(y, errors="coerce").dropna().to_numpy()
    n1 = len(x_vals)
    n2 = len(y_vals)
    if n1 == 0 or n2 == 0:
        return float("nan")
    combined = pd.Series(np.concatenate([x_vals, y_vals]))
    ranks = combined.rank(method="average").to_numpy()
    r1 = float(np.sum(ranks[:n1]))
    u1 = r1 - (n1 * (n1 + 1) / 2.0)
    mean_u = n1 * n2 / 2.0
    _, tie_counts = np.unique(combined.to_numpy(), return_counts=True)
    n = n1 + n2
    tie_term = float(np.sum(tie_counts**3 - tie_counts))
    variance = (n1 * n2 / 12.0) * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if variance <= 0:
        return float("nan")
    # Continuity correction toward the mean.
    diff = u1 - mean_u
    correction = 0.5 if diff > 0 else -0.5 if diff < 0 else 0.0
    z = (diff - correction) / math.sqrt(variance)
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def build_comparison(df: pd.DataFrame) -> pd.DataFrame:
    included = df[df["analysis_inclusion"]].copy()
    dg = included[included["analysis_region"] == "dentate_gyrus"]
    cb = included[included["analysis_region"] == "cerebellum"]
    rows: list[dict[str, Any]] = []
    for metric in SUMMARY_METRICS + ["compact_branch_index"]:
        x = pd.to_numeric(dg[metric], errors="coerce").dropna()
        y = pd.to_numeric(cb[metric], errors="coerce").dropna()
        p_value = mannwhitney_p_approx(x, y)
        rows.append(
            {
                "comparison": "dentate_gyrus_sample_vs_cerebellum_all",
                "metric": metric,
                "dentate_n": int(len(x)),
                "cerebellum_n": int(len(y)),
                "dentate_median": float(x.median()) if len(x) else np.nan,
                "cerebellum_median": float(y.median()) if len(y) else np.nan,
                "median_delta_dentate_minus_cerebellum": float(x.median() - y.median()) if len(x) and len(y) else np.nan,
                "median_ratio_dentate_over_cerebellum": float(x.median() / y.median()) if len(x) and len(y) and y.median() != 0 else np.nan,
                "mannwhitney_p": p_value,
                "cliffs_delta_dentate_vs_cerebellum": cliffs_delta(x, y),
            }
        )
    return pd.DataFrame(rows)


def build_priors(summary: pd.DataFrame, comparison: pd.DataFrame) -> pd.DataFrame:
    def metric_row(region: str, metric: str) -> pd.Series:
        rows = summary[
            (summary["analysis_region"] == region)
            & (summary["species"] == "ALL_INCLUDED")
            & (summary["metric"] == metric)
        ]
        if rows.empty:
            return pd.Series(dtype=object)
        return rows.iloc[0]

    dg_stems = metric_row("dentate_gyrus", "n_stems")
    cb_stems = metric_row("cerebellum", "n_stems")
    dg_branches = metric_row("dentate_gyrus", "n_branch")
    cb_branches = metric_row("cerebellum", "n_branch")
    dg_len = metric_row("dentate_gyrus", "length")
    cb_len = metric_row("cerebellum", "length")

    rows = [
        {
            "aim3_parameter": "input_degree_minimal_dendritic_stems",
            "dentate_empirical_prior": f"{dg_stems.get('median', np.nan):.2f} stems median ({dg_stems.get('q25', np.nan):.2f}-{dg_stems.get('q75', np.nan):.2f} IQR)",
            "cerebellar_empirical_prior": f"{cb_stems.get('median', np.nan):.2f} stems median ({cb_stems.get('q25', np.nan):.2f}-{cb_stems.get('q75', np.nan):.2f} IQR)",
            "recommended_model_use": "Use as a lower-bound anatomical input-sampling prior; it is not equivalent to total synaptic input count.",
            "interpretation": "Cerebellar granule cells show multiple primary dendritic stems, whereas dentate granule cells usually compress input sampling into one primary dendritic tree.",
        },
        {
            "aim3_parameter": "dendritic_branch_complexity",
            "dentate_empirical_prior": f"{dg_branches.get('median', np.nan):.2f} branches median ({dg_branches.get('q25', np.nan):.2f}-{dg_branches.get('q75', np.nan):.2f} IQR)",
            "cerebellar_empirical_prior": f"{cb_branches.get('median', np.nan):.2f} branches median ({cb_branches.get('q25', np.nan):.2f}-{cb_branches.get('q75', np.nan):.2f} IQR)",
            "recommended_model_use": "Use as a morphology-complexity prior separate from sparse input degree.",
            "interpretation": "Branch counts are much closer than raw dendritic length, supporting compact limited-branch architecture rather than identical geometry.",
        },
        {
            "aim3_parameter": "dendritic_field_scale",
            "dentate_empirical_prior": f"{dg_len.get('median', np.nan):.2f} length median ({dg_len.get('q25', np.nan):.2f}-{dg_len.get('q75', np.nan):.2f} IQR)",
            "cerebellar_empirical_prior": f"{cb_len.get('median', np.nan):.2f} length median ({cb_len.get('q25', np.nan):.2f}-{cb_len.get('q75', np.nan):.2f} IQR)",
            "recommended_model_use": "Use as a scale parameter, not as a direct sparse-coding input-degree parameter.",
            "interpretation": "Dentate granule cells have a much larger dendritic field than cerebellar granule cells in this sample.",
        },
        {
            "aim3_parameter": "model_revision",
            "dentate_empirical_prior": "input_degree should not be equated with dendrite number alone",
            "cerebellar_empirical_prior": "input_degree should be separated from branch/field scale",
            "recommended_model_use": "Revise Aim 3 model wording from one morphology parameter to two: primary input-sampling stems/claws and dendritic-field complexity.",
            "interpretation": "The public morphometry supports convergent compact/limited-branch design logic, but also shows strong regional implementation differences.",
        },
    ]
    return pd.DataFrame(rows)


def make_plot(df: pd.DataFrame) -> bool:
    os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    included = df[df["analysis_inclusion"]].copy()
    metrics = [
        ("n_stems", "Primary stems"),
        ("n_branch", "Branches"),
        ("n_bifs", "Bifurcations"),
        ("length", "Dendritic length"),
        ("branch_Order", "Branch order"),
        ("compact_branch_index", "Branches per length"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.ravel()
    colors = {"dentate_gyrus": "#287271", "cerebellum": "#b85634"}
    for ax, (metric, title) in zip(axes, metrics):
        data = [
            pd.to_numeric(included[included["analysis_region"] == region][metric], errors="coerce").dropna()
            for region in ["dentate_gyrus", "cerebellum"]
        ]
        try:
            bp = ax.boxplot(data, tick_labels=["DG", "CB"], patch_artist=True, showfliers=False)
        except TypeError:
            bp = ax.boxplot(data, labels=["DG", "CB"], patch_artist=True, showfliers=False)
        for patch, region in zip(bp["boxes"], ["dentate_gyrus", "cerebellum"]):
            patch.set_facecolor(colors[region])
            patch.set_alpha(0.65)
        for idx, values in enumerate(data, start=1):
            if len(values):
                rng = np.random.default_rng(1234 + idx)
                sample = values.sample(min(len(values), 90), random_state=idx).to_numpy()
                jitter = rng.normal(idx, 0.04, size=len(sample))
                ax.scatter(jitter, sample, s=10, alpha=0.35, color="#222222", linewidths=0)
        if metric in {"length", "surface", "volume", "compact_branch_index"}:
            ax.set_yscale("log")
        ax.set_title(title)
        ax.tick_params(axis="x", length=0)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("NeuroMorpho granule-cell morphometry validation", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)
    return True


def write_report(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    comparison: pd.DataFrame,
    priors: pd.DataFrame,
    plot_built: bool,
) -> None:
    included = df[df["analysis_inclusion"]].copy()
    counts = (
        df.groupby(["analysis_region"], dropna=False)
        .agg(metadata_rows=("neuron_id", "count"), included_rows=("analysis_inclusion", "sum"))
        .reset_index()
    )
    def get_summary(region: str, metric: str) -> pd.Series:
        rows = summary[
            (summary["analysis_region"] == region)
            & (summary["species"] == "ALL_INCLUDED")
            & (summary["metric"] == metric)
        ]
        return rows.iloc[0] if not rows.empty else pd.Series(dtype=object)

    def fmt(region: str, metric: str) -> str:
        row = get_summary(region, metric)
        if row.empty:
            return "NA"
        return f"{row['median']:.2f} ({row['q25']:.2f}-{row['q75']:.2f})"

    lines = [
        "# NeuroMorpho Granule Morphometry Validation",
        "",
        "Date built: 2026-06-23",
        "",
        "## Purpose",
        "",
        "This is the first direct morphology validation layer for Aim 3. It uses all available cerebellar granule-cell hits from the current strict NeuroMorpho query and a reproducible species-stratified dentate granule-cell sample.",
        "",
        "## Sampling",
        "",
        "| Region | Metadata rows | Included rows |",
        "|---|---:|---:|",
    ]
    for _, row in counts.iterrows():
        lines.append(f"| `{row['analysis_region']}` | {int(row['metadata_rows'])} | {int(row['included_rows'])} |")

    lines.extend(
        [
            "",
            "Included rows require dendrites in the reconstruction domain, moderate or complete dendritic integrity, and non-missing `n_stems`, `n_branch`, and `length` morphometry.",
            "",
            "## Main Morphometry Result",
            "",
            "| Metric | Dentate granule median (IQR) | Cerebellar granule median (IQR) | Interpretation |",
            "|---|---:|---:|---|",
            f"| Primary stems | {fmt('dentate_gyrus', 'n_stems')} | {fmt('cerebellum', 'n_stems')} | Different implementation: DG tends toward one primary dendritic tree, CB toward multiple short stems/claws. |",
            f"| Branches | {fmt('dentate_gyrus', 'n_branch')} | {fmt('cerebellum', 'n_branch')} | Branch-count scale is closer than dendritic length. |",
            f"| Bifurcations | {fmt('dentate_gyrus', 'n_bifs')} | {fmt('cerebellum', 'n_bifs')} | Useful proxy for limited branching complexity. |",
            f"| Dendritic length | {fmt('dentate_gyrus', 'length')} | {fmt('cerebellum', 'length')} | DG has a much larger dendritic field in this sample. |",
            f"| Branch order | {fmt('dentate_gyrus', 'branch_Order')} | {fmt('cerebellum', 'branch_Order')} | Branch order separates regional geometry. |",
            "",
            "## Interpretation For The Project",
            "",
            "The morphometry supports a refined version of the morphology hypothesis. Dentate and cerebellar granule cells are not geometrically identical. Instead, both are compact excitatory input-expansion neurons with constrained dendritic branching, but they implement input sampling differently: dentate granule cells through a larger dendritic tree, cerebellar granule cells through several short claw-like dendritic stems.",
            "",
            "This means the Aim 3 computational model should not collapse morphology into one `input_degree` parameter. The better model has at least two anatomical knobs:",
            "",
            "1. `primary_input_sampling_stems_or_claws`: lower-bound input-sampling geometry.",
            "2. `dendritic_field_complexity`: branch count, bifurcation count, and dendritic length/scale.",
            "",
            "## Caveats",
            "",
            "- The dentate set is a species-stratified sample, not the full 9,672-cell NeuroMorpho dentate query.",
            "- The cerebellar set is small and species-heterogeneous under the current strict query.",
            "- Morphometric stems and branches are not direct synaptic input counts.",
            "- Genotype, disease/condition, reconstruction completeness, and archive-specific methods should be filtered more strictly before final manuscript statistics.",
            "",
            "## Outputs",
            "",
            f"- Metadata table: `{rel(OUT_META)}`",
            f"- Morphometry table: `{rel(OUT_MORPH)}`",
            f"- Summary table: `{rel(OUT_SUMMARY)}`",
            f"- DG-vs-CB comparison: `{rel(OUT_COMPARISON)}`",
            f"- Aim 3 empirical priors: `{rel(OUT_PRIORS)}`",
            f"- Plot: `{rel(OUT_PLOT)}`" if plot_built else "- Plot: skipped because Matplotlib was unavailable in the active Python runtime.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    meta = collect_metadata()
    morph = collect_morphometry(meta["neuron_id"].astype(int).tolist())
    morph = numericize(morph, MORPH_FIELDS)
    combined = meta.merge(morph, on="neuron_id", how="left", suffixes=("", "_morph"))
    combined = add_flags(combined)
    summary = build_summary(combined)
    comparison = build_comparison(combined)
    priors = build_priors(summary, comparison)

    meta.to_csv(OUT_META, sep="\t", index=False)
    combined.to_csv(OUT_MORPH, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    comparison.to_csv(OUT_COMPARISON, sep="\t", index=False)
    priors.to_csv(OUT_PRIORS, sep="\t", index=False)
    plot_built = make_plot(combined)
    write_report(combined, summary, comparison, priors, plot_built)

    for path in [OUT_META, OUT_MORPH, OUT_SUMMARY, OUT_COMPARISON, OUT_PRIORS, OUT_PLOT, OUT_MD]:
        if path.exists():
            print(f"Wrote {rel(path)}")
    if not plot_built:
        print("Skipped plot: Matplotlib unavailable in active Python runtime")


if __name__ == "__main__":
    main()
