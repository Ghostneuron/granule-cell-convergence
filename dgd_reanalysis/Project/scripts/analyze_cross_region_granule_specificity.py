#!/usr/bin/env python3
"""Test whether a granule-cell expression configuration transfers across regions.

The analysis trains a regularized classifier on paired library pseudobulks from
one adult Allen brain region and evaluates it, without refitting, in the other
region. Dentate models contrast granule cells with CA1/CA3 pyramidal cells;
cerebellar models contrast granule cells with Purkinje cells. Label
permutations are blocked within library, and matched random-gene panels are
drawn from a deterministic common-gene cache extracted from the same matrices.

This is a molecular specificity audit. It does not establish a causal program
for morphology or generalize beyond the tested adult mouse populations.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

import analyze_allen_consensus_comparators as allen


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

BASE_EXPRESSION = RESULTS / "dgd_allen_consensus_library_gene_expression.tsv.gz"
MODULES = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
TIERS = RESULTS / "primary_core_manuscript_candidate_tiers.tsv"

POOL_EXPRESSION = RESULTS / "dgd_allen_cross_region_feature_pool.tsv.gz"
POOL_MANIFEST = RESULTS / "dgd_allen_cross_region_feature_manifest.tsv"
OUT_METRICS = RESULTS / "dgd_allen_cross_region_transfer_metrics.tsv"
OUT_SUMMARY = RESULTS / "dgd_allen_cross_region_specificity_summary.tsv"
OUT_PERMUTATIONS = RESULTS / "dgd_allen_cross_region_label_permutations.tsv.gz"
OUT_MATCHED_NULL = RESULTS / "dgd_allen_cross_region_matched_gene_null.tsv.gz"
OUT_WEIGHTS = RESULTS / "dgd_allen_cross_region_feature_weights.tsv"
OUT_STAGE = RESULTS / "dgd_allen_cross_region_stage_sensitivity.tsv"
OUT_CONCORDANCE = RESULTS / "dgd_allen_cross_region_contrast_concordance.tsv"
OUT_GENE_CONTRASTS = RESULTS / "dgd_allen_cross_region_gene_contrasts.tsv"
OUT_FIGURE = RESULTS / "dgd_allen_cross_region_transfer.png"
OUT_FIGURE_PDF = RESULTS / "dgd_allen_cross_region_transfer.pdf"
OUT_REPORT = RESULTS / "dgd_allen_cross_region_transfer.md"

RNG_SEED = 205921
N_BACKGROUND_GENES = 2_500
N_LABEL_PERMUTATIONS = 2_000
N_MATCHED_PANELS = 2_000
N_BOOTSTRAP = 5_000
NEAREST_MATCHES = 100

CB_TARGET = "Cerebellar granule"
CB_COMPARATORS = ["Purkinje"]
DG_MATURE = "Dentate granule, mature"
DG_IMMATURE = "Dentate granule, immature"
DG_COMPARATORS = ["CA1/ProS pyramidal", "CA3 pyramidal"]

TIER1_LABEL = "Tier 1 core convergent program"
TIER2_LABEL = "Tier 2 high-confidence wiring/synaptic executor"

CANONICAL_MATRIX = {
    CB_TARGET: "CB",
    "Purkinje": "CB",
    DG_MATURE: "HPF",
    DG_IMMATURE: "HPF",
    "CA1/ProS pyramidal": "HPF",
    "CA3 pyramidal": "HPF",
    "Cortical L4/5 IT excitatory": "Isocortex",
    "Olfactory-bulb GABA, mature proxy": "OLF",
    "Olfactory-bulb GABA, immature proxy": "OLF",
}

FEATURE_LABELS = {
    "tier1": "Tier 1 candidates",
    "tier1_tier2": "Tier 1+2 candidates",
    "downstream_all": "All downstream genes",
    "neurite_morphology": "Neurite/morphology genes",
    "synaptic_excitability": "Synaptic/excitability genes",
    "all_curated": "All curated module genes",
}


def read_matrix_symbols(path: Path) -> set[str]:
    with h5py.File(path, "r") as handle:
        return set(map(str, allen.read_h5ad_column(handle["var"], "gene_symbol")))


def load_feature_sets() -> tuple[dict[str, list[str]], set[str]]:
    modules = pd.read_csv(MODULES, sep="\t")
    tiers = pd.read_csv(TIERS, sep="\t")
    tier1 = tiers.loc[tiers["manuscript_tier"].eq(TIER1_LABEL), "mouse_symbol"].astype(str)
    tier12 = tiers.loc[
        tiers["manuscript_tier"].isin([TIER1_LABEL, TIER2_LABEL]), "mouse_symbol"
    ].astype(str)
    downstream = modules.loc[
        modules["module_family"].eq("downstream_circuit_morphology"),
        "default_mouse_symbol",
    ].astype(str)
    neurite = modules.loc[
        modules["module_id"].eq("downstream_neurite_morphology"),
        "default_mouse_symbol",
    ].astype(str)
    synaptic = modules.loc[
        modules["module_id"].eq("downstream_synaptic_excitability"),
        "default_mouse_symbol",
    ].astype(str)
    feature_sets = {
        "tier1": list(dict.fromkeys(tier1)),
        "tier1_tier2": list(dict.fromkeys(tier12)),
        "downstream_all": list(dict.fromkeys(downstream)),
        "neurite_morphology": list(dict.fromkeys(neurite)),
        "synaptic_excitability": list(dict.fromkeys(synaptic)),
        "all_curated": list(dict.fromkeys(modules["default_mouse_symbol"].astype(str))),
    }
    reserved = set(feature_sets["all_curated"]) | set(feature_sets["tier1_tier2"])
    return feature_sets, reserved


def build_feature_pool(reserved: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if POOL_EXPRESSION.exists() and POOL_MANIFEST.exists():
        return (
            pd.read_csv(POOL_EXPRESSION, sep="\t"),
            pd.read_csv(POOL_MANIFEST, sep="\t"),
        )

    common = set.intersection(*(read_matrix_symbols(path) for path in allen.MATRICES.values()))
    eligible = sorted(common - reserved)
    rng = np.random.default_rng(RNG_SEED)
    n_pool = min(N_BACKGROUND_GENES, len(eligible))
    background = sorted(rng.choice(eligible, size=n_pool, replace=False).tolist())
    selected = sorted(reserved) + background

    taxonomy = allen.build_taxonomy_map()
    alias_to_population = {
        int(alias): population
        for alias, subclass in taxonomy.get("subclass", pd.Series(dtype=object)).items()
        if (population := allen.population_from_subclass(subclass)) is not None
    }
    membership = allen.load_macosko_membership()
    frames: list[pd.DataFrame] = []
    for matrix_name, path in allen.MATRICES.items():
        _, expression = allen.aggregate_matrix(
            matrix_name,
            path,
            membership,
            alias_to_population,
            selected,
        )
        frames.append(expression)
    result = pd.concat(frames, ignore_index=True)
    result.to_csv(POOL_EXPRESSION, sep="\t", index=False, compression="gzip")

    manifest = pd.DataFrame(
        {
            "gene_symbol": selected,
            "feature_role": ["reserved" if gene in reserved else "background_pool" for gene in selected],
            "present_in_all_four_matrices": True,
            "selection_seed": RNG_SEED,
        }
    )
    manifest.to_csv(POOL_MANIFEST, sep="\t", index=False)
    return result, manifest


def canonical_expression(expression: pd.DataFrame) -> pd.DataFrame:
    keep = expression.apply(
        lambda row: CANONICAL_MATRIX.get(str(row["population"])) == str(row["matrix"]), axis=1
    )
    return expression.loc[keep].copy()


def paired_frame(
    expression: pd.DataFrame,
    matrix: str,
    target: str,
    comparators: list[str],
    features: list[str],
) -> pd.DataFrame:
    populations = [target, *comparators]
    subset = expression[
        expression["matrix"].eq(matrix)
        & expression["population"].isin(populations)
        & expression["gene_symbol"].isin(features)
    ]
    wide = subset.pivot_table(
        index=["library", "population"],
        columns="gene_symbol",
        values="mean_log2_expression",
        aggfunc="first",
    )
    available = [gene for gene in features if gene in wide.columns]
    library_sets = [set(wide.xs(pop, level="population").index) for pop in populations]
    paired_libraries = set.intersection(*library_sets)
    wide = wide.loc[wide.index.get_level_values("library").isin(paired_libraries), available]
    population_order = {population: index for index, population in enumerate(populations)}
    order = sorted(
        wide.index,
        key=lambda value: (str(value[0]), population_order[str(value[1])]),
    )
    return wide.loc[order]


def make_model() -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    penalty="l2",
                    C=1.0,
                    class_weight="balanced",
                    solver="liblinear",
                    max_iter=5_000,
                    random_state=RNG_SEED,
                ),
            ),
        ]
    )


def labels_for(frame: pd.DataFrame, target: str) -> np.ndarray:
    return np.asarray(frame.index.get_level_values("population") == target, dtype=int)


def permuted_block_labels(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    labels = np.zeros(len(frame), dtype=int)
    libraries = frame.index.get_level_values("library")
    for library in libraries.unique():
        positions = np.flatnonzero(libraries == library)
        labels[int(rng.choice(positions))] = 1
    return labels


def paired_margins(
    frame: pd.DataFrame,
    scores: np.ndarray,
    target: str,
    comparators: list[str],
) -> pd.Series:
    score_series = pd.Series(scores, index=frame.index, name="score")
    rows: dict[str, float] = {}
    for library in frame.index.get_level_values("library").unique():
        values = score_series.xs(library, level="library")
        rows[str(library)] = float(values.loc[target] - values.reindex(comparators).mean())
    return pd.Series(rows, name="paired_margin")


def bootstrap_median_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    if len(values) < 2:
        return np.nan, np.nan
    draws = rng.choice(values, size=(N_BOOTSTRAP, len(values)), replace=True)
    medians = np.median(draws, axis=1)
    return tuple(map(float, np.quantile(medians, [0.025, 0.975])))


def score_transfer(
    train: pd.DataFrame,
    test: pd.DataFrame,
    train_target: str,
    test_target: str,
    test_comparators: list[str],
) -> tuple[dict[str, float | int], Pipeline, pd.Series]:
    y_train = labels_for(train, train_target)
    y_test = labels_for(test, test_target)
    model = make_model().fit(train.to_numpy(), y_train)
    scores = model.decision_function(test.to_numpy())
    predictions = (scores >= 0).astype(int)
    margins = paired_margins(test, scores, test_target, test_comparators)
    margin_values = margins.to_numpy(dtype=float)
    rng = np.random.default_rng(RNG_SEED + len(train) + len(test) + len(train.columns))
    ci_low, ci_high = bootstrap_median_ci(margin_values, rng)
    nonzero = margin_values[margin_values != 0]
    sign_p = (
        float(stats.binomtest(int(np.sum(nonzero > 0)), len(nonzero), 0.5, alternative="greater").pvalue)
        if len(nonzero)
        else np.nan
    )
    try:
        wilcoxon_p = float(stats.wilcoxon(margin_values, alternative="greater").pvalue)
    except ValueError:
        wilcoxon_p = np.nan
    metrics: dict[str, float | int] = {
        "n_features": len(train.columns),
        "n_train_libraries": train.index.get_level_values("library").nunique(),
        "n_test_libraries": test.index.get_level_values("library").nunique(),
        "roc_auc": float(roc_auc_score(y_test, scores)),
        "balanced_accuracy_at_training_threshold": float(
            balanced_accuracy_score(y_test, predictions)
        ),
        "mean_paired_margin": float(np.mean(margin_values)),
        "median_paired_margin": float(np.median(margin_values)),
        "bootstrap_95ci_low_median_margin": ci_low,
        "bootstrap_95ci_high_median_margin": ci_high,
        "fraction_positive_paired_margins": float(np.mean(margin_values > 0)),
        "paired_sign_p_greater": sign_p,
        "paired_wilcoxon_p_greater": wilcoxon_p,
    }
    return metrics, model, margins


def model_weights(model: Pipeline, features: list[str]) -> pd.DataFrame:
    coefficients = model.named_steps["model"].coef_[0]
    return pd.DataFrame(
        {
            "gene_symbol": features,
            "standardized_logistic_coefficient": coefficients,
            "absolute_coefficient": np.abs(coefficients),
        }
    )


def transfer_specifications(
    expression: pd.DataFrame,
    features: list[str],
    dentate_state: str,
) -> list[dict[str, object]]:
    cb = paired_frame(expression, "CB", CB_TARGET, CB_COMPARATORS, features)
    dg = paired_frame(expression, "HPF", dentate_state, DG_COMPARATORS, features)
    return [
        {
            "direction": "dentate_to_cerebellum",
            "train": dg,
            "test": cb,
            "train_target": dentate_state,
            "test_target": CB_TARGET,
            "test_comparators": CB_COMPARATORS,
        },
        {
            "direction": "cerebellum_to_dentate",
            "train": cb,
            "test": dg,
            "train_target": CB_TARGET,
            "test_target": dentate_state,
            "test_comparators": DG_COMPARATORS,
        },
    ]


def run_observed_and_permutations(
    expression: pd.DataFrame,
    feature_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[tuple[str, str, str], pd.Series]]:
    metric_rows: list[dict[str, object]] = []
    permutation_rows: list[dict[str, object]] = []
    weight_rows: list[pd.DataFrame] = []
    margin_cache: dict[tuple[str, str, str], pd.Series] = {}
    rng = np.random.default_rng(RNG_SEED + 10)

    for state_label, dentate_state in [("mature", DG_MATURE), ("immature", DG_IMMATURE)]:
        for feature_set, requested_features in feature_sets.items():
            available = sorted(set(requested_features) & set(expression["gene_symbol"]))
            for spec in transfer_specifications(expression, available, dentate_state):
                metrics, model, margins = score_transfer(
                    spec["train"],
                    spec["test"],
                    spec["train_target"],
                    spec["test_target"],
                    spec["test_comparators"],
                )
                key = (state_label, feature_set, str(spec["direction"]))
                margin_cache[key] = margins
                row = {
                    "dentate_state": state_label,
                    "feature_set": feature_set,
                    "feature_set_label": FEATURE_LABELS[feature_set],
                    "direction": spec["direction"],
                    **metrics,
                }
                metric_rows.append(row)

                weights = model_weights(model, available)
                weights.insert(0, "direction", spec["direction"])
                weights.insert(0, "feature_set", feature_set)
                weights.insert(0, "dentate_state", state_label)
                weight_rows.append(weights)

                y_test = labels_for(spec["test"], spec["test_target"])
                observed_auc = float(metrics["roc_auc"])
                exceedances = 0
                for permutation in range(N_LABEL_PERMUTATIONS):
                    y_permuted = permuted_block_labels(spec["train"], rng)
                    perm_model = make_model().fit(spec["train"].to_numpy(), y_permuted)
                    perm_scores = perm_model.decision_function(spec["test"].to_numpy())
                    perm_auc = float(roc_auc_score(y_test, perm_scores))
                    exceedances += int(perm_auc >= observed_auc)
                    permutation_rows.append(
                        {
                            "dentate_state": state_label,
                            "feature_set": feature_set,
                            "direction": spec["direction"],
                            "permutation": permutation + 1,
                            "permuted_training_label_auc": perm_auc,
                        }
                    )
                row["training_label_permutation_p"] = float(
                    (1 + exceedances) / (1 + N_LABEL_PERMUTATIONS)
                )
                print(
                    f"Completed label permutations: {state_label} {feature_set} "
                    f"{spec['direction']}",
                    flush=True,
                )

    metrics = pd.DataFrame(metric_rows)
    for state in metrics["dentate_state"].unique():
        mask = metrics["dentate_state"].eq(state)
        metrics.loc[mask, "training_label_permutation_q_bh"] = multipletests(
            metrics.loc[mask, "training_label_permutation_p"], method="fdr_bh"
        )[1]
    return (
        metrics,
        pd.DataFrame(permutation_rows),
        pd.concat(weight_rows, ignore_index=True),
        margin_cache,
    )


def gene_matching_features(expression: pd.DataFrame) -> pd.DataFrame:
    stats_frame = expression.groupby("gene_symbol", as_index=True).agg(
        mean_expression=("mean_log2_expression", "mean"),
        mean_detection=("detection_fraction", "mean"),
    )
    standardized = stats_frame.sub(stats_frame.mean()).div(
        stats_frame.std(ddof=0).replace(0, np.nan)
    )
    return standardized.fillna(0)


def matching_neighborhoods(
    targets: list[str],
    background: list[str],
    matching: pd.DataFrame,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    neighborhoods: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for target in targets:
        distances = ((matching.loc[background] - matching.loc[target]) ** 2).sum(axis=1)
        nearest = distances.nsmallest(min(NEAREST_MATCHES, len(distances)))
        weights = np.exp(-(nearest.to_numpy() - nearest.min()) / 0.35)
        neighborhoods[target] = (
            nearest.index.to_numpy(dtype=object),
            weights / weights.sum(),
        )
    return neighborhoods


def matched_panel(
    targets: list[str],
    background: list[str],
    neighborhoods: dict[str, tuple[np.ndarray, np.ndarray]],
    rng: np.random.Generator,
) -> list[str]:
    target_order = list(targets)
    rng.shuffle(target_order)
    selected: list[str] = []
    available = set(background)
    for target in target_order:
        candidates, base_weights = neighborhoods[target]
        keep = np.asarray([candidate in available for candidate in candidates], dtype=bool)
        if keep.any():
            usable = candidates[keep]
            weights = base_weights[keep]
            weights = weights / weights.sum()
            choice = str(rng.choice(usable, p=weights))
        else:
            choice = str(rng.choice(np.asarray(sorted(available), dtype=object)))
        selected.append(choice)
        available.remove(choice)
    return selected


def prepare_designs(expression: pd.DataFrame) -> dict[str, pd.DataFrame]:
    all_features = sorted(expression["gene_symbol"].unique())
    return {
        "cerebellum": paired_frame(
            expression, "CB", CB_TARGET, CB_COMPARATORS, all_features
        ),
        "dentate_mature": paired_frame(
            expression, "HPF", DG_MATURE, DG_COMPARATORS, all_features
        ),
        "dentate_immature": paired_frame(
            expression, "HPF", DG_IMMATURE, DG_COMPARATORS, all_features
        ),
    }


def bidirectional_auc(
    designs: dict[str, pd.DataFrame],
    features: list[str],
    dentate_state: str = DG_MATURE,
) -> tuple[float, float]:
    dentate_key = "dentate_mature" if dentate_state == DG_MATURE else "dentate_immature"
    cb = designs["cerebellum"].loc[:, features]
    dg = designs[dentate_key].loc[:, features]

    dg_model = make_model().fit(dg.to_numpy(), labels_for(dg, dentate_state))
    cb_scores = dg_model.decision_function(cb.to_numpy())
    left_auc = float(roc_auc_score(labels_for(cb, CB_TARGET), cb_scores))

    cb_model = make_model().fit(cb.to_numpy(), labels_for(cb, CB_TARGET))
    dg_scores = cb_model.decision_function(dg.to_numpy())
    right_auc = float(roc_auc_score(labels_for(dg, dentate_state), dg_scores))
    return left_auc, right_auc


def run_matched_gene_null(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    feature_sets: dict[str, list[str]],
    observed_metrics: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matching = gene_matching_features(expression)
    designs = prepare_designs(expression)
    background = manifest.loc[
        manifest["feature_role"].eq("background_pool"), "gene_symbol"
    ].astype(str)
    background = [gene for gene in background if gene in matching.index]
    rng = np.random.default_rng(RNG_SEED + 20)
    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    mature = observed_metrics[observed_metrics["dentate_state"].eq("mature")]
    for feature_set, requested in feature_sets.items():
        targets = [gene for gene in requested if gene in matching.index]
        neighborhoods = matching_neighborhoods(targets, background, matching)
        observed = mature[mature["feature_set"].eq(feature_set)].set_index("direction")
        observed_left = float(observed.loc["dentate_to_cerebellum", "roc_auc"])
        observed_right = float(observed.loc["cerebellum_to_dentate", "roc_auc"])
        observed_min = min(observed_left, observed_right)
        exceedances = 0
        for replicate in range(N_MATCHED_PANELS):
            genes = matched_panel(targets, background, neighborhoods, rng)
            left_auc, right_auc = bidirectional_auc(designs, genes)
            minimum = min(left_auc, right_auc)
            exceedances += int(minimum >= observed_min)
            rows.append(
                {
                    "feature_set": feature_set,
                    "replicate": replicate + 1,
                    "n_genes": len(genes),
                    "dentate_to_cerebellum_auc": left_auc,
                    "cerebellum_to_dentate_auc": right_auc,
                    "minimum_bidirectional_auc": minimum,
                    "matched_genes": ",".join(genes),
                }
            )
        print(f"Completed matched panels: {feature_set}", flush=True)
        values = np.asarray(
            [row["minimum_bidirectional_auc"] for row in rows if row["feature_set"] == feature_set],
            dtype=float,
        )
        summaries.append(
            {
                "feature_set": feature_set,
                "observed_dentate_to_cerebellum_auc": observed_left,
                "observed_cerebellum_to_dentate_auc": observed_right,
                "observed_minimum_bidirectional_auc": observed_min,
                "matched_null_median_minimum_auc": float(np.median(values)),
                "matched_null_95ci_low_minimum_auc": float(np.quantile(values, 0.025)),
                "matched_null_95ci_high_minimum_auc": float(np.quantile(values, 0.975)),
                "matched_gene_panel_empirical_p": float(
                    (1 + exceedances) / (1 + N_MATCHED_PANELS)
                ),
                "matching_rule": (
                    "same-size unique genes from 2500-gene common-matrix pool; "
                    "nearest overall mean expression and detection"
                ),
            }
        )
    summary = pd.DataFrame(summaries)
    summary["matched_gene_panel_q_bh"] = multipletests(
        summary["matched_gene_panel_empirical_p"], method="fdr_bh"
    )[1]
    return pd.DataFrame(rows), summary


def build_stage_sensitivity(
    margin_cache: dict[tuple[str, str, str], pd.Series],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for feature_set in FEATURE_LABELS:
        mature = margin_cache[("mature", feature_set, "cerebellum_to_dentate")]
        immature = margin_cache[("immature", feature_set, "cerebellum_to_dentate")]
        common = mature.index.intersection(immature.index)
        differences = immature.loc[common] - mature.loc[common]
        try:
            p_value = float(stats.wilcoxon(differences, alternative="two-sided").pvalue)
        except ValueError:
            p_value = np.nan
        rows.append(
            {
                "feature_set": feature_set,
                "feature_set_label": FEATURE_LABELS[feature_set],
                "n_shared_libraries": len(common),
                "median_mature_margin": float(mature.loc[common].median()),
                "median_immature_margin": float(immature.loc[common].median()),
                "median_immature_minus_mature_margin": float(differences.median()),
                "paired_wilcoxon_p_two_sided": p_value,
                "scope": "adult HPF immature-neuron state versus mature DG; not developmental age",
            }
        )
    result = pd.DataFrame(rows)
    result["paired_wilcoxon_q_bh"] = multipletests(
        result["paired_wilcoxon_p_two_sided"].fillna(1), method="fdr_bh"
    )[1]
    return result


def mean_paired_contrast(
    frame: pd.DataFrame,
    target: str,
    comparators: list[str],
) -> pd.Series:
    rows: list[pd.Series] = []
    for library in frame.index.get_level_values("library").unique():
        values = frame.xs(library, level="library")
        rows.append(values.loc[target] - values.reindex(comparators).mean(axis=0))
    return pd.DataFrame(rows).mean(axis=0)


def contrast_statistics(
    genes: list[str],
    cerebellar_delta: pd.Series,
    dentate_delta: pd.Series,
) -> dict[str, float]:
    left = cerebellar_delta.reindex(genes).to_numpy(dtype=float)
    right = dentate_delta.reindex(genes).to_numpy(dtype=float)
    valid = np.isfinite(left) & np.isfinite(right)
    left = left[valid]
    right = right[valid]
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    cosine = float(np.dot(left, right) / denominator) if denominator else np.nan
    spearman = (
        float(stats.spearmanr(left, right).statistic)
        if len(left) >= 3 and np.std(left) > 0 and np.std(right) > 0
        else np.nan
    )
    return {
        "n_genes": len(left),
        "cosine_concordance": cosine,
        "spearman_contrast_correlation": spearman,
        "same_sign_fraction": float(np.mean(np.sign(left) == np.sign(right))),
        "both_positive_fraction": float(np.mean((left > 0) & (right > 0))),
    }


def run_contrast_concordance(
    expression: pd.DataFrame,
    feature_sets: dict[str, list[str]],
    matched_null: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    designs = prepare_designs(expression)
    cb_delta = mean_paired_contrast(designs["cerebellum"], CB_TARGET, CB_COMPARATORS)
    dg_delta = mean_paired_contrast(designs["dentate_mature"], DG_MATURE, DG_COMPARATORS)
    rows: list[dict[str, object]] = []
    metrics = [
        "cosine_concordance",
        "spearman_contrast_correlation",
        "same_sign_fraction",
        "both_positive_fraction",
    ]
    for feature_set, requested in feature_sets.items():
        genes = [gene for gene in requested if gene in cb_delta.index and gene in dg_delta.index]
        observed = contrast_statistics(genes, cb_delta, dg_delta)
        null_rows = matched_null[matched_null["feature_set"].eq(feature_set)]
        null_values: dict[str, list[float]] = defaultdict(list)
        for value in null_rows["matched_genes"]:
            result = contrast_statistics(str(value).split(","), cb_delta, dg_delta)
            for metric in metrics:
                null_values[metric].append(float(result[metric]))
        row: dict[str, object] = {
            "feature_set": feature_set,
            "feature_set_label": FEATURE_LABELS[feature_set],
            **observed,
        }
        for metric in metrics:
            values = np.asarray(null_values[metric], dtype=float)
            values = values[np.isfinite(values)]
            row[f"matched_null_median_{metric}"] = float(np.median(values))
            row[f"matched_null_95ci_low_{metric}"] = float(np.quantile(values, 0.025))
            row[f"matched_null_95ci_high_{metric}"] = float(np.quantile(values, 0.975))
            row[f"matched_null_p_greater_{metric}"] = float(
                (1 + np.sum(values >= float(observed[metric]))) / (1 + len(values))
            )
            row[f"matched_null_p_less_{metric}"] = float(
                (1 + np.sum(values <= float(observed[metric]))) / (1 + len(values))
            )
        rows.append(row)
    summary = pd.DataFrame(rows)
    for metric in metrics:
        summary[f"matched_null_q_bh_greater_{metric}"] = multipletests(
            summary[f"matched_null_p_greater_{metric}"], method="fdr_bh"
        )[1]
        summary[f"matched_null_q_bh_less_{metric}"] = multipletests(
            summary[f"matched_null_p_less_{metric}"], method="fdr_bh"
        )[1]

    genes = sorted(set().union(*map(set, feature_sets.values())))
    gene_rows = pd.DataFrame(
        {
            "gene_symbol": genes,
            "cerebellar_granule_minus_purkinje": cb_delta.reindex(genes).to_numpy(),
            "dentate_granule_minus_mean_ca1_ca3": dg_delta.reindex(genes).to_numpy(),
        }
    )
    gene_rows["same_direction"] = np.sign(
        gene_rows["cerebellar_granule_minus_purkinje"]
    ) == np.sign(gene_rows["dentate_granule_minus_mean_ca1_ca3"])
    gene_rows["both_positive"] = (
        gene_rows["cerebellar_granule_minus_purkinje"].gt(0)
        & gene_rows["dentate_granule_minus_mean_ca1_ca3"].gt(0)
    )
    for feature_set, feature_genes in feature_sets.items():
        gene_rows[f"member_{feature_set}"] = gene_rows["gene_symbol"].isin(feature_genes)
    return summary, gene_rows


def build_summary(metrics: pd.DataFrame, matched_summary: pd.DataFrame) -> pd.DataFrame:
    mature = metrics[metrics["dentate_state"].eq("mature")]
    rows: list[dict[str, object]] = []
    for feature_set in FEATURE_LABELS:
        subset = mature[mature["feature_set"].eq(feature_set)].set_index("direction")
        left = subset.loc["dentate_to_cerebellum"]
        right = subset.loc["cerebellum_to_dentate"]
        iut_p = max(
            float(left["training_label_permutation_p"]),
            float(right["training_label_permutation_p"]),
        )
        left_auc = float(left["roc_auc"])
        right_auc = float(right["roc_auc"])
        if left_auc > 0.5 and right_auc > 0.5 and iut_p < 0.05:
            interpretation = "bidirectional_positive_transfer_before_matched_panel_test"
        elif left_auc > 0.5 and right_auc > 0.5:
            interpretation = "descriptive_positive_ranking_not_label_null_exceeding"
        elif left_auc < 0.5 and right_auc < 0.5:
            interpretation = "bidirectional_reversed_transfer"
        else:
            interpretation = "directionally_inconsistent_or_null"
        rows.append(
            {
                "feature_set": feature_set,
                "feature_set_label": FEATURE_LABELS[feature_set],
                "n_features": int(left["n_features"]),
                "dentate_to_cerebellum_auc": left_auc,
                "cerebellum_to_dentate_auc": right_auc,
                "minimum_bidirectional_auc": min(left_auc, right_auc),
                "intersection_union_label_permutation_p": iut_p,
                "transfer_interpretation": interpretation,
            }
        )
    result = pd.DataFrame(rows).merge(matched_summary, on="feature_set", how="left")
    positive_but_matched_null = (
        result["transfer_interpretation"].str.startswith("bidirectional_positive")
        & result["matched_gene_panel_empirical_p"].ge(0.05)
    )
    result.loc[
        positive_but_matched_null, "transfer_interpretation"
    ] = "descriptive_positive_ranking_not_matched_null_exceeding"
    result["intersection_union_label_permutation_q_bh"] = multipletests(
        result["intersection_union_label_permutation_p"], method="fdr_bh"
    )[1]
    return result


def build_figure(
    metrics: pd.DataFrame,
    summary: pd.DataFrame,
    matched_null: pd.DataFrame,
    stage: pd.DataFrame,
    gene_contrasts: pd.DataFrame,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 11.5), constrained_layout=True)

    ax = axes[0, 0]
    ax.text(-0.08, 1.03, "a", transform=ax.transAxes, fontsize=16, fontweight="bold")
    ax.axis("off")
    boxes = [
        (0.02, 0.65, 0.28, 0.20, "Train in dentate\nDG vs CA1/CA3", "#DCEAF7"),
        (0.70, 0.65, 0.28, 0.20, "Test in cerebellum\nGC vs Purkinje", "#FBE4D5"),
        (0.02, 0.20, 0.28, 0.20, "Train in cerebellum\nGC vs Purkinje", "#FBE4D5"),
        (0.70, 0.20, 0.28, 0.20, "Test in dentate\nDG vs CA1/CA3", "#DCEAF7"),
    ]
    for x, y, width, height, label, color in boxes:
        rect = plt.Rectangle((x, y), width, height, facecolor=color, edgecolor="#374151", lw=1.2)
        ax.add_patch(rect)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center", fontsize=11)
    for y in [0.75, 0.30]:
        ax.annotate("", xy=(0.68, y), xytext=(0.32, y), arrowprops={"arrowstyle": "->", "lw": 1.8})
    ax.text(0.50, 0.90, "No refitting in the test region", ha="center", fontsize=11, fontweight="bold")
    ax.text(0.50, 0.52, "Library-blocked label permutations", ha="center", fontsize=10)

    ax = axes[0, 1]
    ax.text(-0.08, 1.03, "b", transform=ax.transAxes, fontsize=16, fontweight="bold")
    mature = metrics[metrics["dentate_state"].eq("mature")].pivot(
        index="feature_set_label", columns="direction", values="roc_auc"
    )
    order = [FEATURE_LABELS[key] for key in FEATURE_LABELS]
    mature = mature.reindex(order)
    mature = mature.rename(
        columns={
            "dentate_to_cerebellum": "Dentate to cerebellum",
            "cerebellum_to_dentate": "Cerebellum to dentate",
        }
    )
    sns.heatmap(
        mature,
        ax=ax,
        cmap="RdBu_r",
        center=0.5,
        vmin=0,
        vmax=1,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 11},
        cbar_kws={"label": "Cross-region ROC AUC", "shrink": 0.8},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=25)

    ax = axes[1, 0]
    ax.text(-0.08, 1.03, "c", transform=ax.transAxes, fontsize=16, fontweight="bold")
    shown_sets = ["tier1", "tier1_tier2", "downstream_all"]
    colors = {"tier1": "#1B9E77", "tier1_tier2": "#7570B3", "downstream_all": "#D95F02"}
    positions = np.arange(len(shown_sets))
    for position, feature_set in zip(positions, shown_sets, strict=True):
        values = matched_null.loc[
            matched_null["feature_set"].eq(feature_set), "minimum_bidirectional_auc"
        ].to_numpy()
        parts = ax.violinplot(values, positions=[position], widths=0.72, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor("#CBD5E1")
            body.set_edgecolor("#64748B")
            body.set_alpha(0.8)
        observed = float(
            summary.loc[
                summary["feature_set"].eq(feature_set), "minimum_bidirectional_auc"
            ].iloc[0]
        )
        ax.scatter(position, observed, s=80, color=colors[feature_set], edgecolor="white", zorder=3)
    ax.axhline(0.5, color="black", lw=0.9, ls="--")
    ax.set_xticks(positions, ["Tier 1", "Tier 1+2", "All downstream"])
    ax.set_ylabel("Minimum AUC across both directions")
    ax.set_ylim(-0.03, 1.03)

    ax = axes[1, 1]
    ax.text(-0.08, 1.03, "d", transform=ax.transAxes, fontsize=16, fontweight="bold")
    tier1_contrasts = gene_contrasts[gene_contrasts["member_tier1"]].copy()
    ax.scatter(
        tier1_contrasts["cerebellar_granule_minus_purkinje"],
        tier1_contrasts["dentate_granule_minus_mean_ca1_ca3"],
        s=75,
        color="#1B9E77",
        edgecolor="white",
        linewidth=0.7,
    )
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    for row in tier1_contrasts.itertuples(index=False):
        ax.annotate(
            row.gene_symbol,
            (row.cerebellar_granule_minus_purkinje, row.dentate_granule_minus_mean_ca1_ca3),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=9,
        )
    ax.set_xlabel("Cerebellar granule minus Purkinje\nmean log2 expression")
    ax.set_ylabel("Dentate granule minus mean CA1/CA3\nmean log2 expression")

    fig.savefig(OUT_FIGURE, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(
    metrics: pd.DataFrame,
    summary: pd.DataFrame,
    stage: pd.DataFrame,
    concordance: pd.DataFrame,
) -> None:
    indexed = summary.set_index("feature_set")
    tier1 = indexed.loc["tier1"]
    tier12 = indexed.loc["tier1_tier2"]
    downstream = indexed.loc["downstream_all"]
    stage_index = stage.set_index("feature_set")
    concordance_index = concordance.set_index("feature_set")
    tier1_concordance = concordance_index.loc["tier1"]
    downstream_concordance = concordance_index.loc["downstream_all"]
    lines = [
        "# Cross-region granule-cell configuration transfer",
        "",
        "## Question",
        "",
        "Does a molecular pattern that distinguishes dentate granule cells from CA1/CA3 pyramidal cells also distinguish cerebellar granule cells from Purkinje cells, and vice versa?",
        "",
        "## Design",
        "",
        "- The unit of analysis is a library-population pseudobulk, not an individual cell.",
        "- Training used only complete within-library sets: cerebellar granule plus Purkinje or dentate granule plus CA1 plus CA3.",
        "- The fitted classifier was transferred to the other region without refitting or threshold tuning.",
        "- Training-label permutations were blocked within library. Matched random panels were drawn from 2,500 genes present in all four Allen matrices and matched on overall expression and detection.",
        "- The mature-state test used 63 cerebellar library pairs and 17 dentate library triplets. The immature dentate sensitivity test used nine triplets.",
        "",
        "## Main result",
        "",
        f"- Tier 1 ranked the held-out populations in the expected direction (dentate-to-cerebellum AUC {tier1.dentate_to_cerebellum_auc:.3f}; cerebellum-to-dentate AUC {tier1.cerebellum_to_dentate_auc:.3f}; minimum AUC {tier1.minimum_bidirectional_auc:.3f}), but this was not exceptional under blocked training-label permutations (intersection-union P={tier1.intersection_union_label_permutation_p:.4g}) or expression/detection-matched panels (P={tier1.matched_gene_panel_empirical_p:.4g}). The training-derived decision threshold also failed to transfer (balanced accuracies 0.508 and 0.500).",
        f"- Tier 1+2 showed the same descriptive pattern at lower rank performance (AUCs {tier12.dentate_to_cerebellum_auc:.3f} and {tier12.cerebellum_to_dentate_auc:.3f}; minimum {tier12.minimum_bidirectional_auc:.3f}; matched-panel P={tier12.matched_gene_panel_empirical_p:.4g}).",
        f"- The broad downstream panel did not transfer as a common configuration. Both directions were reversed relative to the trained label (AUCs {downstream.dentate_to_cerebellum_auc:.3f} and {downstream.cerebellum_to_dentate_auc:.3f}).",
        f"- Direct contrast concordance separated recurrence from configuration: {tier1_concordance.both_positive_fraction:.3f} of Tier 1 genes were positive in both local contrasts (matched-panel P={tier1_concordance.matched_null_p_greater_both_positive_fraction:.4g}), but their cross-region contrast correlation was only {tier1_concordance.spearman_contrast_correlation:.3f} (matched-panel P={tier1_concordance.matched_null_p_greater_spearman_contrast_correlation:.4g}). The broad downstream contrast vectors were anti-aligned (Spearman rho={downstream_concordance.spearman_contrast_correlation:.3f}).",
        "",
        "## Stage sensitivity",
        "",
        f"- For the cerebellum-trained downstream model, the median paired margin changed from {stage_index.loc['downstream_all', 'median_mature_margin']:.3f} in mature dentate granule cells to {stage_index.loc['downstream_all', 'median_immature_margin']:.3f} in the adult immature-neuron state (immature-minus-mature {stage_index.loc['downstream_all', 'median_immature_minus_mature_margin']:.3f}; paired P={stage_index.loc['downstream_all', 'paired_wilcoxon_p_two_sided']:.4g}).",
        f"- Tier 1 remained positive in both dentate states; its immature-minus-mature paired-margin shift was {stage_index.loc['tier1', 'median_immature_minus_mature_margin']:.3f} (paired P={stage_index.loc['tier1', 'paired_wilcoxon_p_two_sided']:.4g}).",
        "- The immature group is an adult HPF transcriptomic state, not a developmental-age series. Cerebellar developmental-stage transfer cannot be estimated from this adult Allen matrix.",
        "",
        "## Interpretation",
        "",
        "The result separates same-direction candidate recurrence from a shared multigene configuration. Tier 1 is enriched for genes that are positive in both local granule-versus-comparator contrasts, consistent with the existing branch-local analysis. However, neither classifier transfer nor the relative cross-region contrast pattern exceeded the matched configuration null. The full downstream neurite/synaptic panel was anti-aligned across mature regional contrasts. The current data therefore support limited candidate recurrence but do not establish a granule-cell-specific molecular configuration or a causal explanation of morphology.",
        "",
        "## Limits",
        "",
        "- Both directions use one adult mouse platform; anatomical region and comparator identity remain linked.",
        "- Candidate tiers were selected in the primary discovery datasets, so the Allen transfer is external validation of a prespecified panel, not de novo discovery.",
        "- The analysis tests the populations represented here and cannot establish uniqueness relative to every neuronal class.",
        "- Morphology is not measured in the expression matrix.",
        "",
        "## Outputs",
        "",
    ]
    for path in [
        OUT_METRICS,
        OUT_SUMMARY,
        OUT_PERMUTATIONS,
        OUT_MATCHED_NULL,
        OUT_WEIGHTS,
        OUT_STAGE,
        OUT_CONCORDANCE,
        OUT_GENE_CONTRASTS,
        OUT_FIGURE,
        OUT_FIGURE_PDF,
    ]:
        lines.append(f"- `{path.relative_to(ROOT)}`")
    OUT_REPORT.write_text("\n".join(lines) + "\n")


def main() -> None:
    feature_sets, reserved = load_feature_sets()
    base_expression = pd.read_csv(BASE_EXPRESSION, sep="\t")
    pool_expression, manifest = build_feature_pool(reserved)
    expression = pd.concat(
        [
            base_expression[base_expression["gene_symbol"].isin(reserved)],
            pool_expression[pool_expression["gene_symbol"].isin(set(manifest["gene_symbol"]) - reserved)],
        ],
        ignore_index=True,
    ).drop_duplicates(["matrix", "population", "library", "gene_symbol"])
    expression = canonical_expression(expression)

    metrics, permutations, weights, margins = run_observed_and_permutations(
        expression, feature_sets
    )
    metrics.to_csv(OUT_METRICS, sep="\t", index=False)
    permutations.to_csv(OUT_PERMUTATIONS, sep="\t", index=False, compression="gzip")
    weights.to_csv(OUT_WEIGHTS, sep="\t", index=False)
    matched_null, matched_summary = run_matched_gene_null(
        expression, manifest, feature_sets, metrics
    )
    stage = build_stage_sensitivity(margins)
    concordance, gene_contrasts = run_contrast_concordance(
        expression, feature_sets, matched_null
    )
    summary = build_summary(metrics, matched_summary)

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    matched_null.to_csv(OUT_MATCHED_NULL, sep="\t", index=False, compression="gzip")
    stage.to_csv(OUT_STAGE, sep="\t", index=False)
    concordance.to_csv(OUT_CONCORDANCE, sep="\t", index=False)
    gene_contrasts.to_csv(OUT_GENE_CONTRASTS, sep="\t", index=False)
    build_figure(metrics, summary, matched_null, stage, gene_contrasts)
    write_report(metrics, summary, stage, concordance)
    print(f"Wrote {OUT_REPORT}")
    print(f"Wrote {OUT_FIGURE}")


if __name__ == "__main__":
    main()
