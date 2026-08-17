#!/usr/bin/env python3
"""Common-matrix Allen comparison for external candidate testing.

The script reads the official 2025 Consensus-WMB-Macosko-10X log2 matrices
directly with h5py, maps cells to the integrated Allen taxonomy, and retains
only predeclared manuscript module/candidate genes. Cells are summarized by
population and library before any uncertainty calculation. The resulting
analysis is an adult mouse, common-platform specificity test; it is not a
developmental replication of the primary-core analysis.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
ALLEN = Path(
    os.environ.get(
        "ALLEN_CONSENSUS_DIR",
        "/Volumes/VV 2021 backup drive 01/Codex_Project_Archive/"
        "Hippocanpus&Cerebellum/External_Data/Allen_Institute/"
        "Consensus_WMB_Macosko_20251031",
    )
)

MODULES = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
TIERS = RESULTS / "primary_core_manuscript_candidate_tiers.tsv"
MEMBERSHIP = ALLEN / "cell_to_cluster_membership.csv"
TAXONOMY = ALLEN / "cluster_to_annotation.csv"

MATRICES = {
    "CB": ALLEN / "Macosko-10X-CB-log2.h5ad",
    "HPF": ALLEN / "Macosko-10X-HPF-log2.h5ad",
    "OLF": ALLEN / "Macosko-10X-OLF-log2.h5ad",
    "Isocortex": ALLEN / "Macosko-10X-Isocortex-log2.h5ad",
}

OUT_COUNTS = RESULTS / "dgd_allen_consensus_population_counts.tsv"
OUT_LIBRARY_GENE = RESULTS / "dgd_allen_consensus_library_gene_expression.tsv.gz"
OUT_POP_GENE = RESULTS / "dgd_allen_consensus_population_gene_expression.tsv"
OUT_LIBRARY_MODULE = RESULTS / "dgd_allen_consensus_library_module_scores.tsv"
OUT_POP_MODULE = RESULTS / "dgd_allen_consensus_population_module_scores.tsv"
OUT_SIMILARITY = RESULTS / "dgd_allen_consensus_pair_similarity.tsv"
OUT_SPECIFICITY = RESULTS / "dgd_allen_consensus_candidate_specificity.tsv"
OUT_LOCAL_GENE = RESULTS / "dgd_allen_consensus_local_gene_contrasts.tsv"
OUT_LOCAL_MODULE = RESULTS / "dgd_allen_consensus_local_module_contrasts.tsv"
OUT_MATCHED_NULL = RESULTS / "dgd_allen_consensus_candidate_matched_null.tsv"
OUT_FIGURE = RESULTS / "dgd_allen_consensus_comparators.png"
OUT_REPORT = RESULTS / "dgd_allen_consensus_comparators.md"

TIER1 = "Tier 1 core convergent program"
TIER2 = "Tier 2 high-confidence wiring/synaptic executor"
MIN_LIBRARY_CELLS = 50
CHUNK_ROWS = 10_000
N_BOOTSTRAP = 2_000
RNG_SEED = 205920

POPULATION_ORDER = [
    "Cerebellar granule",
    "Dentate granule, mature",
    "Dentate granule, immature",
    "Purkinje",
    "CA1/ProS pyramidal",
    "CA3 pyramidal",
    "Cortical L4/5 IT excitatory",
    "Olfactory-bulb GABA, mature proxy",
    "Olfactory-bulb GABA, immature proxy",
]

DISPLAY_NAMES = {
    "Cerebellar granule": "CB granule",
    "Dentate granule, mature": "DG mature",
    "Dentate granule, immature": "DG immature",
    "Purkinje": "Purkinje",
    "CA1/ProS pyramidal": "CA1/ProS",
    "CA3 pyramidal": "CA3",
    "Cortical L4/5 IT excitatory": "Cortical L4/5 IT",
    "Olfactory-bulb GABA, mature proxy": "OB GABA mature proxy",
    "Olfactory-bulb GABA, immature proxy": "OB GABA immature proxy",
}

PAIR_ORDER = [
    ("Cerebellar granule", "Dentate granule, mature", "target granule pair"),
    ("Cerebellar granule", "Dentate granule, immature", "target/immature bridge"),
    ("Cerebellar granule", "Purkinje", "cerebellar comparator"),
    ("Dentate granule, mature", "CA1/ProS pyramidal", "hippocampal comparator"),
    ("Dentate granule, mature", "CA3 pyramidal", "hippocampal comparator"),
    ("Cerebellar granule", "Cortical L4/5 IT excitatory", "compact excitatory comparator"),
    ("Dentate granule, mature", "Cortical L4/5 IT excitatory", "compact excitatory comparator"),
    ("Cerebellar granule", "Olfactory-bulb GABA, mature proxy", "granule-named proxy"),
    ("Dentate granule, mature", "Olfactory-bulb GABA, mature proxy", "granule-named proxy"),
    ("Cerebellar granule", "Olfactory-bulb GABA, immature proxy", "granule-lineage proxy"),
    ("Dentate granule, mature", "Olfactory-bulb GABA, immature proxy", "granule-lineage proxy"),
]

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values],
        dtype=object,
    )


def read_h5ad_column(group: h5py.Group, name: str) -> np.ndarray:
    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        values = obj[:]
        return decode_strings(values) if values.dtype.kind in {"O", "S", "U"} else values
    categories = decode_strings(obj["categories"][:])
    codes = obj["codes"][:]
    result = np.full(len(codes), None, dtype=object)
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    return result


def build_taxonomy_map() -> pd.DataFrame:
    annotations = pd.read_csv(TAXONOMY)
    pivot = annotations.pivot_table(
        index="cluster_alias",
        columns="cluster_annotation_term_set_name",
        values="cluster_annotation_term_name",
        aggfunc="first",
    )
    pivot.index = pd.to_numeric(pivot.index, errors="coerce").astype("Int64")
    return pivot


def population_from_subclass(subclass: object) -> str | None:
    label = str(subclass)
    if label == "387 CB_Granule_Glut":
        return "Cerebellar granule"
    if label == "386 CBX_Purkinje_Gaba":
        return "Purkinje"
    if label == "043 DG_Glut":
        return "Dentate granule, mature"
    if label == "044 DG_Glut-IMN":
        return "Dentate granule, immature"
    if label == "027 CA1-ProS_Glut":
        return "CA1/ProS pyramidal"
    if label == "029 CA3_Glut":
        return "CA3 pyramidal"
    if label == "003 L4/5-IT-CTX_Glut":
        return "Cortical L4/5 IT excitatory"
    if label.startswith(("046 ", "047 ", "048 ", "049 ", "050 ", "051 ", "052 ")):
        return "Olfactory-bulb GABA, mature proxy"
    if label == "053 OB-STR-CTX_Mex3a_Gaba-IMN":
        return "Olfactory-bulb GABA, immature proxy"
    return None


def load_macosko_membership() -> pd.Series:
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        MEMBERSHIP,
        usecols=["cell_label", "cluster_alias"],
        dtype={"cell_label": "string", "cluster_alias": "Int64"},
        chunksize=500_000,
    ):
        keep = chunk["cell_label"].str.startswith("p", na=False)
        if keep.any():
            pieces.append(chunk.loc[keep])
    membership = pd.concat(pieces, ignore_index=True).drop_duplicates("cell_label")
    return membership.set_index("cell_label")["cluster_alias"]


def load_gene_sets() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    modules = pd.read_csv(MODULES, sep="\t")
    modules["mouse_symbol"] = modules["default_mouse_symbol"].astype(str)
    tiers = pd.read_csv(TIERS, sep="\t")
    tiers = tiers[tiers["manuscript_tier"].isin([TIER1, TIER2])].copy()
    tiers["canonical_gene"] = tiers["gene"].astype(str)
    tier_symbols = tiers["mouse_symbol"].astype(str).tolist()
    symbols = list(dict.fromkeys(modules["mouse_symbol"].tolist() + tier_symbols))
    return modules, tiers, symbols


def aggregate_matrix(
    matrix_name: str,
    path: Path,
    membership: pd.Series,
    alias_to_population: dict[int, str],
    target_symbols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    print(f"Processing {matrix_name}: {path.name}", flush=True)
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        cell_labels = read_h5ad_column(obs, "cell_label")
        libraries = read_h5ad_column(obs, "library_label")
        aliases = membership.reindex(pd.Index(cell_labels, dtype="string")).to_numpy()
        populations = np.asarray(
            [alias_to_population.get(int(x)) if pd.notna(x) else None for x in aliases],
            dtype=object,
        )

        gene_symbols = read_h5ad_column(handle["var"], "gene_symbol")
        symbol_to_index: dict[str, int] = {}
        for index, symbol in enumerate(gene_symbols):
            symbol_to_index.setdefault(str(symbol), index)
        present_symbols = [symbol for symbol in target_symbols if symbol in symbol_to_index]
        gene_indices = np.asarray([symbol_to_index[symbol] for symbol in present_symbols], dtype=int)

        selected = np.asarray([value is not None for value in populations])
        group_keys = np.asarray(
            [f"{populations[i]}\t{libraries[i]}" if selected[i] else "" for i in range(len(selected))],
            dtype=object,
        )
        unique_keys = sorted(set(group_keys[selected]))
        group_to_code = {key: code for code, key in enumerate(unique_keys)}
        group_codes = np.full(len(selected), -1, dtype=np.int32)
        group_codes[selected] = [group_to_code[key] for key in group_keys[selected]]

        sums = np.zeros((len(unique_keys), len(present_symbols)), dtype=np.float64)
        detections = np.zeros_like(sums)
        cell_counts = np.zeros(len(unique_keys), dtype=np.int64)

        matrix = handle["X"]
        indptr_ds = matrix["indptr"]
        data_ds = matrix["data"]
        indices_ds = matrix["indices"]
        n_cells, n_genes = map(int, matrix.attrs["shape"])

        for start in range(0, n_cells, CHUNK_ROWS):
            stop = min(start + CHUNK_ROWS, n_cells)
            local_codes = group_codes[start:stop]
            if np.all(local_codes < 0):
                continue
            indptr = indptr_ds[start : stop + 1]
            data_start, data_stop = int(indptr[0]), int(indptr[-1])
            chunk = sparse.csr_matrix(
                (
                    data_ds[data_start:data_stop],
                    indices_ds[data_start:data_stop],
                    indptr - data_start,
                ),
                shape=(stop - start, n_genes),
            )
            keep_rows = np.flatnonzero(local_codes >= 0)
            selected_chunk = chunk[keep_rows][:, gene_indices]
            kept_codes = local_codes[keep_rows]
            for code in np.unique(kept_codes):
                rows = np.flatnonzero(kept_codes == code)
                block = selected_chunk[rows]
                sums[code] += np.asarray(block.sum(axis=0)).ravel()
                detections[code] += block.getnnz(axis=0)
                cell_counts[code] += len(rows)

        count_rows: list[dict[str, object]] = []
        expression_rows: list[dict[str, object]] = []
        for code, key in enumerate(unique_keys):
            population, library = key.split("\t", 1)
            count_rows.append(
                {
                    "matrix": matrix_name,
                    "population": population,
                    "library": library,
                    "n_cells": int(cell_counts[code]),
                }
            )
            if cell_counts[code] < MIN_LIBRARY_CELLS:
                continue
            mean_values = sums[code] / cell_counts[code]
            detection_values = detections[code] / cell_counts[code]
            for symbol, mean_value, detection_value in zip(
                present_symbols, mean_values, detection_values, strict=True
            ):
                expression_rows.append(
                    {
                        "matrix": matrix_name,
                        "population": population,
                        "library": library,
                        "gene_symbol": symbol,
                        "n_cells": int(cell_counts[code]),
                        "mean_log2_expression": float(mean_value),
                        "detection_fraction": float(detection_value),
                    }
                )
    return pd.DataFrame(count_rows), pd.DataFrame(expression_rows)


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.nan, np.nan
    draws = rng.choice(values, size=(N_BOOTSTRAP, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]).tolist())


def build_module_scores(
    expression: pd.DataFrame, modules: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wide = expression.pivot_table(
        index=["matrix", "population", "library"],
        columns="gene_symbol",
        values="mean_log2_expression",
        aggfunc="first",
    )
    gene_means = wide.mean(axis=0)
    gene_sds = wide.std(axis=0, ddof=0).replace(0, np.nan)
    zwide = (wide - gene_means) / gene_sds

    module_rows: list[dict[str, object]] = []
    for module_id, definition in modules.groupby("module_id", sort=False):
        genes = [gene for gene in definition["mouse_symbol"] if gene in zwide.columns]
        values = zwide[genes].mean(axis=1, skipna=True)
        label = definition["module_label"].iloc[0]
        family = definition["module_family"].iloc[0]
        for index, value in values.items():
            matrix, population, library = index
            module_rows.append(
                {
                    "matrix": matrix,
                    "population": population,
                    "library": library,
                    "module_id": module_id,
                    "module_label": label,
                    "module_family": family,
                    "n_module_genes": len(genes),
                    "module_z_score": float(value),
                }
            )
    library_modules = pd.DataFrame(module_rows)

    rng = np.random.default_rng(RNG_SEED)
    population_module_rows: list[dict[str, object]] = []
    for keys, sub in library_modules.groupby(
        ["population", "module_id", "module_label", "module_family"], sort=False
    ):
        population, module_id, module_label, module_family = keys
        values = sub["module_z_score"].dropna().to_numpy()
        lo, hi = bootstrap_mean_ci(values, rng)
        population_module_rows.append(
            {
                "population": population,
                "module_id": module_id,
                "module_label": module_label,
                "module_family": module_family,
                "n_libraries": len(values),
                "mean_module_z_score": float(np.mean(values)),
                "bootstrap_95ci_low": lo,
                "bootstrap_95ci_high": hi,
            }
        )
    population_modules = pd.DataFrame(population_module_rows)

    pop_gene = (
        expression.groupby(["population", "gene_symbol"], as_index=False)
        .agg(
            n_libraries=("library", "nunique"),
            mean_log2_expression=("mean_log2_expression", "mean"),
            median_log2_expression=("mean_log2_expression", "median"),
            mean_detection_fraction=("detection_fraction", "mean"),
        )
    )
    return library_modules, population_modules, pop_gene, zwide


def safe_spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3 or np.nanstd(left) == 0 or np.nanstd(right) == 0:
        return np.nan
    return float(stats.spearmanr(left, right, nan_policy="omit").statistic)


def pair_similarity(
    zwide: pd.DataFrame,
    modules: pd.DataFrame,
    tiers: pd.DataFrame,
) -> pd.DataFrame:
    by_population: dict[str, pd.DataFrame] = {}
    for population in zwide.index.get_level_values("population").unique():
        by_population[population] = zwide.xs(population, level="population")

    upstream = modules.loc[
        modules["module_family"].ne("downstream_circuit_morphology"), "mouse_symbol"
    ].tolist()
    downstream = modules.loc[
        modules["module_family"].eq("downstream_circuit_morphology"), "mouse_symbol"
    ].tolist()
    tier1 = tiers.loc[tiers["manuscript_tier"].eq(TIER1), "mouse_symbol"].astype(str).tolist()
    tier12 = tiers["mouse_symbol"].astype(str).tolist()
    gene_sets = {
        "all curated modules": modules["mouse_symbol"].tolist(),
        "upstream fate/niche": upstream,
        "downstream assembly": downstream,
        "Tier 1 candidates": tier1,
        "Tier 1+2 candidates": tier12,
    }

    rng = np.random.default_rng(RNG_SEED + 1)
    rows: list[dict[str, object]] = []
    for left_name, right_name, pair_class in PAIR_ORDER:
        if left_name not in by_population or right_name not in by_population:
            continue
        left_df = by_population[left_name]
        right_df = by_population[right_name]
        for set_name, requested_genes in gene_sets.items():
            genes = [gene for gene in dict.fromkeys(requested_genes) if gene in zwide.columns]
            left_vector = left_df[genes].mean(axis=0).to_numpy(dtype=float)
            right_vector = right_df[genes].mean(axis=0).to_numpy(dtype=float)
            point = safe_spearman(left_vector, right_vector)
            boots = np.full(N_BOOTSTRAP, np.nan)
            for index in range(N_BOOTSTRAP):
                left_sample = left_df.iloc[rng.integers(0, len(left_df), len(left_df))]
                right_sample = right_df.iloc[rng.integers(0, len(right_df), len(right_df))]
                gene_positions = rng.integers(0, len(genes), len(genes))
                left_boot = left_sample[genes].mean(axis=0).to_numpy()[gene_positions]
                right_boot = right_sample[genes].mean(axis=0).to_numpy()[gene_positions]
                boots[index] = safe_spearman(left_boot, right_boot)
            valid = boots[np.isfinite(boots)]
            lo, hi = (np.nan, np.nan) if len(valid) == 0 else np.quantile(valid, [0.025, 0.975])
            rows.append(
                {
                    "left_population": left_name,
                    "right_population": right_name,
                    "pair_class": pair_class,
                    "gene_set": set_name,
                    "n_genes": len(genes),
                    "n_left_libraries": len(left_df),
                    "n_right_libraries": len(right_df),
                    "spearman_similarity": point,
                    "nested_bootstrap_95ci_low": float(lo),
                    "nested_bootstrap_95ci_high": float(hi),
                }
            )
    result = pd.DataFrame(rows)
    upstream_rows = result[result["gene_set"].eq("upstream fate/niche")][
        ["left_population", "right_population", "spearman_similarity"]
    ].rename(columns={"spearman_similarity": "upstream_similarity"})
    result = result.merge(upstream_rows, on=["left_population", "right_population"], how="left")
    result["downstream_minus_upstream_similarity"] = np.where(
        result["gene_set"].eq("downstream assembly"),
        result["spearman_similarity"] - result["upstream_similarity"],
        np.nan,
    )
    return result


def candidate_specificity(pop_gene: pd.DataFrame, tiers: pd.DataFrame) -> pd.DataFrame:
    pivot = pop_gene.pivot(index="gene_symbol", columns="population", values="mean_log2_expression")
    tier_rows = tiers[["canonical_gene", "manuscript_tier", "mechanism_class"]].copy()
    tier_rows["gene_symbol"] = tiers["mouse_symbol"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []
    targets = ["Cerebellar granule", "Dentate granule, mature"]
    comparators = [name for name in POPULATION_ORDER if name not in targets]
    for row in tier_rows.itertuples(index=False):
        if row.gene_symbol not in pivot.index:
            continue
        values = pivot.loc[row.gene_symbol]
        target_values = values.reindex(targets).dropna()
        comparator_values = values.reindex(comparators).dropna()
        target_min = float(target_values.min()) if len(target_values) == 2 else np.nan
        comparator_max = float(comparator_values.max()) if len(comparator_values) else np.nan
        rows.append(
            {
                "canonical_gene": row.canonical_gene,
                "gene_symbol": row.gene_symbol,
                "manuscript_tier": row.manuscript_tier,
                "mechanism_class": row.mechanism_class,
                "minimum_target_mean_log2": target_min,
                "maximum_comparator_mean_log2": comparator_max,
                "target_min_minus_comparator_max": target_min - comparator_max,
                "strict_target_pair_specific": bool(target_min > comparator_max),
            }
        )
    return pd.DataFrame(rows)


def local_contrast_table(
    wide: pd.DataFrame,
    value_label: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Contrast each granule population with its own regional comparator."""
    required = [
        "Cerebellar granule",
        "Purkinje",
        "Dentate granule, mature",
        "CA1/ProS pyramidal",
        "CA3 pyramidal",
    ]
    missing = [name for name in required if name not in wide.index.get_level_values("population")]
    if missing:
        raise ValueError("Missing local-comparator populations: " + ", ".join(missing))

    rows: list[dict[str, object]] = []
    feature_level = [name for name in wide.index.names if name != "population"][0]
    features = wide.index.get_level_values(feature_level).unique()
    for feature in features:
        values: dict[str, np.ndarray] = {}
        for population in required:
            series = wide.xs((population, feature), level=("population", feature_level))
            values[population] = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
        if any(len(values[name]) == 0 for name in required):
            continue

        cb_delta = float(values["Cerebellar granule"].mean() - values["Purkinje"].mean())
        dg_reference = (
            values["CA1/ProS pyramidal"].mean() + values["CA3 pyramidal"].mean()
        ) / 2
        dg_delta = float(values["Dentate granule, mature"].mean() - dg_reference)

        cb_boot = rng.choice(
            values["Cerebellar granule"],
            size=(N_BOOTSTRAP, len(values["Cerebellar granule"])),
            replace=True,
        ).mean(axis=1) - rng.choice(
            values["Purkinje"],
            size=(N_BOOTSTRAP, len(values["Purkinje"])),
            replace=True,
        ).mean(axis=1)
        dg_boot = rng.choice(
            values["Dentate granule, mature"],
            size=(N_BOOTSTRAP, len(values["Dentate granule, mature"])),
            replace=True,
        ).mean(axis=1) - (
            rng.choice(
                values["CA1/ProS pyramidal"],
                size=(N_BOOTSTRAP, len(values["CA1/ProS pyramidal"])),
                replace=True,
            ).mean(axis=1)
            + rng.choice(
                values["CA3 pyramidal"],
                size=(N_BOOTSTRAP, len(values["CA3 pyramidal"])),
                replace=True,
            ).mean(axis=1)
        ) / 2
        shared_boot = np.minimum(cb_boot, dg_boot)
        rows.append(
            {
                feature_level: feature,
                "value_scale": value_label,
                "cerebellar_granule_minus_purkinje": cb_delta,
                "cerebellar_bootstrap_95ci_low": float(np.quantile(cb_boot, 0.025)),
                "cerebellar_bootstrap_95ci_high": float(np.quantile(cb_boot, 0.975)),
                "dentate_granule_minus_mean_ca1_ca3": dg_delta,
                "dentate_bootstrap_95ci_low": float(np.quantile(dg_boot, 0.025)),
                "dentate_bootstrap_95ci_high": float(np.quantile(dg_boot, 0.975)),
                "shared_minimum_local_delta": min(cb_delta, dg_delta),
                "shared_minimum_bootstrap_95ci_low": float(np.quantile(shared_boot, 0.025)),
                "shared_minimum_bootstrap_95ci_high": float(np.quantile(shared_boot, 0.975)),
                "bootstrap_fraction_both_positive": float(np.mean((cb_boot > 0) & (dg_boot > 0))),
                "both_local_deltas_positive": bool(cb_delta > 0 and dg_delta > 0),
                "n_cb_granule_libraries": len(values["Cerebellar granule"]),
                "n_purkinje_libraries": len(values["Purkinje"]),
                "n_dg_granule_libraries": len(values["Dentate granule, mature"]),
                "n_ca1_libraries": len(values["CA1/ProS pyramidal"]),
                "n_ca3_libraries": len(values["CA3 pyramidal"]),
            }
        )
    return pd.DataFrame(rows)


def build_local_contrasts(
    expression: pd.DataFrame,
    library_modules: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_wide = expression.set_index(["population", "gene_symbol", "library"])[
        "mean_log2_expression"
    ]
    module_wide = library_modules.set_index(["population", "module_id", "library"])[
        "module_z_score"
    ]
    gene_local = local_contrast_table(
        gene_wide, "mean Allen log2 expression", np.random.default_rng(RNG_SEED + 2)
    )
    module_local = local_contrast_table(
        module_wide, "mean curated-gene z score", np.random.default_rng(RNG_SEED + 3)
    )
    labels = library_modules[["module_id", "module_label", "module_family"]].drop_duplicates()
    module_local = module_local.merge(labels, on="module_id", how="left")
    return gene_local, module_local


def matched_null_analysis(
    expression: pd.DataFrame,
    gene_local: pd.DataFrame,
    tiers: pd.DataFrame,
) -> pd.DataFrame:
    gene_features = expression.groupby("gene_symbol", as_index=True).agg(
        mean_expression=("mean_log2_expression", "mean"),
        mean_detection=("detection_fraction", "mean"),
    )
    feature_z = gene_features.sub(gene_features.mean()).div(
        gene_features.std(ddof=0).replace(0, np.nan)
    ).fillna(0)
    local = gene_local.set_index("gene_symbol")
    all_candidates = set(tiers["mouse_symbol"].astype(str))
    pool = [gene for gene in feature_z.index if gene not in all_candidates and gene in local.index]
    rng = np.random.default_rng(RNG_SEED + 4)
    rows: list[dict[str, object]] = []

    tier_sets = {
        "Tier 1": tiers.loc[tiers["manuscript_tier"].eq(TIER1), "mouse_symbol"].astype(str).tolist(),
        "Tier 1+2": tiers["mouse_symbol"].astype(str).tolist(),
    }
    for label, genes in tier_sets.items():
        genes = [gene for gene in genes if gene in local.index]
        nearest: dict[str, list[str]] = {}
        for gene in genes:
            distances = ((feature_z.loc[pool] - feature_z.loc[gene]) ** 2).sum(axis=1)
            nearest[gene] = distances.nsmallest(min(10, len(distances))).index.tolist()

        observed_fraction = float(local.loc[genes, "both_local_deltas_positive"].mean())
        observed_median = float(local.loc[genes, "shared_minimum_local_delta"].median())
        null_fraction = np.zeros(10_000)
        null_median = np.zeros(10_000)
        for index in range(10_000):
            sampled = [rng.choice(nearest[gene]) for gene in genes]
            null_fraction[index] = local.loc[sampled, "both_local_deltas_positive"].mean()
            null_median[index] = local.loc[sampled, "shared_minimum_local_delta"].median()
        positive_count = int(local.loc[genes, "both_local_deltas_positive"].sum())
        rows.append(
            {
                "candidate_set": label,
                "n_candidates": len(genes),
                "n_both_local_deltas_positive": positive_count,
                "observed_both_positive_fraction": observed_fraction,
                "exact_sign_p_greater": float(
                    stats.binomtest(positive_count, len(genes), 0.5, alternative="greater").pvalue
                ),
                "matched_null_median_both_positive_fraction": float(np.median(null_fraction)),
                "matched_null_95ci_low_both_positive_fraction": float(np.quantile(null_fraction, 0.025)),
                "matched_null_95ci_high_both_positive_fraction": float(np.quantile(null_fraction, 0.975)),
                "empirical_p_both_positive_fraction": float(
                    (1 + np.sum(null_fraction >= observed_fraction)) / (1 + len(null_fraction))
                ),
                "observed_median_shared_minimum_delta": observed_median,
                "matched_null_median_shared_minimum_delta": float(np.median(null_median)),
                "matched_null_95ci_low_shared_minimum_delta": float(np.quantile(null_median, 0.025)),
                "matched_null_95ci_high_shared_minimum_delta": float(np.quantile(null_median, 0.975)),
                "empirical_p_shared_minimum_delta": float(
                    (1 + np.sum(null_median >= observed_median)) / (1 + len(null_median))
                ),
                "null_matching": "10 nearest noncandidate curated genes by overall mean expression and detection",
            }
        )
    return pd.DataFrame(rows)


def build_figure(
    counts: pd.DataFrame,
    pop_modules: pd.DataFrame,
    pop_gene: pd.DataFrame,
    similarity: pd.DataFrame,
    gene_local: pd.DataFrame,
    tiers: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)

    ax = axes[0, 0]
    ax.text(-0.06, 1.03, "a", transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
    count_summary = counts[counts["n_cells"].ge(MIN_LIBRARY_CELLS)].groupby("population").agg(
        n_cells=("n_cells", "sum"), n_libraries=("library", "nunique")
    ).reindex(POPULATION_ORDER).dropna()
    display_index = [DISPLAY_NAMES[name] for name in count_summary.index]
    bars = ax.barh(
        display_index,
        count_summary["n_libraries"],
        color=["#1B9E77" if "granule" in name.lower() else "#5B6C7D" for name in count_summary.index],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Independent libraries with at least 50 cells")
    ax.set_ylabel("")
    for bar, cells in zip(bars, count_summary["n_cells"], strict=True):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f"  {int(cells):,} cells",
            va="center",
            fontsize=9,
        )

    ax = axes[0, 1]
    ax.text(-0.06, 1.03, "b", transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
    module_matrix = pop_modules.pivot(
        index="module_label", columns="population", values="mean_module_z_score"
    ).reindex(columns=POPULATION_ORDER)
    module_order = [
        "Cerebellar fate/rhombic-lip/SHH",
        "Dentate fate/WNT/PROX1",
        "Shared neurogenic niche/progenitor state",
        "Downstream neurite/morphology",
        "Downstream synaptic/excitability",
    ]
    module_matrix = module_matrix.reindex(module_order)
    module_matrix.columns = [DISPLAY_NAMES[name] for name in module_matrix.columns]
    sns.heatmap(
        module_matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 8},
        cbar_kws={"label": "Mean library-level gene z score", "shrink": 0.7},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45, labelsize=8)

    ax = axes[1, 0]
    ax.text(-0.06, 1.03, "c", transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
    tier_plot = tiers[["gene", "mouse_symbol", "manuscript_tier"]].merge(
        gene_local, left_on="mouse_symbol", right_on="gene_symbol", how="left"
    )
    palette = tier_plot["manuscript_tier"].map({TIER1: "#1B9E77", TIER2: "#7570B3"})
    ax.scatter(
        tier_plot["cerebellar_granule_minus_purkinje"],
        tier_plot["dentate_granule_minus_mean_ca1_ca3"],
        c=palette,
        s=58,
        edgecolor="white",
        linewidth=0.6,
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.8)
    annotation_offsets = {
        "GABRA2": (5, 10),
        "GRIN2B": (5, -18),
        "GABRB3": (15, -2),
        "GPM6A": (6, 8),
        "PPP3CA": (4, 8),
        "KCND2": (4, -12),
        "CACNA2D1": (4, -11),
        "NFIB": (4, 7),
    }
    for row in tier_plot.itertuples(index=False):
        offset = annotation_offsets.get(row.gene, (4, 3))
        ax.annotate(
            row.gene,
            (row.cerebellar_granule_minus_purkinje, row.dentate_granule_minus_mean_ca1_ca3),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("Cerebellar granule minus Purkinje\nmean log2 expression")
    ax.set_ylabel("Dentate granule minus mean CA1/CA3\nmean log2 expression")

    ax = axes[1, 1]
    ax.text(-0.06, 1.03, "d", transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
    downstream = similarity[similarity["gene_set"].eq("downstream assembly")].copy()
    downstream["pair"] = (
        downstream["left_population"].map(DISPLAY_NAMES)
        + "\nvs "
        + downstream["right_population"].map(DISPLAY_NAMES)
    )
    downstream = downstream.sort_values("downstream_minus_upstream_similarity", ascending=True)
    colors = ["#D95F02" if value == "target granule pair" else "#6B7280" for value in downstream["pair_class"]]
    ax.barh(downstream["pair"], downstream["downstream_minus_upstream_similarity"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Downstream minus upstream Spearman similarity")
    ax.set_ylabel("")

    fig.savefig(OUT_FIGURE, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(
    counts: pd.DataFrame,
    pop_modules: pd.DataFrame,
    similarity: pd.DataFrame,
    specificity: pd.DataFrame,
    gene_local: pd.DataFrame,
    module_local: pd.DataFrame,
    matched_null: pd.DataFrame,
    tiers: pd.DataFrame,
) -> None:
    target = similarity[
        similarity["pair_class"].eq("target granule pair")
        & similarity["gene_set"].isin(["upstream fate/niche", "downstream assembly"])
    ].set_index("gene_set")
    downstream = float(target.loc["downstream assembly", "spearman_similarity"])
    upstream = float(target.loc["upstream fate/niche", "spearman_similarity"])
    strict_n = int(specificity["strict_target_pair_specific"].sum())
    candidate_symbols = tiers["mouse_symbol"].astype(str)
    candidate_local = gene_local[gene_local["gene_symbol"].isin(candidate_symbols)]
    tier1_symbols = tiers.loc[tiers["manuscript_tier"].eq(TIER1), "mouse_symbol"].astype(str)
    tier1_local = gene_local[gene_local["gene_symbol"].isin(tier1_symbols)]
    tier1_null = matched_null[matched_null["candidate_set"].eq("Tier 1")].iloc[0]
    tier12_null = matched_null[matched_null["candidate_set"].eq("Tier 1+2")].iloc[0]
    local_module_text = "; ".join(
        f"{row.module_label}: CB {row.cerebellar_granule_minus_purkinje:.2f}, DG {row.dentate_granule_minus_mean_ca1_ca3:.2f}"
        for row in module_local.itertuples(index=False)
    )
    library_summary = (
        counts[counts["n_cells"].ge(MIN_LIBRARY_CELLS)]
        .groupby("population")["library"]
        .nunique()
        .reindex(POPULATION_ORDER)
        .dropna()
    )
    lines = [
        "# Allen consensus common-matrix comparator analysis",
        "",
        "## Scope",
        "",
        "This analysis uses adult mouse Consensus-WMB-Macosko-10X log2 expression and the integrated 2025 Allen taxonomy. It is an external common-platform specificity test, not a matched developmental replication.",
        "",
        "## Main observations",
        "",
        f"- Direct adult target-pair similarity did not reproduce the proposed downstream-over-upstream pattern: the cerebellar-granule versus mature-dentate-granule Spearman similarity was {downstream:.3f} across downstream assembly genes and {upstream:.3f} across upstream fate/niche genes (difference {downstream - upstream:.3f}).",
        f"- In the comparator-relative analysis that most closely matches the primary rank-delta design, {int(tier1_local['both_local_deltas_positive'].sum())}/{len(tier1_local)} Tier 1 and {int(candidate_local['both_local_deltas_positive'].sum())}/{len(candidate_local)} Tier 1/2 genes were higher in cerebellar granule than Purkinje cells and in dentate granule than the mean CA1/CA3 reference.",
        f"- Expression/detection-matched external null tests gave empirical p={tier1_null.empirical_p_both_positive_fraction:.4g} for Tier 1 and p={tier12_null.empirical_p_both_positive_fraction:.4g} for Tier 1/2 dual-positive fractions.",
        f"- Only {strict_n}/{len(specificity)} Tier 1/2 candidates had mean expression in both target granule populations above every tested comparator. Candidate tiers therefore remain prioritization sets, not a granule-cell-exclusive marker set.",
        f"- Branch-local curated-module deltas were heterogeneous rather than uniformly convergent ({local_module_text}).",
        "- All uncertainty summaries use library-population pseudobulks. Cell counts are coverage measures, not replicate counts.",
        "- Olfactory-bulb groups are labeled as GABAergic proxies because the integrated transcriptomic taxonomy does not by itself distinguish classical granule morphology from other local interneuron morphologies.",
        "",
        "## Independent library counts",
        "",
    ]
    lines.extend(f"- {population}: {int(value)}" for population, value in library_summary.items())
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The Allen layer changes the interpretation. It independently supports recurrence of most preselected Tier 1/2 candidates when each target is contrasted with its own regional comparator, but it does not support direct adult transcriptomic convergence across the broad curated downstream modules. The defensible claim is therefore comparator-relative reuse of a limited candidate set within otherwise distinct adult molecular configurations. Developmental causality, a unique morphology program and broad adult target-pair convergence are not established.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_COUNTS.relative_to(ROOT)}`",
            f"- `{OUT_LIBRARY_GENE.relative_to(ROOT)}`",
            f"- `{OUT_POP_GENE.relative_to(ROOT)}`",
            f"- `{OUT_LIBRARY_MODULE.relative_to(ROOT)}`",
            f"- `{OUT_POP_MODULE.relative_to(ROOT)}`",
            f"- `{OUT_SIMILARITY.relative_to(ROOT)}`",
            f"- `{OUT_SPECIFICITY.relative_to(ROOT)}`",
            f"- `{OUT_LOCAL_GENE.relative_to(ROOT)}`",
            f"- `{OUT_LOCAL_MODULE.relative_to(ROOT)}`",
            f"- `{OUT_MATCHED_NULL.relative_to(ROOT)}`",
            f"- `{OUT_FIGURE.relative_to(ROOT)}`",
        ]
    )
    OUT_REPORT.write_text("\n".join(lines) + "\n")


def main() -> None:
    missing = [path for path in [MEMBERSHIP, TAXONOMY, *MATRICES.values()] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing Allen files: " + ", ".join(map(str, missing)))

    modules, tiers, target_symbols = load_gene_sets()
    taxonomy = build_taxonomy_map()
    alias_to_population = {
        int(alias): population
        for alias, subclass in taxonomy.get("subclass", pd.Series(dtype=object)).items()
        if (population := population_from_subclass(subclass)) is not None
    }
    if OUT_COUNTS.exists() and OUT_LIBRARY_GENE.exists():
        print("Reusing cached library-level Allen extraction", flush=True)
        counts = pd.read_csv(OUT_COUNTS, sep="\t")
        expression = pd.read_csv(OUT_LIBRARY_GENE, sep="\t")
    else:
        membership = load_macosko_membership()
        print(f"Loaded {len(membership):,} Macosko cell-to-cluster mappings", flush=True)
        count_frames: list[pd.DataFrame] = []
        expression_frames: list[pd.DataFrame] = []
        for matrix_name, path in MATRICES.items():
            counts_part, expression_part = aggregate_matrix(
                matrix_name, path, membership, alias_to_population, target_symbols
            )
            count_frames.append(counts_part)
            expression_frames.append(expression_part)
        counts = pd.concat(count_frames, ignore_index=True)
        expression = pd.concat(expression_frames, ignore_index=True)
        counts.to_csv(OUT_COUNTS, sep="\t", index=False)
        expression.to_csv(OUT_LIBRARY_GENE, sep="\t", index=False, compression="gzip")

    library_modules, pop_modules, pop_gene, zwide = build_module_scores(expression, modules)
    similarity = pair_similarity(zwide, modules, tiers)
    specificity = candidate_specificity(pop_gene, tiers)
    gene_local, module_local = build_local_contrasts(expression, library_modules)
    matched_null = matched_null_analysis(expression, gene_local, tiers)

    library_modules.to_csv(OUT_LIBRARY_MODULE, sep="\t", index=False)
    pop_modules.to_csv(OUT_POP_MODULE, sep="\t", index=False)
    pop_gene.to_csv(OUT_POP_GENE, sep="\t", index=False)
    similarity.to_csv(OUT_SIMILARITY, sep="\t", index=False)
    specificity.to_csv(OUT_SPECIFICITY, sep="\t", index=False)
    gene_local.to_csv(OUT_LOCAL_GENE, sep="\t", index=False)
    module_local.to_csv(OUT_LOCAL_MODULE, sep="\t", index=False)
    matched_null.to_csv(OUT_MATCHED_NULL, sep="\t", index=False)
    build_figure(counts, pop_modules, pop_gene, similarity, gene_local, tiers)
    write_report(
        counts,
        pop_modules,
        similarity,
        specificity,
        gene_local,
        module_local,
        matched_null,
        tiers,
    )
    print(f"Wrote {OUT_REPORT}", flush=True)
    print(f"Wrote {OUT_FIGURE}", flush=True)


if __name__ == "__main__":
    main()
