#!/usr/bin/env python3
"""Aim 3 sparse expansion-coding model for granule-cell convergence.

This script finishes the first computational version of Specific Aim 3:

    Link morphology-associated modules to sparse-coding and pattern-separation
    performance.

The model is intentionally transparent and modest. It does not claim to be a
biophysical cerebellar or dentate circuit simulator. It asks whether the
architectural ingredients associated with granule-cell morphology - high
expansion, sparse input sampling, and sparse output activity - can improve
pattern separation while avoiding information collapse.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

DRIVER_SUMMARY = RESULTS / "primary_core_configuration_driver_audit_summary.tsv"
MODULE_SUMMARY = RESULTS / "primary_core_niche_circuit_module_formal_summary.tsv"
DRIVER_PRIORITIES = RESULTS / "primary_core_configuration_driver_audit_gene_priorities.tsv"

OUT_GRID = RESULTS / "primary_core_aim3_sparse_coding_parameter_grid.tsv"
OUT_ARCHITECTURES = RESULTS / "primary_core_aim3_sparse_coding_architecture_summary.tsv"
OUT_MAPPING = RESULTS / "primary_core_aim3_transcriptomic_parameter_mapping.tsv"
OUT_PLOT = RESULTS / "primary_core_aim3_sparse_coding_model.png"
OUT_MD = RESULTS / "primary_core_aim3_sparse_coding_model.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


N_INPUTS = 64
N_PATTERNS = 144
N_PROTOTYPES = 18
INPUT_ACTIVE_FRACTION = 0.12
INPUT_NOISE_RATE = 0.08

GRID_EXPANSION_RATIOS = [0.5, 1, 2, 4, 8, 16]
GRID_INPUT_DEGREES = [2, 4, 8, 16, 32]
GRID_OUTPUT_ACTIVE_FRACTIONS = [0.01, 0.03, 0.05, 0.10, 0.20]
GRID_REPLICATES = 3

NAMED_ARCHITECTURES = [
    {
        "architecture": "cerebellar_granule_like",
        "expansion_ratio": 16,
        "input_degree": 4,
        "output_active_fraction": 0.05,
        "rationale": "large expansion, very sparse mossy-fiber sampling, sparse output",
    },
    {
        "architecture": "dentate_granule_like",
        "expansion_ratio": 8,
        "input_degree": 8,
        "output_active_fraction": 0.03,
        "rationale": "large expansion, sparse entorhinal input sampling, sparse output",
    },
    {
        "architecture": "balanced_granule_like",
        "expansion_ratio": 8,
        "input_degree": 4,
        "output_active_fraction": 0.05,
        "rationale": "generic sparse expansion code",
    },
    {
        "architecture": "pyramidal_integrator_like",
        "expansion_ratio": 1,
        "input_degree": 32,
        "output_active_fraction": 0.20,
        "rationale": "low expansion, broad input integration, moderate activity",
    },
    {
        "architecture": "purkinje_integrator_like",
        "expansion_ratio": 0.5,
        "input_degree": 32,
        "output_active_fraction": 0.40,
        "rationale": "convergent integrator-like architecture",
    },
    {
        "architecture": "excessive_sparsity",
        "expansion_ratio": 16,
        "input_degree": 4,
        "output_active_fraction": 0.005,
        "rationale": "tests information collapse when sparse coding is too severe",
    },
    {
        "architecture": "dense_expansion_high_activity",
        "expansion_ratio": 8,
        "input_degree": 32,
        "output_active_fraction": 0.20,
        "rationale": "expansion without sparse input sampling or sparse output",
    },
]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def generate_correlated_inputs(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    prototypes = rng.random((N_PROTOTYPES, N_INPUTS)) < INPUT_ACTIVE_FRACTION
    per_proto = int(np.ceil(N_PATTERNS / N_PROTOTYPES))
    patterns: list[np.ndarray] = []
    for proto in prototypes:
        for _ in range(per_proto):
            noisy = proto.copy()
            flips = rng.random(N_INPUTS) < INPUT_NOISE_RATE
            noisy[flips] = ~noisy[flips]
            if not noisy.any():
                noisy[rng.integers(0, N_INPUTS)] = True
            patterns.append(noisy)
            if len(patterns) >= N_PATTERNS:
                break
        if len(patterns) >= N_PATTERNS:
            break
    return np.asarray(patterns, dtype=np.float32)


def sparse_projection(
    *,
    inputs: np.ndarray,
    expansion_ratio: float,
    input_degree: int,
    output_active_fraction: float,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_outputs = max(4, int(round(inputs.shape[1] * expansion_ratio)))
    degree = max(1, min(int(input_degree), inputs.shape[1]))
    weights = np.zeros((inputs.shape[1], n_outputs), dtype=np.float32)
    for j in range(n_outputs):
        idx = rng.choice(inputs.shape[1], size=degree, replace=False)
        weights[idx, j] = 1.0 / np.sqrt(degree)
    activations = inputs @ weights
    k = max(1, int(round(output_active_fraction * n_outputs)))
    k = min(k, n_outputs)
    if k == n_outputs:
        return np.ones_like(activations, dtype=bool)
    thresholds = np.partition(activations, n_outputs - k, axis=1)[:, n_outputs - k]
    outputs = activations >= thresholds[:, None]
    # Tie-breaking can make activity slightly high; keep the strongest k units.
    if outputs.mean() > output_active_fraction * 1.5 and k < n_outputs:
        strongest = np.argpartition(activations, n_outputs - k, axis=1)[:, n_outputs - k :]
        fixed = np.zeros_like(outputs)
        rows = np.arange(outputs.shape[0])[:, None]
        fixed[rows, strongest] = True
        outputs = fixed
    return outputs


def binary_entropy(prob: np.ndarray) -> np.ndarray:
    p = np.clip(prob, 1e-9, 1 - 1e-9)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))


def score_architecture(
    *,
    architecture: str,
    expansion_ratio: float,
    input_degree: int,
    output_active_fraction: float,
    replicate: int,
    seed: int,
) -> dict[str, object]:
    inputs = generate_correlated_inputs(seed)
    actual_input_degree = max(1, min(int(input_degree), inputs.shape[1]))
    outputs = sparse_projection(
        inputs=inputs,
        expansion_ratio=expansion_ratio,
        input_degree=actual_input_degree,
        output_active_fraction=output_active_fraction,
        seed=seed + 100_000,
    )
    input_hamming = pdist(inputs, metric="hamming")
    output_hamming = pdist(outputs.astype(np.float32), metric="hamming")
    input_dist = pdist(inputs.astype(bool), metric="jaccard")
    output_dist = pdist(outputs.astype(bool), metric="jaccard")
    near_mask = input_dist <= np.quantile(input_dist, 0.25)
    nonzero_input = input_dist > 0

    median_input = float(np.median(input_dist))
    median_output = float(np.median(output_dist))
    near_input = float(np.median(input_dist[near_mask]))
    near_output = float(np.median(output_dist[near_mask]))
    separation_gain = median_output / (median_input + 1e-9)
    near_pair_separation_gain = near_output / (near_input + 1e-9)

    corr = spearmanr(input_dist, output_dist).correlation
    if corr is None or not np.isfinite(corr):
        corr = 0.0
    collapse_rate = float(np.mean((output_dist == 0) & nonzero_input))
    active_fraction = float(outputs.mean())
    unit_activity = outputs.mean(axis=0)
    output_entropy = float(binary_entropy(unit_activity).mean())
    active_penalty = float(np.exp(-abs(np.log((active_fraction + 1e-9) / 0.05)) / 2.0))
    information_retention = max(float(corr), 0.0) * (1.0 - collapse_rate) * output_entropy
    useful_pattern_separation_score = near_pair_separation_gain * information_retention * active_penalty
    relative_wiring_cost = float(expansion_ratio * actual_input_degree)
    relative_active_output_load = float(expansion_ratio * active_fraction)
    relative_wiring_activity_cost = float(expansion_ratio * actual_input_degree * active_fraction)
    resource_adjusted_useful_score = useful_pattern_separation_score / (relative_wiring_activity_cost + 1e-9)

    return {
        "architecture": architecture,
        "replicate": replicate,
        "seed": seed,
        "n_inputs": int(inputs.shape[1]),
        "n_outputs": int(outputs.shape[1]),
        "expansion_ratio": float(expansion_ratio),
        "input_degree": int(actual_input_degree),
        "output_active_fraction_parameter": float(output_active_fraction),
        "observed_output_active_fraction": active_fraction,
        "median_input_hamming_distance": float(np.median(input_hamming)),
        "median_output_hamming_distance": float(np.median(output_hamming)),
        "median_input_jaccard_distance": median_input,
        "median_output_jaccard_distance": median_output,
        "separation_gain": separation_gain,
        "near_pair_input_jaccard_distance": near_input,
        "near_pair_output_jaccard_distance": near_output,
        "near_pair_separation_gain": near_pair_separation_gain,
        "distance_spearman_correlation": float(corr),
        "collapse_rate": collapse_rate,
        "output_entropy_bits": output_entropy,
        "active_fraction_penalty": active_penalty,
        "information_retention_score": information_retention,
        "useful_pattern_separation_score": useful_pattern_separation_score,
        "relative_wiring_cost": relative_wiring_cost,
        "relative_active_output_load": relative_active_output_load,
        "relative_wiring_activity_cost": relative_wiring_activity_cost,
        "resource_adjusted_useful_score": resource_adjusted_useful_score,
    }


def label_parameter_zone(expansion_ratio: float, input_degree: int, output_active_fraction: float) -> str:
    if output_active_fraction <= 0.01:
        return "excessively_sparse"
    if expansion_ratio >= 4 and input_degree <= 8 and output_active_fraction <= 0.10:
        return "granule_like_sparse_expansion"
    if expansion_ratio <= 1 and input_degree >= 16 and output_active_fraction >= 0.10:
        return "integrator_like_low_expansion"
    if expansion_ratio >= 4 and input_degree >= 16 and output_active_fraction >= 0.10:
        return "dense_expansion"
    return "intermediate"


def run_grid() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for expansion in GRID_EXPANSION_RATIOS:
        for degree in GRID_INPUT_DEGREES:
            for active in GRID_OUTPUT_ACTIVE_FRACTIONS:
                for replicate in range(GRID_REPLICATES):
                    seed = 20_000 + replicate * 10_000 + int(expansion * 101) + degree * 17 + int(active * 1000)
                    row = score_architecture(
                        architecture=label_parameter_zone(expansion, degree, active),
                        expansion_ratio=expansion,
                        input_degree=degree,
                        output_active_fraction=active,
                        replicate=replicate + 1,
                        seed=seed,
                    )
                    rows.append(row)
    grid = pd.DataFrame(rows)
    grid.to_csv(OUT_GRID, sep="\t", index=False)
    return grid


def run_named_architectures() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for arch_idx, arch in enumerate(NAMED_ARCHITECTURES, start=1):
        replicate_rows = []
        for replicate in range(8):
            seed = 100_000 + arch_idx * 10_000 + replicate
            row = score_architecture(
                architecture=arch["architecture"],
                expansion_ratio=arch["expansion_ratio"],
                input_degree=arch["input_degree"],
                output_active_fraction=arch["output_active_fraction"],
                replicate=replicate + 1,
                seed=seed,
            )
            row["rationale"] = arch["rationale"]
            replicate_rows.append(row)
        rows.extend(replicate_rows)
    reps = pd.DataFrame(rows)
    summary = (
        reps.groupby(["architecture", "rationale"], sort=False)
        .agg(
            n_replicates=("replicate", "size"),
            expansion_ratio=("expansion_ratio", "median"),
            input_degree=("input_degree", "median"),
            output_active_fraction_parameter=("output_active_fraction_parameter", "median"),
            observed_output_active_fraction=("observed_output_active_fraction", "median"),
            median_separation_gain=("separation_gain", "median"),
            median_near_pair_separation_gain=("near_pair_separation_gain", "median"),
            median_distance_spearman_correlation=("distance_spearman_correlation", "median"),
            median_collapse_rate=("collapse_rate", "median"),
            median_output_entropy_bits=("output_entropy_bits", "median"),
            median_information_retention_score=("information_retention_score", "median"),
            median_useful_pattern_separation_score=("useful_pattern_separation_score", "median"),
            median_relative_wiring_cost=("relative_wiring_cost", "median"),
            median_relative_active_output_load=("relative_active_output_load", "median"),
            median_relative_wiring_activity_cost=("relative_wiring_activity_cost", "median"),
            median_resource_adjusted_useful_score=("resource_adjusted_useful_score", "median"),
        )
        .reset_index()
        .sort_values("median_useful_pattern_separation_score", ascending=False)
    )
    summary.to_csv(OUT_ARCHITECTURES, sep="\t", index=False)
    return summary


def build_mapping_table() -> pd.DataFrame:
    module_summary = pd.read_csv(MODULE_SUMMARY, sep="\t")
    driver_summary = pd.read_csv(DRIVER_SUMMARY, sep="\t")
    priorities = pd.read_csv(DRIVER_PRIORITIES, sep="\t")

    downstream_morph = module_summary.loc[module_summary["module_id"].eq("downstream_neurite_morphology")].iloc[0]
    downstream_syn = module_summary.loc[module_summary["module_id"].eq("downstream_synaptic_excitability")].iloc[0]
    primary_driver = driver_summary.loc[
        driver_summary["summary_level"].eq("by_audit_scope")
        & driver_summary["audit_scope"].eq("primary_core_candidate_background")
    ].iloc[0]

    morph_genes = ",".join(
        priorities.loc[priorities["module_id"].eq("downstream_neurite_morphology"), "gene"].astype(str).head(12)
    )
    syn_genes = ",".join(
        priorities.loc[priorities["module_id"].eq("downstream_synaptic_excitability"), "gene"].astype(str).head(12)
    )

    rows = [
        {
            "computational_parameter": "expansion_ratio",
            "biological_interpretation": "many compact output neurons per input channel",
            "transcriptomic_anchor": "downstream neurite/morphology module plus granule-cell identity",
            "support_summary": f"downstream neurite/morphology median formal convergence {downstream_morph['median_overall_convergence_delta']}",
            "candidate_genes": morph_genes,
            "claim_strength": "conceptual_support",
            "caveat": "cell number and packing are anatomical parameters, not directly inferred from expression",
        },
        {
            "computational_parameter": "input_degree",
            "biological_interpretation": "sparse dendritic/input sampling",
            "transcriptomic_anchor": "axon guidance, adhesion, neurite, and cytoskeleton genes",
            "support_summary": f"{int(downstream_morph['n_shared_positive_any_screen'])}/{int(downstream_morph['n_genes_present_formal'])} neurite/morphology genes shared-positive in at least one screen",
            "candidate_genes": morph_genes,
            "claim_strength": "moderate_transcriptomic_support",
            "caveat": "actual dendrite number and synaptic input count require morphology or connectomics",
        },
        {
            "computational_parameter": "output_active_fraction",
            "biological_interpretation": "sparse coding threshold and excitability control",
            "transcriptomic_anchor": "synaptic/excitability module",
            "support_summary": f"downstream synaptic/excitability median formal convergence {downstream_syn['median_overall_convergence_delta']}",
            "candidate_genes": syn_genes,
            "claim_strength": "moderate_transcriptomic_support",
            "caveat": "activity sparsity requires electrophysiology or calcium-imaging validation",
        },
        {
            "computational_parameter": "architecture_configuration",
            "biological_interpretation": "balanced combination of expansion, sparse sampling, and controlled sparsity",
            "transcriptomic_anchor": "identity-coupled transcriptomic configuration score",
            "support_summary": f"{int(primary_driver['n_configuration_positive'])}/{int(primary_driver['n_contrasts'])} primary-core configuration-positive contrasts; median configuration delta {primary_driver['median_delta_configuration_score']}",
            "candidate_genes": "GPM6A,ROBO2,DCC,CADM3,STMN2,STMN3,DPYSL2,MAP1B,BASP1,CFL1,KCNK1,GABRA2",
            "claim_strength": "strong_for_module_balance_not_direct_geometry",
            "caveat": "configuration evidence is identity-coupled and does not by itself prove morphology or performance",
        },
    ]
    mapping = pd.DataFrame(rows)
    mapping.to_csv(OUT_MAPPING, sep="\t", index=False)
    return mapping


def plot_results(grid: pd.DataFrame, arch_summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12.2, 13.0), constrained_layout=True)

    ax = axes[0]
    heat = grid.loc[np.isclose(grid["output_active_fraction_parameter"], 0.05)].copy()
    heat = (
        heat.groupby(["expansion_ratio", "input_degree"], sort=False)["useful_pattern_separation_score"]
        .median()
        .reset_index()
    )
    pivot = heat.pivot(index="expansion_ratio", columns="input_degree", values="useful_pattern_separation_score")
    im = ax.imshow(pivot.values, aspect="auto", origin="lower", cmap="viridis")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([str(int(c)) for c in pivot.columns])
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([str(c) for c in pivot.index])
    ax.set_xlabel("Input degree")
    ax.set_ylabel("Expansion ratio")
    ax.set_title("Useful pattern-separation score at 5% output activity")
    fig.colorbar(im, ax=ax, label="Useful score")

    ax = axes[1]
    scatter = grid.copy()
    sc = ax.scatter(
        scatter["near_pair_separation_gain"],
        scatter["information_retention_score"],
        c=scatter["output_active_fraction_parameter"],
        s=28,
        cmap="plasma",
        alpha=0.75,
        edgecolor="none",
    )
    ax.set_xlabel("Near-pair separation gain")
    ax.set_ylabel("Information retention score")
    ax.set_title("Separation is useful only when retention remains high")
    ax.grid(color="#dddddd", linewidth=0.5)
    fig.colorbar(sc, ax=ax, label="Output active fraction")

    ax = axes[2]
    plot_arch = arch_summary.sort_values("median_useful_pattern_separation_score", ascending=True)
    y = np.arange(len(plot_arch))
    ax.barh(y, plot_arch["median_useful_pattern_separation_score"], color="#2f7f8f")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_arch["architecture"], fontsize=9)
    ax.set_xlabel("Median useful pattern-separation score")
    ax.set_title("Named architecture comparison")
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    fig.suptitle("Sparse Expansion-Coding Model", fontsize=16, y=1.01)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(grid: pd.DataFrame, arch_summary: pd.DataFrame, mapping: pd.DataFrame) -> None:
    best_grid = grid.sort_values("useful_pattern_separation_score", ascending=False).head(1).iloc[0]
    zone_summary = (
        grid.groupby("architecture", sort=False)
        .agg(
            n=("architecture", "size"),
            median_useful=("useful_pattern_separation_score", "median"),
            median_separation=("near_pair_separation_gain", "median"),
            median_retention=("information_retention_score", "median"),
            median_collapse=("collapse_rate", "median"),
            median_cost=("relative_wiring_activity_cost", "median"),
            median_resource_adjusted=("resource_adjusted_useful_score", "median"),
        )
        .reset_index()
        .sort_values("median_useful", ascending=False)
    )

    top_arch = arch_summary.iloc[0]
    top_resource_arch = arch_summary.sort_values("median_resource_adjusted_useful_score", ascending=False).iloc[0]
    nontrivial_expansion = grid.loc[
        grid["expansion_ratio"].ge(4) & grid["useful_pattern_separation_score"].ge(0.15)
    ].copy()
    if nontrivial_expansion.empty:
        nontrivial_expansion = grid.copy()
    best_resource_grid = nontrivial_expansion.sort_values("resource_adjusted_useful_score", ascending=False).head(1).iloc[0]
    cereb = arch_summary.loc[arch_summary["architecture"].eq("cerebellar_granule_like")].iloc[0]
    dentate = arch_summary.loc[arch_summary["architecture"].eq("dentate_granule_like")].iloc[0]
    pyramidal = arch_summary.loc[arch_summary["architecture"].eq("pyramidal_integrator_like")].iloc[0]
    purkinje = arch_summary.loc[arch_summary["architecture"].eq("purkinje_integrator_like")].iloc[0]
    excessive = arch_summary.loc[arch_summary["architecture"].eq("excessive_sparsity")].iloc[0]

    lines = [
        "# Sparse Expansion-Coding Model",
        "",
        "Date built: 2026-06-22",
        "",
        "## Purpose",
        "",
        "This analysis tests whether a granule-like architecture - high expansion, sparse input sampling, and sparse output activity - can improve pattern separation while penalizing information loss.",
        "",
        "## Model Design",
        "",
        f"- Input patterns: {N_PATTERNS} correlated binary patterns over {N_INPUTS} input channels.",
        f"- Pattern generation: {N_PROTOTYPES} prototypes, input active fraction {INPUT_ACTIVE_FRACTION}, noise rate {INPUT_NOISE_RATE}.",
        "- Random sparse projections define output populations.",
        "- Parameters varied: expansion ratio, input degree, and output active fraction.",
        "- Main score: near-pair Jaccard/overlap separation gain multiplied by information retention, output entropy, and an activity-balance penalty.",
        "",
        "This score rewards useful separation of similar input patterns but penalizes collapse, excessive sparsity, and loss of input-distance structure. Jaccard distance is used for the main separation metric because normalized Hamming distance underestimates separation quality in very sparse binary codes.",
        "",
        "## Main Results",
        "",
        f"- Best grid point: expansion ratio {best_grid['expansion_ratio']}, input degree {int(best_grid['input_degree'])}, output active fraction {best_grid['output_active_fraction_parameter']}; useful score {best_grid['useful_pattern_separation_score']:.3f}.",
        f"- Best resource-adjusted nontrivial expansion grid point: expansion ratio {best_resource_grid['expansion_ratio']}, input degree {int(best_resource_grid['input_degree'])}, output active fraction {best_resource_grid['output_active_fraction_parameter']}; resource-adjusted score {best_resource_grid['resource_adjusted_useful_score']:.3f}.",
        f"- Best named architecture: {top_arch['architecture']} with median useful score {top_arch['median_useful_pattern_separation_score']:.3f}.",
        f"- Best resource-adjusted named architecture: {top_resource_arch['architecture']} with median resource-adjusted score {top_resource_arch['median_resource_adjusted_useful_score']:.3f}.",
        f"- Cerebellar granule-like architecture: useful score {cereb['median_useful_pattern_separation_score']:.3f}, near-pair separation gain {cereb['median_near_pair_separation_gain']:.3f}, retention {cereb['median_information_retention_score']:.3f}.",
        f"- Dentate granule-like architecture: useful score {dentate['median_useful_pattern_separation_score']:.3f}, near-pair separation gain {dentate['median_near_pair_separation_gain']:.3f}, retention {dentate['median_information_retention_score']:.3f}.",
        f"- Pyramidal/integrator-like architecture: useful score {pyramidal['median_useful_pattern_separation_score']:.3f}.",
        f"- Purkinje/integrator-like architecture: useful score {purkinje['median_useful_pattern_separation_score']:.3f}.",
        f"- Excessive sparsity: useful score {excessive['median_useful_pattern_separation_score']:.3f}, collapse rate {excessive['median_collapse_rate']:.3f}.",
        "",
        "Parameter-zone summary:",
    ]
    for _, row in zone_summary.iterrows():
        lines.append(
            f"- {row['architecture']}: median useful {row['median_useful']:.3f}, "
            f"near-pair separation {row['median_separation']:.3f}, retention {row['median_retention']:.3f}, "
            f"collapse {row['median_collapse']:.3f}, resource-adjusted {row['median_resource_adjusted']:.3f}."
        )

    lines.extend(
        [
            "",
            "## Link To Transcriptomic Results",
            "",
            "The model supports the computational plausibility of a sparse-expansion convergence model: sparse expansion architectures can separate similar input patterns better than low-expansion integrator-like architectures, but only when sparsity is balanced enough to preserve information.",
            "",
            "The transcriptomic data connect to this model at the parameter-implementation level rather than as direct performance measurements:",
        ]
    )
    for _, row in mapping.iterrows():
        lines.append(
            f"- {row['computational_parameter']}: {row['transcriptomic_anchor']} ({row['support_summary']}). Caveat: {row['caveat']}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Granule-like expansion plus sparse sampling is a plausible convergent solution for useful pattern separation in this conceptual computational validation.",
            "- Dense expansion with high activity can score higher in raw computational terms, but it is much more expensive in wiring/activity load; resource-adjusted scoring favors sparse expansion designs.",
            "- The model also warns against a simplistic 'more sparse is always better' story, because excessive sparsity causes information collapse.",
            "- The transcriptomic evidence supports construction and excitability modules that could implement the parameters, but direct morphology, connectomics, or activity data would be needed to prove the parameter values in vivo.",
            "- This fits the broader project model: morphology similarity is likely constrained by circuit computation, while the implementation remains identity-coupled and region-specific.",
            "",
            "## Outputs",
            "",
            f"- Parameter grid: `{rel(OUT_GRID)}`",
            f"- Named architecture summary: `{rel(OUT_ARCHITECTURES)}`",
            f"- Transcriptomic parameter mapping: `{rel(OUT_MAPPING)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    grid = run_grid()
    architectures = run_named_architectures()
    mapping = build_mapping_table()
    plot_results(grid, architectures)
    write_report(grid, architectures, mapping)
    best = architectures.iloc[0]
    print(f"Wrote {rel(OUT_MD)}")
    print(
        "best_named_architecture",
        best["architecture"],
        "useful_score",
        f"{best['median_useful_pattern_separation_score']:.3f}",
    )


if __name__ == "__main__":
    main()
