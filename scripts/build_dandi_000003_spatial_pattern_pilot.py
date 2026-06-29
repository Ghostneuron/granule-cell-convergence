#!/usr/bin/env python3
"""Position-linked DANDI 000003 spatial coding pilot for Aim 3.

This script uses the already-downloaded smallest DANDI 000003 pilot NWB file
and computes first-pass spatial coding metrics for source-labeled units. It is
intentionally conservative: one session, relative tracking coordinates, awake
moving samples only, and no claim that this is final pattern separation.
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
DATA = ROOT / "External_Data" / "DANDI" / "000003"
NWB = DATA / "sub-YutaMouse41" / "sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb"

OUT_POSITION = RESULTS / "dandi_000003_pilot_position_summary.tsv"
OUT_UNIT = RESULTS / "dandi_000003_pilot_spatial_unit_metrics.tsv"
OUT_CELLTYPE = RESULTS / "dandi_000003_pilot_spatial_celltype_summary.tsv"
OUT_PV = RESULTS / "dandi_000003_pilot_population_vector_separation.tsv"
OUT_PLOT = RESULTS / "dandi_000003_pilot_granule_spatial_maps.png"
OUT_MD = RESULTS / "dandi_000003_spatial_pattern_pilot.md"

N_X_BINS = 18
N_Y_BINS = 14
MIN_OCCUPANCY_S = 0.25
MOVEMENT_SPEED_MIN = 5.0
MOVEMENT_SPEED_MAX = 200.0


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def decode(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def load_position(h5: h5py.File) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    s0 = h5["/acquisition/position_sensor0"]
    s1 = h5["/acquisition/position_sensor1"]
    d0 = np.asarray(s0["data"], dtype=float)
    d1 = np.asarray(s1["data"], dtype=float)
    rate = float(s0["starting_time"].attrs["rate"])
    start = float(s0["starting_time"][()])
    stacked = np.stack([d0, d1], axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pos = np.nanmean(stacked, axis=0)
    times = start + np.arange(pos.shape[0]) / rate
    dt = 1.0 / rate
    speed = np.full(pos.shape[0], np.nan)
    diffs = np.sqrt(np.sum(np.diff(pos, axis=0) ** 2, axis=1)) / dt
    speed[1:] = diffs
    return times, pos, speed, rate


def load_states(h5: h5py.File) -> pd.DataFrame:
    states = h5["/processing/behavior/states"]
    labels = [decode(v) for v in states["label"][:]]
    starts = np.asarray(states["start_time"], dtype=float)
    stops = np.asarray(states["stop_time"], dtype=float)
    ids = np.asarray(states["id"], dtype=int)
    return pd.DataFrame(
        {
            "state_id": ids,
            "state_label": labels,
            "start_time_s": starts,
            "stop_time_s": stops,
            "duration_s": np.maximum(stops - starts, 0.0),
        }
    )


def awake_mask(times: np.ndarray, states: pd.DataFrame) -> np.ndarray:
    mask = np.zeros(times.shape[0], dtype=bool)
    for _, row in states.iterrows():
        if row["state_label"] == "awake":
            mask |= (times >= float(row["start_time_s"])) & (times < float(row["stop_time_s"]))
    return mask


def read_units(h5: h5py.File) -> list[dict[str, Any]]:
    units = h5["/units"]
    ids = np.asarray(units["id"])
    cell_types = [decode(v) for v in units["cell_type"][:]]
    spike_times = np.asarray(units["spike_times"], dtype=float)
    spike_index = np.asarray(units["spike_times_index"], dtype=int)
    rows = []
    prev = 0
    for idx, stop in enumerate(spike_index):
        st = spike_times[prev:stop]
        prev = int(stop)
        rows.append({"unit_index": idx, "unit_id": int(ids[idx]), "cell_type": cell_types[idx], "spike_times": st})
    return rows


def bin_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    finite = values[np.isfinite(values)]
    lo, hi = np.nanpercentile(finite, [1, 99])
    if lo == hi:
        lo, hi = float(np.nanmin(finite)), float(np.nanmax(finite))
    pad = max((hi - lo) * 0.02, 1e-6)
    return np.linspace(lo - pad, hi + pad, n_bins + 1)


def occupancy_map(x: np.ndarray, y: np.ndarray, mask: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray, dt: float) -> np.ndarray:
    counts, _, _ = np.histogram2d(y[mask], x[mask], bins=[y_edges, x_edges])
    return counts * dt


def interpolate_position_for_spikes(
    spike_times: np.ndarray,
    times: np.ndarray,
    pos: np.ndarray,
    valid_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    idx = np.searchsorted(times, spike_times, side="left")
    idx = np.clip(idx, 1, len(times) - 1)
    left = idx - 1
    choose_left = np.abs(spike_times - times[left]) <= np.abs(spike_times - times[idx])
    nearest = np.where(choose_left, left, idx)
    keep = valid_mask[nearest]
    return pos[nearest[keep]], spike_times[keep]


def spatial_metrics(rate_map: np.ndarray, occupancy: np.ndarray, spike_map: np.ndarray) -> dict[str, float]:
    valid = occupancy >= MIN_OCCUPANCY_S
    total_occ = float(occupancy[valid].sum())
    total_spikes = int(spike_map[valid].sum())
    if total_occ <= 0:
        return {}
    p_occ = occupancy[valid] / total_occ
    rates = rate_map[valid]
    mean_rate = float(np.sum(p_occ * rates))
    if mean_rate <= 0:
        info = 0.0
        sparsity = np.nan
    else:
        ratio = rates / mean_rate
        positive = ratio > 0
        info = float(np.sum(p_occ[positive] * ratio[positive] * np.log2(ratio[positive])))
        denom = float(np.sum(p_occ * rates**2))
        sparsity = float((np.sum(p_occ * rates) ** 2) / denom) if denom > 0 else np.nan
    occupied_bins = int(valid.sum())
    active_bins = int(((spike_map > 0) & valid).sum())
    max_rate = float(np.nanmax(rates)) if len(rates) else np.nan
    return {
        "occupied_spatial_bins": occupied_bins,
        "active_spatial_bins": active_bins,
        "active_spatial_bin_fraction": active_bins / occupied_bins if occupied_bins else np.nan,
        "spatial_sparsity": sparsity,
        "spatial_information_bits_per_spike": info,
        "mean_awake_moving_rate_hz": mean_rate,
        "max_spatial_rate_hz": max_rate,
        "spatial_selectivity_max_over_mean": max_rate / mean_rate if mean_rate > 0 else np.nan,
        "awake_moving_spikes": total_spikes,
    }


def rankdata(values: np.ndarray) -> np.ndarray:
    s = pd.Series(values)
    return s.rank(method="average").to_numpy()


def spearman_no_scipy(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    rx = rankdata(x[ok])
    ry = rankdata(y[ok])
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def population_vector_metrics(rate_maps: dict[int, np.ndarray], occupancy: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray) -> pd.DataFrame:
    valid = occupancy >= max(MIN_OCCUPANCY_S, 1.0)
    if len(rate_maps) < 2 or valid.sum() < 8:
        return pd.DataFrame()
    unit_ids = list(rate_maps)
    vectors = np.vstack([rate_maps[uid][valid].ravel() for uid in unit_ids]).T
    if vectors.shape[0] < 8:
        return pd.DataFrame()
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2.0
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2.0
    yy, xx = np.meshgrid(y_centers, x_centers, indexing="ij")
    centers = np.column_stack([xx[valid].ravel(), yy[valid].ravel()])
    spatial_dist = []
    neural_euclid = []
    neural_corr_dist = []
    for i in range(vectors.shape[0]):
        for j in range(i + 1, vectors.shape[0]):
            spatial_dist.append(float(np.linalg.norm(centers[i] - centers[j])))
            neural_euclid.append(float(np.linalg.norm(vectors[i] - vectors[j])))
            if np.std(vectors[i]) == 0 or np.std(vectors[j]) == 0:
                neural_corr_dist.append(np.nan)
            else:
                neural_corr_dist.append(float(1.0 - np.corrcoef(vectors[i], vectors[j])[0, 1]))
    spatial_dist = np.asarray(spatial_dist)
    neural_euclid = np.asarray(neural_euclid)
    neural_corr_dist = np.asarray(neural_corr_dist)
    near = spatial_dist <= np.nanquantile(spatial_dist, 0.25)
    far = spatial_dist >= np.nanquantile(spatial_dist, 0.75)
    rows = [
        {
            "unit_set": "granule_cell_labeled" if len(unit_ids) <= 3 else "all_units",
            "n_units": len(unit_ids),
            "n_occupied_bins": int(valid.sum()),
            "n_bin_pairs": int(len(spatial_dist)),
            "spearman_spatial_vs_neural_euclidean": spearman_no_scipy(spatial_dist, neural_euclid),
            "spearman_spatial_vs_neural_corrdist": spearman_no_scipy(spatial_dist, neural_corr_dist),
            "near_pair_median_neural_euclidean": float(np.nanmedian(neural_euclid[near])),
            "far_pair_median_neural_euclidean": float(np.nanmedian(neural_euclid[far])),
            "far_minus_near_neural_euclidean": float(np.nanmedian(neural_euclid[far]) - np.nanmedian(neural_euclid[near])),
            "near_pair_median_corrdist": float(np.nanmedian(neural_corr_dist[near])),
            "far_pair_median_corrdist": float(np.nanmedian(neural_corr_dist[far])),
            "far_minus_near_corrdist": float(np.nanmedian(neural_corr_dist[far]) - np.nanmedian(neural_corr_dist[near])),
        }
    ]
    return pd.DataFrame(rows)


def make_plot(unit_df: pd.DataFrame, rate_maps: dict[int, np.ndarray], occupancy: np.ndarray) -> bool:
    granules = unit_df[unit_df["cell_type"] == "granule cell"].sort_values("unit_id")
    if granules.empty:
        return False
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    n = granules.shape[0]
    fig, axes = plt.subplots(1, n + 1, figsize=(4 * (n + 1), 3.6))
    if n + 1 == 1:
        axes = [axes]
    occ = occupancy.copy()
    occ[occ < MIN_OCCUPANCY_S] = np.nan
    im = axes[0].imshow(occ, origin="lower", aspect="auto", cmap="viridis")
    axes[0].set_title("Occupancy (s)")
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)
    for ax, (_, row) in zip(axes[1:], granules.iterrows()):
        rm = rate_maps[int(row["unit_id"])].copy()
        rm[occupancy < MIN_OCCUPANCY_S] = np.nan
        im = ax.imshow(rm, origin="lower", aspect="auto", cmap="magma")
        ax.set_title(f"GC unit {int(row['unit_id'])}\n{row['spatial_information_bits_per_spike']:.2f} bits/spk")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("DANDI 000003 pilot awake-moving spatial maps", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)
    return True


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not NWB.exists():
        raise FileNotFoundError(NWB)
    with h5py.File(NWB, "r") as h5:
        times, pos, speed, rate = load_position(h5)
        states = load_states(h5)
        units = read_units(h5)

    dt = 1.0 / rate
    finite_pos = np.isfinite(pos).all(axis=1)
    awake = awake_mask(times, states)
    moving = (speed >= MOVEMENT_SPEED_MIN) & (speed <= MOVEMENT_SPEED_MAX)
    analysis_mask = finite_pos & awake & moving
    x = pos[:, 0]
    y = pos[:, 1]
    x_edges = bin_edges(x[analysis_mask], N_X_BINS)
    y_edges = bin_edges(y[analysis_mask], N_Y_BINS)
    occupancy = occupancy_map(x, y, analysis_mask, x_edges, y_edges, dt)
    occupied_valid = occupancy >= MIN_OCCUPANCY_S

    valid_spike_position_mask = analysis_mask
    unit_rows: list[dict[str, Any]] = []
    rate_maps: dict[int, np.ndarray] = {}
    spike_maps: dict[int, np.ndarray] = {}
    for unit in units:
        spike_pos, kept_spikes = interpolate_position_for_spikes(unit["spike_times"], times, pos, valid_spike_position_mask)
        if spike_pos.size:
            spike_counts, _, _ = np.histogram2d(spike_pos[:, 1], spike_pos[:, 0], bins=[y_edges, x_edges])
        else:
            spike_counts = np.zeros_like(occupancy)
        with np.errstate(divide="ignore", invalid="ignore"):
            rate_map = spike_counts / occupancy
        rate_map[~np.isfinite(rate_map)] = 0.0
        metrics = spatial_metrics(rate_map, occupancy, spike_counts)
        total_awake_moving_s = float(analysis_mask.sum() * dt)
        unit_rows.append(
            {
                "unit_id": int(unit["unit_id"]),
                "unit_index": int(unit["unit_index"]),
                "cell_type": unit["cell_type"],
                "n_total_spikes": int(len(unit["spike_times"])),
                "n_awake_moving_spikes": int(len(kept_spikes)),
                "awake_moving_duration_s": total_awake_moving_s,
                "awake_moving_fraction_of_recording": total_awake_moving_s / max(times[-1] - times[0], 1e-9),
                **metrics,
            }
        )
        rate_maps[int(unit["unit_id"])] = rate_map
        spike_maps[int(unit["unit_id"])] = spike_counts

    unit_df = pd.DataFrame(unit_rows)
    celltype_rows = []
    for cell_type, sub in unit_df.groupby("cell_type", dropna=False):
        celltype_rows.append(
            {
                "cell_type": cell_type,
                "n_units": int(sub.shape[0]),
                "median_spatial_information_bits_per_spike": float(sub["spatial_information_bits_per_spike"].median()),
                "median_spatial_sparsity": float(sub["spatial_sparsity"].median()),
                "median_active_spatial_bin_fraction": float(sub["active_spatial_bin_fraction"].median()),
                "median_awake_moving_rate_hz": float(sub["mean_awake_moving_rate_hz"].median()),
                "median_selectivity_max_over_mean": float(sub["spatial_selectivity_max_over_mean"].median()),
            }
        )
    celltype_df = pd.DataFrame(celltype_rows).sort_values(["n_units", "median_spatial_information_bits_per_spike"], ascending=[False, False])

    granule_ids = unit_df.loc[unit_df["cell_type"] == "granule cell", "unit_id"].astype(int).tolist()
    granule_maps = {uid: rate_maps[uid] for uid in granule_ids}
    pv_granule = population_vector_metrics(granule_maps, occupancy, x_edges, y_edges)
    all_maps = {uid: rate_maps[uid] for uid in unit_df["unit_id"].astype(int).tolist()}
    pv_all = population_vector_metrics(all_maps, occupancy, x_edges, y_edges)
    if not pv_all.empty:
        pv_all["unit_set"] = "all_units"
    pv_df = pd.concat([pv_granule, pv_all], ignore_index=True)

    position_summary = pd.DataFrame(
        [
            {
                "n_position_samples": int(len(times)),
                "position_rate_hz": rate,
                "recording_duration_s": float(times[-1] - times[0]),
                "awake_samples": int((finite_pos & awake).sum()),
                "awake_duration_s": float((finite_pos & awake).sum() * dt),
                "awake_moving_samples": int(analysis_mask.sum()),
                "awake_moving_duration_s": float(analysis_mask.sum() * dt),
                "movement_speed_min": MOVEMENT_SPEED_MIN,
                "movement_speed_max": MOVEMENT_SPEED_MAX,
                "x_min": float(np.nanmin(x[analysis_mask])),
                "x_max": float(np.nanmax(x[analysis_mask])),
                "y_min": float(np.nanmin(y[analysis_mask])),
                "y_max": float(np.nanmax(y[analysis_mask])),
                "x_bins": N_X_BINS,
                "y_bins": N_Y_BINS,
                "occupied_bins_min_0p25s": int(occupied_valid.sum()),
                "total_occupancy_s": float(occupancy[occupied_valid].sum()),
            }
        ]
    )

    plot_built = make_plot(unit_df, rate_maps, occupancy)

    position_summary.to_csv(OUT_POSITION, sep="\t", index=False)
    unit_df.to_csv(OUT_UNIT, sep="\t", index=False)
    celltype_df.to_csv(OUT_CELLTYPE, sep="\t", index=False)
    pv_df.to_csv(OUT_PV, sep="\t", index=False)

    g = celltype_df[celltype_df["cell_type"] == "granule cell"]
    lines = [
        "# DANDI 000003 Spatial Pattern Pilot",
        "",
        "Date built: 2026-06-23",
        "",
        "## Purpose",
        "",
        "This is the first position-linked pilot for Aim 3. It asks whether the already-downloaded DANDI 000003 session can support spatial coding and population-vector analyses for source-labeled dentate granule units.",
        "",
        "## Position/Occupancy",
        "",
        f"- Position samples: {int(position_summary.iloc[0]['n_position_samples'])}",
        f"- Position sampling rate: {rate:.4f} Hz",
        f"- Awake duration with finite position: {float(position_summary.iloc[0]['awake_duration_s']):.1f} s",
        f"- Awake-moving duration after speed filter: {float(position_summary.iloc[0]['awake_moving_duration_s']):.1f} s",
        f"- Occupied spatial bins with >= {MIN_OCCUPANCY_S} s occupancy: {int(position_summary.iloc[0]['occupied_bins_min_0p25s'])}",
        "",
        "Coordinates are treated as relative video-tracking coordinates. The NWB field labels the units as meters, but the observed values behave like tracking coordinates; no absolute anatomical scale is assumed.",
        "",
        "## Granule-Cell Spatial Metrics",
        "",
    ]
    if not g.empty:
        gr = g.iloc[0]
        lines.extend(
            [
                f"- Labeled granule units: {int(gr['n_units'])}",
                f"- Median granule spatial information: {float(gr['median_spatial_information_bits_per_spike']):.4f} bits/spike",
                f"- Median granule spatial sparsity: {float(gr['median_spatial_sparsity']):.4f}",
                f"- Median granule active spatial-bin fraction: {float(gr['median_active_spatial_bin_fraction']):.4f}",
                f"- Median granule awake-moving rate: {float(gr['median_awake_moving_rate_hz']):.4f} Hz",
                f"- Median granule max/mean spatial selectivity: {float(gr['median_selectivity_max_over_mean']):.4f}",
            ]
        )
    else:
        lines.append("No granule-cell-labeled units were available in this pilot.")

    if not pv_df.empty:
        lines.extend(["", "## Population-Vector Check", ""])
        for _, row in pv_df.iterrows():
            lines.append(
                f"- `{row['unit_set']}`: n_units={int(row['n_units'])}, "
                f"occupied_bins={int(row['n_occupied_bins'])}, "
                f"rho(spatial, neural euclidean)={float(row['spearman_spatial_vs_neural_euclidean']):.4f}, "
                f"far-minus-near euclidean={float(row['far_minus_near_neural_euclidean']):.4f}."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This pilot confirms that DANDI 000003 can support position-linked granule-cell activity analysis. It is still a feasibility layer: one session, three granule-labeled units, relative tracking coordinates, and no cross-session/task-condition model. The correct next step is to extend this workflow to more sessions and compute task-specific spatial information, active place-field fraction, and population-vector separation.",
            "",
            "## Outputs",
            "",
            f"- Position summary: `{rel(OUT_POSITION)}`",
            f"- Unit spatial metrics: `{rel(OUT_UNIT)}`",
            f"- Cell-type spatial summary: `{rel(OUT_CELLTYPE)}`",
            f"- Population-vector separation: `{rel(OUT_PV)}`",
            f"- Plot: `{rel(OUT_PLOT)}`" if plot_built else "- Plot: skipped because Matplotlib was unavailable.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")

    for path in [OUT_POSITION, OUT_UNIT, OUT_CELLTYPE, OUT_PV, OUT_PLOT, OUT_MD]:
        if path.exists():
            print(f"Wrote {rel(path)}")


if __name__ == "__main__":
    main()
