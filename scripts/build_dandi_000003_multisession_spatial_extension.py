#!/usr/bin/env python3
"""Multi-session extension for DANDI 000003 spatial coding analysis.

The one-session pilot proved that position-linked analysis is feasible. This
script generalizes that workflow across every DANDI 000003 NWB file currently
available under the legacy internal data folder and the external-drive raw-data
folder. It writes a ranked download plan for the next feasible sessions. It
still treats the result as validation support, not as a final standardized
behavioral pattern-separation assay.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

import build_dandi_000003_spatial_pattern_pilot as spatial


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
LEGACY_DATA = ROOT / "External_Data" / "DANDI" / "000003"
EXTERNAL_DATA = ROOT / "Project" / "external_data" / "raw" / "dandi" / "000003"
DATA_ROOTS = [LEGACY_DATA, EXTERNAL_DATA]
PREFERRED_DOWNLOAD_ROOT = EXTERNAL_DATA if EXTERNAL_DATA.parent.exists() else LEGACY_DATA

MANIFEST = RESULTS / "dandi_000003_asset_manifest.tsv"

OUT_DOWNLOAD_PLAN = RESULTS / "dandi_000003_multisession_download_plan.tsv"
OUT_SESSION = RESULTS / "dandi_000003_multisession_session_summary.tsv"
OUT_UNIT = RESULTS / "dandi_000003_multisession_spatial_unit_metrics.tsv"
OUT_CELLTYPE_SESSION = RESULTS / "dandi_000003_multisession_spatial_celltype_by_session.tsv"
OUT_CELLTYPE_POOLED = RESULTS / "dandi_000003_multisession_spatial_celltype_pooled.tsv"
OUT_PV = RESULTS / "dandi_000003_multisession_population_vector_separation.tsv"
OUT_PLOT = RESULTS / "dandi_000003_multisession_spatial_summary.png"
OUT_MD = RESULTS / "dandi_000003_multisession_spatial_extension.md"
OUT_TARGETED_PRIORITY = RESULTS / "dandi_000003_targeted_download_priority.tsv"

MAX_DOWNLOAD_CANDIDATES = 12

KNOWN_LOW_YIELD_ASSET_IDS = {
    # Previously downloaded and deleted after analysis because it had 0
    # source-labeled granule units under the conservative NWB cell_type rule.
    "7a344712-5124-434b-8a8e-fa36a9213294",
}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def resolve_local_path(relative_path: str) -> Path:
    for root in DATA_ROOTS:
        candidate = root / str(relative_path)
        if candidate.exists():
            return candidate
    return PREFERRED_DOWNLOAD_ROOT / str(relative_path)


def resolve_local_root(local_path: str) -> str:
    path = Path(local_path)
    for root in DATA_ROOTS:
        try:
            path.relative_to(root)
            return str(root)
        except ValueError:
            continue
    return ""


def load_manifest() -> pd.DataFrame:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest["local_path"] = manifest["path"].map(lambda p: str(resolve_local_path(str(p))))
    manifest["local_root"] = manifest["local_path"].map(resolve_local_root)
    manifest["local_exists"] = manifest["local_path"].map(lambda p: Path(p).exists())
    manifest = manifest.sort_values(["size_bytes", "path"]).reset_index(drop=True)
    return manifest


def build_download_plan(manifest: pd.DataFrame) -> pd.DataFrame:
    missing = manifest[
        (~manifest["local_exists"]) & (~manifest["asset_id"].isin(KNOWN_LOW_YIELD_ASSET_IDS))
    ].copy().sort_values(["size_bytes", "path"]).head(MAX_DOWNLOAD_CANDIDATES)
    if missing.empty:
        return missing
    missing["download_rank"] = range(1, len(missing) + 1)
    missing["download_url_api"] = missing["asset_id"].map(
        lambda asset_id: f"https://api.dandiarchive.org/api/assets/{asset_id}/download/"
    )
    missing["local_dir"] = missing["local_path"].map(lambda p: str(Path(p).parent))
    missing["download_command"] = missing.apply(
        lambda row: (
            f"mkdir -p '{row['local_dir']}' && "
            "curl -L --fail --continue-at - "
            f"-o '{row['local_path']}' '{row['download_url_api']}'"
        ),
        axis=1,
    )
    return missing[
        [
            "download_rank",
            "asset_id",
            "path",
            "subject",
            "session",
            "size_bytes",
            "size_gb",
            "local_path",
            "download_url_api",
            "download_command",
        ]
    ]


def session_metadata_from_path(path: str) -> tuple[str, str]:
    parts = Path(path).parts
    subject = parts[0].replace("sub-", "") if parts else ""
    session = ""
    name = Path(path).name
    marker = "_ses-"
    if marker in name:
        session = name.split(marker, 1)[1].split("_", 1)[0]
    return subject, session


def summarize_by_celltype(unit_df: pd.DataFrame, by_session: bool) -> pd.DataFrame:
    if unit_df.empty:
        return pd.DataFrame()
    group_cols = ["cell_type"]
    if by_session:
        group_cols = ["asset_id", "subject", "session", "cell_type"]
    rows: list[dict[str, Any]] = []
    for keys, sub in unit_df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_units": int(sub.shape[0]),
                "median_spatial_information_bits_per_spike": float(sub["spatial_information_bits_per_spike"].median()),
                "median_spatial_sparsity": float(sub["spatial_sparsity"].median()),
                "median_active_spatial_bin_fraction": float(sub["active_spatial_bin_fraction"].median()),
                "median_awake_moving_rate_hz": float(sub["mean_awake_moving_rate_hz"].median()),
                "median_selectivity_max_over_mean": float(sub["spatial_selectivity_max_over_mean"].median()),
            }
        )
        rows.append(row)
    out = pd.DataFrame(rows)
    sort_cols = ["n_units", "median_spatial_information_bits_per_spike"]
    return out.sort_values(sort_cols, ascending=[False, False]).reset_index(drop=True)


def analyze_one_asset(row: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    local_path = Path(str(row["local_path"]))
    subject = str(row.get("subject") or "")
    session = str(row.get("session") or "")
    if not subject or not session:
        subject, session = session_metadata_from_path(str(row["path"]))
    asset_id = str(row["asset_id"])

    with h5py.File(local_path, "r") as h5:
        times, pos, speed, rate = spatial.load_position(h5)
        states = spatial.load_states(h5)
        units = spatial.read_units(h5)

    dt = 1.0 / rate
    finite_pos = np.isfinite(pos).all(axis=1)
    awake = spatial.awake_mask(times, states)
    moving = (speed >= spatial.MOVEMENT_SPEED_MIN) & (speed <= spatial.MOVEMENT_SPEED_MAX)
    analysis_mask = finite_pos & awake & moving
    x = pos[:, 0]
    y = pos[:, 1]
    if analysis_mask.sum() < 100:
        raise ValueError("too few awake-moving position samples")

    x_edges = spatial.bin_edges(x[analysis_mask], spatial.N_X_BINS)
    y_edges = spatial.bin_edges(y[analysis_mask], spatial.N_Y_BINS)
    occupancy = spatial.occupancy_map(x, y, analysis_mask, x_edges, y_edges, dt)
    occupied_valid = occupancy >= spatial.MIN_OCCUPANCY_S

    unit_rows: list[dict[str, Any]] = []
    rate_maps: dict[int, np.ndarray] = {}
    for unit in units:
        spike_pos, kept_spikes = spatial.interpolate_position_for_spikes(
            unit["spike_times"], times, pos, analysis_mask
        )
        if spike_pos.size:
            spike_counts, _, _ = np.histogram2d(spike_pos[:, 1], spike_pos[:, 0], bins=[y_edges, x_edges])
        else:
            spike_counts = np.zeros_like(occupancy)
        with np.errstate(divide="ignore", invalid="ignore"):
            rate_map = spike_counts / occupancy
        rate_map[~np.isfinite(rate_map)] = 0.0
        metrics = spatial.spatial_metrics(rate_map, occupancy, spike_counts)
        total_awake_moving_s = float(analysis_mask.sum() * dt)
        unit_rows.append(
            {
                "asset_id": asset_id,
                "subject": subject,
                "session": session,
                "path": str(row["path"]),
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

    unit_df = pd.DataFrame(unit_rows)
    granule_ids = unit_df.loc[unit_df["cell_type"] == "granule cell", "unit_id"].astype(int).tolist()
    granule_maps = {uid: rate_maps[uid] for uid in granule_ids}
    pv_granule = spatial.population_vector_metrics(granule_maps, occupancy, x_edges, y_edges)
    if not pv_granule.empty:
        pv_granule["unit_set"] = "granule_cell_labeled"
    all_maps = {uid: rate_maps[uid] for uid in unit_df["unit_id"].astype(int).tolist()}
    pv_all = spatial.population_vector_metrics(all_maps, occupancy, x_edges, y_edges)
    if not pv_all.empty:
        pv_all["unit_set"] = "all_units"
    pv_df = pd.concat([pv_granule, pv_all], ignore_index=True)
    if not pv_df.empty:
        pv_df.insert(0, "session", session)
        pv_df.insert(0, "subject", subject)
        pv_df.insert(0, "asset_id", asset_id)
        pv_df.insert(3, "path", str(row["path"]))

    session_df = pd.DataFrame(
        [
            {
                "asset_id": asset_id,
                "subject": subject,
                "session": session,
                "path": str(row["path"]),
                "local_path": str(local_path),
                "size_gb": float(row["size_gb"]),
                "n_units": int(unit_df.shape[0]),
                "n_granule_units": int((unit_df["cell_type"] == "granule cell").sum()),
                "n_position_samples": int(len(times)),
                "position_rate_hz": float(rate),
                "recording_duration_s": float(times[-1] - times[0]),
                "awake_duration_s": float((finite_pos & awake).sum() * dt),
                "awake_moving_duration_s": float(analysis_mask.sum() * dt),
                "occupied_bins_min_0p25s": int(occupied_valid.sum()),
                "total_occupancy_s": float(occupancy[occupied_valid].sum()),
                "x_bins": spatial.N_X_BINS,
                "y_bins": spatial.N_Y_BINS,
                "movement_speed_min": spatial.MOVEMENT_SPEED_MIN,
                "movement_speed_max": spatial.MOVEMENT_SPEED_MAX,
            }
        ]
    )
    return session_df, unit_df, pv_df


def write_report(
    manifest: pd.DataFrame,
    download_plan: pd.DataFrame,
    session_df: pd.DataFrame,
    unit_df: pd.DataFrame,
    pooled_celltype: pd.DataFrame,
    pv_df: pd.DataFrame,
    failures: list[dict[str, str]],
    plot_built: bool,
) -> None:
    local_count = int(manifest["local_exists"].sum())
    local_gb = float(manifest.loc[manifest["local_exists"], "size_gb"].sum())
    granules = unit_df[unit_df["cell_type"] == "granule cell"] if not unit_df.empty else pd.DataFrame()
    g_summary = pooled_celltype[pooled_celltype["cell_type"] == "granule cell"] if not pooled_celltype.empty else pd.DataFrame()
    sessions_with_granules = 0
    zero_granule_sessions: list[str] = []
    if not session_df.empty:
        sessions_with_granules = int((session_df["n_granule_units"] > 0).sum())
        zero_granule_sessions = session_df.loc[session_df["n_granule_units"] == 0, "session"].astype(str).tolist()

    lines = [
        "# DANDI 000003 Multi-Session Spatial Extension",
        "",
        "Date built: 2026-06-24",
        "",
        "## Purpose",
        "",
        "This extends the one-session DANDI 000003 spatial pilot into a reusable multi-session workflow. It analyzes all locally available NWB files from the legacy internal data folder and the external-drive raw-data folder, then writes a ranked plan for downloading additional feasible sessions.",
        "",
        "## Current Local Coverage",
        "",
        f"- Local NWB files detected: {local_count}",
        f"- Local DANDI 000003 size analyzed or available: {local_gb:.2f} GB",
        f"- Sessions successfully analyzed: {session_df.shape[0]}",
        f"- Sessions with labeled granule units: {sessions_with_granules}",
        f"- Units analyzed: {unit_df.shape[0] if not unit_df.empty else 0}",
        f"- Labeled granule units analyzed: {granules.shape[0] if not granules.empty else 0}",
        "",
    ]
    if not g_summary.empty:
        g = g_summary.iloc[0]
        lines.extend(
            [
                "## Pooled Granule Metrics",
                "",
                f"- Granule units: {int(g['n_units'])}",
                f"- Median spatial information: {float(g['median_spatial_information_bits_per_spike']):.4f} bits/spike",
                f"- Median spatial sparsity: {float(g['median_spatial_sparsity']):.4f}",
                f"- Median active spatial-bin fraction: {float(g['median_active_spatial_bin_fraction']):.4f}",
                f"- Median awake-moving firing rate: {float(g['median_awake_moving_rate_hz']):.4f} Hz",
                "",
            ]
        )
    if not pv_df.empty:
        granule_pv = pv_df[pv_df["unit_set"] == "granule_cell_labeled"]
        if not granule_pv.empty:
            lines.extend(["## Granule Population-Vector Checks", ""])
            for _, row in granule_pv.iterrows():
                lines.append(
                    f"- `{row['session']}`: n_units={int(row['n_units'])}, "
                    f"occupied_bins={int(row['n_occupied_bins'])}, "
                    f"rho={float(row['spearman_spatial_vs_neural_euclidean']):.4f}, "
                    f"far-minus-near euclidean={float(row['far_minus_near_neural_euclidean']):.4f}."
                )
            lines.append("")
    if zero_granule_sessions:
        lines.extend(
            [
                "## Granule-Label Coverage Caveat",
                "",
                "The following analyzed sessions have usable units and position data but no source-labeled granule units in the NWB `cell_type` column:",
                "",
            ]
        )
        for session in zero_granule_sessions:
            lines.append(f"- `{session}`")
        lines.extend(
            [
                "",
                "These sessions contribute comparator/all-unit spatial context, but not direct granule-cell spatial evidence under the current conservative label rule.",
                "",
            ]
        )
    if not download_plan.empty:
        next_row = download_plan.iloc[0]
        lines.extend(
            [
                "## Next Size-Ranked Download Candidate",
                "",
                f"- Session: `{next_row['session']}`",
                f"- Subject: `{next_row['subject']}`",
                f"- Size: {float(next_row['size_gb']):.2f} GB",
                f"- Asset: `{next_row['asset_id']}`",
                "",
            ]
        )
    if OUT_TARGETED_PRIORITY.exists():
        try:
            targeted = pd.read_csv(OUT_TARGETED_PRIORITY, sep="\t")
        except Exception:
            targeted = pd.DataFrame()
        if not targeted.empty:
            top = targeted.iloc[0]
            lines.extend(
                [
                    "## Targeted Download Priority",
                    "",
                    "The size-ranked candidate is not necessarily the best biological next step. The targeted priority table favors expected granule-label yield per GB and subject-breadth tradeoffs.",
                    "",
                    f"- Targeted session: `{top['session']}`",
                    f"- Targeted subject: `{top['subject']}`",
                    f"- Targeted track: `{top['recommended_track']}`",
                    f"- Targeted size: {float(top['size_gb']):.2f} GB",
                    f"- Reason: {top['priority_reason']}",
                    f"- Priority table: `{rel(OUT_TARGETED_PRIORITY)}`",
                    "",
                ]
            )
    if failures:
        lines.extend(["## Analysis Failures", ""])
        for item in failures:
            lines.append(f"- `{item['path']}`: {item['error']}")
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "This file is an expansion scaffold. With the current local sessions it provides a first multi-session spatial result; as more NWB files are downloaded, the same outputs will scale directly. A manuscript-level pattern-separation claim still needs source-paper unit validation and task/trajectory-specific comparisons.",
            "",
            "## Outputs",
            "",
            f"- Download plan: `{rel(OUT_DOWNLOAD_PLAN)}`",
            f"- Session summary: `{rel(OUT_SESSION)}`",
            f"- Unit spatial metrics: `{rel(OUT_UNIT)}`",
            f"- Cell-type by session: `{rel(OUT_CELLTYPE_SESSION)}`",
            f"- Pooled cell-type summary: `{rel(OUT_CELLTYPE_POOLED)}`",
            f"- Population-vector separation: `{rel(OUT_PV)}`",
            f"- Plot: `{rel(OUT_PLOT)}`" if plot_built else "- Plot: skipped because Matplotlib was unavailable.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def make_summary_plot(session_df: pd.DataFrame, unit_df: pd.DataFrame, pv_df: pd.DataFrame) -> bool:
    if session_df.empty or unit_df.empty:
        return False
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    granules = unit_df[unit_df["cell_type"] == "granule cell"].copy()
    if granules.empty:
        return False
    sessions = session_df["session"].tolist()
    session_index = {session: idx for idx, session in enumerate(sessions)}
    granules["x"] = granules["session"].map(session_index).astype(float)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2))
    ax = axes[0, 0]
    ax.bar(range(len(sessions)), session_df["n_granule_units"], color="#5b8db8")
    ax.set_title("Granule-labeled units")
    ax.set_ylabel("n units")
    ax.set_xticks(range(len(sessions)), sessions, rotation=25, ha="right")

    ax = axes[0, 1]
    ax.scatter(
        granules["x"],
        granules["spatial_information_bits_per_spike"],
        s=38,
        color="#d1843f",
        edgecolor="white",
        linewidth=0.5,
    )
    med = granules.groupby("session")["spatial_information_bits_per_spike"].median()
    ax.plot([session_index[s] for s in med.index], med.values, color="#6f3f1d", linewidth=1.5)
    ax.set_title("Granule spatial information")
    ax.set_ylabel("bits/spike")
    ax.set_xticks(range(len(sessions)), sessions, rotation=25, ha="right")

    ax = axes[1, 0]
    ax.scatter(
        granules["x"],
        granules["active_spatial_bin_fraction"],
        s=38,
        color="#5f9e6e",
        edgecolor="white",
        linewidth=0.5,
    )
    med = granules.groupby("session")["active_spatial_bin_fraction"].median()
    ax.plot([session_index[s] for s in med.index], med.values, color="#2f6d43", linewidth=1.5)
    ax.set_title("Granule active spatial-bin fraction")
    ax.set_ylabel("fraction")
    ax.set_xticks(range(len(sessions)), sessions, rotation=25, ha="right")

    ax = axes[1, 1]
    granule_pv = pv_df[pv_df["unit_set"] == "granule_cell_labeled"].copy() if not pv_df.empty else pd.DataFrame()
    if not granule_pv.empty:
        x = [session_index[s] for s in granule_pv["session"]]
        ax.bar(x, granule_pv["far_minus_near_neural_euclidean"], color="#7c6bb0")
        ax2 = ax.twinx()
        ax2.plot(x, granule_pv["spearman_spatial_vs_neural_euclidean"], color="#2c2c2c", marker="o")
        ax2.set_ylabel("rho")
    ax.set_title("Population-vector separation")
    ax.set_ylabel("far-near Euclidean")
    ax.set_xticks(range(len(sessions)), sessions, rotation=25, ha="right")

    fig.suptitle("DANDI 000003 multi-session granule spatial pilot", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)
    return True


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest()
    download_plan = build_download_plan(manifest)

    session_tables: list[pd.DataFrame] = []
    unit_tables: list[pd.DataFrame] = []
    pv_tables: list[pd.DataFrame] = []
    failures: list[dict[str, str]] = []
    for _, row in manifest[manifest["local_exists"]].iterrows():
        try:
            session_df, unit_df, pv_df = analyze_one_asset(row)
        except Exception as exc:  # keep expanding even if one large file is odd
            failures.append({"path": str(row["path"]), "error": repr(exc)})
            continue
        session_tables.append(session_df)
        unit_tables.append(unit_df)
        if not pv_df.empty:
            pv_tables.append(pv_df)

    session_df = pd.concat(session_tables, ignore_index=True) if session_tables else pd.DataFrame()
    unit_df = pd.concat(unit_tables, ignore_index=True) if unit_tables else pd.DataFrame()
    pv_df = pd.concat(pv_tables, ignore_index=True) if pv_tables else pd.DataFrame()
    by_session = summarize_by_celltype(unit_df, by_session=True)
    pooled = summarize_by_celltype(unit_df, by_session=False)
    plot_built = make_summary_plot(session_df, unit_df, pv_df)

    download_plan.to_csv(OUT_DOWNLOAD_PLAN, sep="\t", index=False)
    session_df.to_csv(OUT_SESSION, sep="\t", index=False)
    unit_df.to_csv(OUT_UNIT, sep="\t", index=False)
    by_session.to_csv(OUT_CELLTYPE_SESSION, sep="\t", index=False)
    pooled.to_csv(OUT_CELLTYPE_POOLED, sep="\t", index=False)
    pv_df.to_csv(OUT_PV, sep="\t", index=False)
    write_report(manifest, download_plan, session_df, unit_df, pooled, pv_df, failures, plot_built)

    for path in [
        OUT_DOWNLOAD_PLAN,
        OUT_SESSION,
        OUT_UNIT,
        OUT_CELLTYPE_SESSION,
        OUT_CELLTYPE_POOLED,
        OUT_PV,
        OUT_PLOT,
        OUT_MD,
    ]:
        if path.exists():
            print(f"Wrote {rel(path)}")


if __name__ == "__main__":
    main()
