#!/usr/bin/env python3
"""Build a DANDI 000003 dentate activity/sparsity pilot layer.

The full DANDI 000003 archive is multi-terabyte scale, so this script first
builds a complete asset manifest and then analyzes only locally downloaded
pilot NWB files. It uses h5py rather than PyNWB so it can run in the available
workspace runtime.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CACHE = ROOT / "Project" / "dataset_search_cache" / "phys_morph"
DATA = ROOT / "External_Data" / "DANDI" / "000003"

ASSET_PAGES = [
    CACHE / "dandi_000003_assets_page1.json",
    CACHE / "dandi_000003_assets_page2.json",
]
MIN_DETAIL = CACHE / "dandi_000003_asset_min_detail.json"

OUT_MANIFEST = RESULTS / "dandi_000003_asset_manifest.tsv"
OUT_PILOT = RESULTS / "dandi_000003_pilot_asset_plan.tsv"
OUT_STRUCTURE = RESULTS / "dandi_000003_pilot_nwb_structure.tsv"
OUT_UNITS = RESULTS / "dandi_000003_pilot_unit_sparsity.tsv"
OUT_SESSION = RESULTS / "dandi_000003_pilot_session_sparsity_summary.tsv"
OUT_CELLTYPE = RESULTS / "dandi_000003_pilot_celltype_sparsity_summary.tsv"
OUT_STATES = RESULTS / "dandi_000003_pilot_behavior_states.tsv"
OUT_STATE_UNIT = RESULTS / "dandi_000003_pilot_state_unit_firing.tsv"
OUT_STATE_CELLTYPE = RESULTS / "dandi_000003_pilot_state_celltype_summary.tsv"
OUT_PLOT = RESULTS / "dandi_000003_pilot_activity_sparsity.png"
OUT_MD = RESULTS / "dandi_000003_activity_sparsity_pilot.md"

PILOT_ASSET_ID = "70aaf3be-9f95-40af-8abc-a8316cbe50ca"
PILOT_RELATIVE_PATH = "sub-YutaMouse41/sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb"
BIN_SECONDS = 1.0


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    return json.loads(path.read_text())


def asset_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for page in ASSET_PAGES:
        data = load_json(page)
        rows.extend(data.get("results") or [])
    return rows


def parse_subject_session(path: str) -> tuple[str, str]:
    subject = ""
    session = ""
    parts = path.split("/")
    if parts:
        subject = parts[0].replace("sub-", "")
    match = re.search(r"_ses-([^_]+)", path)
    if match:
        session = match.group(1)
    return subject, session


def format_gb(size: int | float) -> float:
    return float(size) / 1e9


def build_manifest() -> pd.DataFrame:
    rows = []
    for item in asset_rows():
        subject, session = parse_subject_session(item.get("path", ""))
        rows.append(
            {
                "asset_id": item.get("asset_id", ""),
                "path": item.get("path", ""),
                "subject": subject,
                "session": session,
                "size_bytes": int(item.get("size") or 0),
                "size_gb": format_gb(int(item.get("size") or 0)),
                "created": item.get("created", ""),
                "modified": item.get("modified", ""),
                "local_path": str(DATA / item.get("path", "")),
                "local_exists": (DATA / item.get("path", "")).exists(),
            }
        )
    manifest = pd.DataFrame(rows).sort_values(["size_bytes", "path"]).reset_index(drop=True)
    return manifest


def pilot_plan(manifest: pd.DataFrame) -> pd.DataFrame:
    detail = load_json(MIN_DETAIL)
    content_urls = detail.get("contentUrl") or []
    direct_url = ""
    api_url = ""
    for url in content_urls:
        if "/download/" in url:
            api_url = url
        elif "dandiarchive.s3.amazonaws.com" in url:
            direct_url = url
    pilot = manifest[manifest["asset_id"] == PILOT_ASSET_ID].copy()
    if pilot.empty:
        pilot = manifest.head(1).copy()
    pilot["pilot_rank"] = range(1, len(pilot) + 1)
    pilot["download_url_api"] = api_url
    pilot["download_url_s3"] = direct_url
    pilot["download_command"] = pilot.apply(
        lambda row: (
            "curl -L --fail --continue-at - "
            f"-o '{row['local_path']}' '{direct_url or api_url}'"
        ),
        axis=1,
    )
    return pilot[
        [
            "pilot_rank",
            "asset_id",
            "path",
            "subject",
            "session",
            "size_bytes",
            "size_gb",
            "local_path",
            "local_exists",
            "download_url_api",
            "download_url_s3",
            "download_command",
        ]
    ]


def h5py_available() -> bool:
    try:
        import h5py  # noqa: F401

        return True
    except ModuleNotFoundError:
        return False


def list_hdf5_structure(local_path: Path, max_rows: int = 300) -> pd.DataFrame:
    import h5py

    rows: list[dict[str, Any]] = []

    def visitor(name: str, obj: Any) -> None:
        if len(rows) >= max_rows:
            return
        kind = "group"
        shape = ""
        dtype = ""
        if hasattr(obj, "shape"):
            kind = "dataset"
            shape = str(obj.shape)
            dtype = str(obj.dtype)
        rows.append(
            {
                "path": "/" + name,
                "kind": kind,
                "shape": shape,
                "dtype": dtype,
            }
        )

    with h5py.File(local_path, "r") as h5:
        h5.visititems(visitor)
    return pd.DataFrame(rows)


def read_dataset(group: Any, key: str) -> np.ndarray | None:
    if key not in group:
        return None
    return np.asarray(group[key])


def decode_if_bytes(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def read_unit_column(units: Any, name: str, n_units: int) -> list[Any]:
    if name not in units:
        return [None] * n_units
    data = np.asarray(units[name])
    if data.ndim == 0:
        return [decode_if_bytes(data.item())] * n_units
    return [decode_if_bytes(v) for v in data[:n_units]]


def extract_units(local_path: Path, subject: str, session: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    import h5py

    with h5py.File(local_path, "r") as h5:
        if "units" not in h5:
            return pd.DataFrame(), pd.DataFrame()
        units = h5["units"]
        ids = read_dataset(units, "id")
        if ids is None:
            n_units = 0
            ids = np.array([])
        else:
            n_units = len(ids)
        spike_times = read_dataset(units, "spike_times")
        spike_index = read_dataset(units, "spike_times_index")
        if spike_times is None or spike_index is None or n_units == 0:
            return pd.DataFrame(), pd.DataFrame()

        electrode_group = read_unit_column(units, "electrode_group", n_units)
        region = read_unit_column(units, "region", n_units)
        cell_type = read_unit_column(units, "cell_type", n_units)
        quality = read_unit_column(units, "quality", n_units)

        global_start = float(np.nanmin(spike_times)) if len(spike_times) else 0.0
        global_stop = float(np.nanmax(spike_times)) if len(spike_times) else 0.0
        duration = max(global_stop - global_start, 1e-9)
        n_bins = max(1, int(math.ceil(duration / BIN_SECONDS)))
        population_active = np.zeros(n_bins, dtype=int)

        unit_rows = []
        previous = 0
        for i, stop in enumerate(spike_index):
            stop = int(stop)
            st = spike_times[previous:stop]
            previous = stop
            spike_count = int(st.size)
            firing_rate_hz = spike_count / duration
            if spike_count:
                bins = np.unique(np.floor((st - global_start) / BIN_SECONDS).astype(int))
                bins = bins[(bins >= 0) & (bins < n_bins)]
                active_bins = int(bins.size)
                population_active[bins] += 1
                isi = np.diff(st)
                mean_isi = float(np.mean(isi)) if isi.size else np.nan
                cv_isi = float(np.std(isi) / mean_isi) if isi.size and mean_isi > 0 else np.nan
            else:
                active_bins = 0
                mean_isi = np.nan
                cv_isi = np.nan
            unit_rows.append(
                {
                    "subject": subject,
                    "session": session,
                    "unit_index": i,
                    "unit_id": int(ids[i]) if i < len(ids) else i,
                    "electrode_group": electrode_group[i],
                    "region": region[i],
                    "cell_type": cell_type[i],
                    "quality": quality[i],
                    "spike_count": spike_count,
                    "recording_duration_s": duration,
                    "firing_rate_hz": firing_rate_hz,
                    "active_bin_fraction_1s": active_bins / n_bins,
                    "mean_isi_s": mean_isi,
                    "cv_isi": cv_isi,
                    "first_spike_s": float(st[0]) if spike_count else np.nan,
                    "last_spike_s": float(st[-1]) if spike_count else np.nan,
                }
            )

        unit_df = pd.DataFrame(unit_rows)
        session_df = pd.DataFrame(
            [
                {
                    "subject": subject,
                    "session": session,
                    "local_path": str(local_path),
                    "n_units": n_units,
                    "total_spikes": int(len(spike_times)),
                    "recording_start_s": global_start,
                    "recording_stop_s": global_stop,
                    "recording_duration_s": duration,
                    "bin_seconds": BIN_SECONDS,
                    "n_bins": n_bins,
                    "mean_population_active_units_per_bin": float(population_active.mean()),
                    "median_population_active_units_per_bin": float(np.median(population_active)),
                    "mean_population_active_fraction_per_bin": float(population_active.mean() / n_units) if n_units else np.nan,
                    "median_population_active_fraction_per_bin": float(np.median(population_active) / n_units) if n_units else np.nan,
                    "fraction_bins_with_any_unit_active": float(np.mean(population_active > 0)),
                    "median_unit_firing_rate_hz": float(unit_df["firing_rate_hz"].median()) if not unit_df.empty else np.nan,
                    "median_unit_active_bin_fraction_1s": float(unit_df["active_bin_fraction_1s"].median()) if not unit_df.empty else np.nan,
                }
            ]
        )
        return unit_df, session_df


def summarize_celltypes(unit_df: pd.DataFrame) -> pd.DataFrame:
    if unit_df.empty or "cell_type" not in unit_df:
        return pd.DataFrame()
    rows = []
    for cell_type, sub in unit_df.groupby("cell_type", dropna=False):
        rows.append(
            {
                "cell_type": cell_type,
                "n_units": int(sub.shape[0]),
                "total_spikes": int(sub["spike_count"].sum()),
                "median_firing_rate_hz": float(sub["firing_rate_hz"].median()),
                "mean_firing_rate_hz": float(sub["firing_rate_hz"].mean()),
                "median_active_bin_fraction_1s": float(sub["active_bin_fraction_1s"].median()),
                "mean_active_bin_fraction_1s": float(sub["active_bin_fraction_1s"].mean()),
                "min_firing_rate_hz": float(sub["firing_rate_hz"].min()),
                "max_firing_rate_hz": float(sub["firing_rate_hz"].max()),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_units", "total_spikes"], ascending=[False, False])


def extract_behavior_states(local_path: Path) -> pd.DataFrame:
    import h5py

    with h5py.File(local_path, "r") as h5:
        if "/processing/behavior/states" not in h5:
            return pd.DataFrame()
        states = h5["/processing/behavior/states"]
        labels = [decode_if_bytes(v) for v in np.asarray(states["label"])]
        starts = np.asarray(states["start_time"], dtype=float)
        stops = np.asarray(states["stop_time"], dtype=float)
        ids = np.asarray(states["id"], dtype=int) if "id" in states else np.arange(len(labels))
        rows = []
        for i, label, start, stop in zip(ids, labels, starts, stops):
            rows.append(
                {
                    "state_id": int(i),
                    "state_label": label,
                    "start_time_s": float(start),
                    "stop_time_s": float(stop),
                    "duration_s": float(max(stop - start, 0.0)),
                }
            )
        return pd.DataFrame(rows)


def extract_state_unit_firing(local_path: Path, subject: str, session: str, states_df: pd.DataFrame) -> pd.DataFrame:
    import h5py

    if states_df.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    with h5py.File(local_path, "r") as h5:
        if "units" not in h5:
            return pd.DataFrame()
        units = h5["units"]
        ids = read_dataset(units, "id")
        spike_times = read_dataset(units, "spike_times")
        spike_index = read_dataset(units, "spike_times_index")
        if ids is None or spike_times is None or spike_index is None:
            return pd.DataFrame()
        n_units = len(ids)
        cell_type = read_unit_column(units, "cell_type", n_units)
        previous = 0
        for unit_i, stop_i in enumerate(spike_index):
            stop_i = int(stop_i)
            st = spike_times[previous:stop_i]
            previous = stop_i
            for _, state in states_df.iterrows():
                start = float(state["start_time_s"])
                stop = float(state["stop_time_s"])
                duration = max(stop - start, 0.0)
                if duration <= 0:
                    count = 0
                    rate = np.nan
                else:
                    count = int(np.searchsorted(st, stop, side="left") - np.searchsorted(st, start, side="left"))
                    rate = count / duration
                rows.append(
                    {
                        "subject": subject,
                        "session": session,
                        "unit_index": unit_i,
                        "unit_id": int(ids[unit_i]),
                        "cell_type": cell_type[unit_i],
                        "state_id": int(state["state_id"]),
                        "state_label": state["state_label"],
                        "state_duration_s": duration,
                        "spike_count": count,
                        "firing_rate_hz": rate,
                        "active_in_state": count > 0,
                    }
                )
    return pd.DataFrame(rows)


def summarize_state_celltypes(state_unit: pd.DataFrame) -> pd.DataFrame:
    if state_unit.empty:
        return pd.DataFrame()
    rows = []
    for keys, sub in state_unit.groupby(["state_label", "cell_type"], dropna=False):
        state_label, cell_type = keys
        n_units = int(sub["unit_id"].nunique())
        state_durations = sub.drop_duplicates("state_id")["state_duration_s"]
        total_state_duration = float(state_durations.sum())
        rows.append(
            {
                "state_label": state_label,
                "cell_type": cell_type,
                "n_units": n_units,
                "n_state_intervals": int(sub["state_id"].nunique()),
                "total_state_duration_s": total_state_duration,
                "total_spikes": int(sub["spike_count"].sum()),
                "pooled_unit_rate_hz": float(sub["spike_count"].sum() / (total_state_duration * n_units)) if total_state_duration and n_units else np.nan,
                "median_unit_state_rate_hz": float(sub["firing_rate_hz"].median()),
                "fraction_unit_intervals_active": float(sub["active_in_state"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["state_label", "n_units", "total_spikes"], ascending=[True, False, False])


def make_plot(celltype_df: pd.DataFrame, state_celltype_df: pd.DataFrame) -> bool:
    if celltype_df.empty:
        return False
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    plot_df = celltype_df.copy()
    plot_df["cell_type_short"] = plot_df["cell_type"].astype(str).str.replace(
        "wide waveform cell \\(narrower, exclude opto tagged SST cell\\)",
        "wide waveform narrower",
        regex=True,
    )
    plot_df["cell_type_short"] = plot_df["cell_type_short"].str.replace("positive negative waveform unit", "pos-neg waveform", regex=False)
    plot_df["cell_type_short"] = plot_df["cell_type_short"].str.replace("positive waveform unit (bursty)", "positive bursty", regex=False)
    plot_df = plot_df.sort_values("median_firing_rate_hz", ascending=True)
    colors = ["#287271" if v == "granule cell" else "#777777" for v in plot_df["cell_type"]]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.4, 1.2, 1.0]})
    axes[0].barh(plot_df["cell_type_short"], plot_df["median_firing_rate_hz"], color=colors, alpha=0.85)
    axes[0].set_xlabel("Median firing rate (Hz)")
    axes[0].set_title("Pilot unit firing by cell type")
    axes[0].grid(axis="x", alpha=0.2)

    axes[1].barh(plot_df["cell_type_short"], plot_df["median_active_bin_fraction_1s"], color=colors, alpha=0.85)
    axes[1].set_xlabel("Median active 1 s-bin fraction")
    axes[1].set_xlim(0, 1)
    axes[1].set_title("Temporal sparsity")
    axes[1].grid(axis="x", alpha=0.2)

    if not state_celltype_df.empty:
        g = state_celltype_df[state_celltype_df["cell_type"] == "granule cell"].copy()
        order = ["awake", "nrem", "rem", "transit"]
        g["state_label"] = pd.Categorical(g["state_label"], categories=order, ordered=True)
        g = g.sort_values("state_label")
        axes[2].bar(g["state_label"].astype(str), g["pooled_unit_rate_hz"], color="#b85634", alpha=0.85)
        axes[2].set_ylabel("Pooled granule rate (Hz)")
        axes[2].set_title("Granule state rates")
        axes[2].grid(axis="y", alpha=0.2)
    else:
        axes[2].axis("off")

    fig.suptitle("DANDI 000003 pilot activity/sparsity", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)
    return True


def write_empty_outputs() -> None:
    pd.DataFrame().to_csv(OUT_STRUCTURE, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_UNITS, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_SESSION, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_CELLTYPE, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_STATES, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_STATE_UNIT, sep="\t", index=False)
    pd.DataFrame().to_csv(OUT_STATE_CELLTYPE, sep="\t", index=False)


def write_report(
    manifest: pd.DataFrame,
    pilot: pd.DataFrame,
    unit_df: pd.DataFrame,
    session_df: pd.DataFrame,
    structure_df: pd.DataFrame,
    celltype_df: pd.DataFrame | None = None,
    states_df: pd.DataFrame | None = None,
    state_celltype_df: pd.DataFrame | None = None,
    plot_built: bool = False,
) -> None:
    total_gb = manifest["size_gb"].sum()
    smallest = manifest.iloc[0]
    local_exists = bool(pilot.iloc[0]["local_exists"]) if not pilot.empty else False
    lines = [
        "# DANDI 000003 Activity/Sparsity Pilot",
        "",
        "Date built: 2026-06-23",
        "",
        "## Purpose",
        "",
        "This layer prepares direct dentate activity validation for Aim 3 using DANDI 000003. The full archive is too large for blind download, so this step builds a complete asset manifest, chooses the smallest pilot NWB file, and analyzes local NWB files only when present.",
        "",
        "## Asset Inventory",
        "",
        f"- Assets listed: {len(manifest)}",
        f"- Total archive size from manifest: {total_gb:.2f} GB",
        f"- Smallest asset: `{smallest['path']}` ({smallest['size_gb']:.2f} GB)",
        f"- Pilot local file present: {local_exists}",
        "",
        "## Pilot Asset",
        "",
        "| Asset | Subject | Session | Size GB | Local file |",
        "|---|---|---|---:|---|",
    ]
    for _, row in pilot.iterrows():
        lines.append(
            f"| `{row['path']}` | `{row['subject']}` | `{row['session']}` | "
            f"{row['size_gb']:.2f} | `{row['local_path']}` |"
        )

    if not local_exists:
        lines.extend(
            [
                "",
                "## Status",
                "",
                "The pilot NWB file has not yet been downloaded, so this run produced a manifest and download/extraction plan but no unit-level firing sparsity values. The exact resumable download command is stored in `Project/results/dandi_000003_pilot_asset_plan.tsv`.",
            ]
        )
    elif unit_df.empty:
        lines.extend(
            [
                "",
                "## Status",
                "",
                "A local pilot file was found, but the script did not extract a standard `/units/spike_times` table. Inspect `Project/results/dandi_000003_pilot_nwb_structure.tsv` to adjust the HDF5 path logic.",
            ]
        )
    else:
        s = session_df.iloc[0]
        granule = pd.DataFrame()
        if celltype_df is not None and not celltype_df.empty:
            granule = celltype_df[celltype_df["cell_type"] == "granule cell"]
        lines.extend(
            [
                "",
                "## Pilot Activity Result",
                "",
                f"- Units: {int(s['n_units'])}",
                f"- Total spikes: {int(s['total_spikes'])}",
                f"- Recording duration: {float(s['recording_duration_s']):.2f} s",
                f"- Median unit firing rate: {float(s['median_unit_firing_rate_hz']):.4f} Hz",
                f"- Median unit active-bin fraction at 1 s bins: {float(s['median_unit_active_bin_fraction_1s']):.4f}",
                f"- Mean population active fraction per 1 s bin: {float(s['mean_population_active_fraction_per_bin']):.4f}",
            ]
        )
        if not granule.empty:
            g = granule.iloc[0]
            lines.extend(
                [
                    "",
                    "Granule-cell-labeled pilot subset:",
                    f"- Labeled granule units: {int(g['n_units'])}",
                    f"- Median granule-cell firing rate: {float(g['median_firing_rate_hz']):.4f} Hz",
                    f"- Median granule-cell active-bin fraction at 1 s bins: {float(g['median_active_bin_fraction_1s']):.4f}",
                ]
            )
        if states_df is not None and not states_df.empty:
            state_summary = (
                states_df.groupby("state_label")
                .agg(n_intervals=("state_id", "count"), total_duration_s=("duration_s", "sum"))
                .reset_index()
            )
            lines.extend(["", "Behavior-state intervals in this pilot:"])
            for _, row in state_summary.iterrows():
                lines.append(
                    f"- `{row['state_label']}`: {int(row['n_intervals'])} intervals, {float(row['total_duration_s']):.1f} s"
                )
        if state_celltype_df is not None and not state_celltype_df.empty:
            granule_state = state_celltype_df[state_celltype_df["cell_type"] == "granule cell"]
            if not granule_state.empty:
                lines.extend(["", "Granule-cell state-dependent pooled rates:"])
                for _, row in granule_state.iterrows():
                    lines.append(
                        f"- `{row['state_label']}`: {float(row['pooled_unit_rate_hz']):.4f} Hz pooled unit rate; "
                        f"{float(row['fraction_unit_intervals_active']):.3f} active unit-interval fraction"
                    )
        lines.extend(
            [
                "",
                "These are pilot values only. The source labels identify three units as granule cells, but broader conclusions require more sessions and source-paper validation of unit identity conventions.",
            ]
        )

    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- DANDI 000003 is excellent for direct activity validation, but the archive is 2.56 TB and the smallest NWB is 4.66 GB.",
            "- Pilot extraction uses generic NWB `/units` HDF5 fields and does not yet classify units as granule cells versus mossy cells or other hippocampal units.",
            "- Firing sparsity here is an electrophysiological activity measure, not a transcriptomic property.",
            "- Pattern-separation metrics require position/task epochs and population-vector analysis after unit identity and behavior timestamps are verified.",
            "",
            "## Outputs",
            "",
            f"- Asset manifest: `{rel(OUT_MANIFEST)}`",
            f"- Pilot download plan: `{rel(OUT_PILOT)}`",
            f"- NWB structure table: `{rel(OUT_STRUCTURE)}`",
            f"- Unit sparsity table: `{rel(OUT_UNITS)}`",
            f"- Session sparsity summary: `{rel(OUT_SESSION)}`",
            f"- Cell-type sparsity summary: `{rel(OUT_CELLTYPE)}`",
            f"- Behavior states: `{rel(OUT_STATES)}`",
            f"- State-unit firing table: `{rel(OUT_STATE_UNIT)}`",
            f"- State-cell-type summary: `{rel(OUT_STATE_CELLTYPE)}`",
            f"- Plot: `{rel(OUT_PLOT)}`" if plot_built else "- Plot: skipped because no plotted cell-type summary was available.",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest()
    pilot = pilot_plan(manifest)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)
    pilot.to_csv(OUT_PILOT, sep="\t", index=False)

    if not h5py_available():
        write_empty_outputs()
        write_report(manifest, pilot, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    else:
        local_path = Path(pilot.iloc[0]["local_path"])
        if not local_path.exists():
            write_empty_outputs()
            write_report(manifest, pilot, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
        else:
            structure = list_hdf5_structure(local_path)
            structure.to_csv(OUT_STRUCTURE, sep="\t", index=False)
            subject = str(pilot.iloc[0]["subject"])
            session = str(pilot.iloc[0]["session"])
            unit_df, session_df = extract_units(local_path, subject, session)
            celltype_df = summarize_celltypes(unit_df)
            states_df = extract_behavior_states(local_path)
            state_unit_df = extract_state_unit_firing(local_path, subject, session, states_df)
            state_celltype_df = summarize_state_celltypes(state_unit_df)
            plot_built = make_plot(celltype_df, state_celltype_df)
            unit_df.to_csv(OUT_UNITS, sep="\t", index=False)
            session_df.to_csv(OUT_SESSION, sep="\t", index=False)
            celltype_df.to_csv(OUT_CELLTYPE, sep="\t", index=False)
            states_df.to_csv(OUT_STATES, sep="\t", index=False)
            state_unit_df.to_csv(OUT_STATE_UNIT, sep="\t", index=False)
            state_celltype_df.to_csv(OUT_STATE_CELLTYPE, sep="\t", index=False)
            write_report(manifest, pilot, unit_df, session_df, structure, celltype_df, states_df, state_celltype_df, plot_built)

    for path in [
        OUT_MANIFEST,
        OUT_PILOT,
        OUT_STRUCTURE,
        OUT_UNITS,
        OUT_SESSION,
        OUT_CELLTYPE,
        OUT_STATES,
        OUT_STATE_UNIT,
        OUT_STATE_CELLTYPE,
        OUT_PLOT,
        OUT_MD,
    ]:
        if path.exists():
            print(f"Wrote {rel(path)}")


if __name__ == "__main__":
    main()
