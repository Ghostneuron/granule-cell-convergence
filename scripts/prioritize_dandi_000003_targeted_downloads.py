#!/usr/bin/env python3
"""Prioritize targeted DANDI 000003 downloads under disk constraints.

The multi-session workflow ranks missing NWB files by size. That is useful for
bootstrap sampling, but after several sessions have been analyzed we can do
better: prefer sessions from animals that already yielded dentate granule units,
especially adjacent recording days, while retaining a separate note about new
subject breadth.
"""

from __future__ import annotations

import math
import re
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
LEGACY_DATA = ROOT / "External_Data" / "DANDI" / "000003"
EXTERNAL_DATA = ROOT / "Project" / "external_data" / "raw" / "dandi" / "000003"
DATA_ROOTS = [LEGACY_DATA, EXTERNAL_DATA]
PREFERRED_DOWNLOAD_ROOT = EXTERNAL_DATA if EXTERNAL_DATA.parent.exists() else LEGACY_DATA

MANIFEST = RESULTS / "dandi_000003_asset_manifest.tsv"
SESSION_SUMMARY = RESULTS / "dandi_000003_multisession_session_summary.tsv"

OUT_TSV = RESULTS / "dandi_000003_targeted_download_priority.tsv"
OUT_MD = RESULTS / "dandi_000003_targeted_download_priority.md"

KNOWN_LOW_YIELD_ASSET_IDS = {
    # Previously downloaded and deleted after analysis because it had 0
    # source-labeled granule units under the conservative NWB cell_type rule.
    "7a344712-5124-434b-8a8e-fa36a9213294",
}


def resolve_local_path(relative_path: str) -> Path:
    for root in DATA_ROOTS:
        candidate = root / str(relative_path)
        if candidate.exists():
            return candidate
    return PREFERRED_DOWNLOAD_ROOT / str(relative_path)


def parse_session_date(session: str) -> datetime | None:
    """Parse Yuta session strings such as YutaMouse55-160908 or YutaMouse51b160516."""
    match = re.search(r"(\d{6})$", str(session))
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%y%m%d")
    except ValueError:
        return None


def day_gap(session: str, local_sessions: list[str]) -> int | None:
    session_date = parse_session_date(session)
    local_dates = [parse_session_date(item) for item in local_sessions]
    local_dates = [item for item in local_dates if item is not None]
    if session_date is None or not local_dates:
        return None
    return min(abs((session_date - local_date).days) for local_date in local_dates)


def adjacency_score(days: int | None) -> float:
    if days is None:
        return 0.0
    if days <= 1:
        return 3.0
    if days <= 3:
        return 2.3
    if days <= 7:
        return 1.5
    if days <= 14:
        return 1.0
    if days <= 30:
        return 0.5
    return 0.0


def score_reason(row: pd.Series) -> str:
    reasons: list[str] = []
    if row["same_subject_with_granule_evidence"]:
        reasons.append(
            f"same subject has {int(row['local_subject_granule_units'])} local granule units"
        )
    elif row["local_subject_sessions"] > 0:
        reasons.append("same subject is already local but has weaker granule-label yield")
    else:
        reasons.append("new subject adds biological breadth")
    if pd.notna(row["nearest_local_day_gap"]):
        reasons.append(f"{int(row['nearest_local_day_gap'])}-day gap from local session")
    reasons.append(f"{row['size_gb']:.2f} GB")
    return "; ".join(reasons)


def build_priority() -> pd.DataFrame:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    sessions = pd.read_csv(SESSION_SUMMARY, sep="\t")

    manifest["local_path"] = manifest["path"].map(lambda p: str(resolve_local_path(str(p))))
    manifest["local_exists"] = manifest["local_path"].map(lambda p: Path(p).exists())
    local = sessions.copy()

    subject_sessions = local.groupby("subject")["session"].apply(list).to_dict()
    subject_summary = (
        local.groupby("subject")
        .agg(
            local_subject_sessions=("session", "count"),
            local_subject_units=("n_units", "sum"),
            local_subject_granule_units=("n_granule_units", "sum"),
            local_subject_max_granule_units=("n_granule_units", "max"),
            local_subject_mean_granule_units=("n_granule_units", "mean"),
        )
        .reset_index()
    )

    missing = manifest[
        (~manifest["local_exists"]) & (~manifest["asset_id"].isin(KNOWN_LOW_YIELD_ASSET_IDS))
    ].copy()
    missing = missing.merge(subject_summary, on="subject", how="left")
    for col in [
        "local_subject_sessions",
        "local_subject_units",
        "local_subject_granule_units",
        "local_subject_max_granule_units",
        "local_subject_mean_granule_units",
    ]:
        missing[col] = missing[col].fillna(0)

    missing["nearest_local_day_gap"] = missing.apply(
        lambda row: day_gap(row["session"], subject_sessions.get(row["subject"], [])),
        axis=1,
    )
    missing["new_subject"] = missing["local_subject_sessions"] == 0
    missing["same_subject_with_granule_evidence"] = missing["local_subject_granule_units"] > 0
    missing["same_subject_no_granule_evidence"] = (
        (missing["local_subject_sessions"] > 0) & (missing["local_subject_granule_units"] == 0)
    )

    min_size = float(missing["size_gb"].min())
    max_size = float(missing["size_gb"].max())
    size_span = max(max_size - min_size, 1e-9)
    missing["small_file_score"] = 2.0 * (max_size - missing["size_gb"]) / size_span
    missing["same_subject_granule_score"] = missing.apply(
        lambda row: (
            5.0
            + min(4.0, float(row["local_subject_mean_granule_units"]))
            + min(2.0, math.log1p(float(row["local_subject_granule_units"])))
            if row["same_subject_with_granule_evidence"]
            else 0.0
        ),
        axis=1,
    )
    missing["adjacent_session_score"] = missing["nearest_local_day_gap"].map(adjacency_score)
    missing["new_subject_breadth_score"] = missing["new_subject"].map(lambda x: 2.0 if x else 0.0)
    missing["same_subject_continuity_score"] = missing.apply(
        lambda row: 1.0 if row["same_subject_no_granule_evidence"] else 0.0,
        axis=1,
    )
    missing["large_file_penalty"] = missing["size_gb"].map(lambda x: max(0.0, (float(x) - 10.0) * 0.20))

    missing["priority_score"] = (
        missing["same_subject_granule_score"]
        + missing["adjacent_session_score"]
        + missing["new_subject_breadth_score"]
        + missing["same_subject_continuity_score"]
        + missing["small_file_score"]
        - missing["large_file_penalty"]
    )
    missing["priority_per_gb"] = missing["priority_score"] / missing["size_gb"]
    missing["priority_reason"] = missing.apply(score_reason, axis=1)

    missing = missing.sort_values(
        ["priority_per_gb", "priority_score", "size_gb"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    missing["targeted_rank"] = range(1, len(missing) + 1)

    missing["download_url_api"] = missing["asset_id"].map(
        lambda asset_id: f"https://api.dandiarchive.org/api/assets/{asset_id}/download/"
    )
    missing["recommended_track"] = missing.apply(
        lambda row: (
            "yield_first"
            if row["same_subject_with_granule_evidence"]
            else ("breadth_first" if row["new_subject"] else "continuity")
        ),
        axis=1,
    )
    missing["recommended_action"] = "hold"
    if not missing.empty:
        missing.loc[0, "recommended_action"] = "download_next_single_file"
    return missing[
        [
            "targeted_rank",
            "recommended_action",
            "recommended_track",
            "asset_id",
            "path",
            "subject",
            "session",
            "size_gb",
            "priority_score",
            "priority_per_gb",
            "local_subject_sessions",
            "local_subject_units",
            "local_subject_granule_units",
            "local_subject_max_granule_units",
            "local_subject_mean_granule_units",
            "nearest_local_day_gap",
            "new_subject",
            "same_subject_with_granule_evidence",
            "same_subject_granule_score",
            "adjacent_session_score",
            "new_subject_breadth_score",
            "small_file_score",
            "large_file_penalty",
            "priority_reason",
            "local_path",
            "download_url_api",
        ]
    ]


def write_report(priority: pd.DataFrame) -> None:
    top = priority.head(12).copy()
    next_row = priority.iloc[0] if not priority.empty else None

    lines = [
        "# DANDI 000003 Targeted Download Priority",
        "",
        "Date built: 2026-06-24",
        "",
        "## Rationale",
        "",
        "The default multi-session download plan ranks missing NWB files by size. This targeted plan instead asks which single additional file is most likely to improve the dentate-granule physiology validation layer per GB.",
        "",
        "Scoring favors: same-subject sessions from animals that already produced labeled granule units, adjacent recording days, small file size, and a smaller bonus for new-subject breadth. The score is heuristic and should guide downloads, not replace source-level unit validation.",
        "",
    ]
    if next_row is not None:
        lines.extend(
            [
                "## Recommended Next File",
                "",
                f"- Session: `{next_row['session']}`",
                f"- Subject: `{next_row['subject']}`",
                f"- Asset: `{next_row['asset_id']}`",
                f"- Size: {float(next_row['size_gb']):.2f} GB",
                f"- Track: `{next_row['recommended_track']}`",
                f"- Reason: {next_row['priority_reason']}",
                "",
                "This is preferred over the smallest missing file because it is an adjacent follow-up from a subject that already yielded labeled granule units locally.",
                "",
            ]
        )

    lines.extend(
        [
            "## Top Candidates",
            "",
            "| Rank | Session | Subject | Track | Size GB | Score/GB | Reason |",
            "|---:|---|---|---|---:|---:|---|",
        ]
    )
    for _, row in top.iterrows():
        lines.append(
            f"| {int(row['targeted_rank'])} | `{row['session']}` | `{row['subject']}` | "
            f"`{row['recommended_track']}` | {float(row['size_gb']):.2f} | "
            f"{float(row['priority_per_gb']):.3f} | {row['priority_reason']} |"
        )

    lines.extend(
        [
            "",
            "## Output",
            "",
            f"- Full priority table: `Project/results/{OUT_TSV.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    priority = build_priority()
    priority.to_csv(OUT_TSV, sep="\t", index=False)
    write_report(priority)
    print(f"Wrote Project/results/{OUT_TSV.name}")
    print(f"Wrote Project/results/{OUT_MD.name}")
    if not priority.empty:
        row = priority.iloc[0]
        print(
            "Next targeted file: "
            f"{row['session']} ({row['asset_id']}), {float(row['size_gb']):.2f} GB"
        )


if __name__ == "__main__":
    main()
