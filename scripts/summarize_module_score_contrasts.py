#!/usr/bin/env python3
"""Summarize per-cell module-score output into identity and structural contrasts."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
SUMMARY_IN = RESULTS / "per_cell_marker_module_score_summary.tsv"
CONTRAST_OUT = RESULTS / "module_score_identity_structural_contrasts.tsv"

STRUCTURAL_PANELS = ["shared_granule_neuronal", "morphogenesis_cytoskeleton", "axon_guidance_synapse"]


def context_note(dataset: str, sample: str, group: str, region: str) -> str:
    if dataset == "GSE104323" and group in {"GC-adult", "GC-juv", "Immature-GC", "Neuroblast"}:
        return "curated_dentate_granule_lineage"
    if dataset in {"GSE214309", "GSE214905", "GSE292261"}:
        return "dentate_granule_enriched"
    if region in {"cerebellum", "cerebellum_spatial"}:
        return "cerebellum_mixed_or_spatial"
    if region == "organoid":
        return "organoid_mixed"
    return "mixed_or_reference"


def dominant_identity(dentate: float, cerebellar: float) -> str:
    if dentate > cerebellar:
        return "dentate_identity_higher"
    if cerebellar > dentate:
        return "cerebellar_identity_higher"
    return "identity_tie"


def main() -> None:
    grouped: dict[tuple[str, str, str, str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    with SUMMARY_IN.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["dataset"], row["sample"], row["group"], row["species"], row["region"], row["platform"])
            grouped[key][row["panel"]] = row

    fields = [
        "dataset",
        "sample",
        "group",
        "species",
        "region",
        "platform",
        "n_cells_or_spots",
        "dentate_identity_log1p",
        "cerebellar_identity_log1p",
        "identity_contrast_dentate_minus_cerebellar",
        "dominant_identity_panel",
        "shared_granule_neuronal_log1p",
        "morphogenesis_cytoskeleton_log1p",
        "axon_guidance_synapse_log1p",
        "structural_program_mean_log1p",
        "context_note",
    ]

    rows = []
    for key, panels in sorted(grouped.items()):
        dataset, sample, group, species, region, platform = key
        dentate = float(panels.get("dentate_identity", {}).get("mean_log1p_expression_panel", 0) or 0)
        cerebellar = float(panels.get("cerebellar_identity", {}).get("mean_log1p_expression_panel", 0) or 0)
        structural_scores = [
            float(panels.get(panel, {}).get("mean_log1p_expression_panel", 0) or 0)
            for panel in STRUCTURAL_PANELS
        ]
        n_cells = next((panel_row["n_cells_or_spots"] for panel_row in panels.values()), "")
        rows.append(
            {
                "dataset": dataset,
                "sample": sample,
                "group": group,
                "species": species,
                "region": region,
                "platform": platform,
                "n_cells_or_spots": n_cells,
                "dentate_identity_log1p": f"{dentate:.6g}",
                "cerebellar_identity_log1p": f"{cerebellar:.6g}",
                "identity_contrast_dentate_minus_cerebellar": f"{dentate - cerebellar:.6g}",
                "dominant_identity_panel": dominant_identity(dentate, cerebellar),
                "shared_granule_neuronal_log1p": panels.get("shared_granule_neuronal", {}).get("mean_log1p_expression_panel", "0"),
                "morphogenesis_cytoskeleton_log1p": panels.get("morphogenesis_cytoskeleton", {}).get("mean_log1p_expression_panel", "0"),
                "axon_guidance_synapse_log1p": panels.get("axon_guidance_synapse", {}).get("mean_log1p_expression_panel", "0"),
                "structural_program_mean_log1p": f"{sum(structural_scores) / len(structural_scores):.6g}",
                "context_note": context_note(dataset, sample, group, region),
            }
        )

    with CONTRAST_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} contrast rows to {CONTRAST_OUT}")


if __name__ == "__main__":
    main()
