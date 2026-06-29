#!/usr/bin/env python3
"""Focused sender-receiver ligand-receptor prediction for Aim 2.

This is the direct niche-cell upgrade to the earlier pathway-readiness audit.
It scores ligand expression in explicit niche sender classes and receptor
expression in granule-lineage receiver classes for:

- GSE122357 mouse cerebellum: Purkinje, astrocyte/Bergmann-proxy, microglia,
  endothelial, and oligodendroglial senders toward granule precursors/cells.
- GSE104323 mouse dentate gyrus/SGZ: astrocyte, endothelial, PVM/macrophage
  proxy, vascular/support, and oligodendroglial senders toward the dentate
  granule lineage.

The output is expression-based LR prediction, not spatial adjacency,
secreted-protein abundance, or perturbation evidence.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
import tarfile
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project/results"

MGI_MAP = RESULTS / "primary_core_mgi_ortholog_meta_model_map.tsv"
AIM2_LR_PAIRS = RESULTS / "primary_core_aim2_ligand_receptor_pairs.tsv"

GSE122357_TAR = EXTERNAL / "GEO/GSE122357/GSE122357_RAW.tar"
GSE122357_LABELS = EXTERNAL / "GEO/GSE122357/GSE122357_cell_number.xlsx"
GSE104323_EXPR = EXTERNAL / "GEO/GSE104323/GSE104323_10X_expression_data_V2.tab.gz"
GSE104323_META = EXTERNAL / "GEO/GSE104323/GSE104323_metadata_barcodes_24185cells.txt.gz"

OUT_PAIRS = RESULTS / "aim2_sender_receiver_lr_pairs.tsv"
OUT_GROUP_EXPR = RESULTS / "aim2_sender_receiver_lr_group_expression.tsv.gz"
OUT_PREDICTIONS = RESULTS / "aim2_sender_receiver_lr_predictions.tsv.gz"
OUT_SUMMARY = RESULTS / "aim2_sender_receiver_lr_summary.tsv"
OUT_TOP = RESULTS / "aim2_sender_receiver_lr_top_predictions.tsv"
OUT_COVERAGE = RESULTS / "aim2_sender_receiver_lr_coverage.tsv"
OUT_PLOT = RESULTS / "aim2_sender_receiver_lr.png"
OUT_MD = RESULTS / "aim2_sender_receiver_lr.md"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CEREBELLAR_SAMPLES = [
    {
        "dataset": "GSE122357",
        "region": "cerebellum",
        "sample": "GSM3464549_P0",
        "stage": "P0",
        "prefix": "p0",
        "tar_member": "GSM3464549_P0.csv.gz",
    },
    {
        "dataset": "GSE122357",
        "region": "cerebellum",
        "sample": "GSM3464550_P8a",
        "stage": "P8a",
        "prefix": "p8a",
        "tar_member": "GSM3464550_P8a.csv.gz",
    },
    {
        "dataset": "GSE122357",
        "region": "cerebellum",
        "sample": "GSM3464551_P8b",
        "stage": "P8b",
        "prefix": "p8b",
        "tar_member": "GSM3464551_P8b.csv.gz",
    },
]

DG_SAMPLE = {
    "dataset": "GSE104323",
    "region": "dentate_sgz",
    "sample": "10X_all_cells",
    "stage": "P120_adult",
    "path": GSE104323_EXPR,
}

CEREBELLAR_CELL_META = {
    "Granule precursor": {
        "niche_role": "cerebellar_granule_precursor_receiver",
        "display": "GC precursor",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "Granule cells": {
        "niche_role": "cerebellar_granule_cell_receiver",
        "display": "GC",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "Purkinje cells": {
        "niche_role": "cerebellar_purkinje_sender",
        "display": "Purkinje",
        "cell_axis": "sender",
        "sender_class": "purkinje",
        "receiver_class": "",
        "core_focus": True,
    },
    "Astrocytes": {
        "niche_role": "cerebellar_astrocyte_bergmann_proxy_sender",
        "display": "Astro/Bergmann proxy",
        "cell_axis": "sender",
        "sender_class": "astroglial_bergmann_proxy",
        "receiver_class": "",
        "core_focus": True,
    },
    "Microglia": {
        "niche_role": "cerebellar_microglia_sender",
        "display": "Microglia",
        "cell_axis": "sender",
        "sender_class": "microglia",
        "receiver_class": "",
        "core_focus": True,
    },
    "Endothelial cells": {
        "niche_role": "cerebellar_endothelial_sender",
        "display": "Endothelial",
        "cell_axis": "sender",
        "sender_class": "endothelial",
        "receiver_class": "",
        "core_focus": True,
    },
    "Oligodendrocytes": {
        "niche_role": "cerebellar_oligodendroglial_sender",
        "display": "Oligodendroglial",
        "cell_axis": "sender",
        "sender_class": "oligodendroglial_support",
        "receiver_class": "",
        "core_focus": False,
    },
    "Interneuron": {
        "niche_role": "cerebellar_interneuron_context",
        "display": "Interneuron context",
        "cell_axis": "context",
        "sender_class": "",
        "receiver_class": "",
        "core_focus": False,
    },
}

DG_CLUSTER_META = {
    "RGL_young": {
        "niche_role": "sgz_rgl_young_receiver",
        "display": "RGL young",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "RGL": {
        "niche_role": "sgz_rgl_receiver",
        "display": "RGL",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "nIPC": {
        "niche_role": "sgz_nipc_receiver",
        "display": "nIPC",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "nIPC-perin": {
        "niche_role": "sgz_nipc_perin_receiver",
        "display": "nIPC-perin",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "Neuroblast": {
        "niche_role": "sgz_neuroblast_receiver",
        "display": "Neuroblast",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "Immature-GC": {
        "niche_role": "sgz_immature_gc_receiver",
        "display": "Immature GC",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "GC-juv": {
        "niche_role": "sgz_juvenile_gc_receiver",
        "display": "Juvenile GC",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "GC-adult": {
        "niche_role": "sgz_adult_gc_receiver",
        "display": "Adult GC",
        "cell_axis": "receiver",
        "sender_class": "",
        "receiver_class": "granule_lineage",
        "core_focus": True,
    },
    "Astro-adult": {
        "niche_role": "sgz_astrocyte_sender",
        "display": "Astro adult",
        "cell_axis": "sender",
        "sender_class": "astrocyte",
        "receiver_class": "",
        "core_focus": True,
    },
    "Astro-juv": {
        "niche_role": "sgz_astrocyte_sender",
        "display": "Astro juvenile",
        "cell_axis": "sender",
        "sender_class": "astrocyte",
        "receiver_class": "",
        "core_focus": True,
    },
    "Immature-Astro": {
        "niche_role": "sgz_astrocyte_sender",
        "display": "Immature astro",
        "cell_axis": "sender",
        "sender_class": "astrocyte",
        "receiver_class": "",
        "core_focus": True,
    },
    "Endothelial": {
        "niche_role": "sgz_endothelial_sender",
        "display": "Endothelial",
        "cell_axis": "sender",
        "sender_class": "endothelial",
        "receiver_class": "",
        "core_focus": True,
    },
    "PVM": {
        "niche_role": "sgz_pvm_microglia_proxy_sender",
        "display": "PVM/microglia proxy",
        "cell_axis": "sender",
        "sender_class": "microglia_macrophage_proxy",
        "receiver_class": "",
        "core_focus": True,
    },
    "VLMC": {
        "niche_role": "sgz_vascular_meningeal_sender",
        "display": "VLMC vascular support",
        "cell_axis": "sender",
        "sender_class": "vascular_support",
        "receiver_class": "",
        "core_focus": False,
    },
    "OPC": {
        "niche_role": "sgz_opc_oligodendroglial_sender",
        "display": "OPC",
        "cell_axis": "sender",
        "sender_class": "oligodendroglial_support",
        "receiver_class": "",
        "core_focus": False,
    },
    "NFOL": {
        "niche_role": "sgz_opc_oligodendroglial_sender",
        "display": "NFOL",
        "cell_axis": "sender",
        "sender_class": "oligodendroglial_support",
        "receiver_class": "",
        "core_focus": False,
    },
    "MOL": {
        "niche_role": "sgz_opc_oligodendroglial_sender",
        "display": "MOL",
        "cell_axis": "sender",
        "sender_class": "oligodendroglial_support",
        "receiver_class": "",
        "core_focus": False,
    },
}

AUGMENTED_PAIRS = [
    ("CXCL12_CXCR4", "cxcl12_cxcr4", "CXCL12/CXCR4", "neurogenic_permissive", "CXCL12 to CXCR4", "CXCL12", "CXCR4"),
    ("CXCL12_ACKR3", "cxcl12_cxcr4", "CXCL12/CXCR4", "neurogenic_permissive", "CXCL12 to ACKR3", "CXCL12", "ACKR3"),
    ("VEGFA_KDR", "vegf_vascular", "VEGF/vascular", "vascular_niche", "VEGFA to KDR", "VEGFA", "KDR"),
    ("VEGFA_FLT1", "vegf_vascular", "VEGF/vascular", "vascular_niche", "VEGFA to FLT1", "VEGFA", "FLT1"),
    ("ANGPT1_TEK", "angiopoietin_vascular", "Angiopoietin/TEK", "vascular_niche", "ANGPT1 to TEK", "ANGPT1", "TEK"),
    ("ANGPT2_TEK", "angiopoietin_vascular", "Angiopoietin/TEK", "vascular_niche", "ANGPT2 to TEK", "ANGPT2", "TEK"),
    ("APOE_LRP1", "apoe_lrp", "APOE/LRP", "glial_lipid_trophic", "APOE to LRP1", "APOE", "LRP1"),
    ("APOE_LDLR", "apoe_lrp", "APOE/LDLR", "glial_lipid_trophic", "APOE to LDLR", "APOE", "LDLR"),
    ("SPP1_ITGAV", "spp1_integrin", "SPP1/integrin", "immune_matrix", "SPP1 to ITGAV", "SPP1", "ITGAV"),
    ("SPP1_CD44", "spp1_cd44", "SPP1/CD44", "immune_matrix", "SPP1 to CD44", "SPP1", "CD44"),
    ("CSF1_CSF1R", "csf1_microglia", "CSF1/CSF1R", "microglia_macrophage", "CSF1 to CSF1R", "CSF1", "CSF1R"),
    ("C1QA_LRP1", "c1q_lrp", "C1Q/LRP", "microglia_macrophage", "C1QA to LRP1", "C1QA", "LRP1"),
    ("CCL2_CCR2", "ccl2_ccr2", "CCL2/CCR2", "immune_chemokine", "CCL2 to CCR2", "CCL2", "CCR2"),
    ("IGF1_IGF1R", "igf_trophic", "IGF1/IGF1R", "trophic_maturation", "IGF1 to IGF1R", "IGF1", "IGF1R"),
    ("CNTF_CNTFR", "cntf_lifr", "CNTF/CNTFR", "astroglial_trophic", "CNTF to CNTFR", "CNTF", "CNTFR"),
    ("LIF_LIFR", "lif_lifr", "LIF/LIFR", "astroglial_trophic", "LIF to LIFR", "LIF", "LIFR"),
    ("GDNF_GFRA1", "gdnf_ret", "GDNF/GFRA1", "trophic_maturation", "GDNF to GFRA1", "GDNF", "GFRA1"),
    ("NRG1_ERBB4", "nrg_erbb", "NRG1/ERBB4", "synaptic_glial", "NRG1 to ERBB4", "NRG1", "ERBB4"),
    ("KITLG_KIT", "kitlg_kit", "KITLG/KIT", "trophic_maturation", "KITLG to KIT", "KITLG", "KIT"),
    ("HBEGF_EGFR", "egf_mapk", "HBEGF/EGFR", "trophic_maturation", "HBEGF to EGFR", "HBEGF", "EGFR"),
    ("EFNB2_EPHB2", "ephrin_guidance", "Ephrin/EPH", "axon_guidance_matrix", "EFNB2 to EPHB2", "EFNB2", "EPHB2"),
    ("SLIT2_ROBO1", "slit_robo", "SLIT/ROBO", "axon_guidance_matrix", "SLIT2 to ROBO1", "SLIT2", "ROBO1"),
    ("SLIT2_ROBO2", "slit_robo", "SLIT/ROBO", "axon_guidance_matrix", "SLIT2 to ROBO2", "SLIT2", "ROBO2"),
    ("THBS1_CD47", "thbs_cd47", "THBS/CD47", "synaptic_matrix", "THBS1 to CD47", "THBS1", "CD47"),
]


def canon(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').strip("'").upper()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_mouse_to_human() -> dict[str, str]:
    table = pd.read_csv(MGI_MAP, sep="\t")
    table = table.loc[table["mgi_one_to_one_human_mouse"].astype(str).str.lower().eq("true")].copy()
    table["mouse_key"] = table["mouse_symbol"].map(canon)
    table["human_key"] = table["human_symbol"].map(canon)
    return dict(zip(table["mouse_key"], table["human_key"]))


def load_pairs() -> pd.DataFrame:
    base = pd.read_csv(AIM2_LR_PAIRS, sep="\t")
    base["pair_source"] = "aim2_curated"
    aug_rows = []
    for idx, (pair_key, pathway_id, pathway_label, family, label, ligand, receptor) in enumerate(AUGMENTED_PAIRS, start=1):
        aug_rows.append(
            {
                "pair_id": f"SRLR{idx:03d}",
                "pathway_id": pathway_id,
                "pathway_label": pathway_label,
                "pathway_family": family,
                "pair_label": label,
                "ligand_gene": ligand,
                "ligand_canonical_gene": canon(ligand),
                "receptor_gene": receptor,
                "receptor_canonical_gene": canon(receptor),
                "pair_source": "niche_sender_receiver_augmented",
            }
        )
    pairs = pd.concat([base, pd.DataFrame(aug_rows)], ignore_index=True)
    pairs["ligand_canonical_gene"] = pairs["ligand_canonical_gene"].map(canon)
    pairs["receptor_canonical_gene"] = pairs["receptor_canonical_gene"].map(canon)
    pairs = pairs.drop_duplicates(["ligand_canonical_gene", "receptor_canonical_gene"]).reset_index(drop=True)
    pairs["sender_receiver_pair_id"] = [f"SRP{i:03d}" for i in range(1, len(pairs) + 1)]
    ordered = [
        "sender_receiver_pair_id",
        "pair_id",
        "pair_source",
        "pathway_id",
        "pathway_label",
        "pathway_family",
        "pair_label",
        "ligand_gene",
        "ligand_canonical_gene",
        "receptor_gene",
        "receptor_canonical_gene",
    ]
    pairs = pairs[ordered]
    pairs.to_csv(OUT_PAIRS, sep="\t", index=False)
    return pairs


def load_cerebellar_label_maps() -> dict[str, dict[str, str]]:
    labels = pd.read_excel(GSE122357_LABELS, sheet_name="Sheet1")
    maps: dict[str, dict[str, str]] = {}
    for sample in CEREBELLAR_SAMPLES:
        prefix = sample["prefix"] + "_"
        label_map: dict[str, str] = {}
        for raw_group in labels.columns:
            if raw_group not in CEREBELLAR_CELL_META:
                continue
            for value in labels[raw_group].dropna().astype(str):
                if not value.startswith(prefix):
                    continue
                label_map[value[len(prefix) :]] = raw_group
        maps[sample["sample"]] = label_map
    return maps


def load_dg_label_map() -> dict[str, str]:
    meta = pd.read_csv(GSE104323_META, sep="\t")
    return dict(
        zip(
            meta["Sample name (24185 single cells)"].astype(str),
            meta["characteristics: cell cluster"].astype(str),
        )
    )


def tar_member_text(tf: tarfile.TarFile, member_name: str):
    raw = tf.extractfile(member_name)
    if raw is None:
        raise FileNotFoundError(member_name)
    if member_name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), newline="")
    return io.TextIOWrapper(raw, newline="")


def group_metadata(region: str, raw_group: str) -> dict[str, object]:
    if region == "cerebellum":
        meta = CEREBELLAR_CELL_META.get(raw_group)
    else:
        meta = DG_CLUSTER_META.get(raw_group)
    if meta is None:
        return {
            "niche_role": f"{region}_other_context",
            "display": raw_group,
            "cell_axis": "context",
            "sender_class": "",
            "receiver_class": "",
            "core_focus": False,
        }
    return meta


def extract_group_expression_from_reader(
    *,
    dataset: str,
    region: str,
    sample: str,
    stage: str,
    reader,
    delimiter: str,
    label_map: dict[str, str],
    wanted_genes: set[str],
    mouse_to_human: dict[str, str],
    source_path: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    csv_reader = csv.reader(reader, delimiter=delimiter)
    header = next(csv_reader)
    cells = [str(cell).strip().strip('"') for cell in header[1:]]
    raw_labels = np.array([label_map.get(cell, "unmapped") for cell in cells], dtype=object)
    group_counts = Counter(raw_labels.tolist())
    group_masks = {group: np.flatnonzero(raw_labels == group) for group in group_counts}

    totals: dict[str, dict[str, float]] = {group: defaultdict(float) for group in group_counts}
    log_totals: dict[str, dict[str, float]] = {group: defaultdict(float) for group in group_counts}
    nonzeros: dict[str, dict[str, int]] = {group: defaultdict(int) for group in group_counts}
    source_gene_symbol: dict[str, str] = {}
    present_genes: set[str] = set()
    rows_scanned = 0
    n_value_cells = len(cells)

    for row in csv_reader:
        if not row:
            continue
        rows_scanned += 1
        raw_gene = str(row[0]).strip().strip('"')
        canonical_gene = mouse_to_human.get(canon(raw_gene), canon(raw_gene))
        if canonical_gene not in wanted_genes:
            continue
        values = row[1:]
        if len(values) < n_value_cells:
            values = [*values, *["0"] * (n_value_cells - len(values))]
        arr = np.fromiter((float(v) if v else 0.0 for v in values[:n_value_cells]), dtype=float, count=n_value_cells)
        present_genes.add(canonical_gene)
        source_gene_symbol[canonical_gene] = raw_gene
        for group, idx in group_masks.items():
            if group == "unmapped" or len(idx) == 0:
                continue
            sub = arr[idx]
            totals[group][canonical_gene] += float(sub.sum())
            log_totals[group][canonical_gene] += float(np.log1p(sub).sum())
            nonzeros[group][canonical_gene] += int(np.count_nonzero(sub > 0))

    rows: list[dict[str, object]] = []
    for group, n_cells in sorted(group_counts.items()):
        if group == "unmapped":
            continue
        meta = group_metadata(region, group)
        for gene in sorted(present_genes):
            nonzero = int(nonzeros[group].get(gene, 0))
            total = float(totals[group].get(gene, 0.0))
            log_total = float(log_totals[group].get(gene, 0.0))
            rows.append(
                {
                    "dataset": dataset,
                    "region": region,
                    "sample": sample,
                    "stage": stage,
                    "source_path": source_path,
                    "raw_group": group,
                    "niche_role": meta["niche_role"],
                    "display_group": meta["display"],
                    "cell_axis": meta["cell_axis"],
                    "sender_class": meta["sender_class"],
                    "receiver_class": meta["receiver_class"],
                    "core_focus": bool(meta["core_focus"]),
                    "n_cells": int(n_cells),
                    "gene": gene,
                    "source_gene_symbol": source_gene_symbol.get(gene, gene),
                    "nonzero_cells": nonzero,
                    "detection_fraction": nonzero / n_cells if n_cells else np.nan,
                    "total_expression": total,
                    "mean_expression": total / n_cells if n_cells else np.nan,
                    "mean_log1p_expression": log_total / n_cells if n_cells else np.nan,
                }
            )

    coverage = {
        "dataset": dataset,
        "region": region,
        "sample": sample,
        "stage": stage,
        "source_path": source_path,
        "n_matrix_cells": len(cells),
        "n_labeled_cells": int(sum(v for k, v in group_counts.items() if k != "unmapped")),
        "n_unmapped_cells": int(group_counts.get("unmapped", 0)),
        "n_gene_rows_scanned": rows_scanned,
        "n_lr_genes_requested": len(wanted_genes),
        "n_lr_genes_present": len(present_genes),
        "lr_genes_present": ",".join(sorted(present_genes)),
        "lr_genes_missing": ",".join(sorted(wanted_genes - present_genes)),
        "raw_group_counts": ";".join(f"{k}:{v}" for k, v in sorted(group_counts.items())),
    }
    return rows, coverage


def build_group_expression(pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted = set(pairs["ligand_canonical_gene"]) | set(pairs["receptor_canonical_gene"])
    mouse_to_human = load_mouse_to_human()
    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    cerebellar_maps = load_cerebellar_label_maps()
    with tarfile.open(GSE122357_TAR) as tf:
        for sample in CEREBELLAR_SAMPLES:
            with tar_member_text(tf, sample["tar_member"]) as fh:
                sample_rows, coverage = extract_group_expression_from_reader(
                    dataset=sample["dataset"],
                    region=sample["region"],
                    sample=sample["sample"],
                    stage=sample["stage"],
                    reader=fh,
                    delimiter=",",
                    label_map=cerebellar_maps[sample["sample"]],
                    wanted_genes=wanted,
                    mouse_to_human=mouse_to_human,
                    source_path=f"{rel(GSE122357_TAR)}:{sample['tar_member']}",
                )
            rows.extend(sample_rows)
            coverage_rows.append(coverage)

    dg_label_map = load_dg_label_map()
    with gzip.open(GSE104323_EXPR, "rt", newline="") as fh:
        sample_rows, coverage = extract_group_expression_from_reader(
            dataset=DG_SAMPLE["dataset"],
            region=DG_SAMPLE["region"],
            sample=DG_SAMPLE["sample"],
            stage=DG_SAMPLE["stage"],
            reader=fh,
            delimiter="\t",
            label_map=dg_label_map,
            wanted_genes=wanted,
            mouse_to_human=mouse_to_human,
            source_path=rel(GSE104323_EXPR),
        )
    rows.extend(sample_rows)
    coverage_rows.append(coverage)

    expr = pd.DataFrame(rows)
    if expr.empty:
        raise RuntimeError("No LR genes were extracted from focused sender-receiver sources")
    expr["sample_gene_max_log1p"] = expr.groupby(["dataset", "sample", "gene"])["mean_log1p_expression"].transform("max")
    expr["relative_expression"] = np.where(
        expr["sample_gene_max_log1p"].gt(0),
        expr["mean_log1p_expression"] / expr["sample_gene_max_log1p"],
        0.0,
    )
    expr["gene_detection_rank"] = expr.groupby(["dataset", "sample", "gene"])["detection_fraction"].rank(
        pct=True, method="average"
    )
    expr.to_csv(OUT_GROUP_EXPR, sep="\t", index=False, compression="gzip")
    coverage_df = pd.DataFrame(coverage_rows)
    coverage_df.to_csv(OUT_COVERAGE, sep="\t", index=False)
    return expr, coverage_df


def expression_lookup(expr: pd.DataFrame) -> dict[tuple[str, str, str, str], dict[str, object]]:
    return {
        (row.dataset, row.sample, row.raw_group, row.gene): row._asdict()
        for row in expr.itertuples(index=False)
    }


def support_from_scores(score: float, ligand_det: float, receptor_det: float) -> str:
    if score >= 0.25 and ligand_det >= 0.05 and receptor_det >= 0.05:
        return "high_expression_support"
    if score >= 0.10 and ligand_det >= 0.02 and receptor_det >= 0.02:
        return "moderate_expression_support"
    if score > 0 and ligand_det > 0 and receptor_det > 0:
        return "weak_detected"
    return "not_detected_or_one_sided"


def build_predictions(expr: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    lookup = expression_lookup(expr)
    rows: list[dict[str, object]] = []
    group_meta_cols = [
        "dataset",
        "region",
        "sample",
        "stage",
        "raw_group",
        "niche_role",
        "display_group",
        "cell_axis",
        "sender_class",
        "receiver_class",
        "core_focus",
        "n_cells",
    ]
    groups = expr[group_meta_cols].drop_duplicates()
    for (dataset, sample), sub in groups.groupby(["dataset", "sample"], sort=False):
        senders = sub.loc[sub["cell_axis"].eq("sender")].copy()
        receivers = sub.loc[sub["cell_axis"].eq("receiver")].copy()
        for sender in senders.itertuples(index=False):
            for receiver in receivers.itertuples(index=False):
                for pair in pairs.itertuples(index=False):
                    ligand = lookup.get((dataset, sample, sender.raw_group, pair.ligand_canonical_gene))
                    receptor = lookup.get((dataset, sample, receiver.raw_group, pair.receptor_canonical_gene))
                    ligand_rel = float(ligand["relative_expression"]) if ligand else 0.0
                    receptor_rel = float(receptor["relative_expression"]) if receptor else 0.0
                    ligand_det = float(ligand["detection_fraction"]) if ligand else 0.0
                    receptor_det = float(receptor["detection_fraction"]) if receptor else 0.0
                    lr_score = math.sqrt(max(ligand_rel * receptor_rel, 0.0)) * math.sqrt(
                        max(ligand_det * receptor_det, 0.0)
                    )
                    support = support_from_scores(lr_score, ligand_det, receptor_det)
                    rows.append(
                        {
                            "dataset": dataset,
                            "region": sender.region,
                            "sample": sample,
                            "stage": sender.stage,
                            "sender_raw_group": sender.raw_group,
                            "sender_display": sender.display_group,
                            "sender_niche_role": sender.niche_role,
                            "sender_class": sender.sender_class,
                            "sender_core_focus": bool(sender.core_focus),
                            "sender_n_cells": int(sender.n_cells),
                            "receiver_raw_group": receiver.raw_group,
                            "receiver_display": receiver.display_group,
                            "receiver_niche_role": receiver.niche_role,
                            "receiver_class": receiver.receiver_class,
                            "receiver_n_cells": int(receiver.n_cells),
                            "sender_receiver_pair_id": pair.sender_receiver_pair_id,
                            "source_pair_id": pair.pair_id,
                            "pair_source": pair.pair_source,
                            "pathway_id": pair.pathway_id,
                            "pathway_label": pair.pathway_label,
                            "pathway_family": pair.pathway_family,
                            "pair_label": pair.pair_label,
                            "ligand_gene": pair.ligand_canonical_gene,
                            "receptor_gene": pair.receptor_canonical_gene,
                            "ligand_sender_mean_log1p": float(ligand["mean_log1p_expression"]) if ligand else 0.0,
                            "ligand_sender_detection_fraction": ligand_det,
                            "ligand_sender_relative_expression": ligand_rel,
                            "receptor_receiver_mean_log1p": float(receptor["mean_log1p_expression"]) if receptor else 0.0,
                            "receptor_receiver_detection_fraction": receptor_det,
                            "receptor_receiver_relative_expression": receptor_rel,
                            "lr_expression_score": lr_score,
                            "support_class": support,
                            "supported": support in {"high_expression_support", "moderate_expression_support"},
                            "core_sender_receiver_focus": bool(sender.core_focus),
                        }
                    )
    pred = pd.DataFrame(rows)
    pred.to_csv(OUT_PREDICTIONS, sep="\t", index=False, compression="gzip")
    return pred


def summarize_predictions(pred: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    groupings = [
        ("region", ["region"]),
        ("region_sender_class", ["region", "sender_class"]),
        ("region_sender_pathway", ["region", "sender_class", "pathway_label"]),
        ("region_pathway", ["region", "pathway_label"]),
        ("region_pair", ["region", "pair_label"]),
    ]
    focus = pred.loc[pred["core_sender_receiver_focus"]].copy()
    for level, cols in groupings:
        for keys, sub in focus.groupby(cols, sort=False, dropna=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            row: dict[str, object] = {"summary_level": level}
            for col, key in zip(cols, keys):
                row[col] = key
            row.update(
                {
                    "n_predictions": int(len(sub)),
                    "n_supported": int(sub["supported"].sum()),
                    "fraction_supported": float(sub["supported"].mean()) if len(sub) else np.nan,
                    "median_lr_expression_score": float(sub["lr_expression_score"].median()),
                    "mean_lr_expression_score": float(sub["lr_expression_score"].mean()),
                    "max_lr_expression_score": float(sub["lr_expression_score"].max()),
                    "n_high_support": int(sub["support_class"].eq("high_expression_support").sum()),
                    "n_moderate_support": int(sub["support_class"].eq("moderate_expression_support").sum()),
                }
            )
            summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    top = (
        focus.loc[focus["supported"]]
        .sort_values(
            [
                "lr_expression_score",
                "ligand_sender_detection_fraction",
                "receptor_receiver_detection_fraction",
            ],
            ascending=False,
        )
        .head(250)
    )
    top.to_csv(OUT_TOP, sep="\t", index=False)
    return summary, top


def plot_results(pred: pd.DataFrame, summary: pd.DataFrame, coverage: pd.DataFrame) -> None:
    focus = pred.loc[pred["core_sender_receiver_focus"]].copy()
    sender_order = [
        "purkinje",
        "astroglial_bergmann_proxy",
        "microglia",
        "endothelial",
        "astrocyte",
        "microglia_macrophage_proxy",
        "vascular_support",
    ]
    heat = (
        focus.groupby(["region", "sender_class", "pathway_label"])["lr_expression_score"]
        .median()
        .reset_index()
    )
    heat["sender_context"] = heat["region"].map({"cerebellum": "CB", "dentate_sgz": "DG"}) + " " + heat[
        "sender_class"
    ].astype(str)
    heat["sender_context"] = pd.Categorical(
        heat["sender_context"],
        categories=[
            "CB purkinje",
            "CB astroglial_bergmann_proxy",
            "CB microglia",
            "CB endothelial",
            "DG astrocyte",
            "DG microglia_macrophage_proxy",
            "DG vascular_support",
            "DG endothelial",
        ],
        ordered=True,
    )
    heat_pivot = heat.pivot_table(
        index="pathway_label",
        columns="sender_context",
        values="lr_expression_score",
        aggfunc="median",
        observed=False,
    )
    heat_pivot = heat_pivot.dropna(axis=1, how="all").fillna(0.0)
    pathway_order = heat_pivot.mean(axis=1).sort_values(ascending=False).index.tolist()
    heat_pivot = heat_pivot.loc[pathway_order]

    top_pairs = (
        focus.loc[focus["supported"]]
        .groupby(["region", "sender_class", "pair_label"])["lr_expression_score"]
        .median()
        .reset_index()
        .sort_values("lr_expression_score", ascending=False)
        .head(14)
    )
    top_pairs["label"] = (
        top_pairs["region"].map({"cerebellum": "CB", "dentate_sgz": "DG"})
        + " "
        + top_pairs["sender_class"].str.replace("_", " ")
        + "\n"
        + top_pairs["pair_label"]
    )

    counts_rows = []
    for _, cov in coverage.iterrows():
        for part in str(cov["raw_group_counts"]).split(";"):
            if not part or ":" not in part:
                continue
            raw, count = part.rsplit(":", 1)
            if raw == "unmapped":
                continue
            meta = group_metadata(cov["region"], raw)
            if meta["cell_axis"] not in {"sender", "receiver"}:
                continue
            counts_rows.append(
                {
                    "region": cov["region"],
                    "sample": cov["sample"],
                    "raw_group": raw,
                    "display": meta["display"],
                    "cell_axis": meta["cell_axis"],
                    "core_focus": bool(meta["core_focus"]),
                    "n_cells": int(count),
                }
            )
    counts = pd.DataFrame(counts_rows)
    counts = counts.loc[counts["core_focus"]].copy()
    count_plot = counts.groupby(["region", "cell_axis"])["n_cells"].sum().reset_index()

    fig = plt.figure(figsize=(15, 9), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.45], height_ratios=[1.0, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    xlabels = [f"{r.region}\n{r.cell_axis}" for r in count_plot.itertuples()]
    ax.bar(np.arange(len(count_plot)), count_plot["n_cells"], color=["#7f4e8a", "#2f7f8f", "#7f4e8a", "#2f7f8f"][: len(count_plot)])
    ax.set_xticks(np.arange(len(count_plot)))
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_ylabel("Cells contributing to LR screen")
    ax.set_title("Sender/receiver coverage")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)

    ax = fig.add_subplot(gs[:, 1])
    im = ax.imshow(heat_pivot.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(0.35, float(heat_pivot.max().max())))
    ax.set_xticks(np.arange(heat_pivot.shape[1]))
    ax.set_xticklabels([str(c).replace("_", " ") for c in heat_pivot.columns], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(heat_pivot.shape[0]))
    ax.set_yticklabels(heat_pivot.index, fontsize=8)
    ax.set_title("Median sender->granule-lineage LR score")
    for i in range(heat_pivot.shape[0]):
        for j in range(heat_pivot.shape[1]):
            val = heat_pivot.iat[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.5, color="#1f2933")
    cbar = fig.colorbar(im, ax=ax, shrink=0.75)
    cbar.set_label("LR expression score")

    ax = fig.add_subplot(gs[1, 0])
    if not top_pairs.empty:
        y = np.arange(len(top_pairs))
        colors = top_pairs["region"].map({"cerebellum": "#7f4e8a", "dentate_sgz": "#2f7f8f"}).fillna("#777777")
        ax.barh(y, top_pairs["lr_expression_score"], color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels(top_pairs["label"], fontsize=7)
        ax.invert_yaxis()
    ax.set_xlabel("Median supported LR score")
    ax.set_title("Top supported sender/pair summaries")
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)

    fig.suptitle("Focused Niche Sender-Receiver Ligand-Receptor Prediction", fontsize=15)
    fig.savefig(OUT_PLOT, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(pred: pd.DataFrame, summary: pd.DataFrame, top: pd.DataFrame, coverage: pd.DataFrame, pairs: pd.DataFrame) -> None:
    focus = pred.loc[pred["core_sender_receiver_focus"]].copy()
    region_summary = summary.loc[summary["summary_level"].eq("region")].copy()
    sender_summary = summary.loc[summary["summary_level"].eq("region_sender_class")].copy()
    pathway_summary = summary.loc[summary["summary_level"].eq("region_pathway")].copy()

    lines = [
        "# Aim 2 Focused Sender-Receiver Ligand-Receptor Prediction",
        "",
        f"Date built: {date.today().isoformat()}",
        "",
        "## Purpose",
        "",
        "This analysis upgrades the earlier Aim 2 pathway-readiness audit by scoring ligand expression in niche sender classes and receptor expression in granule-lineage receiver classes.",
        "",
        "## Scope",
        "",
        "- Cerebellum: `GSE122357` P0, P8a, and P8b mouse cerebellum; senders are Purkinje cells, astrocytes as a Bergmann/astroglial proxy, microglia, endothelial cells, and supporting oligodendroglia; receivers are granule precursors and granule cells.",
        "- Dentate SGZ: `GSE104323` adult mouse dentate gyrus; senders are astrocyte states, endothelial cells, PVM/macrophage as a microglia-like proxy, and vascular/support classes; receivers are RGL, nIPC, neuroblast, immature GC, juvenile GC, and adult GC states.",
        "- LR database: previous Aim 2 curated pairs plus niche/glial/vascular pairs such as CXCL12/CXCR4, VEGF/FLT1/KDR, APOE/LRP, SPP1, CSF1, IGF1, SLIT/ROBO, and Ephrin/EPH.",
        "",
        "## Method",
        "",
        "- For each sample, raw matrices were reduced to LR genes and grouped by sender/receiver labels.",
        "- Ligand support is sender ligand relative expression multiplied by detection support.",
        "- Receptor support is receiver receptor relative expression multiplied by detection support.",
        "- The final LR expression score is the geometric support product, so a pair scores well only when the ligand is present in the sender and the receptor is present in the receiver.",
        "",
        "## Scale",
        "",
        f"- LR pairs scored: {pairs.shape[0]:,}.",
        f"- Sender-receiver-pair predictions: {focus.shape[0]:,} core-focus rows.",
        f"- Supported predictions: {int(focus['supported'].sum()):,} moderate/high expression-supported rows.",
        f"- Source samples: {coverage.shape[0]:,}.",
        "",
        "Region-level support:",
    ]
    for _, row in region_summary.iterrows():
        lines.append(
            f"- {row['region']}: {int(row['n_supported'])}/{int(row['n_predictions'])} supported, "
            f"median score {row['median_lr_expression_score']:.3f}, max {row['max_lr_expression_score']:.3f}."
        )

    lines.extend(["", "Sender-class highlights:"])
    for _, row in sender_summary.sort_values(["region", "median_lr_expression_score"], ascending=[True, False]).iterrows():
        lines.append(
            f"- {row['region']} / {row['sender_class']}: {int(row['n_supported'])}/{int(row['n_predictions'])} supported, "
            f"median score {row['median_lr_expression_score']:.3f}."
        )

    lines.extend(["", "Top pathway summaries:"])
    top_path = pathway_summary.sort_values("median_lr_expression_score", ascending=False).head(12)
    for _, row in top_path.iterrows():
        lines.append(
            f"- {row['region']} / {row['pathway_label']}: median score {row['median_lr_expression_score']:.3f}; "
            f"{int(row['n_supported'])}/{int(row['n_predictions'])} supported."
        )

    lines.extend(["", "Top individual predictions:"])
    display_top = top.head(12)
    for _, row in display_top.iterrows():
        lines.append(
            f"- {row['region']} {row['sender_display']} -> {row['receiver_display']}: "
            f"{row['pair_label']} ({row['ligand_gene']}->{row['receptor_gene']}), score {row['lr_expression_score']:.3f}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This gives us the first focused sender-receiver LR layer for Aim 2.",
            "- The cerebellar analysis directly includes Purkinje, microglia, and endothelial senders; Bergmann glia are not separately annotated in `GSE122357`, so the astrocyte class should be described as a Bergmann/astroglial proxy.",
            "- The dentate analysis directly includes astrocyte and endothelial senders; `GSE104323` has PVM rather than a clean microglia class, so SGZ microglia claims should be phrased as PVM/microglia-macrophage proxy unless we add a human or mouse dataset with explicit microglia labels.",
            "- The result is stronger than the previous pathway-readiness audit for niche directionality, but still not proof of spatial contact or protein secretion.",
            "",
            "## Outputs",
            "",
            f"- LR pair table: `{rel(OUT_PAIRS)}`",
            f"- Group expression table: `{rel(OUT_GROUP_EXPR)}`",
            f"- Prediction table: `{rel(OUT_PREDICTIONS)}`",
            f"- Summary table: `{rel(OUT_SUMMARY)}`",
            f"- Top predictions: `{rel(OUT_TOP)}`",
            f"- Coverage table: `{rel(OUT_COVERAGE)}`",
            f"- Plot: `{rel(OUT_PLOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    pairs = load_pairs()
    expr, coverage = build_group_expression(pairs)
    pred = build_predictions(expr, pairs)
    summary, top = summarize_predictions(pred)
    plot_results(pred, summary, coverage)
    write_report(pred, summary, top, coverage, pairs)
    print(f"Wrote {rel(OUT_MD)}")
    print(summary.loc[summary["summary_level"].isin(["region", "region_sender_class"])].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
