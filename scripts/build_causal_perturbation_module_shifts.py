#!/usr/bin/env python3
"""Build module-shift signatures from public perturbation datasets.

This script turns the perturbation-resource triage into a preliminary
gene-module evidence layer. It intentionally scores curated module genes rather
than attempting a full causal model: the public datasets differ in modality,
lineage specificity, and annotation depth, so the most reproducible first pass
is to ask whether each perturbation shifts the same fate, niche, pathway, and
construction modules used elsewhere in the manuscript.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "Project/dataset_search_cache/causal_perturbation"
RESULTS = ROOT / "Project/results"
PROCESSED = ROOT / "Project/processed/causal_perturbation"

NICHE_MODULES = RESULTS / "primary_core_niche_circuit_module_gene_sets.tsv"
PATHWAY_MODULES = RESULTS / "primary_core_aim2_pathway_gene_sets.tsv"
TRIAGE = RESULTS / "causal_perturbation_dataset_triage.tsv"

GSE84786_CSV = CACHE / "GSE84786_gene_exp.csv"
GSE71916_RBFOX_DIFF = CACHE / "GSE71916_Expression_Cufflinks_siRbfox1_3_EGFP_vs_siNT_EGFP.txt.gz"
GSE242199_COUNTS = CACHE / "GSE242199_Raw_counts.txt.gz"
GSE81962_MATRIX = CACHE / "GSE81962_series_matrix.txt.gz"
GPL6887_ANNOT = CACHE / "GPL6887.annot.gz"
GSE107252_MATRIX = CACHE / "GSE107252_series_matrix.txt.gz"

GSE268609_SELECTED = RESULTS / "gse268609_selected_gene_presence.tsv"
GSE268609_PEAKS = RESULTS / "gse268609_epigenomic_peak_gene_summary.tsv"
GSE322785_H5 = RESULTS / "gse322785_human_h5_epigenomic_gene_summary.tsv"

OUT_MODULE_CATALOG = RESULTS / "causal_perturbation_module_catalog.tsv"
OUT_GENE_EFFECTS = RESULTS / "causal_perturbation_module_shift_gene_effects.tsv"
OUT_SUMMARY = RESULTS / "causal_perturbation_module_shift_summary.tsv"
OUT_NODE_SUMMARY = RESULTS / "causal_perturbation_module_shift_node_summary.tsv"
OUT_STATUS = RESULTS / "causal_perturbation_processing_status.tsv"
OUT_MD = RESULTS / "causal_perturbation_module_shift.md"
OUT_PNG = RESULTS / "causal_perturbation_module_shift_heatmap.png"


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_symbol(value: Any) -> str:
    """Return an uppercase single-gene symbol for cross-species module scoring."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip().strip('"').strip()
    if not text or text.lower() == "nan":
        return ""
    text = re.split(r"\s*///\s*|[;,]\s*", text)[0].strip()
    text = re.sub(r"\s+", "", text)
    return text.upper()


def safe_float(value: Any, default: float = np.nan) -> float:
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def effect_call(value: float, threshold: float = 0.10) -> str:
    if not np.isfinite(value):
        return "not_scored"
    if value >= threshold:
        return "up_shift"
    if value <= -threshold:
        return "down_shift"
    return "near_zero"


def support_class(abs_score: float, coverage: float, n_genes: int) -> str:
    if n_genes < 2 or coverage < 0.15:
        return "low_coverage"
    if abs_score >= 0.35 and coverage >= 0.35:
        return "strong_shift"
    if abs_score >= 0.20 and coverage >= 0.25:
        return "moderate_shift"
    if abs_score >= 0.10:
        return "weak_shift"
    return "near_zero_or_mixed"


def context_match_for(node: str, module_id: str) -> str:
    """Flag whether a module is biologically matched to the perturbation context."""
    rbfox_modules = {
        "downstream_neurite_morphology",
        "downstream_synaptic_excitability",
        "semaphorin_guidance",
        "reelin_migration_stop",
        "bdnf_trkb_mapk",
        "fgf_mapk",
    }
    bdnf_modules = {
        "bdnf_trkb_mapk",
        "downstream_neurite_morphology",
        "downstream_synaptic_excitability",
        "fgf_mapk",
        "tgf_beta_smad",
    }
    shh_modules = {
        "cerebellar_fate_rhombic_lip_shh",
        "shh_granule_expansion",
        "shared_neurogenic_niche_state",
        "downstream_neurite_morphology",
        "downstream_synaptic_excitability",
    }
    if node.startswith("RBFOX3"):
        return "matched_context" if module_id in rbfox_modules else "off_context"
    if node.startswith("BDNF"):
        return "matched_context" if module_id in bdnf_modules else "off_context"
    if node.startswith("SHH"):
        return "matched_context" if module_id in shh_modules else "off_context"
    return "unclassified_context"


def build_module_catalog() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    niche = read_tsv(NICHE_MODULES)
    for _, row in niche.iterrows():
        rows.append(
            {
                "module_id": row["module_id"],
                "module_label": row["module_label"],
                "module_family": row["module_family"],
                "hypothesis_role": row["hypothesis_role"],
                "gene_symbol": clean_symbol(row["canonical_gene"]),
                "source_table": NICHE_MODULES.name,
            }
        )
    pathways = read_tsv(PATHWAY_MODULES)
    for _, row in pathways.iterrows():
        rows.append(
            {
                "module_id": row["pathway_id"],
                "module_label": row["pathway_label"],
                "module_family": row["pathway_family"],
                "hypothesis_role": row["hypothesis_role"],
                "gene_symbol": clean_symbol(row["canonical_gene"]),
                "source_table": PATHWAY_MODULES.name,
            }
        )
    catalog = pd.DataFrame(rows)
    catalog = catalog.loc[catalog["gene_symbol"].ne("")].drop_duplicates()
    return catalog


def module_gene_sets(catalog: pd.DataFrame) -> dict[str, set[str]]:
    return {
        module_id: set(sub["gene_symbol"])
        for module_id, sub in catalog.groupby("module_id", sort=False)
    }


def module_metadata(catalog: pd.DataFrame) -> dict[str, dict[str, str]]:
    meta = {}
    for module_id, sub in catalog.groupby("module_id", sort=False):
        row = sub.iloc[0]
        meta[module_id] = {
            "module_label": str(row["module_label"]),
            "module_family": str(row["module_family"]),
            "hypothesis_role": str(row["hypothesis_role"]),
        }
    return meta


def build_ensembl_symbol_map() -> dict[str, str]:
    """Use project-local selected-feature maps before any external annotation."""
    ens_to_symbol: dict[str, str] = {}

    if GSE268609_SELECTED.exists():
        df = read_tsv(GSE268609_SELECTED)
        for _, row in df.iterrows():
            gene = clean_symbol(row.get("gene"))
            ens = str(row.get("gse268609_feature_id", "")).split(".")[0]
            if gene and ens.startswith("ENSG"):
                ens_to_symbol.setdefault(ens, gene)

    for path in [GSE268609_PEAKS, GSE322785_H5]:
        if not path.exists():
            continue
        df = read_tsv(path)
        for _, row in df.iterrows():
            gene = clean_symbol(row.get("gene"))
            ens = str(row.get("gene_feature_id", "")).split(".")[0]
            if gene and ens.startswith("ENSG"):
                ens_to_symbol.setdefault(ens, gene)

    return ens_to_symbol


def add_effect(
    rows: list[dict[str, Any]],
    *,
    dataset: str,
    accession: str,
    node: str,
    contrast_id: str,
    contrast_label: str,
    species: str,
    model_or_tissue: str,
    gene_symbol: str,
    effect_value: float,
    effect_metric: str,
    source_file: str,
    p_value: float = np.nan,
    q_value: float = np.nan,
    source_gene_id: str = "",
    notes: str = "",
) -> None:
    gene_symbol = clean_symbol(gene_symbol)
    if not gene_symbol or not np.isfinite(effect_value):
        return
    rows.append(
        {
            "dataset": dataset,
            "accession": accession,
            "node": node,
            "contrast_id": contrast_id,
            "contrast_label": contrast_label,
            "species": species,
            "model_or_tissue": model_or_tissue,
            "gene_symbol": gene_symbol,
            "source_gene_id": source_gene_id,
            "effect_metric": effect_metric,
            "effect_value": effect_value,
            "effect_direction": effect_call(effect_value),
            "p_value": p_value,
            "q_value": q_value,
            "source_file": source_file,
            "notes": notes,
        }
    )


def process_gse84786(module_union: set[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not GSE84786_CSV.exists():
        return rows, {"accession": "GSE84786", "status": "missing", "n_gene_effects": 0, "notes": "CSV converted from xls not found"}

    df = pd.read_csv(GSE84786_CSV)
    for _, row in df.iterrows():
        gene = clean_symbol(row.get("gene"))
        if gene not in module_union:
            continue
        wt = safe_float(row.get("value_1"))
        ko = safe_float(row.get("value_2"))
        if not (np.isfinite(wt) and np.isfinite(ko)):
            continue
        effect = float(np.log2((ko + 0.25) / (wt + 0.25)))
        add_effect(
            rows,
            dataset="Rbfox3 KO hippocampus",
            accession="GSE84786",
            node="RBFOX3",
            contrast_id="GSE84786_Rbfox3_KO_vs_WT",
            contrast_label="Rbfox3 homozygous KO hippocampus vs WT",
            species="mouse",
            model_or_tissue="hippocampus; dentate phenotype reported",
            gene_symbol=gene,
            source_gene_id=str(row.get("gene_id", "")),
            effect_value=effect,
            effect_metric="log2((KO_FPKM+0.25)/(WT_FPKM+0.25))",
            p_value=safe_float(row.get("p_value")),
            q_value=safe_float(row.get("q_value")),
            source_file=GSE84786_CSV.name,
            notes="GEO processed Cuffdiff-style gene table; low RNA-seq sample count in source study.",
        )
    return rows, {"accession": "GSE84786", "status": "processed", "n_gene_effects": len(rows), "notes": "direct RBFOX3 hippocampal KO table"}


def process_gse71916(module_union: set[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not GSE71916_RBFOX_DIFF.exists():
        return rows, {"accession": "GSE71916", "status": "missing", "n_gene_effects": 0, "notes": "differential table not found"}

    df = pd.read_csv(GSE71916_RBFOX_DIFF, sep="\t")
    for _, row in df.iterrows():
        gene = clean_symbol(row.get("gene"))
        if gene not in module_union:
            continue
        effect = safe_float(row.get("log2.fold_change"))
        add_effect(
            rows,
            dataset="Rbfox1/3 knockdown hippocampal neurons",
            accession="GSE71916",
            node="RBFOX3/RBFOX1",
            contrast_id="GSE71916_siRbfox1_3_vs_siNT",
            contrast_label="siRbfox1/3 hippocampal neurons vs non-targeting control",
            species="mouse",
            model_or_tissue="cultured hippocampal neurons",
            gene_symbol=gene,
            source_gene_id=str(row.get("gene_id", "")),
            effect_value=effect,
            effect_metric="reported_log2_fold_change",
            p_value=safe_float(row.get("p_value")),
            q_value=safe_float(row.get("q_value")),
            source_file=GSE71916_RBFOX_DIFF.name,
            notes="Double Rbfox1/Rbfox3 knockdown; not RBFOX3-specific.",
        )
    return rows, {"accession": "GSE71916", "status": "processed", "n_gene_effects": len(rows), "notes": "direct RBFOX-family neuronal knockdown table"}


def process_gse242199(module_union: set[str], ens_to_symbol: dict[str, str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not GSE242199_COUNTS.exists():
        return rows, {"accession": "GSE242199", "status": "missing", "n_gene_effects": 0, "notes": "count matrix not found"}

    df = pd.read_csv(GSE242199_COUNTS, sep="\t")
    df.columns = [str(c).strip('"') for c in df.columns]
    df["ensembl_id"] = df["Gene_id"].astype(str).str.strip('"').str.split(".").str[0]
    df["gene_symbol"] = df["ensembl_id"].map(ens_to_symbol).fillna("")
    df = df.loc[df["gene_symbol"].isin(module_union)].copy()
    if df.empty:
        return rows, {"accession": "GSE242199", "status": "processed_low_mapping", "n_gene_effects": 0, "notes": "no module genes mapped from local Ensembl-symbol map"}

    sample_cols = [c for c in df.columns if c not in {"Gene_id", "ensembl_id", "gene_symbol"}]
    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib = counts.sum(axis=0).replace(0, np.nan)
    log_cpm = np.log1p(counts.div(lib, axis=1) * 1_000_000.0)
    log_cpm.index = df.index

    contrasts = [
        ("GSE242199_WT_3h_BDNF_vs_vehicle", "WT 3h BDNF vs WT 3h vehicle", "WT_3h_BDNF", "WT_3h_vehicle"),
        ("GSE242199_WT_24h_BDNF_vs_vehicle", "WT 24h BDNF vs WT 24h vehicle", "WT_24h_BDNF", "WT_24h_vehicle"),
        ("GSE242199_NTRK2ko_3h_BDNF_vs_vehicle", "NTRK2 KO 3h BDNF vs NTRK2 KO 3h vehicle", "NTRK2ko_3h_BDNF", "NTRK2ko_3h_vehicle"),
        ("GSE242199_NTRK2ko_24h_BDNF_vs_vehicle", "NTRK2 KO 24h BDNF vs NTRK2 KO 24h vehicle", "NTRK2ko_24h_BDNF", "NTRK2ko_24h_vehicle"),
        ("GSE242199_NTRK2ko_0h_vs_WT_0h", "NTRK2 KO baseline vs WT baseline", "NTRK2ko_0h", "WT_0h"),
    ]
    for contrast_id, label, target_prefix, control_prefix in contrasts:
        target_cols = [c for c in sample_cols if c.startswith(target_prefix)]
        control_cols = [c for c in sample_cols if c.startswith(control_prefix)]
        if not target_cols or not control_cols:
            continue
        effect = log_cpm[target_cols].mean(axis=1) - log_cpm[control_cols].mean(axis=1)
        for idx, value in effect.items():
            add_effect(
                rows,
                dataset="NTRK2 KO and BDNF response ReNcell VM",
                accession="GSE242199",
                node="BDNF/NTRK2",
                contrast_id=contrast_id,
                contrast_label=label,
                species="human",
                model_or_tissue="ReNcell VM neural progenitor/stem-cell differentiation model",
                gene_symbol=df.at[idx, "gene_symbol"],
                source_gene_id=df.at[idx, "ensembl_id"],
                effect_value=float(value),
                effect_metric="mean_log1p_CPM_difference",
                source_file=GSE242199_COUNTS.name,
                notes="Targeted module scoring using project-local Ensembl-symbol maps.",
            )
    return rows, {"accession": "GSE242199", "status": "processed", "n_gene_effects": len(rows), "notes": f"mapped {df['gene_symbol'].nunique()} module genes"}


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    metadata: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t"))
                metadata[parts[0]] = [p.strip('"') for p in parts[1:]]
    if not table_lines:
        raise ValueError(f"No matrix table found in {path}")
    matrix = pd.read_csv(pd.io.common.StringIO("\n".join(table_lines)), sep="\t")
    matrix.columns = [str(c).strip('"') for c in matrix.columns]
    matrix["ID_REF"] = matrix["ID_REF"].astype(str).str.strip('"')
    return matrix, metadata


def load_gpl6887_annotation() -> pd.DataFrame:
    rows: list[str] = []
    in_table = False
    with gzip.open(GPL6887_ANNOT, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                rows.append(line)
    annot = pd.read_csv(pd.io.common.StringIO("\n".join(rows)), sep="\t", dtype=str)
    annot = annot[["ID", "Gene symbol"]].copy()
    annot["gene_symbol"] = annot["Gene symbol"].map(clean_symbol)
    annot = annot.loc[annot["gene_symbol"].ne("")]
    return annot[["ID", "gene_symbol"]].drop_duplicates()


def gse81962_group_from_title(title: str) -> str:
    text = title.lower()
    if "medullob" in text:
        return "tumor_excluded"
    if "ndp ko ptch het" in text or "ndp knockout ptch" in text:
        return "Ndp_KO_Ptch_het_GNP"
    if "ptch heterozygous" in text:
        return "Ptch_het_GNP"
    if "ndp knockout" in text:
        return "Ndp_KO_GNP"
    if text.startswith("wt "):
        return "WT_GNP"
    return "unassigned"


def process_gse81962(module_union: set[str]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not (GSE81962_MATRIX.exists() and GPL6887_ANNOT.exists()):
        return rows, {"accession": "GSE81962", "status": "missing", "n_gene_effects": 0, "notes": "matrix or GPL6887 annotation missing"}

    matrix, metadata = parse_series_matrix(GSE81962_MATRIX)
    sample_ids = [c for c in matrix.columns if c != "ID_REF"]
    titles = metadata.get("!Sample_title", [])
    groups = {sample: gse81962_group_from_title(title) for sample, title in zip(sample_ids, titles)}
    keep_samples = [s for s in sample_ids if groups.get(s) not in {"tumor_excluded", "unassigned"}]

    values = matrix[["ID_REF", *keep_samples]].copy()
    for col in keep_samples:
        values[col] = pd.to_numeric(values[col], errors="coerce")
    annot = load_gpl6887_annotation()
    values = values.merge(annot, left_on="ID_REF", right_on="ID", how="inner")
    values = values.loc[values["gene_symbol"].isin(module_union)]
    if values.empty:
        return rows, {"accession": "GSE81962", "status": "processed_low_mapping", "n_gene_effects": 0, "notes": "no module genes mapped through GPL6887"}

    gene_expr = values.groupby("gene_symbol", as_index=True)[keep_samples].median()
    group_to_cols: dict[str, list[str]] = defaultdict(list)
    for sample in keep_samples:
        group_to_cols[groups[sample]].append(sample)

    contrasts = [
        ("GSE81962_Ptch_het_GNP_vs_WT_GNP", "P6 Ptch heterozygous GNP vs WT GNP", "Ptch_het_GNP", "WT_GNP"),
        ("GSE81962_Ndp_KO_GNP_vs_WT_GNP", "P6 Ndp KO GNP vs WT GNP", "Ndp_KO_GNP", "WT_GNP"),
        ("GSE81962_Ndp_KO_Ptch_het_GNP_vs_Ptch_het_GNP", "P6 Ndp KO/Ptch het GNP vs Ptch het GNP", "Ndp_KO_Ptch_het_GNP", "Ptch_het_GNP"),
        ("GSE81962_Ndp_KO_Ptch_het_GNP_vs_WT_GNP", "P6 Ndp KO/Ptch het GNP vs WT GNP", "Ndp_KO_Ptch_het_GNP", "WT_GNP"),
    ]
    for contrast_id, label, target_group, control_group in contrasts:
        target_cols = group_to_cols.get(target_group, [])
        control_cols = group_to_cols.get(control_group, [])
        if not target_cols or not control_cols:
            continue
        effect = gene_expr[target_cols].mean(axis=1) - gene_expr[control_cols].mean(axis=1)
        for gene, value in effect.items():
            add_effect(
                rows,
                dataset="P6 cerebellar GNP Ptch/Ndp array",
                accession="GSE81962",
                node="SHH/PTCH/Norrin",
                contrast_id=contrast_id,
                contrast_label=label,
                species="mouse",
                model_or_tissue="P6 cerebellar granule neuron progenitors",
                gene_symbol=gene,
                effect_value=float(value),
                effect_metric="mean_log_array_expression_difference",
                source_file=f"{GSE81962_MATRIX.name};{GPL6887_ANNOT.name}",
                notes="Tumor samples excluded; probe-level values collapsed to gene medians.",
            )
    return rows, {"accession": "GSE81962", "status": "processed", "n_gene_effects": len(rows), "notes": f"kept {len(keep_samples)} GNP samples; mapped {gene_expr.shape[0]} module genes"}


def not_processed_status() -> list[dict[str, Any]]:
    rows = []
    if GSE107252_MATRIX.exists():
        rows.append(
            {
                "accession": "GSE107252",
                "status": "deferred_annotation",
                "n_gene_effects": 0,
                "notes": "series matrix downloaded, but GPL16570 platform annotation is too large for this first pass; retained for later targeted annotation or external chip map",
            }
        )
    return rows


def score_modules(gene_effects: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    gene_sets = module_gene_sets(catalog)
    meta = module_metadata(catalog)
    rows: list[dict[str, Any]] = []

    contrast_cols = [
        "dataset",
        "accession",
        "node",
        "contrast_id",
        "contrast_label",
        "species",
        "model_or_tissue",
    ]
    for contrast_id, sub in gene_effects.groupby("contrast_id", sort=False):
        info = sub.iloc[0][contrast_cols].to_dict()
        effect_by_gene = sub.groupby("gene_symbol")["effect_value"].median()
        for module_id, genes in gene_sets.items():
            vals = effect_by_gene.reindex(sorted(genes)).dropna()
            n_total = len(genes)
            n_avail = int(vals.shape[0])
            coverage = n_avail / n_total if n_total else np.nan
            if n_avail:
                median_effect = float(vals.median())
                mean_effect = float(vals.mean())
                pos_fraction = float((vals > 0).mean())
                neg_fraction = float((vals < 0).mean())
                abs_median = abs(median_effect)
                score = median_effect * math.sqrt(max(coverage, 0.0))
                top_genes = ";".join(vals.abs().sort_values(ascending=False).head(5).index.tolist())
            else:
                median_effect = np.nan
                mean_effect = np.nan
                pos_fraction = np.nan
                neg_fraction = np.nan
                abs_median = np.nan
                score = np.nan
                top_genes = ""
            rows.append(
                {
                    **info,
                    "module_id": module_id,
                    "module_label": meta[module_id]["module_label"],
                    "module_family": meta[module_id]["module_family"],
                    "hypothesis_role": meta[module_id]["hypothesis_role"],
                    "n_module_genes": n_total,
                    "n_genes_available": n_avail,
                    "coverage_fraction": coverage,
                    "median_effect": median_effect,
                    "mean_effect": mean_effect,
                    "positive_fraction": pos_fraction,
                    "negative_fraction": neg_fraction,
                    "signed_module_shift_score": score,
                    "direction_call": effect_call(median_effect),
                    "support_class": support_class(abs_median, coverage, n_avail),
                    "context_match": context_match_for(info["node"], module_id),
                    "top_abs_effect_genes": top_genes,
                }
            )
    return pd.DataFrame(rows)


def summarize_by_node(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (node, module_id), sub in summary.groupby(["node", "module_id"], sort=False):
        valid = sub.dropna(subset=["signed_module_shift_score"])
        if valid.empty:
            continue
        meta = valid.iloc[0]
        score = float(valid["signed_module_shift_score"].median())
        abs_score = abs(score)
        rows.append(
            {
                "node": node,
                "module_id": module_id,
                "module_label": meta["module_label"],
                "module_family": meta["module_family"],
                "context_match": "matched_context" if (valid["context_match"] == "matched_context").any() else valid["context_match"].iloc[0],
                "n_contrasts": valid["contrast_id"].nunique(),
                "median_signed_module_shift_score": score,
                "median_abs_module_shift_score": abs_score,
                "median_coverage_fraction": float(valid["coverage_fraction"].median()),
                "dominant_direction": effect_call(score),
                "node_module_support_class": support_class(abs_score, float(valid["coverage_fraction"].median()), int(valid["n_genes_available"].median())),
                "strongest_contrast": valid.loc[valid["signed_module_shift_score"].abs().idxmax(), "contrast_id"],
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["context_sort"] = out["context_match"].map({"matched_context": 0, "off_context": 1}).fillna(2).astype(int)
    out = out.sort_values(["context_sort", "median_abs_module_shift_score", "median_coverage_fraction"], ascending=[True, False, False])
    return out.drop(columns=["context_sort"])


def plot_heatmap(summary: pd.DataFrame) -> None:
    plot_df = summary.pivot_table(
        index="module_label",
        columns="contrast_id",
        values="signed_module_shift_score",
        aggfunc="median",
    )
    if plot_df.empty:
        return

    preferred_modules = [
        "Cerebellar fate/rhombic-lip/SHH",
        "SHH/PTCH/GLI",
        "Shared neurogenic niche/progenitor state",
        "BDNF/TrkB/MAPK",
        "TGF-beta/SMAD",
        "BMP/SMAD",
        "WNT/beta-catenin",
        "Downstream neurite/morphology",
        "Downstream synaptic/excitability",
        "Semaphorin/plexin guidance",
        "Reelin/migration-stop",
    ]
    plot_df = plot_df.reindex([m for m in preferred_modules if m in plot_df.index])
    vmax = float(np.nanmax(np.abs(plot_df.values))) if np.isfinite(plot_df.values).any() else 1.0
    vmax = max(vmax, 0.25)

    fig_w = max(10.5, 0.72 * plot_df.shape[1] + 4.5)
    fig_h = max(6.0, 0.42 * plot_df.shape[0] + 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(plot_df.values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(np.arange(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(plot_df.shape[0]))
    ax.set_yticklabels(plot_df.index, fontsize=9)
    ax.set_title("Public perturbation module-shift signatures", fontsize=13, pad=12)

    for i in range(plot_df.shape[0]):
        for j in range(plot_df.shape[1]):
            val = plot_df.iat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="black")

    ax.set_xlabel("Perturbation contrast")
    ax.set_ylabel("Granule convergence module")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Signed module-shift score")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300)
    plt.close(fig)


def write_markdown(summary: pd.DataFrame, node_summary: pd.DataFrame, status: pd.DataFrame) -> None:
    processed = status.loc[status["status"].str.startswith("processed", na=False)]
    lines = [
        "# Causal Perturbation Module-Shift Layer",
        "",
        "This analysis converts public perturbation resources into preliminary module-shift signatures. The scored modules are the same curated fate, niche, pathway, neurite, synaptic, and excitability modules used in the primary manuscript analyses.",
        "",
        "## Processing Status",
        "",
        f"- Processed accessions: {', '.join(processed['accession'].astype(str)) if len(processed) else 'none'}",
        f"- Deferred accessions: {', '.join(status.loc[~status['status'].str.startswith('processed', na=False), 'accession'].astype(str)) if len(status) > len(processed) else 'none'}",
        f"- Gene-level module effects: {len(pd.read_csv(OUT_GENE_EFFECTS, sep='\t')) if OUT_GENE_EFFECTS.exists() else 0}",
        f"- Contrast-module summaries: {len(summary)}",
        "",
        "## Strongest Matched-Context Node-Module Shifts",
        "",
    ]
    if node_summary.empty:
        lines.append("No node-module shifts were scored.")
    else:
        lines.extend(
            [
                "| Node | Module | Direction | Median signed score | Contrasts | Context | Strongest contrast |",
                "|---|---|---|---:|---:|---|---|",
            ]
        )
        display = node_summary.loc[node_summary["context_match"].eq("matched_context")].head(16)
        if display.empty:
            display = node_summary.head(16)
        for _, row in display.iterrows():
            lines.append(
                f"| {row['node']} | {row['module_label']} | {row['dominant_direction']} | "
                f"{row['median_signed_module_shift_score']:.3f} | {int(row['n_contrasts'])} | "
                f"{row['context_match']} | {row['strongest_contrast']} |"
            )

    lines.extend(
        [
            "",
            "Off-context shifts are retained in the tables as exploratory signals, but the manuscript interpretation should emphasize matched-context shifts.",
            "",
            "## Strongest Matched-Context Contrast-Level Shifts",
            "",
            "| Contrast | Node | Module | Direction | Signed score | Available genes | Top shifted genes |",
            "|---|---|---|---|---:|---:|---|",
        ]
    )
    contrast_display = summary.loc[summary["context_match"].eq("matched_context")].copy()
    if contrast_display.empty:
        lines.append("| Not scored | Not scored | Not scored | Not scored | 0.000 | 0 |  |")
    else:
        contrast_display["abs_score"] = contrast_display["signed_module_shift_score"].abs()
        for _, row in contrast_display.sort_values("abs_score", ascending=False).head(16).iterrows():
            lines.append(
                f"| {row['contrast_id']} | {row['node']} | {row['module_label']} | {row['direction_call']} | "
                f"{row['signed_module_shift_score']:.3f} | {int(row['n_genes_available'])} | {row['top_abs_effect_genes']} |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The output should be read as perturbation-sensitive module evidence, not as a causal mixed-effects model. It is strongest when a contrast is lineage-relevant, has adequate module-gene coverage, and shifts biologically matched modules in the expected branch. The current layer is therefore useful for prioritizing follow-up perturbations and for supporting the hierarchical integrative model, while still requiring matched single-cell or perturb-seq data for formal causality.",
            "",
            "## Output Files",
            "",
            f"- Module catalog: `{OUT_MODULE_CATALOG.name}`",
            f"- Gene-level effects: `{OUT_GENE_EFFECTS.name}`",
            f"- Contrast-module summary: `{OUT_SUMMARY.name}`",
            f"- Node-module summary: `{OUT_NODE_SUMMARY.name}`",
            f"- Processing status: `{OUT_STATUS.name}`",
            f"- Heatmap: `{OUT_PNG.name}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    catalog = build_module_catalog()
    catalog.to_csv(OUT_MODULE_CATALOG, sep="\t", index=False)
    module_union = set(catalog["gene_symbol"])
    ens_to_symbol = build_ensembl_symbol_map()

    effect_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []

    for fn in [
        process_gse84786,
        process_gse71916,
    ]:
        rows, status = fn(module_union)
        effect_rows.extend(rows)
        status_rows.append(status)

    rows, status = process_gse242199(module_union, ens_to_symbol)
    effect_rows.extend(rows)
    status_rows.append(status)

    rows, status = process_gse81962(module_union)
    effect_rows.extend(rows)
    status_rows.append(status)
    status_rows.extend(not_processed_status())

    gene_effects = pd.DataFrame(effect_rows)
    if gene_effects.empty:
        raise RuntimeError("No gene effects were generated.")
    gene_effects = gene_effects.sort_values(["accession", "contrast_id", "gene_symbol"])
    gene_effects.to_csv(OUT_GENE_EFFECTS, sep="\t", index=False)

    summary = score_modules(gene_effects, catalog)
    summary = summary.sort_values(["accession", "contrast_id", "module_id"])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    node_summary = summarize_by_node(summary)
    node_summary.to_csv(OUT_NODE_SUMMARY, sep="\t", index=False)

    status_df = pd.DataFrame(status_rows)
    status_df.to_csv(OUT_STATUS, sep="\t", index=False)

    plot_heatmap(summary)
    write_markdown(summary, node_summary, status_df)

    print(f"Wrote {OUT_MODULE_CATALOG.relative_to(ROOT)}")
    print(f"Wrote {OUT_GENE_EFFECTS.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_NODE_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_STATUS.relative_to(ROOT)}")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")
    print(f"Wrote {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
