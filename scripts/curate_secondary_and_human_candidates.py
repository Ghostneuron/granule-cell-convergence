#!/usr/bin/env python3
"""Curate promotion status for secondary dentate datasets and human DG candidates."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "External_Data"
RESULTS = ROOT / "Project" / "results"

PROMOTED_METADATA_OUT = RESULTS / "promoted_secondary_dataset_metadata.tsv"
READINESS_OUT = RESULTS / "secondary_dataset_promotion_readiness.tsv"
HUMAN_CANDIDATES_OUT = RESULTS / "human_dentate_hippocampus_candidate_datasets.tsv"
INTERPRETATION_OUT = RESULTS / "secondary_promotion_and_human_dentate_interpretation.md"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip().strip('"').strip("'")


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Parse GEO series-matrix sample-level fields and characteristics."""
    rows: list[dict[str, str]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            field = parts[0].replace("!Sample_", "")
            values = [clean(value) for value in parts[1:]]
            while len(rows) < len(values):
                rows.append({})
            if field == "characteristics_ch1":
                for i, value in enumerate(values):
                    if ": " in value:
                        key, val = value.split(": ", 1)
                        rows[i][key.lower().replace(" ", "_")] = val
                    else:
                        rows[i].setdefault("characteristics_unparsed", value)
                continue
            for i, value in enumerate(values):
                rows[i][field] = value
    return pd.DataFrame(rows)


def count_header_columns(path: Path, delimiter: str, header_has_gene_col: bool = True) -> int:
    with open_text(path) as fh:
        header = next(csv.reader(fh, delimiter=delimiter))
    return len(header) - 1 if header_has_gene_col else len(header)


def header_values(path: Path, delimiter: str, skip_first: bool = True) -> set[str]:
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader, [])
    values = header[1:] if skip_first else header
    return {clean(value) for value in values}


def split_gse214309_title(title: str) -> dict[str, str]:
    parts = [part.strip() for part in title.split(",")]
    raw_state = parts[0].lower().replace(" ", "_") if parts else ""
    timepoint = parts[1].strip() if len(parts) > 1 else ""
    cell_id = parts[2].strip() if len(parts) > 2 else ""
    maturation = "immature" if raw_state.startswith("immature") else "mature" if raw_state.startswith("mature") else ""
    activity = "active" if raw_state.endswith("active") else "inactive_or_unsorted"
    curated_group = "_".join(part for part in [maturation, activity, timepoint] if part)
    return {
        "cell_id": cell_id,
        "maturation_state": maturation,
        "activity_state": activity,
        "timepoint": timepoint,
        "curated_group": curated_group,
    }


def curate_gse292261() -> pd.DataFrame:
    path = EXTERNAL / "GEO/GSE292261/GSE292261_sample_data_SS2_filtered.csv.gz"
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "cell_id"})
    df["dataset"] = "GSE292261"
    df["species"] = "mouse"
    df["region"] = "dentate_gyrus"
    df["platform"] = "Smart-seq2"
    df["developmental_stage"] = df["Sample"].fillna("").astype(str)
    df["postnatal_day"] = df["developmental_stage"].str.extract(r"P(\d+)", expand=False).fillna("")
    df["curated_group"] = df["developmental_stage"] + "_Leiden" + df["Leiden"].astype(str)
    df["count_matrix_present"] = True
    df["gene_identifier_status"] = "gene symbols already present in count matrix header"
    df["promotion_status"] = "promote_to_primary_dentate_developmental_validation"
    df["source_note"] = "filtered Smart-seq2 metadata; postnatal DG samples with QC metrics and clustering labels"
    keep = [
        "dataset",
        "cell_id",
        "species",
        "region",
        "platform",
        "curated_group",
        "developmental_stage",
        "postnatal_day",
        "Leiden",
        "louvain",
        "nGene",
        "nCounts",
        "percent_mito",
        "count_matrix_present",
        "gene_identifier_status",
        "promotion_status",
        "source_note",
    ]
    return df[keep]


def curate_gse214309() -> pd.DataFrame:
    meta = parse_series_matrix(EXTERNAL / "GEO/GSE214309/GSE214309_series_matrix.txt.gz")
    title_parts = meta["title"].map(split_gse214309_title).apply(pd.Series)
    meta = pd.concat([meta, title_parts], axis=1)
    matrix_cells = set()
    counts = EXTERNAL / "GEO/GSE214309/GSE214309_counts.txt.gz"
    with open_text(counts) as fh:
        header = next(csv.reader(fh, delimiter=","))
        matrix_cells = {clean(cell) for cell in header}
    meta["dataset"] = "GSE214309"
    meta["species"] = "mouse"
    meta["region"] = "dentate_gyrus"
    meta["platform"] = "snRNA-seq"
    meta["count_matrix_present"] = meta["cell_id"].isin(matrix_cells)
    meta["gene_identifier_status"] = "Ensembl gene IDs in rows; marker genes currently resolved by feature-alias lookup; full object analysis needs complete Ensembl-to-symbol map"
    meta["promotion_status"] = "promote_to_primary_dentate_activity_maturation_validation"
    meta["source_note"] = "PROX1+CTIP2+ dentate granule neurons split by mature/immature, activity state, and 1hr/4hr timepoint"
    keep = [
        "dataset",
        "cell_id",
        "geo_accession",
        "title",
        "species",
        "region",
        "platform",
        "curated_group",
        "maturation_state",
        "activity_state",
        "timepoint",
        "cell_state",
        "cell_type",
        "sex",
        "mouseid",
        "count_matrix_present",
        "gene_identifier_status",
        "promotion_status",
        "source_note",
    ]
    return meta[keep]


def curate_gse214905() -> pd.DataFrame:
    meta = parse_series_matrix(EXTERNAL / "GEO/GSE214905/GSE214905_series_matrix.txt.gz")
    count_cells = header_values(EXTERNAL / "GEO/GSE214905/GSE214905_Data-counts.tsv.gz", "\t")
    meta["dataset"] = "GSE214905"
    meta["cell_id"] = meta["title"]
    meta["species"] = "mouse"
    meta["region"] = "dentate_gyrus"
    meta["platform"] = "patch-seq"
    meta["projection_label"] = meta.get("cell_type", "").astype(str).str.replace(" cells", "", regex=False)
    meta["curated_group"] = (
        meta.get("treatment", "").astype(str).str.replace(" ", "_", regex=False)
        + "_"
        + meta["projection_label"].str.replace("+", "plus", regex=False).str.replace("-", "minus", regex=False)
        + "_"
        + meta.get("quality_control", "").astype(str)
    )
    meta["count_matrix_present"] = meta["cell_id"].isin(count_cells)
    meta["gene_identifier_status"] = "gene symbols already present in GENESYMBOL column"
    meta["promotion_status"] = "supporting_targeted_patch_seq_validation_not_primary_discovery"
    meta["source_note"] = "transcription linked to EGFP+/EGFP- projection physiology; small targeted patch-seq design"
    keep = [
        "dataset",
        "cell_id",
        "geo_accession",
        "species",
        "region",
        "platform",
        "curated_group",
        "treatment",
        "projection_label",
        "quality_control",
        "count_matrix_present",
        "gene_identifier_status",
        "promotion_status",
        "source_note",
    ]
    return meta[keep]


def summarize_candidate_calls() -> pd.DataFrame:
    path = RESULTS / "refined_candidate_granule_cell_calls.tsv.gz"
    calls = pd.read_csv(path, sep="\t")
    rows = []
    for dataset, sub in calls.groupby("dataset"):
        counts = sub["candidate_call"].value_counts().to_dict()
        rows.append(
            {
                "dataset": dataset,
                "local_observations": len(sub),
                "dentate_candidate_calls": int(counts.get("candidate_dentate_granule", 0)),
                "cerebellar_candidate_calls": int(counts.get("candidate_cerebellar_granule", 0)),
                "cerebellum_warning_calls": int(counts.get("cerebellum_dentate_panel_high_warning", 0)),
            }
        )
    return pd.DataFrame(rows)


def write_readiness() -> pd.DataFrame:
    call_summary = summarize_candidate_calls()
    manual = pd.DataFrame(
        [
            {
                "dataset": "GSE292261",
                "previous_tier": "priority_secondary",
                "revised_tier": "primary_validation",
                "can_be_primary": "yes",
                "primary_role": "mouse postnatal dentate developmental validation",
                "cleanup_done_now": "curated postnatal day, Leiden/louvain groups, QC metrics, gene-symbol status",
                "remaining_cleanup": "build object-level normalized matrix and align clusters to dentate lineage states",
                "reason": "The limitation was mostly annotation organization, not biology or unusable identifiers.",
            },
            {
                "dataset": "GSE214309",
                "previous_tier": "priority_secondary",
                "revised_tier": "primary_validation",
                "can_be_primary": "yes",
                "primary_role": "mouse adult dentate maturation/activity-state validation",
                "cleanup_done_now": "split maturation, activity, timepoint, sex, mouse ID, and matrix presence",
                "remaining_cleanup": "add full Ensembl-to-symbol mapping for whole-transcriptome object analysis",
                "reason": "The metadata are strong and every local count column is annotated; gene IDs only need systematic symbol harmonization.",
            },
            {
                "dataset": "GSE214905",
                "previous_tier": "supporting",
                "revised_tier": "supporting_validation",
                "can_be_primary": "no",
                "primary_role": "targeted physiology/projection validation",
                "cleanup_done_now": "curated treatment, EGFP projection label, QC status, matrix presence, gene-symbol status",
                "remaining_cleanup": "use passed-QC cells only and treat as targeted validation",
                "reason": "It is secondary because it is small patch-seq with a targeted design, not because annotation/gene symbols are poor.",
            },
            {
                "dataset": "GSE242688",
                "previous_tier": "supporting",
                "revised_tier": "supporting_spatial_context",
                "can_be_primary": "no",
                "primary_role": "cerebellar spatial/proteomics-linked context",
                "cleanup_done_now": "not applicable to cell-level dentate promotion",
                "remaining_cleanup": "keep as spot-level spatial validation with separate statistical treatment",
                "reason": "It is secondary because Visium spots are not single cells and should not be mixed into cell-level discovery statistics.",
            },
        ]
    )
    out = manual.merge(call_summary, on="dataset", how="left")
    out.to_csv(READINESS_OUT, sep="\t", index=False)
    return out


def write_human_candidates() -> pd.DataFrame:
    rows = [
        {
            "dataset": "GSE325391",
            "organism": "human",
            "public_date": "2026-04-20",
            "title": "Unique transcriptional profiles of adult human immature neurons in healthy aging, Alzheimer's disease, and cognitive resilience",
            "region_or_sample": "adult dentate gyrus punch plus fetal hippocampal sections",
            "assay": "single-cell/single-nucleus expression object",
            "local_status": "adult RDS downloaded; GEO metadata curated; RDS inspected; selected sparse bridge built; human-core label projection complete",
            "recommended_tier": "core build-first human dentate primary",
            "use_in_project": "primary adult human DG anchor for mature and differentiating dentate granule comparison",
            "download_note": "adultgc RDS is local; 59,075 nuclei mapped into the human-core selected feature space; fetal RDS remains deferred",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325391",
        },
        {
            "dataset": "GSE268609",
            "organism": "human",
            "public_date": "2025-12-15",
            "title": "A roadmap to human hippocampal neurogenesis in adulthood, aging and AD",
            "region_or_sample": "postmortem human hippocampus across young adult, aging, MCI, and AD",
            "assay": "single-nucleus RNA plus ATAC multiome",
            "local_status": "RNA matrix/barcodes/features downloaded; GEO metadata curated; selected sparse bridge built; human-core label projection complete; ATAC/full Seurat RDS deferred",
            "recommended_tier": "core build-first human hippocampus multiome",
            "use_in_project": "strong human hippocampal neurogenesis reference; useful for dentate/immature-neuron programs and disease sensitivity",
            "download_note": "39 RNA samples and 366,175 nuclei/barcodes are local in a 2,169 selected-gene bridge; projected labels are available, but source cell-type annotations from the full Seurat object remain deferred",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268609",
        },
        {
            "dataset": "GSE186538",
            "organism": "human; macaque; pig",
            "public_date": "2021-11-18",
            "title": "Transcriptomic Taxonomy and Neurogenic Trajectories of Adult Human, Macaque and Pig Hippocampal and Entorhinal Cells",
            "region_or_sample": "human hippocampal-entorhinal subregions including DG samples",
            "assay": "single-nucleus RNA-seq",
            "local_status": "human files downloaded; sparse subset built; QC harmonized; marker validated; GEO metadata curated; normalized reduced object built; labels tuned; dataset-aware module tests complete",
            "recommended_tier": "core build-first human DG taxonomy reference",
            "use_in_project": "important tuned-label DG anchor for human DG identity and cross-species limits of neurogenesis signatures",
            "download_note": "32,067 DG GC sparse subset built, marker-validated, donor-curated, included in the normalized reduced human-core object, and used as the curated DG anchor in module tests",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538",
        },
        {
            "dataset": "GSE185277",
            "organism": "human",
            "public_date": "2022-04-19",
            "title": "Molecular landscape of immature neurons in the human hippocampus across the lifespan",
            "region_or_sample": "postnatal human hippocampus across lifespan",
            "assay": "single-nucleus RNA-seq",
            "local_status": "raw archive downloaded; sparse objects built; QC harmonized; marker validated; GEO metadata curated; normalized reduced object built; labels tuned; dataset-aware module tests complete",
            "recommended_tier": "core build-first human dentate/imGC reference",
            "use_in_project": "useful for human immature dentate granule markers and lifespan framing; now provides tuned high-confidence and immature/neurogenic candidate labels before larger human primary datasets",
            "download_note": "linked to Nature 2022 PMCID PMC9316413; sparse scaffold built, marker-scored, metadata-curated, included in the normalized reduced human-core object, label-tuned, and module-tested",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185277",
        },
        {
            "dataset": "GSE185553",
            "organism": "human",
            "public_date": "2022-04-19",
            "title": "Dissecting the transcriptome landscape of the human hippocampus",
            "region_or_sample": "human hippocampus across lifespan",
            "assay": "single-nucleus RNA-seq",
            "local_status": "raw archive downloaded; sparse objects built; QC harmonized; marker validated; GEO metadata curated; normalized reduced object built; labels tuned; dataset-aware module tests complete",
            "recommended_tier": "core build-first human hippocampus reference",
            "use_in_project": "companion hippocampal context for distinguishing tuned DG/imGC candidates from broader hippocampal neuronal and background states",
            "download_note": "linked to PMID 35794479; sparse companion reference built, marker-scored, metadata-curated, included in the normalized reduced human-core object, label-tuned, and module-tested",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185553",
        },
        {
            "dataset": "GSE198323",
            "organism": "human",
            "public_date": "2022-04-19",
            "title": "Molecular landscape of immature neurons in the human hippocampus in Alzheimer's disease",
            "region_or_sample": "human hippocampus in Alzheimer's disease",
            "assay": "expression profiling by high throughput sequencing",
            "local_status": "source listing saved; not downloaded",
            "recommended_tier": "human disease reference",
            "use_in_project": "useful for AD-related immature-neuron changes; keep separate from healthy/reference analyses",
            "download_note": "linked to Nature 2022 PMCID PMC9316413; disease-context layer",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198323",
        },
        {
            "dataset": "GSE216877",
            "organism": "human",
            "public_date": "2022-11-03",
            "title": "Multi-modal characterization and simulation of human epileptic circuitry",
            "region_or_sample": "human temporal-lobe epilepsy hippocampus/granule-cell context",
            "assay": "multi-modal expression/electrophysiology-oriented dataset",
            "local_status": "not downloaded",
            "recommended_tier": "disease-context validation",
            "use_in_project": "useful after healthy/reference datasets; avoid as the only human primary because epilepsy can reshape activity and morphology programs",
            "download_note": "smaller disease-focused validation candidate",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216877",
        },
        {
            "dataset": "GSE317381",
            "organism": "human",
            "public_date": "2026-01-30",
            "title": "spaTransfer [DG]",
            "region_or_sample": "dentate gyrus spatial/annotation-transfer resource",
            "assay": "spatial/annotation transfer",
            "local_status": "not downloaded",
            "recommended_tier": "human DG spatial context",
            "use_in_project": "supportive spatial context, not primary cell-level discovery",
            "download_note": "treat separately from single-cell/nucleus datasets",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317381",
        },
        {
            "dataset": "spatial_DG_lifespan",
            "organism": "human",
            "public_date": "2025-02-25",
            "title": "Spatiotemporal analysis of gene expression in the human dentate gyrus reveals age-associated changes in cellular maturation and neuroinflammation",
            "region_or_sample": "human dentate gyrus Visium spatial transcriptomics across infant, teen, adult, and elderly donors",
            "assay": "spatial transcriptomics",
            "local_status": "article PDF present in Literature; data not downloaded",
            "recommended_tier": "human DG spatial lifespan context",
            "use_in_project": "use as human spatial/lifespan validation for granule-cell maturation, extracellular matrix, and neuroinflammation modules",
            "download_note": "processed SpatialExperiment objects are on Zenodo; keep separate from single-cell/nucleus primary discovery statistics",
            "source_url": "https://doi.org/10.5281/zenodo.10126687",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(HUMAN_CANDIDATES_OUT, sep="\t", index=False)
    return out


def write_interpretation(readiness: pd.DataFrame, human: pd.DataFrame) -> None:
    promoted = readiness.loc[readiness["can_be_primary"] == "yes", "dataset"].tolist()
    text = f"""# Secondary dataset promotion and human dentate options

## Answer

The secondary label was not one thing. For `GSE292261` and `GSE214309`, it mostly reflected unfinished annotation/gene-identifier curation, so these two can be promoted to primary validation datasets. For `GSE214905` and `GSE242688`, the secondary label is intrinsic to the experimental design: `GSE214905` is small targeted patch-seq, and `GSE242688` is spatial/proteomics-linked spot-level data rather than cell-level single-cell/nucleus data.

Promoted local datasets: {", ".join(promoted)}.

## Local promotion decisions

- `GSE292261`: promote to primary dentate developmental validation. The count matrix already uses gene symbols, and the metadata now records postnatal stage, Leiden/louvain group, QC metrics, and matrix presence.
- `GSE214309`: promote to primary dentate maturation/activity validation. The sample metadata are strong and now split maturation state, activity state, timepoint, sex, and mouse ID. The remaining technical requirement is a full Ensembl-to-symbol map for whole-transcriptome object-level analysis.
- `GSE214905`: keep as supporting targeted validation. Gene symbols and metadata are usable, but the sample size and patch-seq/projection design make it a validation layer rather than a discovery backbone.
- `GSE242688`: keep as supporting context. Spatial spots should not be combined with cell-level observations in primary cell-level statistics.

## Human dentate/hippocampal core construction candidates

- `GSE185277`: first sparse human hippocampal/dentate imGC scaffold is built, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
- `GSE185553`: broader human hippocampal companion sparse reference is built, marker-scored, GEO specimen/age curated, included in the normalized reduced object, label-tuned, and module-tested.
- `GSE186538`: human DG GC sparse subset is built from 32,067 candidate cells, marker-validated, donor metadata curated, and included as the tuned DG anchor in the normalized reduced object.
- `GSE325391`: primary modern adult human dentate acquisition is now downloaded, inspected, converted to a selected sparse bridge, and projected into the tuned human-core label convention.
- `GSE268609`: RNA matrix branch is now downloaded, selected-gene bridged, and projected into human-core labels as the broader human hippocampal aging/AD multiome expansion; defer ATAC/full Seurat object initially.
- `GSE198323`: keep as disease context after the healthy/reference human imGC branch is built.

`GSE216877` and `GSE317381` should be kept as disease/spatial context rather than the first human primary datasets.

## Practical next step

Build primary object-level analyses in two tiers:

1. Human dentate/hippocampal construction: use `GSE325391` as the primary adult human DG anchor alongside the tuned `GSE185277`, `GSE185553`, and `GSE186538` scaffold, with `GSE268609` as the broader aging/AD hippocampal RNA expansion.
2. Main comparative rerun: combine the new human dentate/hippocampal branch with the local mouse dentate and cerebellar backbone (`GSE104323`, `GSE95752`, `GSE122357`, `GSE165657`, `GSE312658`, plus promoted `GSE292261` and `GSE214309`).
"""
    INTERPRETATION_OUT.write_text(text)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    metadata = pd.concat([curate_gse292261(), curate_gse214309(), curate_gse214905()], ignore_index=True, sort=False)
    metadata.to_csv(PROMOTED_METADATA_OUT, sep="\t", index=False)
    readiness = write_readiness()
    human = write_human_candidates()
    write_interpretation(readiness, human)
    print(f"Wrote {PROMOTED_METADATA_OUT}")
    print(f"Wrote {READINESS_OUT}")
    print(f"Wrote {HUMAN_CANDIDATES_OUT}")
    print(f"Wrote {INTERPRETATION_OUT}")


if __name__ == "__main__":
    main()
