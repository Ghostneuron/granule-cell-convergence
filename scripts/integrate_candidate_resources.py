#!/usr/bin/env python3
"""Write a unified inventory for candidate resources not yet in the analysis matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"

OUT_TSV = RESULTS / "candidate_resource_inventory.tsv"
OUT_MD = RESULTS / "candidate_resource_integration.md"


FIELDNAMES = [
    "resource_group",
    "component_id",
    "source_type",
    "species",
    "region",
    "modality",
    "sample_scope",
    "local_status",
    "recommended_tier",
    "integration_role",
    "acquisition_note",
    "source_url",
]


ROWS = [
    {
        "resource_group": "human_dentate_core_construction",
        "component_id": "GSE185277",
        "source_type": "GEO",
        "species": "Homo sapiens",
        "region": "hippocampus;dentate_granule_immature_neuron_focus",
        "modality": "snRNA_seq",
        "sample_scope": "5 GEO samples across postnatal lifespan",
        "local_status": "raw_downloaded;sparse_objects_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests",
        "recommended_tier": "core_build_first_human_dentate_reference",
        "integration_role": "human imGC marker scaffold and lifespan neurogenesis framing",
        "acquisition_note": "7 sparse libraries built, QC-harmonized, marker-scored, GEO metadata curated, included in the normalized reduced human-core object, label-tuned, and module-tested",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185277",
    },
    {
        "resource_group": "human_dentate_core_construction",
        "component_id": "GSE185553",
        "source_type": "GEO",
        "species": "Homo sapiens",
        "region": "hippocampus",
        "modality": "snRNA_seq",
        "sample_scope": "5 GEO samples, broader human hippocampus across lifespan",
        "local_status": "raw_downloaded;sparse_objects_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests",
        "recommended_tier": "core_build_first_human_hippocampus_reference",
        "integration_role": "broader human hippocampal context for GSE185277",
        "acquisition_note": "27 sparse libraries built, QC-harmonized, marker-scored, GEO metadata curated, included in the normalized reduced human-core object, label-tuned, and module-tested",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185553",
    },
    {
        "resource_group": "human_dentate_core_construction",
        "component_id": "GSE186538",
        "source_type": "GEO",
        "species": "Homo sapiens;Macaca mulatta;Sus scrofa",
        "region": "dentate_gyrus;hippocampal_entorhinal",
        "modality": "snRNA_seq",
        "sample_scope": "human, macaque, and pig hippocampal-entorhinal nuclei; prioritize human DG-labelled material",
        "local_status": "human_files_downloaded;sparse_subset_built;qc_harmonized;marker_validated;geo_metadata_curated;normalized_reduced_object_built;labels_tuned;dataset_aware_module_tests",
        "recommended_tier": "core_build_first_human_dg_taxonomy",
        "integration_role": "human DG taxonomy and cross-species reference",
        "acquisition_note": "32,067 DG GC sparse subset built, QC-harmonized, marker-validated, donor-curated, included in the normalized reduced human-core object, and used as the tuned-label DG anchor",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538",
    },
    {
        "resource_group": "human_dentate_core_construction",
        "component_id": "GSE325391",
        "source_type": "GEO",
        "species": "Homo sapiens",
        "region": "dentate_gyrus",
        "modality": "single_nucleus_or_single_cell_expression",
        "sample_scope": "adult human dentate gyrus punch and fetal hippocampal material",
        "local_status": "adult_rds_downloaded;geo_metadata_curated;rds_inspected;selected_sparse_bridge_built;human_core_label_projected",
        "recommended_tier": "core_build_first_human_dentate_primary",
        "integration_role": "primary modern adult human dentate comparator",
        "acquisition_note": "adultgc RDS downloaded and inspected; 59,075 adult DG nuclei mapped into the human-core selected feature space and label convention; fetal RDS remains deferred",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325391",
    },
    {
        "resource_group": "human_dentate_core_construction",
        "component_id": "GSE268609",
        "source_type": "GEO",
        "species": "Homo sapiens",
        "region": "hippocampus",
        "modality": "snRNA_seq;snATAC_seq_multiome",
        "sample_scope": "human hippocampal neurogenesis across aging, MCI, and AD",
        "local_status": "source_listing_saved;geo_metadata_curated;rna_matrix_barcodes_features_downloaded;rna_selected_npz_built;human_core_label_projected;atac_and_full_rds_deferred",
        "recommended_tier": "core_build_first_human_hippocampus_multiome",
        "integration_role": "broader human hippocampal RNA/multiome expansion",
        "acquisition_note": "RNA matrix/barcodes/features downloaded; 39 RNA samples and 366,175 nuclei/barcodes streamed into a 2,169 selected-gene human-core object; ATAC fragments and full Seurat object remain deferred",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268609",
    },
    {
        "resource_group": "nature_neuroscience_2025_mouse_dg_aging",
        "component_id": "GSE233363",
        "source_type": "GEO",
        "species": "Mus musculus",
        "region": "dentate_gyrus",
        "modality": "scRNA_seq;spatial_transcriptomics",
        "sample_scope": "20 GEO samples: 3 scRNA-seq samples plus young/middle-age/old spatial transcriptomics samples",
        "local_status": "source_listing_saved;not_downloaded",
        "recommended_tier": "core_validation_after_human_dentate",
        "integration_role": "mouse DG aging, neurogenic-lineage maturation, niche inflammation, and spatial validation",
        "acquisition_note": "download after human dentate/hippocampal branch is constructed",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233363",
    },
    {
        "resource_group": "cell_reports_2025_human_dg_spatial_lifespan",
        "component_id": "spatial_DG_lifespan",
        "source_type": "Zenodo;local_literature_pdf",
        "species": "Homo sapiens",
        "region": "dentate_gyrus",
        "modality": "Visium_spatial_transcriptomics",
        "sample_scope": "16 analyzed human donors across infant, teen, adult, and elderly age groups",
        "local_status": "article_pdf_present;zenodo_metadata_saved;data_not_downloaded",
        "recommended_tier": "supporting_spatial_context",
        "integration_role": "human DG spatial/lifespan validation for granule-cell maturation, ECM, and neuroinflammation modules",
        "acquisition_note": "processed/raw SpatialExperiment archive is about 38 GB; use only after deciding to add a spatial layer",
        "source_url": "https://doi.org/10.5281/zenodo.10126687",
    },
    {
        "resource_group": "nature_2022_human_hippocampal_immature_neurons",
        "component_id": "GSE198323",
        "source_type": "GEO",
        "species": "Homo sapiens",
        "region": "hippocampus",
        "modality": "high_throughput_expression",
        "sample_scope": "6 GEO samples, Alzheimer's disease context",
        "local_status": "source_listing_saved;not_downloaded",
        "recommended_tier": "supporting_disease_reference",
        "integration_role": "AD-related human immature-neuron disease context; keep separate from healthy/reference analyses",
        "acquisition_note": "download only after healthy/reference human imGC analysis is defined",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198323",
    },
]


def write_tsv() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(ROWS)


def write_md() -> None:
    groups = {}
    for row in ROWS:
        groups.setdefault(row["resource_group"], []).append(row)

    lines = [
        "# Candidate resource integration",
        "",
        "Date checked: 2026-06-21",
        "",
        "## Decision",
        "",
        "The newly inspected resources are useful, and the human dentate/hippocampal core expansion now includes a direct adult dentate anchor plus a broader human hippocampal aging/AD RNA branch from `GSE268609`.",
        "",
        "## Resource groups",
        "",
    ]
    for group, rows in groups.items():
        lines.append(f"### {group}")
        for row in rows:
            lines.append(
                f"- `{row['component_id']}`: {row['integration_role']} "
                f"({row['recommended_tier']}; {row['local_status']})."
            )
        lines.append("")

    lines.extend(
        [
            "## Practical use",
            "",
            "1. Use the tuned human-core labels and dataset-aware module tests as the convention for the human dentate/hippocampal branch.",
            "2. Treat the `GSE186538` DG GC subset as the marker-validated human DG taxonomy anchor.",
            "3. Treat `GSE325391` as the primary modern adult human dentate anchor and `GSE268609` as the broader human hippocampal aging/AD RNA expansion with projected, not source-taxonomy, labels.",
            "4. Next, integrate these human dentate/hippocampal layers with the mouse dentate and cerebellar backbone, while keeping `GSE233363`, `spatial_DG_lifespan`, and `GSE198323` as validation/context layers.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def main() -> None:
    write_tsv()
    write_md()
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
