#!/usr/bin/env python3
"""Build the secondary physiology/morphology validation layer.

This layer answers a narrow follow-up question: which public resources can
help parameterize or validate the Aim 3 sparse-expansion model beyond the
primary transcriptomic core?

It intentionally keeps Allen Cell Types as calibration/comparator evidence,
because API probes and Allen's own documentation indicate that the Cell Types
Database is focused on cortex/thalamus rather than dentate or cerebellar
granule cells.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project" / "results"
CACHE = ROOT / "Project" / "dataset_search_cache" / "phys_morph"

OUT_RESOURCE_TSV = RESULTS / "secondary_phys_morph_candidate_resources.tsv"
OUT_PARAMETER_TSV = RESULTS / "secondary_phys_morph_parameter_map.tsv"
OUT_ALLEN_TSV = RESULTS / "secondary_phys_morph_allen_celltypes_probe.tsv"
OUT_NEUROMORPHO_TSV = RESULTS / "secondary_phys_morph_neuromorpho_probe.tsv"
OUT_MD = RESULTS / "secondary_phys_morph_validation_plan.md"


def load_json(name: str) -> dict[str, Any]:
    path = CACHE / name
    if not path.exists() or path.stat().st_size == 0:
        return {}
    return json.loads(path.read_text())


def total_rows(name: str) -> int | None:
    data = load_json(name)
    value = data.get("total_rows")
    return int(value) if value is not None else None


def total_elements(name: str) -> int | None:
    data = load_json(name)
    if data.get("status") == 404:
        return 0
    page = data.get("page") or {}
    value = page.get("totalElements")
    return int(value) if value is not None else None


def first_neuron(name: str) -> dict[str, Any]:
    data = load_json(name)
    embedded = data.get("_embedded") or {}
    resources = embedded.get("neuronResources") or []
    return resources[0] if resources else {}


def morphometry(name: str) -> dict[str, Any]:
    return load_json(name)


def dandi_summary(name: str) -> dict[str, Any]:
    data = load_json(name)
    assets = data.get("assetsSummary") or {}
    return {
        "name": data.get("name", ""),
        "doi": data.get("doi", ""),
        "url": data.get("url", ""),
        "files": assets.get("numberOfFiles"),
        "subjects": assets.get("numberOfSubjects"),
        "bytes": assets.get("numberOfBytes"),
        "variables": ",".join(assets.get("variableMeasured") or []),
        "description": data.get("description", ""),
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_allen_probe() -> list[dict[str, Any]]:
    return [
        {
            "probe": "all_celltypes_api_rows",
            "query_scope": "ApiCellTypesSpecimenDetail; one-row fetch with total_rows",
            "matched_rows": total_rows("allen_celltypes_all_num1.json"),
            "interpretation": "Allen Cell Types has a large ephys/morphology/model table suitable for comparator calibration.",
            "cache_file": "Project/dataset_search_cache/phys_morph/allen_celltypes_all_num1.json",
        },
        {
            "probe": "dentate_structure_name",
            "query_scope": "structure__name contains dentate",
            "matched_rows": total_rows("allen_celltypes_dentate_probe.json"),
            "interpretation": "No direct dentate Cell Types rows in this API probe.",
            "cache_file": "Project/dataset_search_cache/phys_morph/allen_celltypes_dentate_probe.json",
        },
        {
            "probe": "cerebell_structure_name",
            "query_scope": "structure__name contains cerebell",
            "matched_rows": total_rows("allen_celltypes_cerebell_probe.json"),
            "interpretation": "No direct cerebellar Cell Types rows in this API probe.",
            "cache_file": "Project/dataset_search_cache/phys_morph/allen_celltypes_cerebell_probe.json",
        },
        {
            "probe": "DG_structure_acronym",
            "query_scope": "structure__acronym contains DG",
            "matched_rows": total_rows("allen_celltypes_dg_acronym_probe.json"),
            "interpretation": "No dentate-gyrus acronym matches in this API probe.",
            "cache_file": "Project/dataset_search_cache/phys_morph/allen_celltypes_dg_acronym_probe.json",
        },
    ]


def build_neuromorpho_probe() -> list[dict[str, Any]]:
    dg = first_neuron("neuromorpho_dg_granule_count.json")
    cb = first_neuron("neuromorpho_cb_granule_count.json")
    dg_m = morphometry("neuromorpho_dg_102085_morphometry.json")
    cb_m = morphometry("neuromorpho_cb_32462_morphometry.json")

    return [
        {
            "probe": "dentate_gyrus_granule_all_species",
            "matched_neurons": total_elements("neuromorpho_dg_granule_count.json"),
            "mouse_matches": total_elements("neuromorpho_dg_mouse_count.json"),
            "rat_matches": total_elements("neuromorpho_dg_rat_count.json"),
            "human_matches": total_elements("neuromorpho_dg_human_count.json"),
            "example_neuron_id": dg.get("neuron_id", ""),
            "example_species": dg.get("species", ""),
            "example_domain": dg.get("domain", ""),
            "example_n_stems": dg_m.get("n_stems", ""),
            "example_n_bifs": dg_m.get("n_bifs", ""),
            "example_n_branch": dg_m.get("n_branch", ""),
            "example_length": dg_m.get("length", ""),
            "interpretation": "Strong public morphology source for dentate granule dendrite and branch metrics.",
        },
        {
            "probe": "cerebellum_granule_all_species",
            "matched_neurons": total_elements("neuromorpho_cb_granule_count.json"),
            "mouse_matches": total_elements("neuromorpho_cb_mouse_count.json"),
            "rat_matches": total_elements("neuromorpho_cb_rat_count.json"),
            "human_matches": total_elements("neuromorpho_cb_human_count.json"),
            "example_neuron_id": cb.get("neuron_id", ""),
            "example_species": cb.get("species", ""),
            "example_domain": cb.get("domain", ""),
            "example_n_stems": cb_m.get("n_stems", ""),
            "example_n_bifs": cb_m.get("n_bifs", ""),
            "example_n_branch": cb_m.get("n_branch", ""),
            "example_length": cb_m.get("length", ""),
            "interpretation": "Useful but smaller cerebellar granule morphology source; strict mouse/rat filters are sparse or absent in current query.",
        },
    ]


def build_resource_rows() -> list[dict[str, Any]]:
    d000003 = dandi_summary("dandi_000003_version.json")
    d000165 = dandi_summary("dandi_000165_version.json")
    dg_n = total_elements("neuromorpho_dg_granule_count.json")
    cb_n = total_elements("neuromorpho_cb_granule_count.json")
    allen_n = total_rows("allen_celltypes_all_num1.json")
    allen_dg = total_rows("allen_celltypes_dentate_probe.json")
    allen_cb = total_rows("allen_celltypes_cerebell_probe.json")

    return [
        {
            "resource_id": "Allen_Cell_Types_ephys",
            "source_type": "Allen Institute Cell Types Database/API/SDK",
            "species": "Homo sapiens;Mus musculus",
            "region": "cortex;thalamus",
            "modality": "whole_cell_current_clamp_electrophysiology;computed_ephys_features;NWB",
            "evidence_tier": "secondary_comparator_calibration",
            "usable_parameters": "intrinsic_excitability;firing_rate;input_resistance;tau;rheobase;FI_slope;spike_waveform",
            "local_status": f"API_probe_cached;all_rows={allen_n};dentate_rows={allen_dg};cerebellar_rows={allen_cb}",
            "recommended_use": "Use to calibrate transcriptomic excitability modules against measured neuronal intrinsic physiology and to benchmark non-granule comparator physiology.",
            "main_caveat": "Not a direct dentate or cerebellar granule-cell primary dataset in the current API probe.",
            "source_url": "https://celltypes.brain-map.org/",
        },
        {
            "resource_id": "Allen_Cell_Types_morphology_models",
            "source_type": "Allen Institute Cell Types Database/API/SDK",
            "species": "Homo sapiens;Mus musculus",
            "region": "cortex;thalamus",
            "modality": "morphology_reconstruction;SWC;GLIF_and_biophysical_models",
            "evidence_tier": "secondary_comparator_calibration",
            "usable_parameters": "morphology_feature_scaling;single_cell_model_calibration;intrinsic_response_simulation",
            "local_status": "documentation_verified;not_downloaded",
            "recommended_use": "Use only for comparator-cell feature ranges and model-fitting conventions, not for granule-cell morphology claims.",
            "main_caveat": "Allen's public overview states data generation has focused on selected cortex and thalamic neurons.",
            "source_url": "https://allensdk.readthedocs.io/en/latest/cell_types.html",
        },
        {
            "resource_id": "DANDI_000003_DG_granule_mossy_activity",
            "source_type": "DANDI Archive",
            "species": "Mus musculus",
            "region": "hippocampus;dentate_gyrus_context",
            "modality": "extracellular_electrophysiology;behavior;NWB;Units;LFP;Position",
            "evidence_tier": "primary_activity_validation_for_DG",
            "usable_parameters": "firing_sparsity;active_fraction;place_field_activity;spatial_information;pattern_separation_behavior_proxy",
            "local_status": f"metadata_cached;files={d000003['files']};subjects={d000003['subjects']};variables={d000003['variables']}",
            "recommended_use": "Use as the main public dentate activity dataset for deriving DG granule-cell firing sparsity and behavior-linked coding metrics.",
            "main_caveat": "Large NWB archive; cell identity and unit selection must follow the source paper/metadata carefully.",
            "source_url": d000003["url"],
        },
        {
            "resource_id": "DANDI_000165_DG_CA3_network_units",
            "source_type": "DANDI Archive",
            "species": "Mus musculus",
            "region": "DG;CA3;CA1",
            "modality": "extracellular_electrophysiology;behavior;NWB;Units;LFP;Position",
            "evidence_tier": "supporting_network_activity_validation",
            "usable_parameters": "network_state;DG_CA3_drive;LFP_context;unit_activity_context",
            "local_status": f"metadata_cached;files={d000165['files']};subjects={d000165['subjects']};variables={d000165['variables']}",
            "recommended_use": "Use as supporting hippocampal network physiology, especially for DG/CA3 state context and interneuron modulation.",
            "main_caveat": "Not a pure dentate granule-cell morphology/ephys primary dataset.",
            "source_url": d000165["url"],
        },
        {
            "resource_id": "NeuroMorpho_DG_granule_morphometry",
            "source_type": "NeuroMorpho.Org API",
            "species": "Homo sapiens;Mus musculus;Rattus norvegicus;other",
            "region": "dentate_gyrus",
            "modality": "neuron_reconstruction;morphometry;SWC_metadata",
            "evidence_tier": "primary_morphology_validation_for_DG",
            "usable_parameters": "dendrite_stem_count;branch_count;bifurcations;dendritic_length;surface;volume",
            "local_status": f"API_probe_cached;matched_neurons={dg_n}",
            "recommended_use": "Use to estimate dentate granule dendritic stem/branch distributions and validate morphology claims against public reconstructions.",
            "main_caveat": "Metadata heterogeneity, disease/genotype conditions, and reconstruction completeness need filtering before quantitative fitting.",
            "source_url": "https://neuromorpho.org/apiReference.html",
        },
        {
            "resource_id": "NeuroMorpho_cerebellar_granule_morphometry",
            "source_type": "NeuroMorpho.Org API",
            "species": "Homo sapiens;other;limited_mouse_rat_under_strict_query",
            "region": "cerebellum;cerebellar_cortex;granular_layer",
            "modality": "neuron_reconstruction;morphometry;SWC_metadata",
            "evidence_tier": "supporting_morphology_validation_for_cerebellum",
            "usable_parameters": "dendrite_stem_count;branch_count;bifurcations;dendritic_length;surface;volume",
            "local_status": f"API_probe_cached;matched_neurons={cb_n}",
            "recommended_use": "Use for cerebellar granule morphology validation, but treat as a small and species-heterogeneous evidence set.",
            "main_caveat": "The current strict query returns far fewer cerebellar than dentate granule morphologies.",
            "source_url": "https://neuromorpho.org/apiReference.html",
        },
        {
            "resource_id": "GSE214905_patch_seq_DG",
            "source_type": "GEO",
            "species": "Mus musculus",
            "region": "dentate_gyrus",
            "modality": "patch_seq;targeted_expression;physiology_linked",
            "evidence_tier": "supporting_expression_physiology_bridge",
            "usable_parameters": "intrinsic_excitability_module_bridge;DG_maturation_activity_context",
            "local_status": "downloaded;already_listed_as_supporting_validation",
            "recommended_use": "Use for targeted checks linking DG expression state to physiology, not broad discovery.",
            "main_caveat": "Small sample and targeted patch-seq design.",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214905",
        },
        {
            "resource_id": "cerebellar_granule_synaptic_input_literature",
            "source_type": "curated_literature",
            "species": "multiple",
            "region": "cerebellum",
            "modality": "reviewed_anatomy;synaptic_input_count",
            "evidence_tier": "manual_parameter_prior",
            "usable_parameters": "input_degree;mossy_fiber_rosette_convergence;dendritic_claw_count",
            "local_status": "not_yet_curated",
            "recommended_use": "Use curated numerical priors for cerebellar input-degree ranges until a direct public raw connectomics/ephys source is verified.",
            "main_caveat": "Should be cited as literature-derived prior, not public dataset-derived measurement.",
            "source_url": "",
        },
    ]


def build_parameter_rows() -> list[dict[str, Any]]:
    return [
        {
            "model_parameter": "expansion_ratio",
            "biological_quantity": "relative number of output granule cells versus input channels/projection neurons",
            "best_public_resource": "curated_anatomical_literature;Allen_WMB_taxonomy_for_cell_labels",
            "extraction_plan": "Use literature/atlas estimates as priors; do not infer directly from scRNA-seq abundance without sampling correction.",
            "current_evidence_grade": "partial_prior_only",
            "project_use": "Set plausible ranges for the sparse-expansion simulation.",
            "caveat": "Single-cell datasets distort absolute cell abundance and cannot by themselves estimate expansion ratio.",
        },
        {
            "model_parameter": "input_degree",
            "biological_quantity": "number of dendritic stems/claws and effective synaptic input streams per granule cell",
            "best_public_resource": "NeuroMorpho_DG_granule_morphometry;NeuroMorpho_cerebellar_granule_morphometry;curated_cerebellar_synapse_literature",
            "extraction_plan": "Filter reconstructions by region, cell type, species, condition, and completeness; summarize n_stems, n_branch, and length; combine with literature synaptic counts.",
            "current_evidence_grade": "direct_morphology_partial_synapse",
            "project_use": "Replace the arbitrary Aim 3 input-degree grid with empirically plausible DG/CB ranges.",
            "caveat": "Morphological stems are a proxy; true synaptic input count needs EM/connectomics or focused literature.",
        },
        {
            "model_parameter": "output_active_fraction",
            "biological_quantity": "fraction of units active during behavior or stimulus epochs",
            "best_public_resource": "DANDI_000003_DG_granule_mossy_activity;Allen_Cell_Types_ephys_as_calibration",
            "extraction_plan": "From DANDI units, compute active fraction across behavior bins, firing-rate sparsity, and place-field sparsity; use Allen only for intrinsic excitability feature calibration.",
            "current_evidence_grade": "direct_DG_activity_comparator_Allen",
            "project_use": "Empirically constrain the output sparsity term in the Aim 3 model.",
            "caveat": "Cerebellar granule-cell activity remains a public-data gap in the current resource layer.",
        },
        {
            "model_parameter": "pattern_separation_behavior",
            "biological_quantity": "separation of neural representations or behaviorally discriminated similar contexts/positions",
            "best_public_resource": "DANDI_000003_DG_granule_mossy_activity",
            "extraction_plan": "Use spike-position data to compute population-vector decorrelation, spatial information, pairwise representational distance, and within/between trajectory separation.",
            "current_evidence_grade": "derivable_for_DG",
            "project_use": "Test whether fitted sparse parameters improve separation of similar behavioral states.",
            "caveat": "This is a derived behavior-linked proxy, not a direct standardized pattern-separation behavioral assay.",
        },
        {
            "model_parameter": "intrinsic_excitability",
            "biological_quantity": "membrane response and firing properties",
            "best_public_resource": "Allen_Cell_Types_ephys;GSE214905_patch_seq_DG",
            "extraction_plan": "Use Allen features for general neuronal calibration and patch-seq for DG targeted expression-physiology bridge.",
            "current_evidence_grade": "strong_comparator_weak_DG_direct",
            "project_use": "Connect synaptic/excitability transcriptomic modules to measured electrophysiology feature families.",
            "caveat": "Allen intrinsic ephys is not direct DG/CB granule-cell evidence.",
        },
        {
            "model_parameter": "morphology_similarity",
            "biological_quantity": "stem count, branch count, dendritic length, compactness, and branch order",
            "best_public_resource": "NeuroMorpho.Org API",
            "extraction_plan": "Build matched DG and CB granule morphometry tables, harmonize species/age/condition filters, then compare distributions and effect sizes.",
            "current_evidence_grade": "direct_but_unbalanced",
            "project_use": "Quantify the morphology part of the central hypothesis instead of relying on qualitative similarity.",
            "caveat": "DG data are abundant; CB data are much smaller and need careful species filtering.",
        },
        {
            "model_parameter": "secreted_stop_or_maturation_inputs",
            "biological_quantity": "extrinsic ligand environment affecting proliferation and maturation",
            "best_public_resource": "existing_project_secretome_screen;GSE242688_spatial_proteomics;future_validation_assays",
            "extraction_plan": "Use current secretome rankings as RNA-level candidates; test protein support with proteomics or targeted assays.",
            "current_evidence_grade": "RNA_candidate_not_bioactivity",
            "project_use": "Keep Aim 2 mechanistically connected to the 2005 conditioned-medium question.",
            "caveat": "Secreted protein abundance and activity cannot be proven from transcriptomics alone.",
        },
    ]


def write_markdown(
    resource_rows: list[dict[str, Any]],
    allen_rows: list[dict[str, Any]],
    morph_rows: list[dict[str, Any]],
    parameter_rows: list[dict[str, Any]],
) -> None:
    d000003 = dandi_summary("dandi_000003_version.json")
    lines = [
        "# Secondary Physiology And Morphology Validation Plan",
        "",
        "Date built: 2026-06-23",
        "",
        "## Decision",
        "",
        "Yes, public physiology and morphology data can strengthen Aim 3, but they should be added as a secondary validation layer rather than merged into the strict 10-dataset transcriptomic core. The clean division is:",
        "",
        "- Use `NeuroMorpho.Org` for direct dendrite/stem/branch morphometry.",
        "- Use `DANDI:000003` as the main public dentate granule/mossy-cell activity and behavior-linked physiology source.",
        "- Use Allen Cell Types for intrinsic electrophysiology calibration and comparator-cell feature ranges, not as primary dentate/cerebellar granule evidence.",
        "- Use curated literature, and possibly later connectomics, for synaptic input counts where public raw datasets remain incomplete.",
        "",
        "## Why Allen Helps But Is Not A Primary Granule Dataset",
        "",
        "Allen Cell Types is highly useful because it provides whole-cell current-clamp recordings, computed ephys features, morphology reconstructions, NWB files, and neuronal models. However, the public overview states that data generation has focused on selected cortical and thalamic neurons, and the local API probes found no direct dentate or cerebellar structure matches.",
        "",
        "| Probe | Matched rows | Interpretation |",
        "|---|---:|---|",
    ]
    for row in allen_rows:
        lines.append(
            f"| `{row['probe']}` | {row['matched_rows']} | {row['interpretation']} |"
        )

    lines.extend(
        [
            "",
            "Therefore, Allen should be framed as a comparator/calibration resource for `intrinsic_excitability`, not evidence that cerebellar and dentate granule cells have matched physiology.",
            "",
            "## Direct Morphology Evidence",
            "",
            "| Resource | Matched neurons | Human | Mouse | Rat | Example stems | Example branches | Interpretation |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in morph_rows:
        lines.append(
            "| "
            f"`{row['probe']}` | {row['matched_neurons']} | {row['human_matches']} | "
            f"{row['mouse_matches']} | {row['rat_matches']} | {row['example_n_stems']} | "
            f"{row['example_n_branch']} | {row['interpretation']} |"
        )

    lines.extend(
        [
            "",
            "These values are only probes, not final filtered estimates. The next quantitative step is to download matched NeuroMorpho metadata/morphometry tables, exclude disease/genotype or incomplete reconstructions when appropriate, and compare distributions of `n_stems`, `n_branch`, `n_bifs`, length, and branch order.",
            "",
            "## Direct Activity Evidence",
            "",
            f"`DANDI:000003` is the strongest current public activity source for the dentate branch: {d000003['files']} NWB files, {d000003['subjects']} subjects, and measured variables `{d000003['variables']}`. It can support firing sparsity, active fraction, spatial coding, and behavior-linked population-separation analyses.",
            "",
            "`DANDI:000165` is useful as supporting DG/CA3/CA1 network physiology, especially for LFP/unit state context, but it is not a pure dentate granule-cell dataset.",
            "",
            "## Parameter Map For The Sparse-Coding Model",
            "",
            "| Model parameter | Best public resource | Evidence grade | Use in model |",
            "|---|---|---|---|",
        ]
    )
    for row in parameter_rows:
        lines.append(
            f"| `{row['model_parameter']}` | {row['best_public_resource']} | "
            f"{row['current_evidence_grade']} | {row['project_use']} |"
        )

    lines.extend(
        [
            "",
            "## Recommended Next Analysis",
            "",
            "1. Build a filtered NeuroMorpho morphometry table for DG and cerebellar granule cells.",
            "2. Refit the Aim 3 model using empirical `input_degree` priors from stem/branch/claw literature rather than only arbitrary grid points.",
            "3. Download a small DANDI:000003 subset first, verify NWB structure and unit identity conventions, then compute firing sparsity and behavior-linked separation metrics.",
            "4. Add Allen Cell Types ephys as an intrinsic-excitability calibration table for comparator neurons and for mapping transcriptomic excitability modules to measured ephys feature families.",
            "5. Keep cerebellar granule-cell activity and synaptic input count as explicit evidence gaps unless a direct public cerebellar granule ephys/connectomics dataset is verified.",
            "",
            "## Outputs",
            "",
            f"- Resource table: `{OUT_RESOURCE_TSV.relative_to(ROOT)}`",
            f"- Parameter map: `{OUT_PARAMETER_TSV.relative_to(ROOT)}`",
            f"- Allen API probe: `{OUT_ALLEN_TSV.relative_to(ROOT)}`",
            f"- NeuroMorpho API probe: `{OUT_NEUROMORPHO_TSV.relative_to(ROOT)}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    allen_rows = build_allen_probe()
    morph_rows = build_neuromorpho_probe()
    resource_rows = build_resource_rows()
    parameter_rows = build_parameter_rows()

    write_tsv(
        OUT_ALLEN_TSV,
        allen_rows,
        ["probe", "query_scope", "matched_rows", "interpretation", "cache_file"],
    )
    write_tsv(
        OUT_NEUROMORPHO_TSV,
        morph_rows,
        [
            "probe",
            "matched_neurons",
            "mouse_matches",
            "rat_matches",
            "human_matches",
            "example_neuron_id",
            "example_species",
            "example_domain",
            "example_n_stems",
            "example_n_bifs",
            "example_n_branch",
            "example_length",
            "interpretation",
        ],
    )
    write_tsv(
        OUT_RESOURCE_TSV,
        resource_rows,
        [
            "resource_id",
            "source_type",
            "species",
            "region",
            "modality",
            "evidence_tier",
            "usable_parameters",
            "local_status",
            "recommended_use",
            "main_caveat",
            "source_url",
        ],
    )
    write_tsv(
        OUT_PARAMETER_TSV,
        parameter_rows,
        [
            "model_parameter",
            "biological_quantity",
            "best_public_resource",
            "extraction_plan",
            "current_evidence_grade",
            "project_use",
            "caveat",
        ],
    )
    write_markdown(resource_rows, allen_rows, morph_rows, parameter_rows)

    for path in [
        OUT_RESOURCE_TSV,
        OUT_PARAMETER_TSV,
        OUT_ALLEN_TSV,
        OUT_NEUROMORPHO_TSV,
        OUT_MD,
    ]:
        print(f"Wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
