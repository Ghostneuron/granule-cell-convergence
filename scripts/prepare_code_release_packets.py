#!/usr/bin/env python3
"""Prepare GitHub-ready and Zenodo-ready release packets.

The release packets intentionally exclude raw downloads, processed sparse
matrices, render-QA files and local caches. They include scripts, configuration,
manuscript-facing figures/tables, compact result summaries, and reproducibility
documentation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import textwrap
import zipfile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "Project"
RELEASE_ROOT = PROJECT / "release"
VERSION = "v0.9-submission-prep"
SLUG = "granule-cell-convergence"
GITHUB_URL = "https://github.com/Ghostneuron/granule-cell-convergence"
ZENODO_DOI = "10.5281/zenodo.21018501"
ZENODO_URL = f"https://doi.org/{ZENODO_DOI}"
ORCID = "https://orcid.org/0000-0001-6843-9720"
GITHUB_DIR = RELEASE_ROOT / "github_packet" / SLUG
ZENODO_DIR = RELEASE_ROOT / "zenodo_packet" / f"{SLUG}_{VERSION}"
ZENODO_ZIP = RELEASE_ROOT / "zenodo_packet" / f"{SLUG}_{VERSION}.zip"


INCLUDE_RESULT_PATTERNS = {
    ".md",
    ".tsv",
    ".tsv.gz",
    ".png",
}

RESULT_EXCLUDE_PREFIXES = (
    "human_core_",
    "gse268609_human_core_label_projection",
    "gse325391_human_core_label_projection",
    "primary_core_expanded_gene_pseudobulk_expression",
    "primary_core_genomewide_symbol_pseudobulk_expression",
    "primary_core_mgi_ortholog_full_matrix_expression",
)


def clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_tree_files(src_dir: Path, dst_dir: Path, include_suffixes: set[str] | None = None) -> None:
    for src in sorted(src_dir.rglob("*")):
        if not src.is_file():
            continue
        if src.name == ".DS_Store" or "__pycache__" in src.parts:
            continue
        if include_suffixes and not any(src.name.endswith(s) for s in include_suffixes):
            continue
        rel = src.relative_to(src_dir)
        copy_file(src, dst_dir / rel)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).lstrip(), encoding="utf-8")


def write_manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if ".git" in path.parts:
            continue
        if path.name == "included_file_manifest.tsv":
            continue
        rows.append(
            {
                "relative_path": str(path.relative_to(root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    out = root / "included_file_manifest.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["relative_path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_excluded_manifest(root: Path) -> None:
    rows = []
    for base in [PROJECT / "processed", PROJECT / "dataset_search_cache", PROJECT / "manuscript" / "render_QA"]:
        if base.exists():
            rows.append(
                {
                    "path": str(base.relative_to(ROOT)),
                    "reason": "Excluded from code release: raw/processed/cache/render-QA material; regenerate or download via manifests.",
                }
            )
    for path in sorted(PROJECT.rglob("*")):
        if not path.is_file():
            continue
        if RELEASE_ROOT in path.parents:
            continue
        if PROJECT / "manuscript" / "Supplementary tables" in path.parents:
            continue
        if PROJECT / "manuscript" / "Final figures" in path.parents:
            continue
        if PROJECT / "manuscript" / "Supplementary figures" in path.parents:
            continue
        if path.stat().st_size >= 10 * 1024 * 1024:
            rows.append(
                {
                    "path": str(path.relative_to(ROOT)),
                    "reason": "Large local file excluded unless represented by supplementary table archive or figure output.",
                }
            )
    out = root / "excluded_large_or_local_files.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["path", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def copy_core_content(root: Path) -> None:
    copy_tree_files(PROJECT / "scripts", root / "scripts", {".py", ".R"})
    copy_tree_files(PROJECT / "config", root / "config")

    docs_dir = root / "docs"
    copy_file(PROJECT / "downloaded_external_data_manifest.tsv", docs_dir / "downloaded_external_data_manifest.tsv")
    copy_file(PROJECT / "external_dataset_inventory.md", docs_dir / "external_dataset_inventory.md")
    copy_file(PROJECT / "manuscript" / "granule_cell_convergence_Development_formatted_main.md", docs_dir / "manuscript_main.md")
    copy_file(PROJECT / "manuscript" / "granule_cell_convergence_Supplementary_information.md", docs_dir / "supplementary_information.md")
    copy_file(PROJECT / "manuscript" / "development_methods_clean.md", docs_dir / "methods_clean.md")
    if (PROJECT / "results" / "manuscript_figure_plan.tsv").exists():
        copy_file(PROJECT / "results" / "manuscript_figure_plan.tsv", docs_dir / "manuscript_figure_plan.tsv")

    copy_tree_files(PROJECT / "manuscript" / "Final figures", root / "manuscript_outputs" / "final_figures")
    copy_tree_files(PROJECT / "manuscript" / "Supplementary figures", root / "manuscript_outputs" / "supplementary_figures")
    copy_tree_files(PROJECT / "manuscript" / "Supplementary tables", root / "manuscript_outputs" / "supplementary_tables")

    results_out = root / "results_summary"
    results_out.mkdir(parents=True, exist_ok=True)
    for src in sorted((PROJECT / "results").iterdir()):
        if not src.is_file() or src.name == ".DS_Store":
            continue
        if not any(src.name.endswith(s) for s in INCLUDE_RESULT_PATTERNS):
            continue
        if src.stat().st_size >= 10 * 1024 * 1024:
            continue
        if src.name.startswith(RESULT_EXCLUDE_PREFIXES):
            continue
        copy_file(src, results_out / src.name)


def write_release_docs(root: Path, packet_kind: str) -> None:
    today = date.today().isoformat()
    write_text(
        root / "README.md",
        f"""
        # Granule-cell convergence analysis code release ({VERSION})

        This packet accompanies the manuscript **"Distinct Dentate and Cerebellar
        Granule-Cell Lineages Converge through Niche and Circuit Constraints"**.

        Packet kind: `{packet_kind}`

        Prepared: {today}

        ## Contents

        - `scripts/`: analysis, curation, figure and manuscript-generation scripts.
        - `config/`: marker panels and small configuration tables.
        - `docs/`: manuscript-facing documentation, data-source inventory, figure plan and clean methods text.
        - `results_summary/`: compact result summaries, plots and machine-readable outputs used to assemble figures/tables.
        - `manuscript_outputs/final_figures/`: final main figure image files.
        - `manuscript_outputs/supplementary_figures/`: final supplementary figure image files.
        - `manuscript_outputs/supplementary_tables/`: ordered supplementary table packet and table archives.
        - `included_file_manifest.tsv`: SHA-256 checksums for included files.
        - `excluded_large_or_local_files.tsv`: local files deliberately excluded from the release.

        ## What is not included

        Raw public datasets, large sparse matrices, DANDI NWB files, GEO H5/RDS
        downloads, local render-QA files and project caches are not redistributed.
        Public accessions and download sources are listed in `docs/downloaded_external_data_manifest.tsv`,
        `docs/external_dataset_inventory.md`, and the supplementary table packet.

        ## Reproducibility note

        The scripts were written as project-level workflows and assume public raw
        data are downloaded or reconstructed according to the dataset manifests.
        For manuscript review, the recommended reproducibility path is to inspect
        the included scripts together with `manuscript_outputs/supplementary_tables/`
        and `results_summary/`, which preserve the analysis products used for the
        manuscript figures and tables.

        ## License status

        A final license has not yet been selected. Before public GitHub/Zenodo
        release, choose a code license (for example MIT/BSD-3-Clause/GPL) and a
        data/documentation license if desired.

        ## Citation

        Use the metadata in `CITATION.cff`. The current public records are:
        GitHub `{GITHUB_URL}` and Zenodo `{ZENODO_URL}`.
        """,
    )

    write_text(
        root / "requirements.txt",
        """
        numpy
        pandas
        scipy
        scikit-learn
        matplotlib
        h5py
        python-docx
        openpyxl
        requests
        """,
    )

    write_text(
        root / ".gitignore",
        """
        .DS_Store
        __pycache__/
        *.pyc
        .ipynb_checkpoints/
        .venv/
        venv/
        data/
        raw_data/
        processed/
        cache/
        *.h5
        *.h5ad
        *.rds
        *.RDS
        *.nwb
        *.npz
        *.h5seurat
        *.loom
        """,
    )

    write_text(
        root / "LICENSE_PENDING.txt",
        """
        No public reuse license has been selected yet.

        Before uploading this packet to GitHub or Zenodo, choose and add a real
        license file. Common choices for analysis code include MIT, BSD-3-Clause,
        Apache-2.0 or GPL-3.0. Data/table reuse can be licensed separately if desired.
        """,
    )

    write_text(
        root / "CITATION.cff",
        f"""
        cff-version: 1.2.0
        title: "Granule-cell convergence analysis code"
        message: "If you use this code or release packet, please cite the manuscript and archived software release."
        type: software
        authors:
          - family-names: Lu
            given-names: Jie
            orcid: "{ORCID}"
        version: "{VERSION}"
        date-released: "{today}"
        repository-code: "{GITHUB_URL}"
        doi: "{ZENODO_DOI}"
        keywords:
          - granule cell
          - dentate gyrus
          - cerebellum
          - single-cell transcriptomics
          - neuronal morphology
          - sparse coding
        """,
    )

    zenodo = {
        "title": "Granule-cell convergence analysis code and curated manuscript outputs",
        "upload_type": "software",
        "description": (
            "Analysis scripts, manuscript-facing summary outputs, figures and supplementary table packet "
            "for a computational study of dentate and cerebellar granule-cell convergence."
        ),
        "creators": [{"name": "Lu, Jie", "orcid": ORCID}],
        "keywords": [
            "granule cell",
            "dentate gyrus",
            "cerebellum",
            "single-cell transcriptomics",
            "neuronal morphology",
            "sparse coding",
        ],
        "version": VERSION,
        "related_identifiers": [
            {
                "identifier": GITHUB_URL,
                "relation": "isSupplementTo",
                "scheme": "url",
            }
        ],
        "notes": "License metadata should be finalized before public reuse. Raw public datasets and large intermediate matrices are not redistributed.",
    }
    (root / ".zenodo.json").write_text(json.dumps(zenodo, indent=2) + "\n", encoding="utf-8")

    write_text(
        root / "docs" / "run_order.md",
        """
        # Suggested Run Order

        This project was assembled as script-level workflows. Exact reruns require
        public raw datasets or local selected matrices reconstructed from GEO,
        Allen-related resources, NeuroMorpho and DANDI as documented in the
        manifests.

        1. Dataset discovery and metadata curation:
           `curate_secondary_and_human_candidates.py`,
           `integrate_candidate_resources.py`,
           `prioritize_datasets_for_next_phase.py`.
        2. Human bridge/core construction and label validation:
           `inspect_human_seed_archives.py`,
           `curate_human_core_geo_sample_metadata.py`,
           `build_human_seed_sparse_objects.py`,
           `qc_harmonize_human_core_sparse_objects.py`,
           `build_human_core_normalized_reduced_object.py`,
           `tune_human_core_labels_and_test_modules.py`,
           `validate_human_core_marker_programs.py`.
        3. Primary-core annotation, pseudobulk and ortholog rank-meta analyses:
           `classify_candidate_granule_cells.py`,
           `build_primary_core_*pseudobulk.py`,
           `build_primary_core_mgi_ortholog_*meta_model.py`,
           `build_primary_core_mgi_ortholog_formal_rank_model.py`.
        4. Candidate tiering, comparator and configuration analyses:
           `build_primary_core_manuscript_candidate_packet.py`,
           `build_primary_core_granule_specificity_named_comparators.py`,
           `build_primary_core_transcriptomic_configuration_model.py`,
           `build_primary_core_transcriptomic_configuration_primary_validation.py`,
           `build_primary_core_configuration_driver_audit.py`.
        5. Stage, pathway, ligand-receptor and conditioned-medium analyses:
           `build_primary_core_aim2_niche_pathway_model.py`,
           `build_primary_core_aim2b_stage_resolved_tgf_bdnf.py`,
           `fit_aim2_stage_pseudotime_model.py`,
           `build_aim2_sender_receiver_ligand_receptor.py`,
           `build_cerebellar_conditioned_medium_secretome_candidates.py`.
        6. External validation and sparse-coding model:
           `build_neuromorpho_granule_morphometry_validation.py`,
           `prioritize_dandi_000003_targeted_downloads.py`,
           `build_dandi_000003_*`,
           `build_primary_core_aim3_sparse_coding_model.py`,
           `calibrate_aim3_empirical_sparse_model.py`.
        7. Regulatory, hierarchical and hypothesis-support summaries:
           `build_epigenomic_extension_targets.py`,
           `build_gse322785_*`,
           `build_hierarchical_integrative_granule_model.py`,
           `build_hypothesis_support_score_matrix.py`,
           `build_fig7_hypothesis_support_figure.py`.
        8. Manuscript table/figure assembly:
           `organize_supplementary_tables.py`,
           `assemble_manuscript_figures.py`,
           `build_development_submission_manuscript.py`,
           `build_development_formatted_submission_package.py`.
        """,
    )

    write_text(
        root / "docs" / "release_checklist.md",
        """
        # Release Checklist Before Public Upload

        - [ ] Replace author affiliation and email placeholders.
        - [ ] Choose a license and replace `LICENSE_PENDING.txt` with `LICENSE`.
        - [ ] Confirm `CITATION.cff` and `.zenodo.json` carry the final GitHub URL, Zenodo DOI and ORCID.
        - [ ] Add license and related manuscript DOI metadata if available.
        - [ ] Confirm no raw data or private files are present.
        - [ ] Confirm all included files are listed in `included_file_manifest.tsv`.
        - [ ] Tag the GitHub repository, then archive the release on Zenodo.
        """,
    )


def zip_dir(src_dir: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(src_dir.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(src_dir.parent))


def build_packet(root: Path, packet_kind: str) -> None:
    clean_dir(root)
    copy_core_content(root)
    write_release_docs(root, packet_kind)
    write_excluded_manifest(root)
    write_manifest(root)


def main() -> int:
    clean_dir(RELEASE_ROOT / "github_packet")
    clean_dir(RELEASE_ROOT / "zenodo_packet")
    build_packet(GITHUB_DIR, "github-ready repository directory")
    build_packet(ZENODO_DIR, "zenodo-ready archive snapshot")
    zip_dir(ZENODO_DIR, ZENODO_ZIP)
    print(f"Wrote GitHub packet: {GITHUB_DIR}")
    print(f"Wrote Zenodo packet directory: {ZENODO_DIR}")
    print(f"Wrote Zenodo zip: {ZENODO_ZIP}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
