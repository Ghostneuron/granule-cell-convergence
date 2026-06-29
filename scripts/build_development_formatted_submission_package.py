#!/usr/bin/env python3
"""Build a Development-formatted main manuscript and Supplementary Information.

This keeps the existing Development draft untouched and produces a formatted
copy with Development-style section order, a shorter title, <=10 key words,
Materials and Methods heading, and supplementary figure legends moved into a
separate Supplementary Information DOCX.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import build_development_submission_manuscript as dev


ROOT = Path(__file__).resolve().parents[2]
SRC_MD = ROOT / "Project" / "manuscript" / "granule_cell_convergence_Development_submission_draft.md"
OUT_MAIN_MD = ROOT / "Project" / "manuscript" / "granule_cell_convergence_Development_formatted_main.md"
OUT_MAIN_DOCX = ROOT / "Project" / "manuscript" / "granule_cell_convergence_Development_formatted_main.docx"
OUT_SUPP_MD = ROOT / "Project" / "manuscript" / "granule_cell_convergence_Supplementary_information.md"
OUT_SUPP_DOCX = ROOT / "Project" / "manuscript" / "granule_cell_convergence_Supplementary_information.docx"

TITLE = "Distinct Dentate and Cerebellar Granule-Cell Lineages Converge through Niche and Circuit Constraints"
KEY_WORDS = (
    "granule cell; dentate gyrus; cerebellum; developmental convergence; "
    "postmitotic assembly; neuronal morphology; single-cell transcriptomics; "
    "TGF-beta; BDNF; sparse coding/pattern separation"
)

AUTHOR_INFORMATION = """Jie Lu

Affiliation: [to be added]

Correspondence: [email to be added]

ORCID: https://orcid.org/0000-0001-6843-9720"""


def section_body(markdown: str, heading: str) -> str:
    marker = f"## {heading}"
    start = markdown.index(marker)
    next_match = re.search(r"^## ", markdown[start + len(marker) :], flags=re.M)
    end = start + len(marker) + next_match.start() if next_match else len(markdown)
    return markdown[start + len(marker) : end].strip()


def build_markdown() -> tuple[str, str]:
    md = SRC_MD.read_text(encoding="utf-8")

    short_title = section_body(md, "Short Title")
    summary = section_body(md, "Summary Statement")
    abstract = section_body(md, "Abstract")
    introduction = section_body(md, "Introduction")
    results = section_body(md, "Results")
    discussion = section_body(md, "Discussion")
    methods = section_body(md, "Methods")
    data_availability = section_body(md, "Data Availability")
    code_availability = section_body(md, "Code Availability")
    author_contributions = section_body(md, "Author Contributions")
    competing_interests = section_body(md, "Competing Interests")
    funding = section_body(md, "Funding")
    references = section_body(md, "References")
    figure_legends = section_body(md, "Figure Legends")

    supp_marker = "### Supplementary Figure S1."
    if supp_marker not in figure_legends:
        raise RuntimeError("Could not find Supplementary Figure S1 marker in figure legends.")
    split = figure_legends.index(supp_marker)
    main_legends = figure_legends[:split].rstrip()
    supp_legends = figure_legends[split:].strip()

    main_parts = [
        f"# {TITLE}",
        "## Author information",
        AUTHOR_INFORMATION,
        "## Short Title",
        short_title,
        "## Summary Statement",
        summary,
        "## Abstract",
        abstract,
        "## Key words",
        KEY_WORDS,
        "## Introduction",
        introduction,
        "## Results",
        results,
        "## Discussion",
        discussion,
        "## Materials and Methods",
        methods,
        "## Acknowledgements",
        "Not applicable.",
        "## Competing interests",
        competing_interests,
        "## Author contributions",
        author_contributions,
        "## Funding",
        funding,
        "## Data availability",
        data_availability,
        "## Code availability",
        code_availability,
        "## References",
        references,
        "## Figure Legends",
        main_legends,
    ]
    main_md = "\n\n".join(part.strip() for part in main_parts if part.strip()) + "\n"

    supp_parts = [
        "# Supplementary Information",
        "## Article",
        TITLE,
        "## Author",
        "Jie Lu",
        "## Supplementary Figure Legends",
        supp_legends,
    ]
    supp_md = "\n\n".join(part.strip() for part in supp_parts if part.strip()) + "\n"
    return main_md, supp_md


def write_docx(markdown_path: Path, docx_path: Path) -> None:
    dev.build_reference_docx()
    subprocess.run(
        [
            "/opt/anaconda3/bin/pandoc",
            str(markdown_path),
            "--from",
            "markdown+tex_math_dollars+tex_math_single_backslash",
            "--to",
            "docx",
            "--reference-doc",
            str(dev.REFERENCE_DOCX),
            "--output",
            str(docx_path),
        ],
        check=True,
    )
    dev.scrub_metadata_and_enable_line_numbers(docx_path)


def main() -> int:
    main_md, supp_md = build_markdown()
    OUT_MAIN_MD.write_text(main_md, encoding="utf-8")
    OUT_SUPP_MD.write_text(supp_md, encoding="utf-8")
    write_docx(OUT_MAIN_MD, OUT_MAIN_DOCX)
    write_docx(OUT_SUPP_MD, OUT_SUPP_DOCX)
    print(f"Wrote main markdown: {OUT_MAIN_MD}")
    print(f"Wrote main DOCX: {OUT_MAIN_DOCX}")
    print(f"Wrote supplementary markdown: {OUT_SUPP_MD}")
    print(f"Wrote supplementary DOCX: {OUT_SUPP_DOCX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
