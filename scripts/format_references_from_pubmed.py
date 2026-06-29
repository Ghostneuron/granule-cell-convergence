#!/usr/bin/env python3
"""Build formatted reference lists from PMID anchors in the manuscript.

The script extracts PMIDs in first-appearance order, fetches PubMed XML
metadata, and writes a master table plus several common reference styles.
"""

from __future__ import annotations

import html
import re
import sys
import textwrap
import unicodedata
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = ROOT / "Project" / "manuscript" / "granule_cell_convergence_manuscript_draft.md"
OUTDIR = ROOT / "Project" / "manuscript" / "References"
XML_CACHE = OUTDIR / "pubmed_efetch.xml"


UNICODE_REPLACEMENTS = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "κ": "kappa",
    "–": "-",
    "—": "-",
    "−": "-",
    "“": '"',
    "”": '"',
    "‘": "'",
    "’": "'",
}


def clean_text(value: str | None) -> str:
    if value is None:
        return ""
    value = html.unescape(value)
    for source, target in UNICODE_REPLACEMENTS.items():
        value = value.replace(source, target)
    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii")
    value = re.sub(r"\s+", " ", value).strip()
    return value


def node_text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return clean_text("".join(node.itertext()))


def find_text(parent: ET.Element | None, path: str) -> str:
    if parent is None:
        return ""
    return clean_text(parent.findtext(path))


def extract_pmids(markdown: str) -> list[str]:
    seen: set[str] = set()
    pmids: list[str] = []
    for match in re.finditer(r"PMID:\s*([0-9]+)", markdown):
        pmid = match.group(1)
        if pmid not in seen:
            pmids.append(pmid)
            seen.add(pmid)
    return pmids


def fetch_pubmed_xml(pmids: list[str]) -> str:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ids = ",".join(pmids)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={ids}&retmode=xml"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "granule-cell-reference-formatter/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            data = response.read().decode("utf-8")
    except urllib.error.URLError as exc:
        if XML_CACHE.exists():
            return XML_CACHE.read_text(encoding="utf-8")
        raise RuntimeError(f"Could not fetch PubMed metadata and no cache exists: {exc}") from exc
    XML_CACHE.write_text(data, encoding="utf-8")
    return data


def initials_with_periods(initials: str) -> str:
    initials = clean_text(initials).replace(".", "").replace(" ", "")
    if not initials:
        return ""
    return " ".join(f"{char}." for char in initials if char.isalpha())


def author_records(article: ET.Element) -> list[dict[str, str]]:
    authors: list[dict[str, str]] = []
    for author in article.findall("./AuthorList/Author"):
        collective = find_text(author, "CollectiveName")
        if collective:
            authors.append({"collective": collective})
            continue
        last = find_text(author, "LastName")
        fore = find_text(author, "ForeName")
        initials = find_text(author, "Initials")
        if last:
            authors.append({"last": last, "fore": fore, "initials": initials})
    return authors


def vancouver_author(author: dict[str, str]) -> str:
    if "collective" in author:
        return author["collective"]
    return clean_text(f"{author.get('last', '')} {author.get('initials', '')}")


def period_author(author: dict[str, str]) -> str:
    if "collective" in author:
        return author["collective"]
    initials = initials_with_periods(author.get("initials", ""))
    return clean_text(f"{author.get('last', '')}, {initials}")


def bibtex_author(author: dict[str, str]) -> str:
    if "collective" in author:
        return "{" + author["collective"] + "}"
    fore = author.get("fore") or author.get("initials", "")
    return clean_text(f"{author.get('last', '')}, {fore}")


def join_with_ampersand(names: list[str]) -> str:
    names = [name for name in names if name]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"{names[0]} & {names[1]}"
    return f"{', '.join(names[:-1])}, & {names[-1]}"


def truncated(names: list[str], limit: int) -> list[str]:
    if len(names) > limit:
        return names[:limit] + ["et al"]
    return names


def close_sentence(text: str) -> str:
    text = clean_text(text)
    if not text:
        return ""
    return text if text.endswith(".") else text + "."


def parse_pubmed(xml_text: str) -> dict[str, dict[str, object]]:
    root = ET.fromstring(xml_text)
    parsed: dict[str, dict[str, object]] = {}
    for pubmed_article in root.findall(".//PubmedArticle"):
        medline = pubmed_article.find("./MedlineCitation")
        article = medline.find("./Article") if medline is not None else None
        if medline is None or article is None:
            continue
        pmid = find_text(medline, "PMID")
        journal = article.find("./Journal")
        issue = journal.find("./JournalIssue") if journal is not None else None
        pubdate = issue.find("./PubDate") if issue is not None else None
        medline_date = find_text(pubdate, "MedlineDate")
        year = find_text(pubdate, "Year") or find_text(article, "./ArticleDate/Year")
        if not year and medline_date:
            year_match = re.search(r"([12][0-9]{3})", medline_date)
            year = year_match.group(1) if year_match else ""
        doi = ""
        for article_id in pubmed_article.findall("./PubmedData/ArticleIdList/ArticleId"):
            if article_id.attrib.get("IdType") == "doi":
                doi = clean_text(article_id.text)
                break
        if not doi:
            for eloc in article.findall("./ELocationID"):
                if eloc.attrib.get("EIdType") == "doi":
                    doi = clean_text(eloc.text)
                    break
        pages = find_text(article, "./Pagination/MedlinePgn")
        if not pages:
            for eloc in article.findall("./ELocationID"):
                if eloc.attrib.get("EIdType") != "doi":
                    pages = clean_text(eloc.text)
                    break
        authors = author_records(article)
        parsed[pmid] = {
            "pmid": pmid,
            "authors": authors,
            "title": node_text(article.find("./ArticleTitle")).rstrip("."),
            "journal": find_text(journal, "Title"),
            "journal_iso": find_text(journal, "ISOAbbreviation") or find_text(journal, "Title"),
            "year": year,
            "volume": find_text(issue, "Volume"),
            "issue": find_text(issue, "Issue"),
            "pages": pages,
            "doi": doi,
        }
    return parsed


def format_vancouver(record: dict[str, object]) -> str:
    authors = [vancouver_author(a) for a in record["authors"]]  # type: ignore[index]
    author_text = ", ".join(truncated(authors, 6))
    journal = record["journal_iso"] or record["journal"]
    vol = record["volume"]
    issue = f"({record['issue']})" if record["issue"] else ""
    pages = f":{record['pages']}" if record["pages"] else ""
    doi = f" doi: {record['doi']}." if record["doi"] else ""
    return (
        f"{author_text}. {record['title']}. {journal}. "
        f"{record['year']};{vol}{issue}{pages}.{doi} PMID: {record['pmid']}."
    ).replace(";.", ".")


def format_nature(record: dict[str, object]) -> str:
    authors = [period_author(a) for a in record["authors"]]  # type: ignore[index]
    author_text = ", ".join(truncated(authors, 5))
    journal = record["journal_iso"] or record["journal"]
    vol = f" {record['volume']}" if record["volume"] else ""
    pages = f", {record['pages']}" if record["pages"] else ""
    doi = f" https://doi.org/{record['doi']}." if record["doi"] else ""
    return (
        f"{close_sentence(author_text)} {record['title']}. {journal}{vol}{pages} "
        f"({record['year']}).{doi} PMID: {record['pmid']}."
    )


def format_apa(record: dict[str, object]) -> str:
    authors = [period_author(a) for a in record["authors"]]  # type: ignore[index]
    author_text = join_with_ampersand(truncated(authors, 20))
    journal = record["journal"]
    volume = record["volume"]
    issue = f"({record['issue']})" if record["issue"] else ""
    pages = f", {record['pages']}" if record["pages"] else ""
    doi = f" https://doi.org/{record['doi']}" if record["doi"] else ""
    return (
        f"{close_sentence(author_text)} ({record['year']}). {record['title']}. "
        f"{journal}, {volume}{issue}{pages}.{doi} PMID: {record['pmid']}."
    ).replace(", .", ".")


def bibtex_key(record: dict[str, object]) -> str:
    authors = record["authors"]  # type: ignore[assignment]
    first = "Reference"
    if authors:
        first_author = authors[0]
        first = re.sub(r"[^A-Za-z0-9]", "", first_author.get("last", first_author.get("collective", "Reference")))
    return f"{first}{record['year']}PMID{record['pmid']}"


def format_bibtex(record: dict[str, object]) -> str:
    authors = " and ".join(bibtex_author(a) for a in record["authors"])  # type: ignore[index]
    fields = {
        "author": authors,
        "title": "{" + str(record["title"]) + "}",
        "journal": str(record["journal"]),
        "year": str(record["year"]),
        "volume": str(record["volume"]),
        "number": str(record["issue"]),
        "pages": str(record["pages"]),
        "doi": str(record["doi"]),
        "pmid": str(record["pmid"]),
    }
    lines = [f"@article{{{bibtex_key(record)},"]
    for key, value in fields.items():
        if value:
            lines.append(f"  {key} = {{{value}}},")
    lines[-1] = lines[-1].rstrip(",")
    lines.append("}")
    return "\n".join(lines)


def write_outputs(pmids: list[str], records: dict[str, dict[str, object]]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    missing = [pmid for pmid in pmids if pmid not in records]
    ordered = [records[pmid] for pmid in pmids if pmid in records]
    alpha = sorted(
        ordered,
        key=lambda rec: (
            str(rec["authors"][0].get("last", rec["authors"][0].get("collective", ""))) if rec["authors"] else "",
            str(rec["year"]),
            str(rec["title"]),
        ),
    )

    master_header = [
        "reference_number_first_appearance",
        "pmid",
        "first_author",
        "year",
        "title",
        "journal",
        "journal_iso",
        "volume",
        "issue",
        "pages",
        "doi",
        "authors_vancouver",
    ]
    master_lines = ["\t".join(master_header)]
    for idx, rec in enumerate(ordered, 1):
        authors = rec["authors"]  # type: ignore[assignment]
        first = ""
        if authors:
            first = authors[0].get("last", authors[0].get("collective", ""))
        row = [
            str(idx),
            str(rec["pmid"]),
            str(first),
            str(rec["year"]),
            str(rec["title"]),
            str(rec["journal"]),
            str(rec["journal_iso"]),
            str(rec["volume"]),
            str(rec["issue"]),
            str(rec["pages"]),
            str(rec["doi"]),
            "; ".join(vancouver_author(a) for a in authors),
        ]
        master_lines.append("\t".join(clean_text(item) for item in row))
    (OUTDIR / "references_master_pubmed.tsv").write_text("\n".join(master_lines) + "\n", encoding="utf-8")

    map_lines = ["reference_number\tpmid\tfirst_author\tyear\ttitle"]
    for idx, rec in enumerate(ordered, 1):
        authors = rec["authors"]  # type: ignore[assignment]
        first = authors[0].get("last", authors[0].get("collective", "")) if authors else ""
        map_lines.append("\t".join([str(idx), str(rec["pmid"]), str(first), str(rec["year"]), str(rec["title"])]))
    (OUTDIR / "references_numbered_citation_map.tsv").write_text("\n".join(map_lines) + "\n", encoding="utf-8")

    vancouver_entries: list[str] = []
    for idx, rec in enumerate(ordered, 1):
        vancouver_entries.append(f"{idx}. {format_vancouver(rec)}")
    (OUTDIR / "references_vancouver_numbered_by_first_appearance.md").write_text(
        "# Vancouver Numbered References\n\n" + "\n\n".join(vancouver_entries) + "\n",
        encoding="utf-8",
    )

    nature_entries: list[str] = []
    for idx, rec in enumerate(ordered, 1):
        nature_entries.append(f"{idx}. {format_nature(rec)}")
    (OUTDIR / "references_nature_numbered_by_first_appearance.md").write_text(
        "# Nature-Style Numbered References\n\n" + "\n\n".join(nature_entries) + "\n",
        encoding="utf-8",
    )

    apa_entries: list[str] = []
    for rec in alpha:
        apa_entries.append(format_apa(rec))
    (OUTDIR / "references_apa_author_year_alphabetical.md").write_text(
        "# APA/Author-Year References\n\n" + "\n\n".join(apa_entries) + "\n",
        encoding="utf-8",
    )

    bibtex = "\n\n".join(format_bibtex(rec) for rec in alpha)
    (OUTDIR / "references_pubmed.bib").write_text(bibtex + "\n", encoding="utf-8")

    readme = f"""# Reference Formatting Outputs

Generated from PMID anchors in:

`Project/manuscript/granule_cell_convergence_manuscript_draft.md`

Files:

- `references_master_pubmed.tsv`: complete PubMed-derived metadata table.
- `references_numbered_citation_map.tsv`: PMID-to-number map in first-appearance order.
- `references_vancouver_numbered_by_first_appearance.md`: biomedical/Vancouver numbered style.
- `references_nature_numbered_by_first_appearance.md`: Nature-like numbered style.
- `references_apa_author_year_alphabetical.md`: author-year alphabetical style.
- `references_pubmed.bib`: BibTeX export for journal systems or later conversion.
- `pubmed_efetch.xml`: cached PubMed XML used for generation.

The numbered styles are ordered by first PMID appearance in the manuscript. If the manuscript is later converted to a numbered citation style, replace in-text author-year citations using `references_numbered_citation_map.tsv`.

PMIDs found: {len(pmids)}
PubMed records retrieved: {len(ordered)}
Missing records: {", ".join(missing) if missing else "none"}
"""
    (OUTDIR / "README.md").write_text(textwrap.dedent(readme), encoding="utf-8")


def main() -> int:
    markdown = MANUSCRIPT.read_text(encoding="utf-8")
    pmids = extract_pmids(markdown)
    if not pmids:
        print("No PMID anchors found.", file=sys.stderr)
        return 1
    xml_text = fetch_pubmed_xml(pmids)
    records = parse_pubmed(xml_text)
    write_outputs(pmids, records)
    missing = [pmid for pmid in pmids if pmid not in records]
    print(f"PMIDs found: {len(pmids)}")
    print(f"PubMed records retrieved: {len(pmids) - len(missing)}")
    if missing:
        print("Missing PMIDs: " + ", ".join(missing))
    print(f"Wrote references to: {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
