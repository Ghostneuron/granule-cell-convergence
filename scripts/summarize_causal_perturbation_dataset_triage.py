#!/usr/bin/env python3
"""Summarize curated public perturbation datasets for causal follow-up planning."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"
TRIAGE = RESULTS / "causal_perturbation_dataset_triage.tsv"
OUT_SUMMARY = RESULTS / "causal_perturbation_node_summary.tsv"
OUT_MD = RESULTS / "causal_perturbation_dataset_triage.md"


RELEVANCE_RANK = {
    "high": 4,
    "medium_high": 3,
    "medium": 2,
    "low_indirect": 1,
    "low": 1,
}

NODE_CALLS = {
    "NFIA": {
        "priority_call": "strong regulatory follow-up candidate",
        "model_role": "cerebellar GNP regulatory target prior plus hippocampal astrocyte/niche sensitivity",
        "main_caveat": "direct normal dentate granule-cell NFIA perturbation remains missing",
    },
    "BDNF": {
        "priority_call": "usable pathway perturbation candidate",
        "model_role": "BDNF/TrkB response signature for neural differentiation and hippocampal plasticity modules",
        "main_caveat": "available datasets are neural progenitor, hippocampal injury/plasticity, or indirect cerebellar chromatin contexts",
    },
    "TGFBeta": {
        "priority_call": "supporting secreted-factor/niche candidate",
        "model_role": "TGF-beta2 secreted-response or SMAD-patterning signatures for conditioned-medium logic",
        "main_caveat": "few direct granule-cell TGF-beta perturbation transcriptomes were found",
    },
    "TGFBeta_SMAD": {
        "priority_call": "supporting hindbrain-patterning candidate",
        "model_role": "SMAD/Wnt patterning context for hindbrain-like neural stem-cell states",
        "main_caveat": "dual-SMAD patterning is not equivalent to TGF-beta ligand response",
    },
    "SHH": {
        "priority_call": "best causal extension candidate",
        "model_role": "cerebellar GNP proliferation/fate perturbation and SHH-responsive progenitor control",
        "main_caveat": "many datasets are medulloblastoma or tumor-susceptibility models rather than normal granule morphology",
    },
    "SHH_Niche": {
        "priority_call": "supporting niche-response candidate",
        "model_role": "SHH/SAG-responsive astrocyte or glial niche signature",
        "main_caveat": "not granule lineage",
    },
    "RBFOX3": {
        "priority_call": "strong postmitotic functional follow-up candidate",
        "model_role": "dentate/hippocampal synaptic maturation, splicing, and plasticity perturbation signature",
        "main_caveat": "best Rbfox3 KO RNA-seq has low sample count; Rbfox1/3 knockdown datasets are not RBFOX3-specific",
    },
    "HMGN2": {
        "priority_call": "weak public perturbation support",
        "model_role": "exploratory chromatin/RNA-regulatory competence screen only",
        "main_caveat": "no clean neural or granule-lineage HMGN2-focused perturbation dataset was found",
    },
}


def normalized_node(node: str) -> str:
    if node.startswith("BDNF"):
        return "BDNF"
    return node


def main() -> None:
    triage = pd.read_csv(TRIAGE, sep="\t")
    triage["node_group"] = triage["node"].map(normalized_node)
    triage["relevance_rank"] = triage["relevance_to_project"].map(RELEVANCE_RANK).fillna(0).astype(int)

    rows = []
    for node, sub in triage.groupby("node_group", sort=False):
        sub = sub.sort_values(["relevance_rank", "accession"], ascending=[False, True])
        high = sub.loc[sub["relevance_to_project"].eq("high")]
        med_high = sub.loc[sub["relevance_to_project"].eq("medium_high")]
        usable = sub.loc[sub["relevance_rank"].ge(2)]
        top = sub.iloc[0]
        call = NODE_CALLS.get(
            node,
            {
                "priority_call": "candidate",
                "model_role": "candidate perturbation signature",
                "main_caveat": "requires manual review",
            },
        )
        rows.append(
            {
                "node": node,
                "n_candidate_datasets": len(sub),
                "n_high": len(high),
                "n_medium_high": len(med_high),
                "n_usable_medium_or_higher": len(usable),
                "top_accession": top["accession"],
                "top_relevance": top["relevance_to_project"],
                "all_accessions": ";".join(sub["accession"].astype(str)),
                "priority_call": call["priority_call"],
                "recommended_model_role": call["model_role"],
                "main_caveat": call["main_caveat"],
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    lines = [
        "# Public Perturbation Dataset Triage",
        "",
        "This curation asks whether public datasets can move the project from a hierarchical evidence model toward causal or quasi-causal mixed-effects tests. The answer is uneven: `SHH`, `RBFOX3`, and `NFIA` have the most useful perturbation resources; `BDNF` and TGF-beta/SMAD have usable pathway-response resources; `HMGN2` currently has weak public perturbation support.",
        "",
        "## Node Summary",
        "",
        "| Node | Candidate datasets | Top accession | Priority call | Main caveat |",
        "|---|---:|---|---|---|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['node']} | {row['n_candidate_datasets']} | {row['top_accession']} | "
            f"{row['priority_call']} | {row['main_caveat']} |"
        )
    lines.extend(
        [
            "",
            "## Recommended Use",
            "",
            "1. Use `SHH`/cerebellar GNP perturbation resources first for the clearest causal extension of the cerebellar branch.",
            "2. Use `RBFOX3` datasets as postmitotic dentate/hippocampal functional perturbation support, especially synaptic and plasticity modules.",
            "3. Use `NFIA` datasets as regulatory-target evidence, strongest in cerebellar GNPs and weaker but useful in hippocampal astrocyte/niche context.",
            "4. Treat `BDNF` and TGF-beta resources as pathway-response signatures rather than direct granule morphology perturbations.",
            "5. Keep `HMGN2` as a hypothesis-generating chromatin-competence candidate until focused neural perturbation data are found or generated.",
            "",
            f"Detailed dataset table: `{TRIAGE.name}`",
            f"Node summary table: `{OUT_SUMMARY.name}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))

    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
