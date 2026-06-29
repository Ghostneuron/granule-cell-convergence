#!/usr/bin/env python3
"""Score evidence alignment for three granule-cell convergence hypotheses.

This script adds an interpretable hypothesis-comparison layer on top of the
existing hierarchical evidence model. It does not estimate causal effects or
posterior model probabilities. Instead, each hypothesis receives transparent
prediction coefficients for each evidence term, and observed term scores are
combined into a bounded alignment index.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "Project/results"

TERMS = RESULTS / "hierarchical_integrative_model_terms.tsv"
LAYER_SUMMARY = RESULTS / "hierarchical_integrative_model_layer_summary.tsv"

OUT_MATRIX = RESULTS / "hypothesis_support_score_matrix.tsv"
OUT_SCORES = RESULTS / "hypothesis_support_scores.tsv"
OUT_NOTE = RESULTS / "hypothesis_support_score_model.md"
OUT_FIGURE = RESULTS / "hypothesis_support_score_matrix.png"


HYPOTHESES = {
    "H1_shared_granule_fate": {
        "label": "H1: hidden shared granule-cell fate",
        "summary": (
            "Dentate and cerebellar granule cells share a broad hidden fate or "
            "recent lineage-like molecular identity."
        ),
        "prediction": {
            "Y": 0.5,
            "F": -1.0,
            "C": 0.5,
            "I": -1.0,
            "T": -0.5,
            "N": -0.5,
            "E": 0.3,
            "M": 0.0,
            "A": 0.0,
            "R": 0.0,
        },
    },
    "H2_identity_coupled_assembly": {
        "label": "H2: identity-coupled assembly convergence",
        "summary": (
            "Distinct regional lineages preserve branch identity while reusing "
            "overlapping postmitotic construction modules."
        ),
        "prediction": {
            "Y": 1.0,
            "F": 1.0,
            "C": 1.0,
            "I": 1.0,
            "T": 0.6,
            "N": 0.4,
            "E": 0.7,
            "M": 0.2,
            "A": 0.2,
            "R": 0.2,
        },
    },
    "H3_niche_circuit_constraint": {
        "label": "H3: niche/circuit constraint convergence",
        "summary": (
            "Stage-dependent niche signals and sparse-expansion circuit constraints "
            "favor compact granule-like designs."
        ),
        "prediction": {
            "Y": 0.3,
            "F": 0.3,
            "C": 0.2,
            "I": 0.2,
            "T": 1.0,
            "N": 1.0,
            "E": 0.3,
            "M": 1.0,
            "A": 1.0,
            "R": 1.0,
        },
    },
    "H2_H3_integrated_model": {
        "label": "Integrated H2+H3 developmental-convergence model",
        "summary": (
            "Distinct regional lineages enter related assembly states that are "
            "modulated by stage/niche context and constrained by sparse-expansion "
            "circuit design."
        ),
        "prediction": {
            "Y": 1.0,
            "F": 1.0,
            "C": 1.0,
            "I": 1.0,
            "T": 1.0,
            "N": 1.0,
            "E": 0.7,
            "M": 0.7,
            "A": 0.7,
            "R": 0.8,
        },
    },
}


INTERPRETIVE_DOMAIN = {
    "Y": "direct_transcriptomic_configuration",
    "F": "direct_transcriptomic_configuration",
    "C": "direct_transcriptomic_configuration",
    "I": "direct_transcriptomic_configuration",
    "T": "stage_niche_regulatory",
    "N": "stage_niche_regulatory",
    "E": "stage_niche_regulatory",
    "M": "external_morphology_activity_circuit",
    "A": "external_morphology_activity_circuit",
    "R": "external_morphology_activity_circuit",
}


def classify(score: float) -> str:
    if score >= 75:
        return "strong_alignment"
    if score >= 60:
        return "moderate_alignment"
    if score >= 45:
        return "mixed_or_weak_alignment"
    return "low_alignment_or_opposition"


def markdown_score_table(scores: pd.DataFrame) -> str:
    rows = [
        "| Hypothesis | Support index | Support class | Signed alignment |",
        "|---|---:|---|---:|",
    ]
    for _, row in scores.iterrows():
        rows.append(
            "| {label} | {support:.3f} | {support_class} | {alignment:.3f} |".format(
                label=row["hypothesis_label"],
                support=row["support_index_0_to100"],
                support_class=row["support_class"],
                alignment=row["weighted_signed_alignment_minus1_to1"],
            )
        )
    return "\n".join(rows)


def main() -> None:
    terms = pd.read_csv(TERMS, sep="\t")
    layer = pd.read_csv(LAYER_SUMMARY, sep="\t")

    observed = layer.set_index("term_symbol")["weighted_mean_score"].to_dict()
    weights = terms.set_index("term_symbol")["default_weight"].to_dict()
    names = terms.set_index("term_symbol")["term"].to_dict()

    matrix_rows = []
    for term_symbol in names:
        row = {
            "term_symbol": term_symbol,
            "term": names[term_symbol],
            "interpretive_domain": INTERPRETIVE_DOMAIN.get(term_symbol, "other"),
            "observed_term_score": observed.get(term_symbol),
            "evidence_weight": weights.get(term_symbol),
        }
        for hypothesis_id, spec in HYPOTHESES.items():
            coeff = spec["prediction"].get(term_symbol, 0.0)
            row[f"{hypothesis_id}_coefficient"] = coeff
            row[f"{hypothesis_id}_weighted_alignment"] = (
                observed.get(term_symbol, 0.0) * weights.get(term_symbol, 0.0) * coeff
            )
            row[f"{hypothesis_id}_weighted_relevance"] = (
                weights.get(term_symbol, 0.0) * abs(coeff)
            )
        matrix_rows.append(row)

    matrix = pd.DataFrame(matrix_rows)
    matrix.to_csv(OUT_MATRIX, sep="\t", index=False)

    score_rows = []
    for hypothesis_id, spec in HYPOTHESES.items():
        alignment_col = f"{hypothesis_id}_weighted_alignment"
        relevance_col = f"{hypothesis_id}_weighted_relevance"
        denominator = matrix[relevance_col].sum()
        signed_alignment = matrix[alignment_col].sum() / denominator
        support_0_100 = 50.0 * (signed_alignment + 1.0)
        row = {
            "hypothesis_id": hypothesis_id,
            "hypothesis_label": spec["label"],
            "hypothesis_summary": spec["summary"],
            "weighted_signed_alignment_minus1_to1": signed_alignment,
            "support_index_0_to100": support_0_100,
            "support_class": classify(support_0_100),
            "total_weighted_relevance": denominator,
        }
        for domain in sorted(set(INTERPRETIVE_DOMAIN.values())):
            subset = matrix[matrix["interpretive_domain"] == domain]
            domain_denominator = subset[relevance_col].sum()
            if domain_denominator > 0:
                domain_alignment = subset[alignment_col].sum() / domain_denominator
                row[f"{domain}_support_index_0_to100"] = 50.0 * (domain_alignment + 1.0)
            else:
                row[f"{domain}_support_index_0_to100"] = pd.NA
        score_rows.append(row)

    scores = pd.DataFrame(score_rows).sort_values("support_index_0_to100", ascending=False)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    figure_written = False
    try:
        plot_summary(matrix, scores)
        figure_written = True
    except ModuleNotFoundError:
        figure_written = False

    lines = [
        "# Hypothesis Support Score Matrix",
        "",
        "This analysis compares the three working hypotheses using the existing",
        "hierarchical integrative evidence model. It is an evidence-alignment",
        "index, not a causal test, posterior probability, or replacement for",
        "experimental validation.",
        "",
        "For each evidence term \\(t\\), the observed term score \\(S_t\\) comes",
        "from the hierarchical model layer summary and lies on a bounded",
        "\\([-1,1]\\) scale. The evidence weight \\(w_t\\) comes from the",
        "term-specification table. Each hypothesis \\(h\\) is assigned a",
        "transparent prediction coefficient \\(a_{h,t}\\) from -1 to 1, where",
        "positive values indicate that the term is expected under the hypothesis,",
        "negative values indicate that the term argues against that hypothesis,",
        "and zero means the term is not used to discriminate that hypothesis.",
        "",
        "The score is:",
        "",
        "\\[",
        "A_h = \\frac{\\sum_t w_t a_{h,t} S_t}{\\sum_t w_t |a_{h,t}|}",
        "\\]",
        "",
        "\\[",
        "\\mathrm{SupportIndex}_h = 50(1 + A_h)",
        "\\]",
        "",
        "Thus 50 is neutral, values above 50 indicate alignment, and values below",
        "50 indicate opposition or mismatch.",
        "",
        "## Current scores",
        "",
        markdown_score_table(scores),
        "",
        "## Interpretation",
        "",
        "The score matrix argues against a single hidden shared-fate explanation",
        "because branch-matched fate polarity, identity-coupled configuration,",
        "stage-window behavior and branch-specific niche signals are all positive",
        "in the observed data. The niche/circuit hypothesis receives the strongest",
        "standalone numerical alignment because the morphology, activity and",
        "resource-constraint layers are strong. However, those layers do not",
        "replace the direct transcriptomic evidence for identity-coupled assembly.",
        "The manuscript conclusion should therefore remain an integrated",
        "developmental-convergence model: distinct regional lineages reuse related",
        "assembly machinery under stage, niche and sparse-expansion circuit",
        "constraints.",
        "",
        "## Outputs",
        "",
        f"- `{OUT_MATRIX.name}`: term-by-hypothesis coefficient matrix.",
        f"- `{OUT_SCORES.name}`: hypothesis-level support scores and domain scores.",
    ]
    if figure_written:
        lines.append(f"- `{OUT_FIGURE.name}`: visual support-index and coefficient summary.")
    else:
        lines.append("- Figure generation was skipped because matplotlib was not available in the active Python runtime.")
    OUT_NOTE.write_text("\n".join(lines) + "\n")


def plot_summary(matrix: pd.DataFrame, scores: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    score_plot = scores.sort_values("support_index_0_to100")
    hypothesis_ids = [h for h in HYPOTHESES if h in scores["hypothesis_id"].values]
    term_symbols = matrix["term_symbol"].tolist()
    coeff = np.array(
        [
            matrix[f"{hypothesis_id}_coefficient"].to_numpy(dtype=float)
            for hypothesis_id in hypothesis_ids
        ]
    )

    fig = plt.figure(figsize=(12, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[0.85, 1.15])
    ax_bar = fig.add_subplot(grid[0])
    ax_heat = fig.add_subplot(grid[1])

    colors = [
        "#4f8fba" if value >= 75 else "#7aa974" if value >= 60 else "#c98c42" if value >= 45 else "#b85c5c"
        for value in score_plot["support_index_0_to100"]
    ]
    ax_bar.barh(score_plot["hypothesis_label"], score_plot["support_index_0_to100"], color=colors)
    ax_bar.axvline(50, color="#555555", linewidth=1, linestyle="--")
    ax_bar.set_xlim(0, 100)
    ax_bar.set_xlabel("Support index (0-100; 50 = neutral)")
    ax_bar.set_title("Hypothesis support from hierarchical evidence synthesis", loc="left", fontweight="bold")
    for y_pos, value in enumerate(score_plot["support_index_0_to100"]):
        ax_bar.text(min(value + 1.2, 98), y_pos, f"{value:.1f}", va="center", fontsize=9)
    ax_bar.spines[["top", "right"]].set_visible(False)

    im = ax_heat.imshow(coeff, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax_heat.set_xticks(range(len(term_symbols)))
    ax_heat.set_xticklabels(term_symbols)
    ax_heat.set_yticks(range(len(hypothesis_ids)))
    ax_heat.set_yticklabels([HYPOTHESES[h]["label"].replace(": ", ":\n") for h in hypothesis_ids])
    ax_heat.set_xlabel("Evidence term")
    ax_heat.set_title("Prediction coefficients used in the support-index formula", loc="left", fontweight="bold")
    for i in range(coeff.shape[0]):
        for j in range(coeff.shape[1]):
            ax_heat.text(j, i, f"{coeff[i, j]:.1f}", ha="center", va="center", fontsize=8)
    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.8, pad=0.01)
    cbar.set_label("Coefficient")

    fig.savefig(OUT_FIGURE, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
