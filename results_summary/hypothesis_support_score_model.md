# Hypothesis Support Score Matrix

This analysis compares the three working hypotheses using the existing
hierarchical integrative evidence model. It is an evidence-alignment
index, not a causal test, posterior probability, or replacement for
experimental validation.

For each evidence term \(t\), the observed term score \(S_t\) comes
from the hierarchical model layer summary and lies on a bounded
\([-1,1]\) scale. The evidence weight \(w_t\) comes from the
term-specification table. Each hypothesis \(h\) is assigned a
transparent prediction coefficient \(a_{h,t}\) from -1 to 1, where
positive values indicate that the term is expected under the hypothesis,
negative values indicate that the term argues against that hypothesis,
and zero means the term is not used to discriminate that hypothesis.

The score is:

\[
A_h = \frac{\sum_t w_t a_{h,t} S_t}{\sum_t w_t |a_{h,t}|}
\]

\[
\mathrm{SupportIndex}_h = 50(1 + A_h)
\]

Thus 50 is neutral, values above 50 indicate alignment, and values below
50 indicate opposition or mismatch.

## Current scores

| Hypothesis | Support index | Support class | Signed alignment |
|---|---:|---|---:|
| H3: niche/circuit constraint convergence | 76.887 | strong_alignment | 0.538 |
| Integrated H2+H3 developmental-convergence model | 71.794 | moderate_alignment | 0.436 |
| H2: identity-coupled assembly convergence | 69.096 | moderate_alignment | 0.382 |
| H1: hidden shared granule-cell fate | 39.328 | low_alignment_or_opposition | -0.213 |

## Interpretation

The score matrix argues against a single hidden shared-fate explanation
because branch-matched fate polarity, identity-coupled configuration,
stage-window behavior and branch-specific niche signals are all positive
in the observed data. The niche/circuit hypothesis receives the strongest
standalone numerical alignment because the morphology, activity and
resource-constraint layers are strong. However, those layers do not
replace the direct transcriptomic evidence for identity-coupled assembly.
The manuscript conclusion should therefore remain an integrated
developmental-convergence model: distinct regional lineages reuse related
assembly machinery under stage, niche and sparse-expansion circuit
constraints.

## Outputs

- `hypothesis_support_score_matrix.tsv`: term-by-hypothesis coefficient matrix.
- `hypothesis_support_scores.tsv`: hypothesis-level support scores and domain scores.
- `hypothesis_support_score_matrix.png`: visual support-index and coefficient summary.
