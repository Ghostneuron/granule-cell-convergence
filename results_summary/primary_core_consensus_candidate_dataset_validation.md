# Consensus Candidate Dataset-Aware Validation

Date built: 2026-06-22

## Purpose

This validation asks whether the 24 cross-screen consensus candidates remain positive across individual datasets/samples, rather than only in pooled pseudobulk summaries.

## Summary

- Consensus candidates tested: 24.
- Dataset/sample/gene branch-delta rows: 1,608.
- Genes robust across all available screen/branch tests: 6.

## Robust Genes

- `GABRA2`, `GPM6A`, `KCNK1`, `NFIA`, `NFIB`, `RFX3`

## Top Validation Summary

- `GABRA2`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.80.
- `GPM6A`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.92.
- `KCNK1`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.88.
- `NFIA`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 1.00.
- `NFIB`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.86.
- `RFX3`: 4/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.96.
- `CACNA2D1`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.55.
- `GABRB3`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.55.
- `GRIN2B`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.51.
- `KCND2`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.59.
- `KCNJ3`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.73.
- `KCNJ6`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.65.
- `PPP3CA`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.65.
- `STXBP5L`: 3/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.47.
- `CALM2`: 2/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `CAMTA1`: 2/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `MAP3K4`: 2/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `MAPK1`: 2/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `SYNPR`: 2/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `ADD2`: 1/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.43.
- `CACNA1E`: 1/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.63.
- `KCND3`: 1/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.63.
- `STXBP1`: 1/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.67.
- `STXBP5`: 1/4 screen-branch tests robust; median delta 0.500; minimum positive-unit fraction 0.53.

## Robustness Rule

A screen/branch is robust if it has at least 2 dataset/sample units, >=75% positive deltas, and median rank delta > 0.

## Outputs

- Dataset deltas: `Project/results/primary_core_consensus_candidate_dataset_deltas.tsv`
- Validation summary: `Project/results/primary_core_consensus_candidate_dataset_validation.tsv`
- Heatmap: `Project/results/primary_core_consensus_candidate_dataset_validation_heatmap.png`
