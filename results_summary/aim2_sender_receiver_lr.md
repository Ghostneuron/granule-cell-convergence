# Aim 2 Focused Sender-Receiver Ligand-Receptor Prediction

Date built: 2026-06-25

## Purpose

This analysis upgrades the earlier Aim 2 pathway-readiness audit by scoring ligand expression in niche sender classes and receptor expression in granule-lineage receiver classes.

## Scope

- Cerebellum: `GSE122357` P0, P8a, and P8b mouse cerebellum; senders are Purkinje cells, astrocytes as a Bergmann/astroglial proxy, microglia, endothelial cells, and supporting oligodendroglia; receivers are granule precursors and granule cells.
- Dentate SGZ: `GSE104323` adult mouse dentate gyrus; senders are astrocyte states, endothelial cells, PVM/macrophage as a microglia-like proxy, and vascular/support classes; receivers are RGL, nIPC, neuroblast, immature GC, juvenile GC, and adult GC states.
- LR database: previous Aim 2 curated pairs plus niche/glial/vascular pairs such as CXCL12/CXCR4, VEGF/FLT1/KDR, APOE/LRP, SPP1, CSF1, IGF1, SLIT/ROBO, and Ephrin/EPH.

## Method

- For each sample, raw matrices were reduced to LR genes and grouped by sender/receiver labels.
- Ligand support is sender ligand relative expression multiplied by detection support.
- Receptor support is receiver receptor relative expression multiplied by detection support.
- The final LR expression score is the geometric support product, so a pair scores well only when the ligand is present in the sender and the receptor is present in the receiver.

## Scale

- LR pairs scored: 47.
- Sender-receiver-pair predictions: 3,008 core-focus rows.
- Supported predictions: 208 moderate/high expression-supported rows.
- Source samples: 4.

Region-level support:
- cerebellum: 97/1128 supported, median score 0.002, max 0.502.
- dentate_sgz: 111/1880 supported, median score 0.001, max 0.557.

Sender-class highlights:
- cerebellum / astroglial_bergmann_proxy: 15/282 supported, median score 0.003.
- cerebellum / purkinje: 43/282 supported, median score 0.003.
- cerebellum / microglia: 20/282 supported, median score 0.001.
- cerebellum / endothelial: 19/282 supported, median score 0.000.
- dentate_sgz / endothelial: 24/376 supported, median score 0.002.
- dentate_sgz / astrocyte: 64/1128 supported, median score 0.001.
- dentate_sgz / microglia_macrophage_proxy: 23/376 supported, median score 0.000.

Top pathway summaries:
- dentate_sgz / APOE/LRP: median score 0.171; 33/40 supported.
- cerebellum / APOE/LDLR: median score 0.102; 12/24 supported.
- cerebellum / Ephrin/EPH: median score 0.089; 9/24 supported.
- cerebellum / HBEGF/EGFR: median score 0.042; 9/24 supported.
- dentate_sgz / APOE/LDLR: median score 0.030; 5/40 supported.
- cerebellum / IGF1/IGF1R: median score 0.029; 10/24 supported.
- cerebellum / APOE/LRP: median score 0.020; 1/24 supported.
- dentate_sgz / HBEGF/EGFR: median score 0.015; 4/40 supported.
- cerebellum / KITLG/KIT: median score 0.015; 2/24 supported.
- cerebellum / THBS/CD47: median score 0.013; 4/24 supported.
- cerebellum / SLIT/ROBO: median score 0.013; 8/48 supported.
- cerebellum / SPP1/integrin: median score 0.011; 3/24 supported.

Top individual predictions:
- dentate_sgz Immature astro -> Juvenile GC: SEMA6A to PLXNA4 (SEMA6A->PLXNA4), score 0.557.
- dentate_sgz Immature astro -> Immature GC: SEMA6A to PLXNA2 (SEMA6A->PLXNA2), score 0.521.
- cerebellum Purkinje -> GC precursor: IGF1 to IGF1R (IGF1->IGF1R), score 0.502.
- cerebellum Purkinje -> GC precursor: IGF1 to IGF1R (IGF1->IGF1R), score 0.489.
- cerebellum Microglia -> GC precursor: IGF1 to IGF1R (IGF1->IGF1R), score 0.448.
- dentate_sgz Immature astro -> Immature GC: SEMA6A to PLXNA4 (SEMA6A->PLXNA4), score 0.448.
- cerebellum Purkinje -> GC precursor: SHH to SMO readiness (SHH->SMO), score 0.447.
- cerebellum Purkinje -> GC precursor: SHH to PTCH1 (SHH->PTCH1), score 0.418.
- cerebellum Purkinje -> GC precursor: SHH to SMO readiness (SHH->SMO), score 0.409.
- cerebellum Microglia -> GC precursor: IGF1 to IGF1R (IGF1->IGF1R), score 0.391.
- cerebellum Purkinje -> GC precursor: SEMA7A to PLXNC1 (SEMA7A->PLXNC1), score 0.374.
- cerebellum Purkinje -> GC precursor: SLIT2 to ROBO2 (SLIT2->ROBO2), score 0.365.

## Interpretation

- This gives us the first focused sender-receiver LR layer for Aim 2.
- The cerebellar analysis directly includes Purkinje, microglia, and endothelial senders; Bergmann glia are not separately annotated in `GSE122357`, so the astrocyte class should be described as a Bergmann/astroglial proxy.
- The dentate analysis directly includes astrocyte and endothelial senders; `GSE104323` has PVM rather than a clean microglia class, so SGZ microglia claims should be phrased as PVM/microglia-macrophage proxy unless we add a human or mouse dataset with explicit microglia labels.
- The result is stronger than the previous pathway-readiness audit for niche directionality, but still not proof of spatial contact or protein secretion.

## Outputs

- LR pair table: `Project/results/aim2_sender_receiver_lr_pairs.tsv`
- Group expression table: `Project/results/aim2_sender_receiver_lr_group_expression.tsv.gz`
- Prediction table: `Project/results/aim2_sender_receiver_lr_predictions.tsv.gz`
- Summary table: `Project/results/aim2_sender_receiver_lr_summary.tsv`
- Top predictions: `Project/results/aim2_sender_receiver_lr_top_predictions.tsv`
- Coverage table: `Project/results/aim2_sender_receiver_lr_coverage.tsv`
- Plot: `Project/results/aim2_sender_receiver_lr.png`
