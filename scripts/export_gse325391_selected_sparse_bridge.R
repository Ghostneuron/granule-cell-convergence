#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_arg <- cmd_args[grep("^--file=", cmd_args)][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)

base <- file.path(root, "External_Data/GEO/GSE325391")
out_dir <- file.path(root, "Project/processed/gse325391_adult_dg_selected")
results <- file.path(root, "Project/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(base, "GSE325391_adultgc_filtered.RDS")
selected_var_path <- file.path(root, "Project/processed/human_core_normalized_reduced_object/var.tsv")

out_i <- file.path(out_dir, "selected_present_genes_csc_i_int32.bin")
out_p <- file.path(out_dir, "selected_present_genes_csc_p_int32.bin")
out_x <- file.path(out_dir, "selected_present_genes_csc_x_float64.bin")
out_shape <- file.path(out_dir, "selected_present_genes_csc_shape.tsv")
out_features <- file.path(out_dir, "selected_features.tsv")
out_cells <- file.path(out_dir, "cell_metadata.tsv.gz")
out_genes <- file.path(out_dir, "gene_metadata.tsv.gz")
out_umap <- file.path(out_dir, "harmony_umap.tsv.gz")
out_readme <- file.path(out_dir, "README.md")
out_summary <- file.path(results, "gse325391_selected_sparse_bridge_summary.tsv")
out_md <- file.path(results, "gse325391_selected_sparse_bridge_summary.md")

norm_gene <- function(x) toupper(trimws(as.character(x)))

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

cat("Reading", rds_path, "\n")
obj <- readRDS(rds_path)
counts <- obj[["RNA"]]@counts
selected <- read.delim(selected_var_path, stringsAsFactors = FALSE)

if (!("gene" %in% colnames(selected))) {
  stop("Selected var table must contain a gene column")
}

feature_norm <- norm_gene(rownames(counts))
selected_norm <- norm_gene(selected$gene)
match_idx <- match(selected_norm, feature_norm)
present <- !is.na(match_idx)
present_source_rows <- match_idx[present]

cat("Selected genes:", nrow(selected), "present:", sum(present), "missing:", sum(!present), "\n")
selected_present_counts <- counts[present_source_rows, , drop = FALSE]
selected_present_counts <- as(selected_present_counts, "dgCMatrix")

writeBin(as.integer(selected_present_counts@i), out_i, size = 4, endian = "little")
writeBin(as.integer(selected_present_counts@p), out_p, size = 4, endian = "little")
writeBin(as.numeric(selected_present_counts@x), out_x, size = 8, endian = "little")

shape <- data.frame(
  dataset = "GSE325391",
  assay = "RNA",
  orientation = "selected_present_genes_by_cells_csc",
  n_selected_genes = nrow(selected),
  n_present_selected_genes = sum(present),
  n_missing_selected_genes = sum(!present),
  n_cells = ncol(selected_present_counts),
  nnz_present_selected = length(selected_present_counts@x),
  i_dtype = "int32_little_endian",
  p_dtype = "int32_little_endian",
  x_dtype = "float64_little_endian",
  stringsAsFactors = FALSE
)
write.table(shape, out_shape, sep = "\t", row.names = FALSE, quote = FALSE)

selected_out <- selected
selected_out$selected_feature_index <- seq_len(nrow(selected_out)) - 1L
selected_out$gene_norm <- selected_norm
selected_out$present_in_gse325391 <- present
selected_out$source_feature <- ifelse(present, rownames(counts)[match_idx], "")
selected_out$source_row_in_rds_1based <- ifelse(present, match_idx, NA_integer_)
present_rank <- rep(NA_integer_, length(present))
present_rank[present] <- seq_len(sum(present)) - 1L
selected_out$present_matrix_row_index <- present_rank
write.table(selected_out, out_features, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

metadata <- obj[[]]
metadata$cell_name <- colnames(obj)
metadata$cell_id <- paste("GSE325391", "adultgc", metadata$cell_name, sep = ":")
metadata <- metadata[, c("cell_id", "cell_name", setdiff(colnames(metadata), c("cell_id", "cell_name")))]
write_tsv_gz(metadata, out_cells)

gene_metadata <- data.frame(
  gene = rownames(counts),
  gene_norm = feature_norm,
  n_counts = as.numeric(Matrix::rowSums(counts)),
  n_cells = as.integer(Matrix::rowSums(counts != 0)),
  in_human_core_selected_features = feature_norm %in% selected_norm,
  stringsAsFactors = FALSE
)
write_tsv_gz(gene_metadata, out_genes)

if ("harmony.umap" %in% Reductions(obj)) {
  umap <- as.data.frame(Embeddings(obj, "harmony.umap"))
  umap$cell_name <- rownames(umap)
  umap$cell_id <- paste("GSE325391", "adultgc", umap$cell_name, sep = ":")
  umap <- umap[, c("cell_id", "cell_name", setdiff(colnames(umap), c("cell_id", "cell_name")))]
  write_tsv_gz(umap, out_umap)
}

summary <- data.frame(
  dataset = "GSE325391",
  source_rds = file.path("External_Data/GEO/GSE325391", basename(rds_path)),
  output_dir = sub(paste0("^", root, "/"), "", out_dir),
  n_cells = ncol(counts),
  n_source_genes = nrow(counts),
  source_nnz = length(counts@x),
  selected_genes = nrow(selected),
  selected_genes_present = sum(present),
  selected_genes_missing = sum(!present),
  selected_present_nnz = length(selected_present_counts@x),
  stringsAsFactors = FALSE
)
write.table(summary, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

readme <- c(
  "# GSE325391 Adult DG Selected Sparse Bridge",
  "",
  "This folder contains a compact export from the downloaded `GSE325391_adultgc_filtered.RDS` Seurat object.",
  "The exported count matrix uses the same selected feature list as the current normalized human-core object.",
  "",
  "## Matrix",
  "",
  "- Orientation before Python conversion: selected-present genes by cells, CSC arrays.",
  "- Missing selected genes are recorded in `selected_features.tsv` and should be restored as zero rows during Python conversion.",
  "- The source full RNA count matrix remains in the downloaded Seurat RDS.",
  "",
  "## Files",
  "",
  "- `selected_present_genes_csc_i_int32.bin`",
  "- `selected_present_genes_csc_p_int32.bin`",
  "- `selected_present_genes_csc_x_float64.bin`",
  "- `selected_present_genes_csc_shape.tsv`",
  "- `selected_features.tsv`",
  "- `cell_metadata.tsv.gz`",
  "- `gene_metadata.tsv.gz`",
  "- `harmony_umap.tsv.gz`",
  ""
)
writeLines(readme, out_readme)

md <- c(
  "# GSE325391 Selected Sparse Bridge",
  "",
  "Date built: 2026-06-21",
  "",
  "## Summary",
  "",
  paste0("- Cells/nuclei: ", ncol(counts)),
  paste0("- Source genes: ", nrow(counts)),
  paste0("- Source RNA count nnz: ", length(counts@x)),
  paste0("- Human-core selected genes: ", nrow(selected)),
  paste0("- Selected genes present: ", sum(present)),
  paste0("- Selected genes missing: ", sum(!present)),
  paste0("- Selected-present count nnz: ", length(selected_present_counts@x)),
  "",
  "## Interpretation",
  "",
  "This is the bridge export for immediate label mapping into the current human-core feature space. The full downloaded RDS remains the primary raw object.",
  "",
  "## Outputs",
  "",
  paste0("- Bridge directory: `", sub(paste0("^", root, "/"), "", out_dir), "`"),
  paste0("- Summary TSV: `", sub(paste0("^", root, "/"), "", out_summary), "`"),
  ""
)
writeLines(md, out_md)

cat("Wrote", out_dir, "\n")
cat("Wrote", out_summary, "\n")
cat("Wrote", out_md, "\n")
cat("selected_present_nnz=", length(selected_present_counts@x), "\n", sep = "")
