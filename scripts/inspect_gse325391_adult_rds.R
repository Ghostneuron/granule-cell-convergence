#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_arg <- cmd_args[grep("^--file=", cmd_args)][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
base <- file.path(root, "External_Data/GEO/GSE325391")
results <- file.path(root, "Project/results")
processed <- file.path(root, "Project/processed/gse325391_adult_dg")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
dir.create(processed, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(base, "GSE325391_adultgc_filtered.RDS")
out_summary <- file.path(results, "gse325391_adult_rds_summary.tsv")
out_meta_profile <- file.path(results, "gse325391_adult_metadata_profile.tsv")
out_group_summary <- file.path(results, "gse325391_adult_metadata_group_summary.tsv")
out_gene_probe <- file.path(results, "gse325391_adult_marker_gene_presence.tsv")
out_md <- file.path(results, "gse325391_adult_rds_inspection.md")

norm_gene <- function(x) toupper(trimws(as.character(x)))

marker_genes <- unique(norm_gene(c(
  "PROX1", "DCX", "NEUROD1", "NEUROD2", "SOX11", "SOX4", "MEX3A", "CALB2",
  "STMN2", "TUBB3", "BCL11B", "CALB1", "SLC17A7", "GRIN1", "GRIN2A",
  "CAMK2A", "RIMS1", "SYT1", "SNAP25", "RBFOX3", "ATOH1", "PAX6",
  "ZIC1", "ZIC2", "ZIC3", "BARHL1", "PDE1C", "GABRA6", "GABRA1",
  "GRID2", "AQP4", "GFAP", "ALDH1L1", "SLC1A2", "MBP", "PLP1",
  "MOG", "MOBP", "PDGFRA", "CSPG4", "P2RY12", "AIF1", "C1QA",
  "TYROBP", "CLDN5", "PECAM1", "VWF", "TTR", "AQP1"
)))

cat("Reading", rds_path, "\n")
obj <- readRDS(rds_path)
obj_class <- paste(class(obj), collapse = ";")

is_seurat <- inherits(obj, "Seurat")
metadata <- data.frame()
assays <- character()
reductions <- character()
default_assay <- ""
n_features <- NA_integer_
n_cells <- NA_integer_
feature_names <- character()
cell_names <- character()
active_ident_counts <- data.frame()

if (is_seurat) {
  suppressPackageStartupMessages({
    library(SeuratObject)
  })
  assays <- Assays(obj)
  reductions <- Reductions(obj)
  default_assay <- DefaultAssay(obj)
  metadata <- obj[[]]
  n_cells <- ncol(obj)
  n_features <- nrow(obj)
  feature_names <- rownames(obj)
  cell_names <- colnames(obj)
  active_ident_counts <- as.data.frame(table(as.character(Idents(obj))), stringsAsFactors = FALSE)
  colnames(active_ident_counts) <- c("identity", "n_cells")
} else if (inherits(obj, "SingleCellExperiment")) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
  })
  assays <- assayNames(obj)
  reductions <- reducedDimNames(obj)
  metadata <- as.data.frame(colData(obj))
  n_cells <- ncol(obj)
  n_features <- nrow(obj)
  feature_names <- rownames(obj)
  cell_names <- colnames(obj)
} else if (inherits(obj, "SummarizedExperiment")) {
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
  })
  assays <- assayNames(obj)
  reductions <- character()
  metadata <- as.data.frame(colData(obj))
  n_cells <- ncol(obj)
  n_features <- nrow(obj)
  feature_names <- rownames(obj)
  cell_names <- colnames(obj)
} else if (is.list(obj)) {
  metadata <- data.frame()
  assays <- names(obj)
  n_cells <- NA_integer_
  n_features <- NA_integer_
}

summary_rows <- data.frame(
  dataset = "GSE325391",
  object_class = obj_class,
  file = file.path("External_Data/GEO/GSE325391", basename(rds_path)),
  file_size_bytes = file.info(rds_path)$size,
  is_seurat = is_seurat,
  assays = paste(assays, collapse = ";"),
  default_assay = default_assay,
  reductions = paste(reductions, collapse = ";"),
  n_features = n_features,
  n_cells = n_cells,
  n_metadata_columns = ncol(metadata),
  stringsAsFactors = FALSE
)
write.table(summary_rows, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

if (ncol(metadata) > 0) {
  profile <- do.call(rbind, lapply(colnames(metadata), function(col) {
    x <- metadata[[col]]
    x_chr <- as.character(x)
    top <- sort(table(x_chr, useNA = "ifany"), decreasing = TRUE)
    data.frame(
      column = col,
      class = paste(class(x), collapse = ";"),
      n_unique = length(unique(x_chr[!is.na(x_chr)])),
      n_missing = sum(is.na(x)),
      top_values = paste(utils::head(paste(names(top), as.integer(top), sep = "="), 12), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  write.table(profile, out_meta_profile, sep = "\t", row.names = FALSE, quote = FALSE)

  candidate_cols <- grep("sample|orig|group|condition|diagn|disease|batch|run|cell|type|cluster|annotation|ident|age|sex", colnames(metadata), ignore.case = TRUE, value = TRUE)
  group_rows <- list()
  for (col in candidate_cols) {
    values <- as.character(metadata[[col]])
    tab <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
    keep <- utils::head(tab, 100)
    group_rows[[col]] <- data.frame(
      column = col,
      value = names(keep),
      n_cells = as.integer(keep),
      stringsAsFactors = FALSE
    )
  }
  group_df <- if (length(group_rows)) do.call(rbind, group_rows) else data.frame()
  write.table(group_df, out_group_summary, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  write.table(data.frame(), out_meta_profile, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(), out_group_summary, sep = "\t", row.names = FALSE, quote = FALSE)
}

feature_norm <- norm_gene(feature_names)
gene_probe <- data.frame(
  gene = marker_genes,
  present = marker_genes %in% feature_norm,
  matched_feature = vapply(marker_genes, function(g) {
    idx <- which(feature_norm == g)
    if (length(idx)) feature_names[idx[1]] else ""
  }, character(1)),
  stringsAsFactors = FALSE
)
write.table(gene_probe, out_gene_probe, sep = "\t", row.names = FALSE, quote = FALSE)

lines <- c(
  "# GSE325391 Adult RDS Inspection",
  "",
  "Date inspected: 2026-06-21",
  "",
  "## Object",
  "",
  paste0("- Class: `", obj_class, "`"),
  paste0("- Cells/nuclei: ", n_cells),
  paste0("- Features: ", n_features),
  paste0("- Assays: `", paste(assays, collapse = ";"), "`"),
  paste0("- Default assay: `", default_assay, "`"),
  paste0("- Reductions: `", paste(reductions, collapse = ";"), "`"),
  paste0("- Metadata columns: ", ncol(metadata)),
  "",
  "## Marker Probe",
  "",
  paste0("- Marker panel genes present: ", sum(gene_probe$present), " / ", nrow(gene_probe)),
  "",
  "## Outputs",
  "",
  paste0("- Summary: `", sub(paste0("^", root, "/"), "", out_summary), "`"),
  paste0("- Metadata profile: `", sub(paste0("^", root, "/"), "", out_meta_profile), "`"),
  paste0("- Metadata group summary: `", sub(paste0("^", root, "/"), "", out_group_summary), "`"),
  paste0("- Marker gene presence: `", sub(paste0("^", root, "/"), "", out_gene_probe), "`"),
  ""
)
writeLines(lines, out_md)

cat("Wrote", out_summary, "\n")
cat("Wrote", out_meta_profile, "\n")
cat("Wrote", out_group_summary, "\n")
cat("Wrote", out_gene_probe, "\n")
cat("Wrote", out_md, "\n")
cat("cells=", n_cells, " features=", n_features, " metadata_cols=", ncol(metadata), "\n", sep = "")
