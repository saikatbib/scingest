# SPDX-License-Identifier: GPL-3.0-only
#
# scIngest Seurat bridge.
#
# This script is licensed GPL-3.0 (see ./LICENSE) because it calls the R
# `Matrix` package, which is GPL (>= 2). The rest of the scIngest repository
# is MIT licensed.
#
# DESIGN CONTRACT — do not violate without re-reading docs 12.2:
#   * This script is invoked as a SUBPROCESS. Arm's-length interface only:
#     a staging directory in, an exit code and one status line out.
#   * It must contain NO scIngest business logic. All naming, uniquification,
#     sanitisation and metadata joining happens in Python, once, so that the
#     .h5ad and .rds outputs are guaranteed identical.
#   * Never replace this with rpy2 or any in-process embedding.
#
# Usage: Rscript build_seurat.R <staging_dir> <output_rds>

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: build_seurat.R <staging_dir> <output_rds>")
staging <- args[[1]]
out     <- args[[2]]

mf <- jsonlite::fromJSON(file.path(staging, "manifest.json"))

# Emit a legacy (v3/v4) Assay rather than Assay5 when requested.
if (identical(mf$object_version, "v4")) {
  options(Seurat.object.assay.version = "v3")
}

counts   <- Matrix::readMM(file.path(staging, "matrix.mtx"))
features <- utils::read.csv(file.path(staging, "features.csv"), stringsAsFactors = FALSE)
barcodes <- utils::read.csv(file.path(staging, "barcodes.csv"), stringsAsFactors = FALSE)

# Names arrive already unique and already sanitised. Do not transform them here:
# any change made on this side would desynchronise the .rds from the .h5ad.
rownames(counts) <- features$name
colnames(counts) <- barcodes$cell_id
counts <- as(counts, "CsparseMatrix")

meta <- utils::read.csv(file.path(staging, "metadata.csv"), row.names = 1,
                        stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(rownames(meta), colnames(counts)))

# min.cells / min.features are pinned to 0. scIngest does not filter (spec L-2).
obj <- CreateSeuratObject(
  counts       = counts,
  meta.data    = meta,
  project      = mf$project,
  min.cells    = 0,
  min.features = 0
)

# Additional modalities (ADT / HTO / CRISPR) as extra assays.
if (!is.null(mf$modalities) && length(mf$modalities) > 0L) {
  for (m in seq_len(nrow(mf$modalities))) {
    key <- mf$modalities$key[[m]]
    dir <- file.path(staging, "modalities", key)
    mc  <- Matrix::readMM(file.path(dir, "matrix.mtx"))
    mfeat <- utils::read.csv(file.path(dir, "features.csv"), stringsAsFactors = FALSE)
    rownames(mc) <- mfeat$name
    colnames(mc) <- barcodes$cell_id
    obj[[toupper(key)]] <- CreateAssayObject(counts = as(mc, "CsparseMatrix"))
  }
}

obj@misc$scingest <- mf$provenance

saveRDS(obj, file = out, compress = TRUE)

# Status line consumed by the Python side for the parity check.
cat(sprintf("OK\t%d\t%d\n", nrow(obj), ncol(obj)))
