#!/usr/bin/env Rscript
# Score stem-like vs terminal CD8 states per cell and export the lab's cells.json + patients.json.
#
# Inputs
#   --seurat   RDS of a Seurat object (v4 or v5). Needs per-cell metadata columns named by
#              --patient-col (default "patient") and, optionally, --sample-col to keep only
#              infusion-product cells (e.g. --sample-col sample_type --sample-value product).
#   --meta     CSV with one row per patient: patient,response[,persistence30]
#              response must be "Responder" or "Non-responder" (map CR/PR vs SD/PD yourself).
#   --signatures  directory of one-gene-per-line .txt files (default analysis/signatures)
#   --out      output directory (default site/data/lab)
#
# Usage
#   Rscript analysis/01_score_signatures.R --seurat products.rds --meta clinical.csv \
#     --sample-col sample_type --sample-value product
#
# What it does
#   1. keeps CD8 T cells (CD8A/CD8B > 0, CD4 == 0 unless a `celltype` column says otherwise)
#   2. UCell scores for every gene set in --signatures (stem_like.txt -> "stem", terminal.txt -> "term")
#   3. re-runs PCA/UMAP/clustering on the CD8 subset and names clusters by their top scores
#   4. per patient: share of cells in the stem-like cluster(s) and the terminal cluster(s)
#   5. writes compact columnar JSON (<= --max-cells cells, sampled evenly across patients)
suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(UCell); library(jsonlite); library(dplyr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--seurat"), make_option("--meta"),
  make_option("--signatures", default = "analysis/signatures"),
  make_option("--out", default = "site/data/lab"),
  make_option("--patient-col", dest = "patient_col", default = "patient"),
  make_option("--sample-col", dest = "sample_col", default = NA),
  make_option("--sample-value", dest = "sample_value", default = NA),
  make_option("--max-cells", dest = "max_cells", type = "integer", default = 20000),
  make_option("--resolution", type = "double", default = 0.6),
  make_option("--seed", type = "integer", default = 2020)
)))
stopifnot(!is.null(opt$seurat), !is.null(opt$meta))
set.seed(opt$seed)
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

obj  <- readRDS(opt$seurat)
clin <- read.csv(opt$meta, stringsAsFactors = FALSE)
stopifnot(all(c("patient", "response") %in% names(clin)))
obj$patient <- as.character(obj[[opt$patient_col]][, 1])

# ---- 1. subset: infusion product, CD8 T cells
if (!is.na(opt$sample_col)) obj <- obj[, obj[[opt$sample_col]][, 1] == opt$sample_value]
if ("celltype" %in% colnames(obj[[]])) {
  obj <- obj[, grepl("CD8", obj$celltype)]
} else {
  expr <- FetchData(obj, vars = intersect(c("CD8A", "CD8B", "CD4"), rownames(obj)))
  keep <- (rowSums(expr[, intersect(c("CD8A", "CD8B"), colnames(expr)), drop = FALSE]) > 0) &
          (if ("CD4" %in% colnames(expr)) expr$CD4 == 0 else TRUE)
  obj <- obj[, keep]
}
obj <- obj[, obj$patient %in% clin$patient]
message(sprintf("CD8 cells kept: %d from %d patients", ncol(obj), length(unique(obj$patient))))

# ---- 2. signature scores (UCell is rank-based, so scores don't depend on the other cells)
read_set <- function(f) { g <- trimws(readLines(f)); g[nzchar(g) & !startsWith(g, "#")] }
sets <- lapply(list.files(opt$signatures, "\\.txt$", full.names = TRUE), read_set)
names(sets) <- sub("\\.txt$", "", basename(list.files(opt$signatures, "\\.txt$")))
missing <- lapply(sets, setdiff, rownames(obj))
for (s in names(missing)) if (length(missing[[s]])) message("  ", s, ": not in data -> ", paste(missing[[s]], collapse = ", "))
obj <- AddModuleScore_UCell(obj, features = sets, name = "")
stem_col <- "stem_like"; term_col <- "terminal"
stopifnot(stem_col %in% colnames(obj[[]]), term_col %in% colnames(obj[[]]))

# ---- 3. embed + cluster the CD8 subset
obj <- NormalizeData(obj, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
# keep TCR genes out of the embedding so clones don't drive clusters
VariableFeatures(obj) <- grep("^TR[ABDG][VJC]", VariableFeatures(obj), value = TRUE, invert = TRUE)
obj <- ScaleData(obj, verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE) |>
  FindNeighbors(dims = 1:20, verbose = FALSE) |> FindClusters(resolution = opt$resolution, verbose = FALSE) |>
  RunUMAP(dims = 1:20, seed.use = opt$seed, verbose = FALSE)

# name clusters by mean signature score: highest stem -> stem-like, highest terminal -> terminal
cl <- obj[[]] |> group_by(seurat_clusters) |>
  summarise(stem = mean(.data[[stem_col]]), term = mean(.data[[term_col]]), n = n(), .groups = "drop") |>
  mutate(label = case_when(stem == max(stem) ~ "Stem-like (TCF7+ DN)",
                           term == max(term) ~ "Terminal exhausted (CD39+ DP)",
                           TRUE ~ paste("Cluster", seurat_clusters)))
print(cl)
obj$label <- cl$label[match(obj$seurat_clusters, cl$seurat_clusters)]

# ---- 4. per-patient summaries (the Krishna comparison)
pat <- obj[[]] |> group_by(patient) |>
  summarise(nCells = n(), stemFrac = mean(label == "Stem-like (TCF7+ DN)"),
            termFrac = mean(label == "Terminal exhausted (CD39+ DP)"), .groups = "drop") |>
  left_join(clin, by = "patient")
if (!"persistence30" %in% names(pat)) pat$persistence30 <- NA_real_
print(pat)

# ---- 5. export (columnar JSON keeps the file small)
idx <- unlist(lapply(split(seq_len(ncol(obj)), obj$patient), function(i)
  if (length(i) > opt$max_cells / length(unique(obj$patient))) sample(i, floor(opt$max_cells / length(unique(obj$patient)))) else i))
um <- Embeddings(obj, "umap")[idx, ]
labels <- sort(unique(obj$label))
meta <- list(synthetic = FALSE, source = basename(opt$seurat), sample = "infusion product, CD8+ T cells",
             clusters = labels, signatures = list(stem = "signatures/stem_like.txt", term = "signatures/terminal.txt"),
             generated = format(Sys.Date()))
cells <- list(meta = meta, patients = sort(unique(obj$patient)),
              responses = as.list(setNames(clin$response, clin$patient)),
              x = round(um[, 1], 3), y = round(um[, 2], 3),
              stem = round(obj[[stem_col]][idx, 1], 3), term = round(obj[[term_col]][idx, 1], 3),
              cluster = match(obj$label[idx], labels) - 1L, patient = obj$patient[idx])
write_json(cells, file.path(opt$out, "cells.json"), auto_unbox = TRUE, digits = NA)
write_json(list(meta = meta, patients = pat), file.path(opt$out, "patients.json"), auto_unbox = TRUE, digits = 4, pretty = TRUE)
saveRDS(obj, sub("\\.rds$", ".scored.rds", opt$seurat))
message("wrote ", paste(file.path(opt$out, c("cells.json", "patients.json")), collapse = ", "))
