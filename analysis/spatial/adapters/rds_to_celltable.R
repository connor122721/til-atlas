#!/usr/bin/env Rscript
# R objects (Seurat, SingleCellExperiment, SpatialExperiment, or a data.frame) -> standard cell table.
# Used for deposits that ship an analyzed object rather than raw instrument output, e.g.
#   cSCC neoadjuvant PD-1 (Dryad doi:10.5061/dryad.s4mw6m9jh): Seurat_Object_cSCC.RDS
#   melanoma IMC anti-PD-1 (Zenodo 6838169): objects inside Melanoma_code_data.zip
#
# Step 1, look inside (column names differ for every deposit):
#   Rscript analysis/spatial/adapters/rds_to_celltable.R --inspect obj.rds
# Step 2, export with the columns you found:
#   Rscript analysis/spatial/adapters/rds_to_celltable.R obj.rds out.csv.gz \
#     --sample-col fov --celltype-col cell_type --patient-col patient \
#     --response-col response --responders "MPR,pCR" [--x-col x --y-col y] [--um-per-unit 1]
# Then run analysis/spatial/export_spatial.py on out.csv.gz. Map your cell-type names so the CD8
# population is called exactly "CD8 T" and tumor cells "Tumor" (--rename "CD8+ T=CD8 T,Melanoma=Tumor"),
# or pass --focal/--target to export_spatial.py.
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--inspect", action = "store_true", default = FALSE),
  make_option("--sample-col", dest = "sample"), make_option("--celltype-col", dest = "ctype"),
  make_option("--patient-col", dest = "patient", default = NA), make_option("--response-col", dest = "resp", default = NA),
  make_option("--responders", default = NA, help = "comma-separated values of --response-col that mean Responder"),
  make_option("--x-col", dest = "x", default = NA), make_option("--y-col", dest = "y", default = NA),
  make_option("--um-per-unit", dest = "scale", type = "double", default = 1),
  make_option("--rename", default = NA, help = "old=new,old=new for cell types")
)), positional_arguments = c(1, 2))
obj <- readRDS(opt$args[1]); o <- opt$options

meta_of <- function(x) {
  if (inherits(x, "Seurat")) return(as.data.frame(x[[]]))
  if (inherits(x, "SingleCellExperiment")) return(as.data.frame(SummarizedExperiment::colData(x)))
  as.data.frame(x)
}
coords_of <- function(x, md) {
  if (!is.na(o$x)) return(as.matrix(md[, c(o$x, o$y)]))
  if (inherits(x, "SpatialExperiment")) return(SpatialExperiment::spatialCoords(x))
  if (inherits(x, "Seurat") && length(SeuratObject::Images(x))) {
    cc <- do.call(rbind, lapply(SeuratObject::Images(x), function(i) SeuratObject::GetTissueCoordinates(x, image = i)))
    rownames(cc) <- if ("cell" %in% names(cc)) cc$cell else rownames(cc)
    return(as.matrix(cc[colnames(x), sapply(cc, is.numeric)][, 1:2]))
  }
  guess <- grep("^(x|y|center_?x|center_?y|x_centroid|y_centroid|pos_?x|pos_?y)", names(md), ignore.case = TRUE, value = TRUE)
  if (length(guess) >= 2) return(as.matrix(md[, guess[1:2]]))
  stop("no coordinates found: pass --x-col/--y-col (see --inspect)")
}

md <- meta_of(obj)
if (o$inspect) {
  cat("class:", class(obj)[1], " cells:", nrow(md), "\n\nmetadata columns (first values):\n")
  for (n in names(md)) cat(sprintf("  %-30s %s\n", n, paste(head(unique(md[[n]]), 6), collapse = " | ")))
  if (inherits(obj, "Seurat")) cat("\nimages:", paste(SeuratObject::Images(obj), collapse = ", "), "\n")
  if (inherits(obj, "SpatialExperiment")) cat("\nspatialCoords columns:", paste(colnames(SpatialExperiment::spatialCoords(obj)), collapse = ", "), "\n")
  quit(status = 0)
}
stopifnot(length(opt$args) == 2, !is.null(o$sample), !is.null(o$ctype))
xy <- coords_of(obj, md) * o$scale
ct <- as.character(md[[o$ctype]])
if (!is.na(o$rename)) for (p in strsplit(strsplit(o$rename, ",")[[1]], "=")) ct[ct == p[1]] <- p[2]
resp <- if (is.na(o$resp)) NA_character_ else {
  r <- as.character(md[[o$resp]]); ok <- strsplit(o$responders, ",")[[1]]
  ifelse(is.na(r) | r == "", NA_character_, ifelse(r %in% ok, "Responder", "Non-responder"))
}
out <- data.table(cell_id = rownames(md), sample = as.character(md[[o$sample]]), x = xy[, 1], y = xy[, 2], cell_type = ct,
                  patient = if (is.na(o$patient)) as.character(md[[o$sample]]) else as.character(md[[o$patient]]), response = resp)
fwrite(out, opt$args[2])
message(sprintf("wrote %d cells, %d samples, %d cell types -> %s", nrow(out), uniqueN(out$sample), uniqueN(out$cell_type), opt$args[2]))
