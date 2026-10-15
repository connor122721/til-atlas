# Helpers for turning GEO supplementary files into Seurat objects and AIRR tables.
# source("analysis/lib/geo_10x.R")
#
# GEO submitters package 10x data in a few common ways; these helpers find all of them:
#   <prefix>_matrix.mtx.gz + <prefix>_barcodes.tsv.gz + <prefix>_features.tsv.gz (or genes.tsv.gz)
#   <prefix>_filtered_feature_bc_matrix.h5 (or any *.h5)
#   <prefix>_filtered_contig_annotations.csv.gz   (10x V(D)J, for TCR)
# Workflow: find_10x_samples() -> write_sample_template() -> you fill patient/compartment/response
#           -> build_seurat(sheet) and contigs_to_airr(sheet)
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(Matrix) })

TRIPLET_RX <- "^(.*?)[_.-]?(matrix\\.mtx|barcodes\\.tsv|features\\.tsv|genes\\.tsv)(\\.gz)?$"

find_10x_samples <- function(dir) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  out <- list()
  # triplets
  trip <- f[grepl(TRIPLET_RX, basename(f))]
  if (length(trip)) {
    key <- file.path(dirname(trip), sub(TRIPLET_RX, "\\1", basename(trip)))
    for (k in unique(key)) {
      g <- trip[key == k]
      get <- function(rx) g[grepl(rx, basename(g))][1]
      out[[length(out) + 1]] <- data.table(prefix = basename(k), type = "mtx",
        matrix = get("matrix\\.mtx"), barcodes = get("barcodes\\.tsv"), features = get("(features|genes)\\.tsv"))
    }
  }
  for (h in f[grepl("\\.h5$", f)])
    out[[length(out) + 1]] <- data.table(prefix = sub("(_filtered_feature_bc_matrix|_raw_feature_bc_matrix)?\\.h5$", "", basename(h)),
                                         type = "h5", matrix = h, barcodes = NA, features = NA)
  res <- rbindlist(out, fill = TRUE)
  if (!nrow(res)) stop("No 10x matrices found under ", dir, " (untar *_RAW.tar first?)")
  contigs <- f[grepl("contig_annotations\\.csv(\\.gz)?$", f)]
  res[, tcr := vapply(prefix, function(p) { m <- contigs[startsWith(basename(contigs), p)]; if (length(m)) m[1] else NA_character_ }, "")]
  res[]
}

#' Write a sample sheet for you to complete. Columns patient/compartment/response are blank.
write_sample_template <- function(samples, path) {
  if (file.exists(path)) { message("keeping existing ", path); return(invisible(fread(path))) }
  tmpl <- samples[, .(prefix, type, has_tcr = !is.na(tcr), patient = "", compartment = "", response = "")]
  fwrite(tmpl, path)
  message("\nWrote ", path, " with ", nrow(tmpl), " samples.\n",
          "Fill in patient, compartment (tumor_pre / product / blood_d30 / tumor_post) and response\n",
          "(Responder / Non-responder) from the GEO sample titles and the paper's tables, then re-run.\n")
  invisible(tmpl)
}

read_sheet <- function(path) {
  s <- fread(path, colClasses = "character")
  todo <- s[patient == "" | compartment == "" | response == ""]
  if (nrow(todo)) stop(nrow(todo), " rows of ", path, " still need patient/compartment/response: ",
                       paste(head(todo$prefix, 5), collapse = ", "))
  bad <- setdiff(unique(s$response), c("Responder", "Non-responder"))
  if (length(bad)) stop("response must be Responder or Non-responder; found: ", paste(bad, collapse = ", "))
  s
}

read_10x_sample <- function(row) {
  m <- if (row$type == "h5") Read10X_h5(row$matrix) else {
    d <- tempfile("tenx_"); dir.create(d)
    file.symlink(normalizePath(row$matrix), file.path(d, paste0("matrix.mtx", if (grepl("gz$", row$matrix)) ".gz" else "")))
    file.symlink(normalizePath(row$barcodes), file.path(d, paste0("barcodes.tsv", if (grepl("gz$", row$barcodes)) ".gz" else "")))
    feat <- if (grepl("genes\\.tsv", row$features)) "genes.tsv" else "features.tsv"
    file.symlink(normalizePath(row$features), file.path(d, paste0(feat, if (grepl("gz$", row$features)) ".gz" else "")))
    Read10X(d)
  }
  if (is.list(m)) m <- m[["Gene Expression"]] %||% m[[1]]   # multi-modal (e.g. CITE-seq): keep RNA
  m
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Merge every sample in the sheet (optionally only some compartments) into one Seurat object.
build_seurat <- function(sheet, samples, compartments = NULL, min_features = 200, max_mt = 20) {
  rows <- merge(samples, sheet[, .(prefix, patient, compartment, response)], by = "prefix")
  if (!is.null(compartments)) rows <- rows[compartment %in% compartments]
  objs <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i]; m <- read_10x_sample(r)
    o <- CreateSeuratObject(m, project = r$prefix, min.cells = 3, min.features = min_features)
    o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
    o <- subset(o, percent.mt < max_mt)
    o$sample <- r$prefix; o$patient <- r$patient; o$compartment <- r$compartment; o$response <- r$response
    message(sprintf("  %-40s %6d cells  (%s, %s)", r$prefix, ncol(o), r$patient, r$compartment))
    RenameCells(o, add.cell.id = r$prefix)
  })
  obj <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
  obj
}

#' 10x filtered_contig_annotations -> AIRR-style TSV (one row per cell's productive TRB chain).
contigs_to_airr <- function(sheet, samples, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  rows <- merge(samples[!is.na(tcr)], sheet[, .(prefix, patient, compartment)], by = "prefix")
  ss <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i]; d <- fread(r$tcr)
    d <- d[chain == "TRB" & tolower(as.character(productive)) %in% c("true", "t") & tolower(as.character(is_cell)) %in% c("true", "t")]
    d <- d[, .SD[which.max(umis)], by = barcode]                 # one TRB per cell
    airr <- d[, .(junction_aa = cdr3, v_call = v_gene, duplicate_count = 1L, productive = TRUE, locus = "TRB")]
    f <- file.path(outdir, paste0(r$prefix, ".airr.tsv")); fwrite(airr, f, sep = "\t")
    data.table(file = f, patient = r$patient, compartment = r$compartment)
  }))
  fwrite(ss, file.path(outdir, "tcr_samples.csv"))
  message("wrote ", nrow(ss), " AIRR tables and ", file.path(outdir, "tcr_samples.csv"))
  ss
}
