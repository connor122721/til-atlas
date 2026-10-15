#!/usr/bin/env Rscript
# Sade-Feldman et al. 2018 (GSE120575): Smart-seq2 TPM of 16,291 CD45+ cells from 48 melanoma
# biopsies taken before ("Pre") or during ("Post") checkpoint blockade. The per-cell table carries
# patient/timepoint ("Pre_P1") and response, so no manual sample sheet is needed.
#
#   Rscript analysis/datasets/sadefeldman2018/prepare.R data/raw/GSE120575 data/processed/sadefeldman2018 [Pre|Post|all]
#
# Memory: the TPM matrix is read densely before sparsifying; ask for ~64 GB.
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(Matrix) })
args <- commandArgs(TRUE); raw <- args[1]; work <- args[2]; tp <- if (length(args) > 2) args[3] else "Pre"
dir.create(work, showWarnings = FALSE, recursive = TRUE)
find1 <- function(rx) { f <- list.files(raw, rx, recursive = TRUE, full.names = TRUE); stopifnot(length(f) == 1); f }

# ---- per-cell annotation: GEO-style table with a few description lines above the header
ann_lines <- readLines(gzfile(find1("patient_ID_single_cells\\.txt\\.gz$")))
hdr <- grep("response", ann_lines, ignore.case = TRUE)[1]
ann <- fread(text = ann_lines[hdr:length(ann_lines)], sep = "\t", header = TRUE, fill = TRUE)
ann <- ann[, which(nzchar(names(ann)) & !grepl("^V\\d+$", names(ann))), with = FALSE]
col <- function(rx) { n <- grep(rx, names(ann), ignore.case = TRUE, value = TRUE); if (!length(n)) stop("no column matching ", rx, " in: ", paste(names(ann), collapse = " | ")); n[1] }
ann <- ann[, .(cell = get(col("^title$")), sample = get(col("pati?n?e?n?t")), response = get(col("response")))]
ann <- ann[nzchar(cell)]
ann[, timepoint := sub("_.*", "", sample)]
ann[, response := fifelse(grepl("^non", response, ignore.case = TRUE), "Non-responder", "Responder")]
print(ann[, .N, by = .(timepoint, response)])

# ---- TPM matrix: first line = cell IDs; an optional second line repeats sample labels
mf <- find1("single_cells_TPM_GEO\\.txt\\.gz$")
head2 <- strsplit(readLines(gzfile(mf), n = 2), "\t")
ids <- head2[[1]]; second <- head2[[2]]
skip <- if (any(grepl("^(Pre|Post)_", second[-1]))) 2 else 1
tpm <- fread(mf, skip = skip, header = FALSE)
genes <- tpm[[1]]; tpm[, 1 := NULL]
if (length(ids) == ncol(tpm) + 1) ids <- ids[-1]          # header row has a blank over the gene column
stopifnot(length(ids) == ncol(tpm))
keep_rows <- rowSums(tpm > 0) > 0 & !duplicated(genes) & nzchar(genes)
m <- as(as.matrix(tpm[keep_rows]), "CsparseMatrix"); rm(tpm); gc()
dimnames(m) <- list(genes[keep_rows], ids)

# ---- subset timepoint, build Seurat. TPM goes in as "counts"; 01_score_signatures.R log-normalizes it.
ann <- ann[cell %in% colnames(m)]
if (tp != "all") ann <- ann[timepoint == tp]
obj <- CreateSeuratObject(m[, ann$cell], project = "GSE120575", min.features = 200)
obj$patient <- ann$sample[match(colnames(obj), ann$cell)]   # unit = one biopsy (e.g. Pre_P1)
obj$timepoint <- ann$timepoint[match(colnames(obj), ann$cell)]
saveRDS(obj, file.path(work, "scrna.rds"))
fwrite(unique(ann[, .(patient = sample, response)]), file.path(work, "clinical.csv"))
message(sprintf("%d cells from %d biopsies (%s) -> %s", ncol(obj), uniqueN(obj$patient), tp, work))
