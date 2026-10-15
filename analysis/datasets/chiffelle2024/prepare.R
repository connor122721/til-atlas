#!/usr/bin/env Rscript
# Build scrna.rds, clinical.csv and AIRR tables for the Lausanne TIL-ACT cohort.
#   Rscript analysis/datasets/chiffelle2024/prepare.R data/raw data/processed/chiffelle2024
args <- commandArgs(TRUE); raw <- args[1]; work <- args[2]
source("analysis/lib/geo_10x.R")

samples <- rbind(find_10x_samples(file.path(raw, "GSE229861")),
                 tryCatch(find_10x_samples(file.path(raw, "GSE222750")), error = function(e) NULL), fill = TRUE)
print(samples[, .(prefix, type, tcr = !is.na(tcr))])

sheet_path <- file.path(work, "samples.csv")
if (!file.exists(sheet_path)) { write_sample_template(samples, sheet_path); quit(status = 0) }
sheet <- read_sheet(sheet_path)

obj <- build_seurat(sheet, samples)
saveRDS(obj, file.path(work, "scrna.rds"))
fwrite(unique(sheet[, .(patient, response)]), file.path(work, "clinical.csv"))
if (any(!is.na(samples$tcr))) contigs_to_airr(sheet, samples, file.path(work, "airr"))
