#!/usr/bin/env Rscript
# Track TCRβ clonotypes across samples and export the lab's clones.json.
#
# Input: AIRR-format rearrangement TSVs (MiXCR `exportAirr`, 10x `airr_rearrangement.tsv`,
# Adaptive exports converted with immunarch, ...) with at least
#   junction_aa, v_call, duplicate_count (or umi_count / consensus_count)
# plus a sample sheet mapping file -> patient, compartment:
#   file,patient,compartment
#   tcr/P01_tumor.tsv,P01,tumor_pre
#   tcr/P01_TIL.tsv,P01,product
#   tcr/P01_d30.tsv,P01,blood_d30
#   tcr/P01_post.tsv,P01,tumor_post
# Compartments must use the four names above (missing ones are fine).
#
#   Rscript analysis/02_clonotypes.R --samples samples.csv --clinical clinical.csv \
#     [--phenotype product_phenotype.csv] [--top 25]
#
# --phenotype (optional): clone_id,phenotype,reactive from your scRNA/TCR join
#   (e.g. majority label of the clone's cells in the product: "stem-like" / "terminal").
suppressPackageStartupMessages({ library(optparse); library(data.table); library(jsonlite) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--samples"), make_option("--clinical"), make_option("--phenotype", default = NA),
  make_option("--top", type = "integer", default = 25), make_option("--out", default = "site/data/lab")
)))
COMP <- c("tumor_pre", "product", "blood_d30", "tumor_post")
ss   <- fread(opt$samples); clin <- fread(opt$clinical)
stopifnot(all(ss$compartment %in% COMP))

read_airr <- function(f) {
  d <- fread(f, select = function(n) n %in% c("junction_aa", "v_call", "duplicate_count", "umi_count", "consensus_count", "productive", "locus"))
  if ("locus" %in% names(d)) d <- d[locus %in% c("TRB", "")]
  if ("productive" %in% names(d)) d <- d[productive %in% c(TRUE, "T", "TRUE")]
  cnt <- intersect(c("duplicate_count", "umi_count", "consensus_count"), names(d))[1]
  d[, count := as.numeric(get(cnt))]
  d[, v := sub("\\*.*", "", v_call)]                  # drop allele: TRBV5-1*01 -> TRBV5-1
  d[nzchar(junction_aa), .(count = sum(count)), by = .(clone_id = paste(v, junction_aa, sep = "|"), cdr3 = junction_aa)]
}
rep <- rbindlist(lapply(seq_len(nrow(ss)), function(i) {
  d <- read_airr(ss$file[i]); d[, `:=`(patient = ss$patient[i], compartment = ss$compartment[i])]
}))
rep[, freq := count / sum(count), by = .(patient, compartment)]

# ---- repertoire metrics per sample
metrics <- rep[, .(clones = .N, reads = sum(count),
                   clonality = 1 - (-sum(freq * log(freq))) / log(.N)), by = .(patient, compartment)]
# persistence: share of day-30 blood made of clones present in the product
pers <- rep[, {
  prod <- clone_id[compartment == "product"]
  .(persistence30 = sum(freq[compartment == "blood_d30" & clone_id %in% prod]))
}, by = patient]
# Morisita-Horn overlap between product and day-30 blood
mh <- function(x, y) { X <- sum(x); Y <- sum(y); if (!X || !Y) return(NA_real_)
  2 * sum(x * y) / ((sum(x^2) / X^2 + sum(y^2) / Y^2) * X * Y) }
overlap <- rep[compartment %in% c("product", "blood_d30"),
  { w <- dcast(.SD, clone_id ~ compartment, value.var = "count", fill = 0)
    .(morisitaHorn = if (all(c("product", "blood_d30") %in% names(w))) mh(w$product, w$blood_d30) else NA_real_) },
  by = patient]
summ <- Reduce(function(a, b) merge(a, b, by = "patient", all = TRUE), list(clin, pers, overlap))
print(metrics); print(summ)
fwrite(metrics, file.path(opt$out, "repertoire_metrics.csv"))

# ---- top clones per patient (ranked by product frequency), with frequencies in every compartment
ph <- if (!is.na(opt$phenotype)) fread(opt$phenotype) else data.table(clone_id = character(), phenotype = character(), reactive = logical())
out <- list()
for (p in unique(rep$patient)) {
  w <- dcast(rep[patient == p], clone_id + cdr3 ~ compartment, value.var = "freq", fill = NA)
  for (c in setdiff(COMP, names(w))) w[, (c) := NA_real_]
  rank_col <- if ("product" %in% names(w) && any(!is.na(w$product))) "product" else "tumor_pre"
  w <- w[order(-get(rank_col))][seq_len(min(.N, opt$top))]
  w <- merge(w, ph, by = "clone_id", all.x = TRUE, sort = FALSE)
  out[[p]] <- list(response = clin[patient == p, response][1],
    clones = lapply(seq_len(nrow(w)), function(i) list(
      cdr3 = w$cdr3[i], reactive = isTRUE(w$reactive[i]),
      phenotype = ifelse(is.na(w$phenotype[i]), "terminal", w$phenotype[i]),
      freq = setNames(lapply(COMP, function(c) { v <- w[[c]][i]; if (is.na(v) || v <= 0) NULL else round(log10(v), 3) }), COMP))))
}
write_json(list(meta = list(synthetic = FALSE, compartments = COMP, freqScale = "log10 clone frequency",
                            generated = format(Sys.Date())), patients = out),
           file.path(opt$out, "clones.json"), auto_unbox = TRUE, null = "null", digits = NA)
message("wrote ", file.path(opt$out, "clones.json"))
