#!/usr/bin/env bash
# Cervical cancer TIL phase I (Huang et al., JCI 2022): scRNA-seq of 8 TIL infusion products.
# Response per patient is in the paper's tables (not in GEO): fill it into the sample sheet.
#   bash analysis/datasets/huang2022/run.sh
set -euo pipefail
ID=huang2022; GSE=GSE190075
RAW=data/raw; WORK=data/processed/$ID; mkdir -p "$WORK"
[ -d "$RAW/$GSE" ] || OUT=$RAW bash analysis/00_download_geo.sh "$GSE"
for t in "$RAW/$GSE"/*_RAW.tar "$RAW/$GSE"/suppl/*_RAW.tar; do [ -f "$t" ] && tar -xf "$t" -C "$(dirname "$t")" && rm "$t" || true; done

Rscript - "$RAW/$GSE" "$WORK" <<'EOF'
args <- commandArgs(TRUE); source("analysis/lib/geo_10x.R")
samples <- find_10x_samples(args[1]); print(samples[, .(prefix, type)])
sheet_path <- file.path(args[2], "samples.csv")
if (!file.exists(sheet_path)) { write_sample_template(samples, sheet_path); quit(status = 0) }
sheet <- read_sheet(sheet_path)
saveRDS(build_seurat(sheet, samples), file.path(args[2], "scrna.rds"))
fwrite(unique(sheet[, .(patient, response)]), file.path(args[2], "clinical.csv"))
EOF
[ -s "$WORK/clinical.csv" ] || { echo "Fill $WORK/samples.csv (compartment = product), then re-run."; exit 0; }

Rscript analysis/01_score_signatures.R --seurat "$WORK/scrna.rds" --meta "$WORK/clinical.csv" \
  --dataset-id $ID --label "Cervical TIL · products" --study "Huang et al. 2022" \
  --therapy TIL-ACT --sample-label "infusion product"
echo "Done: site/data/lab/$ID/"
