#!/usr/bin/env bash
# Lausanne TIL-ACT trial (NCT03475134): Chiffelle et al. 2024 + Barras et al. 2024.
# Products, baseline and post-ACT tumors with scRNA/scTCR; bulk RNA from Barras.
#
#   bash analysis/datasets/chiffelle2024/run.sh            # from the repo root, on the cluster
#
# Step 1 downloads and writes a sample sheet; you fill it in from the GEO sample titles
# and the papers' patient tables; step 2 (re-run) builds the lab exports.
set -euo pipefail
ID=chiffelle2024
RAW=data/raw; WORK=data/processed/$ID; mkdir -p "$WORK"

# GSE229861 was listed in the preprint as "public at publication" — check it's released first.
for GSE in GSE229861 GSE222750; do
  [ -d "$RAW/$GSE" ] || OUT=$RAW bash analysis/00_download_geo.sh "$GSE"
  for t in "$RAW/$GSE"/*_RAW.tar "$RAW/$GSE"/suppl/*_RAW.tar; do
    [ -f "$t" ] && tar -xf "$t" -C "$(dirname "$t")" && rm "$t" || true
  done
done

Rscript analysis/datasets/$ID/prepare.R "$RAW" "$WORK"
[ -s "$WORK/scrna.rds" ] || { echo "Fill $WORK/samples.csv, then re-run this script."; exit 0; }

# products only: the Krishna comparison is about what was infused
Rscript analysis/01_score_signatures.R --seurat "$WORK/scrna.rds" --meta "$WORK/clinical.csv" \
  --sample-col compartment --sample-value product \
  --dataset-id $ID --label "Lausanne TIL-ACT · products" --study "Chiffelle et al. 2024" \
  --therapy TIL-ACT --sample-label "infusion product"

if [ -s "$WORK/airr/tcr_samples.csv" ]; then
  Rscript analysis/02_clonotypes.R --samples "$WORK/airr/tcr_samples.csv" --clinical "$WORK/clinical.csv" --dataset-id $ID
fi
echo "Done: site/data/lab/$ID/ — commit it and push."
