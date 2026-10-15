#!/usr/bin/env bash
# Sade-Feldman et al. 2018, melanoma under checkpoint blockade (GSE120575, public, no sample sheet needed).
# Same question as Krishna, different therapy: do TCF7+ stem-like CD8 cells in the tumor track with response?
#   sbatch --mem=64G --wrap "bash analysis/datasets/sadefeldman2018/run.sh"
set -euo pipefail
ID=sadefeldman2018; GSE=GSE120575
RAW=data/raw; WORK=data/processed/$ID
[ -d "$RAW/$GSE" ] || OUT=$RAW bash analysis/00_download_geo.sh "$GSE"
Rscript analysis/datasets/$ID/prepare.R "$RAW/$GSE" "$WORK" Pre
Rscript analysis/01_score_signatures.R --seurat "$WORK/scrna.rds" --meta "$WORK/clinical.csv" \
  --dataset-id $ID --label "Melanoma ICB · baseline tumors" --study "Sade-Feldman et al. 2018" \
  --therapy "Checkpoint blockade" --sample-label "tumor biopsy before treatment"
echo "Done: site/data/lab/$ID/"
