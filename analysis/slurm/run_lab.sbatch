#!/bin/bash
#SBATCH --job-name=til-atlas-lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
# Submit from the repo root:  mkdir -p logs && sbatch analysis/slurm/run_lab.sbatch
# Edit the input paths below for your cohort.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate til-atlas

bash analysis/00_download_geo.sh GSE222750 GSE221253

Rscript analysis/01_score_signatures.R \
  --seurat data/processed/products_cd8.rds \
  --meta   data/processed/clinical.csv \
  --sample-col sample_type --sample-value product

Rscript analysis/02_clonotypes.R \
  --samples data/processed/tcr_samples.csv \
  --clinical data/processed/clinical.csv \
  --phenotype data/processed/product_clone_phenotype.csv

echo "Done. Commit site/data/lab/*.json and push."
