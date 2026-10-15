#!/usr/bin/env bash
# Download GEO supplementary files for one or more series into data/raw/<GSE>/.
#   bash analysis/00_download_geo.sh GSE222750 GSE221253
#
# GEO lays series out as /geo/series/GSE<prefix>nnn/<GSE>/suppl/, where <prefix> is the
# accession with its last three digits replaced by "nnn" (GSE222750 -> GSE222nnn).
# Afterwards, LOOK at what arrived (file names, formats) before running 01/02: every
# submitter packages things differently, and the readers in the R scripts assume you'll
# point them at the right files.
set -euo pipefail
OUT=${OUT:-data/raw}

for GSE in "$@"; do
  prefix="${GSE:0:${#GSE}-3}nnn"
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${GSE}"
  mkdir -p "${OUT}/${GSE}"
  echo ">> ${GSE}: ${url}"
  # series matrix carries sample-level metadata (characteristics_ch1 = response, timepoint, ...)
  wget -q -c -P "${OUT}/${GSE}" "${url}/matrix/${GSE}_series_matrix.txt.gz" || echo "   (no series matrix)"
  # mirror everything under suppl/
  wget -q -c -r -np -nH --cut-dirs=4 -R "index.html*" -P "${OUT}/${GSE}" "${url}/suppl/"
  echo "   files:"; find "${OUT}/${GSE}" -type f -printf "   %8k KB  %p\n" | sort -k3
done

cat <<'EOF'

Next:
  - untar any *_RAW.tar:   tar -xf data/raw/GSEXXXX/suppl/GSEXXXX_RAW.tar -C data/raw/GSEXXXX/suppl/
  - build/obtain a Seurat object of CD8 T cells with per-cell `patient` and a clinical table with `response`
  - see analysis/README.md for the expected inputs of 01_score_signatures.R and 02_clonotypes.R
EOF
