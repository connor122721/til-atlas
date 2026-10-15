# TIL Atlas

A static GitHub Pages site about tumor-infiltrating lymphocytes (TILs) and TIL therapy, built entirely from **real, public data**.
Nothing on the site is simulated data; the one model (the clone simulator) starts from measured repertoires and says so.

| Page | What it shows | Real data behind it |
|---|---|---|
| **Trial tracker** (`trials.html`) | Every TIL / neoantigen TCR-T trial, tagged, mapped, filterable, CSV export, weekly changelog | ClinicalTrials.gov API v2, re-pulled every Monday by GitHub Actions |
| **Single-cell TILs** (`lab.html`) | CD8⁺ T cell states, clones followed across tumor / normal tissue / blood, clone-expansion patterns, tumor B cell isotypes and clonality | Wu et al. 2020 (141,623 T cells, scRNA + TCR); Maynard et al. 2020 (1,506 BCRs, 3,505 TCRs) |
| **Bulk RNA** (`bulk.html`) | CD8, stem-like, tumor-reactive, TLS, B cell and plasma cell signatures vs response; forest plot, Kaplan–Meier, AUC heatmap | 9 immunotherapy trials, 918 patients (CRI iAtlas via cBioPortal datahub, ODbL) |
| **Clone simulator** (`simulator.html`) | Infuse a real tumor repertoire into a stochastic model; virtual trial over all 14 repertoires next to published TIL trial results | Starting repertoires from Wu et al. 2020; benchmarks from Rohaan 2022 and C-144-01 |
| **Spatial** (`spatial.html`) | Tissue maps, neighborhood enrichment, CD8→tumor distance by response | Empty until a public spatial cohort is converted on your machine (recipes in `analysis/spatial/README.md`) |
| **Datasets** (`datasets.html`) | Curated catalog of public TIL cohorts (single-cell, TCR/BCR, bulk, spatial) + weekly "new on GEO" watch | `site/data/catalog.json`, `site/data/geo_watch.json` |
| **Learn** (`learn.html`) | Exercises: APIs, signature scoring, clone metrics, bulk signatures, spatial statistics, model fitting | — |

No build step, no server, no framework: plain HTML + D3. Heavy lifting happens in `analysis/` and lands in `site/data/` as small JSON.

## Put it on GitHub (first time)

You need a GitHub account (yours: `connor122721`) and `git`. Two ways:

**A. With the GitHub CLI** (fastest)

```bash
unzip til-atlas.zip && cd til-atlas
gh auth login                                   # once per machine; choose GitHub.com → HTTPS → browser
gh repo create til-atlas --public --source . --push
```

**B. With the website + git**

1. Go to <https://github.com/new>, name it `til-atlas`, keep it **Public**, and **don't** add a README, .gitignore or licence (the repo already has them). Click *Create repository*.
2. In a terminal:
   ```bash
   unzip til-atlas.zip && cd til-atlas
   git remote add origin https://github.com/connor122721/til-atlas.git
   git push -u origin main
   ```
   When asked for a password, paste a personal access token (github.com → Settings → Developer settings → Personal access tokens → Fine-grained, with *Contents: read and write* on this repo).

**Then, either way:**

3. Repo → **Settings → Pages** → *Source*: **GitHub Actions**.
4. Repo → **Settings → Actions → General** → *Workflow permissions*: **Read and write** (lets the weekly job commit new trial data).
5. Repo → **Actions** → **Scrape ClinicalTrials.gov** → **Run workflow**. About a minute later the site is live at `https://connor122721.github.io/til-atlas/`, with real trial data.
6. Optional: Settings → Secrets and variables → Actions → add `NCBI_EMAIL` (and `NCBI_API_KEY`) so the weekly GEO search identifies you to NCBI.

To move it into the lab's organisation later: repo → Settings → *Transfer ownership* → `SmithViralLabLouisville`.

Preview locally any time: `python -m http.server -d site 8000` and open <http://localhost:8000>.

## Where the real data come from

| Data | Source | How it gets into `site/data/` |
|---|---|---|
| Trials | ClinicalTrials.gov API v2 | `scraper/fetch_trials.py`, weekly Action |
| New GEO series | NCBI E-utilities | `scraper/fetch_geo.py`, weekly Action |
| Wu 2020 scRNA + TCR (GSE139555) | processed file on scirpy's GitHub release `d0.1.0` | `python analysis/datasets/wu2020/prepare.py` |
| Maynard 2020 BCR/TCR (PRJNA591860) | processed file on scirpy's GitHub release `d0.1.0` | `python analysis/datasets/maynard2020/prepare.py` |
| 9 checkpoint trials, bulk RNA-seq | cBioPortal datahub (Git LFS), iAtlas-harmonized | `python analysis/bulk/fetch_iatlas.py` |

All three analysis scripts download from GitHub, so they run anywhere (laptop, cluster, Actions) with `numpy pandas scipy h5py`.

Cohorts that live on GEO / Zenodo / Dryad need a machine with general internet access (your cluster), then commit the output:

```bash
python analysis/bulk/geo_cohort.py --expr <GSE100797 matrix> --clinical lauss2017.csv --id lauss2017 \
    --label "Lauss 2017" --cancer Melanoma --therapy "TIL ACT"          # TIL-therapy bulk RNA-seq
sbatch --mem=64G --wrap "bash analysis/datasets/sadefeldman2018/run.sh"  # melanoma ICB scRNA-seq (R)
bash analysis/datasets/chiffelle2024/run.sh                              # Lausanne TIL-ACT (fill the sample sheet)
# spatial: analysis/spatial/README.md
git add site/data && git commit -m "data: add <cohort>" && git push
```

Krishna et al. 2020 has no public deposit, so it is listed in the catalog for reference only.

## Repository layout

```
scraper/                fetch_trials.py, fetch_geo.py, config.json, tags.json, tests/
analysis/
  datasets/wu2020/      prepare.py  → Single-cell tab + simulator repertoires
  datasets/maynard2020/ prepare.py  → B/T cell receptor panel
  bulk/                 fetch_iatlas.py, geo_cohort.py, signatures.json, tests/
  datasets/<cohort>/    run.sh for GEO cohorts (Sade-Feldman, Chiffelle/Barras, Huang) — R/Seurat
  spatial/              export_spatial.py, adapters (Xenium, CosMx, R objects), tests, README
  lib/                  h5ad_lite.py (h5py-only AnnData reader), geo_10x.R
  signatures/           stem_like.txt, terminal.txt (starter panels)
site/                   the website (deployed as-is); data/ holds the real-data JSON
.github/workflows/      scrape.yml (weekly), pages.yml (deploy)
```

## Tests

```bash
python -m unittest discover -s scraper/tests            # trials + GEO watch (stdlib only)
python -m unittest discover -s analysis/bulk/tests      # AUC + cohort scoring
python -m unittest discover -s analysis/spatial/tests   # spatial adapters + exporter (needs h5py)
```

## Licences and attribution

Site code: MIT. Derived data from the cBioPortal datahub / CRI iAtlas is shared under ODbL 1.0 with attribution, as that licence requires.
Wu 2020 and Maynard 2020 derived summaries cite the original GEO/SRA accessions; please cite the original papers when reusing them.
