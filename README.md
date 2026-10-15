# TIL Atlas

A static GitHub Pages site about tumor-infiltrating lymphocyte (TIL) therapy, built to learn from public data:

| Page | What it does | Data |
|---|---|---|
| **Trial tracker** (`trials.html`) | Every TIL / neoantigen TCR-T trial on ClinicalTrials.gov, auto-tagged, mapped, filterable, CSV export, weekly changelog | `site/data/trials.json`, refreshed weekly by GitHub Actions |
| **Clone simulator** (`simulator.html`) | Stochastic model of 80 infused clones (stem-like → effector → exhausted) vs a tumor; 40-patient virtual trials with waterfall plots | runs in the browser (`site/assets/sim-model.js`) |
| **Replication lab** (`lab.html`) | Krishna et al. 2020 stem-like CD39⁻CD69⁻ comparison per dataset and across datasets, UMAP explorer, clone tracking tumor → product → blood | `site/data/lab/<dataset>/`, exported by `analysis/datasets/*/run.sh` |
| **Spatial** (`spatial.html`) | Tissue maps, neighborhood enrichment (label-permutation z-scores), CD8→tumor distance, comparison by response | `site/data/spatial/<dataset>/`, exported by `analysis/spatial/` |
| **Datasets** (`datasets.html`) | Curated catalog of public TIL / checkpoint single-cell and spatial cohorts + weekly "new on GEO" watch | `site/data/catalog.json` (hand-curated), `site/data/geo_watch.json` (weekly) |
| **Learn** (`learn.html`) | Exercises: APIs, signature scoring, clonotype metrics, spatial statistics, model fitting | — |

No build step, no server, no framework: plain HTML + D3 from a CDN. Everything heavy happens on your cluster and lands here as small JSON.

## Quick start

```bash
# 1. put it on GitHub (repo name is up to you)
gh repo create til-atlas --public --source . --push

# 2. Settings → Pages → Source: "GitHub Actions"
# 3. Actions → "Scrape ClinicalTrials.gov" → Run workflow
#    (first real data pull; Pages redeploys automatically when it finishes)
```

Preview locally any time:

```bash
python -m http.server -d site 8000      # open http://localhost:8000
```

Until the first scrape runs, the tracker and GEO feed show **synthetic fixtures** (IDs `NCT9…` / `GSE9000…`, titles start with `[FIXTURE]`) with a banner saying so. The lab and spatial tabs ship **synthetic demo datasets**, flagged on the page, until you export real ones; real exports become the default view automatically.

Optional: add repository secrets `NCBI_EMAIL` (and `NCBI_API_KEY`) so the GEO watch identifies itself to NCBI.

## Repository layout

```
scraper/
  fetch_trials.py      API v2 client: paginate → normalize → tag → diff → write JSON (stdlib only)
  config.json          search topics (Essie query syntax) — add a topic here
  tags.json            regex rules for indication / combination / engineering tags
  tests/               offline unittest suite + synthetic fixture generator
scraper/fetch_geo.py   weekly NCBI E-utilities search for new GEO series (queries in config.json → "geo")
analysis/
  datasets/<id>/run.sh  one script per public cohort: download → prepare → export (sadefeldman2018, chiffelle2024, huang2022)
  lib/geo_10x.R        find 10x matrices/V(D)J files in GEO downloads, sample-sheet template, Seurat + AIRR builders
  spatial/             export_spatial.py, adapters (Xenium, CosMx, R objects), synthetic tissue, tests, README with recipes
  00_download_geo.sh   pull GEO supplementary files (e.g. GSE222750, GSE221253)
  01_score_signatures.R  CD8 subset → UCell scores → UMAP → cells.json, patients.json
  02_clonotypes.R      AIRR TCR tables → frequencies, persistence, Morisita-Horn → clones.json
  signatures/          starter gene panels (replace with the published DN/DP lists)
  make_synthetic_demo.py  demo lab data in the same schema
  environment.yml, slurm/run_lab.sbatch
site/                  the website (deployed as-is)
  data/                trials.json, history.json, changelog.json (+ demo/, lab/)
  assets/              style.css, common.js, charts.js, per-page scripts, sim-model.js
.github/workflows/
  scrape.yml           weekly cron (Mon 11:17 UTC) + manual trigger; commits data if changed
  pages.yml            deploys site/ on push and after each scrape
```

## How the weekly scrape works

1. `scrape.yml` runs the offline tests, then `python scraper/fetch_trials.py`.
2. For each topic in `config.json` the script calls `GET https://clinicaltrials.gov/api/v2/studies` with `pageSize=200`, following `nextPageToken` until it's gone (backoff on 429/5xx).
3. Each study is flattened (`normalize`) to ~20 fields and tagged with the regexes in `tags.json`. Trials found by several topics are merged.
4. The new set is diffed against last week's `trials.json`: **new**, **status changes**, **removed** → prepended to `changelog.json` (two years kept). Totals by status go to `history.json`.
5. If anything changed, the bot commits; `pages.yml` sees the completed run and redeploys.

Storage stays small: a few hundred trials is well under 1 MB of JSON, and git history holds every past week.

Useful commands:

```bash
python scraper/fetch_trials.py --dry-run                 # live fetch, print counts, write nothing
python -m unittest discover -s scraper/tests -v          # offline tests (trials + GEO watch)
python -m unittest discover -s analysis/spatial/tests -v # spatial adapters + exporter (needs h5py)
python scraper/tests/make_fixture.py /tmp/fx && \
  python scraper/fetch_trials.py --fixture /tmp/fx --out /tmp/out   # full offline run
```

## Loading real data

Single-cell (replication lab):

```bash
conda env create -f analysis/environment.yml && conda activate til-atlas
sbatch --mem=64G --wrap "bash analysis/datasets/sadefeldman2018/run.sh"   # public, no manual step
bash analysis/datasets/chiffelle2024/run.sh   # 1st run writes data/processed/chiffelle2024/samples.csv
#   fill patient / compartment / response, then run again
bash analysis/datasets/huang2022/run.sh       # same pattern
git add site/data/lab && git commit -m "lab: real exports" && git push
```

Spatial: follow `analysis/spatial/README.md` (Xenium melanoma, melanoma IMC anti-PD-1, cSCC and HNSCC spatial cohorts).

| Dataset | Access | Status of the pipeline |
|---|---|---|
| Sade-Feldman 2018, GSE120575 | open | ready; reads GEO's per-cell patient/response table |
| Chiffelle 2024 / Barras 2024, GSE229861 + GSE222750 | check GSE229861 is released | ready; needs the sample sheet filled |
| Huang 2022 cervical TIL, GSE190075 | check release | ready; needs the sample sheet filled |
| Xenium melanoma (10x) | open, CC BY 4.0 | ready; tested on mock bundles |
| Melanoma IMC anti-PD-1, Zenodo 6838169 | open | inspect the object, then export (column names unknown in advance) |
| cSCC spatial, Dryad s4mw6m9jh | open | inspect the Seurat object, then export |
| HNSCC CosMx, Dryad n2z34tn7x | open | needs a FOV → core/response map from the deposit README |
| Krishna 2020, Moffitt GeoMx | accession not confirmed | catalog only |

Inputs the R scripts expect (details at the top of each script):

* `01_score_signatures.R`: a Seurat object with per-cell `patient` (and optionally a sample-type column to keep only infusion-product cells) plus `clinical.csv` with `patient,response[,persistence30]`.
* `02_clonotypes.R`: AIRR TSVs and a sample sheet `file,patient,compartment` with compartments `tumor_pre`, `product`, `blood_d30`, `tumor_post`.

Each submitter packages GEO files differently, so expect to write a few lines that turn the downloaded files into those inputs. Check every dataset's data-use terms before publishing derived data; never commit patient-level data from a restricted cohort.

## The simulator model

Per clone, three compartments (stem-like S, effector E, exhausted X), tau-leaped Poisson events every 0.25 days. Stem-like cells self-renew under homeostatic/IL-2 drive and divide asymmetrically under antigen (S → S + E); effectors kill but die within weeks; exhaustion pulls effectors into X. The tumor grows logistically and is killed by avidity-weighted effectors. Equations are on the page and in the header of `sim-model.js`. Parameters are hand-set for qualitative behavior, not fitted; the module runs unchanged in Node:

```js
const { simulate } = require("./site/assets/sim-model.js");
const r = simulate({ stemFrac: 0.05, reactiveFrac: 0.4, seed: 1 });
console.log(r.readouts.best, r.readouts.final, r.readouts.relapse);
```

## Ideas to extend

* Add CAR-T solid-tumor trials as a topic and an antigen tag group (`config.json`, `tags.json`).
* Fit simulator parameters to real persistence curves with approximate Bayesian computation.
* A spatial tab: neighborhood-enrichment heatmaps from Xenium/CODEX/mIF summaries.
* A publications feed: weekly PubMed E-utilities query for TIL papers, same diff-and-commit pattern.
