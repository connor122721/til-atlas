# Spatial pipeline

```
platform files ──adapter──▶ cells.csv.gz ──export_spatial.py──▶ site/data/spatial/<id>/
                            cell_id, sample, x, y (µm),        samples.json  composition, density, CD8→tumor distance
                            cell_type, patient, response        nhood.json    neighborhood z-scores per sample
                                                                cells/*.json  downsampled maps
```

Every dataset ends up as one cell table; `export_spatial.py` does the rest and registers the dataset in
`site/data/spatial/index.json`. Name the CD8 population `CD8 T` and tumor cells `Tumor`, or pass
`--focal` / `--target` with your names.

Tests (no real data needed): `python -m unittest discover -s analysis/spatial/tests`
Python deps: numpy, scipy, pandas (+ h5py for Xenium). R deps for objects: data.table, optparse, Seurat or SpatialExperiment.

## Recipes for the public datasets in the catalog

Download links are on each repository page; very large files are best fetched on the cluster with `wget -c`.

### Xenium melanoma (10x Genomics, CC BY 4.0) — no treatment, a good first map

1. From the [10x dataset page](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin) download the **Xenium Output Bundle**
   (or the tarballs in [Zenodo 17429142](https://zenodo.org/records/17429142)) and unpack into `data/raw/xenium/`.
2. Cell table (marker-based typing from `adapters/markers.json`):
   ```bash
   python analysis/spatial/adapters/xenium.py data/processed/xenium_cells.csv.gz \
       data/raw/xenium/<bundle_dir>:prime5k
   python analysis/spatial/export_spatial.py data/processed/xenium_cells.csv.gz \
       --dataset-id xenium_melanoma --label "Xenium melanoma (10x)" --platform "Xenium Prime 5K" --study "10x Genomics public data"
   ```
   One 112k-cell section exports in well under a minute; the page draws 25,000 of its cells.

### Melanoma IMC, anti-PD-1 (Zenodo 6838169) — 26 patients with response

1. `wget -c "https://zenodo.org/record/6838169/files/Melanoma_code_data.zip?download=1" -O data/raw/imc/Melanoma_code_data.zip && unzip -d data/raw/imc data/raw/imc/Melanoma_code_data.zip`
2. Find the single-cell object(s) (`find data/raw/imc -name "*.rds"`) and the clinical table (`clidata` folder), then look inside:
   ```bash
   Rscript analysis/spatial/adapters/rds_to_celltable.R --inspect data/raw/imc/<object>.rds
   ```
3. Export with the column names you found. If response lives in the clinical table rather than the object,
   join it in R first (one `merge()` on patient ID) and save a data.frame RDS.
   ```bash
   Rscript analysis/spatial/adapters/rds_to_celltable.R data/raw/imc/<object>.rds data/processed/imc_cells.csv.gz \
     --sample-col <image_id> --celltype-col <celltype> --patient-col <patient> \
     --response-col <response> --responders "<values meaning response>" --x-col <x> --y-col <y> \
     --rename "<their CD8 label>=CD8 T,<their tumor label>=Tumor"
   python analysis/spatial/export_spatial.py data/processed/imc_cells.csv.gz --dataset-id imc_melanoma_pd1 \
     --label "Melanoma IMC · anti-PD-1" --platform "Imaging mass cytometry" --therapy "Checkpoint blockade" --radius 20
   ```
   IMC pixels are 1 µm, so coordinates are already in µm (`--um-per-unit 1`). A 20 µm radius suits the smaller IMC images.

### cSCC, neoadjuvant PD-1 (Dryad doi:10.5061/dryad.s4mw6m9jh) — Seurat object

1. Download `Seurat_Object_cSCC.RDS` (2.5 GB) from the Dryad page into `data/raw/cscc/`.
2. `Rscript analysis/spatial/adapters/rds_to_celltable.R --inspect data/raw/cscc/Seurat_Object_cSCC.RDS`
   — this tells you whether response (major pathological response) and coordinates are in the object.
3. Export as above with `--responders` set to the values meaning major pathological response, then `export_spatial.py`
   with `--dataset-id cscc_spatial_pd1`. Keep only one timepoint per patient (pre-treatment) if the object mixes them:
   filter the CSV with one `awk`/pandas line before exporting.

### HNSCC CosMx, checkpoint blockade (Dryad doi:10.5061/dryad.n2z34tn7x) — TMA cores

1. Download and unzip `CosMx_RNA_TMA_Data_Expr.zip`; read the deposit's `README.md` for core annotations and pixel size.
2. Build `fov_map.csv` (`fov,sample,patient,response`) from those annotations.
3. ```bash
   python analysis/spatial/adapters/cosmx.py data/processed/hnscc_cells.csv.gz <exprMat>.csv <metadata>.csv fov_map.csv --um-per-px 0.12
   python analysis/spatial/export_spatial.py data/processed/hnscc_cells.csv.gz --dataset-id hnscc_cosmx_icb \
     --label "HNSCC CosMx · checkpoint blockade" --platform "CosMx 1000-plex" --therapy "Checkpoint blockade"
   ```

## Choosing the radius

`--radius` is the contact distance. 15–20 µm ≈ touching neighbors (good for IMC and segmented mIF);
30 µm (default) ≈ one cell apart (Xenium/CosMx, where segmentation is looser). Re-run with a different radius
before trusting any single z-score; a real interaction holds up across a range.

## Demo data

`make_synthetic_tissue.py` writes 12 invented cores where the CD8 pattern is planted by hand. Regenerate with
```bash
python analysis/spatial/make_synthetic_tissue.py /tmp/syn.csv.gz
python analysis/spatial/export_spatial.py /tmp/syn.csv.gz --dataset-id demo-cores --label "Demo · tissue cores" \
  --platform synthetic --study "Synthetic demo" --synthetic
```
