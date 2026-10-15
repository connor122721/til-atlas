#!/usr/bin/env python3
"""Xenium output bundle(s) -> standard cell table (cell_id, sample, x, y, cell_type[, patient, response]).

Each bundle directory must contain cells.csv.gz (x_centroid, y_centroid in µm) and
cell_feature_matrix.h5 (10x HDF5). Cell types come from marker modules (markers.json).

  python analysis/spatial/adapters/xenium.py OUT.csv.gz BUNDLE_DIR[:sample[:patient[:response]]] ...

Example, the public 10x melanoma sections (no treatment, so no response):
  python analysis/spatial/adapters/xenium.py data/processed/xenium_cells.csv.gz \\
      data/raw/xenium/Xenium_Prime_Human_Skin_FFPE:prime5k \\
      data/raw/xenium/Xenium_V1_hSkin_Melanoma_Base_FFPE:v1_base
  python analysis/spatial/export_spatial.py data/processed/xenium_cells.csv.gz \\
      --dataset-id xenium_melanoma --label "Xenium melanoma (10x)" --platform Xenium --study "10x Genomics public data"
Requires h5py (pip install h5py).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from celltyping import type_cells  # noqa: E402


def read_10x_h5(path: Path):
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = tuple(m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape).T.tocsr()  # cells x features
        names = [x.decode() for x in m["features"]["name"][:]]
        ftype = [x.decode() for x in m["features"]["feature_type"][:]]
        barcodes = [x.decode() for x in m["barcodes"][:]]
    keep = np.array([t == "Gene Expression" for t in ftype])               # drop negative-control probes
    return X[:, keep], [n for n, k in zip(names, keep) if k], barcodes


def one(spec: str) -> pd.DataFrame:
    parts = spec.split(":")
    d = Path(parts[0]); sample = parts[1] if len(parts) > 1 else d.name
    X, genes, bcs = read_10x_h5(d / "cell_feature_matrix.h5")
    cells = pd.read_csv(d / "cells.csv.gz").set_index("cell_id")
    xy = cells.loc[bcs, ["x_centroid", "y_centroid"]].to_numpy()
    labels, rep = type_cells(X, genes)
    print(f"{sample}: {len(bcs):,} cells; {json.dumps(rep['counts'])}", file=sys.stderr)
    return pd.DataFrame({"cell_id": bcs, "sample": sample, "x": xy[:, 0], "y": xy[:, 1], "cell_type": labels,
                         "patient": parts[2] if len(parts) > 2 else sample,
                         "response": parts[3] if len(parts) > 3 else None})


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    pd.concat([one(s) for s in sys.argv[2:]]).to_csv(sys.argv[1], index=False)
    print(f"wrote {sys.argv[1]}", file=sys.stderr)
