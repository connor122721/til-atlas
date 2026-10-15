#!/usr/bin/env python3
"""CosMx SMI flat files -> standard cell table.

CosMx exports per slide: *_exprMat_file.csv(.gz) (fov, cell_ID, one column per gene) and
*_metadata_file.csv(.gz) (fov, cell_ID, CenterX_global_px, CenterY_global_px, ...). On tissue
microarrays each core spans one or more FOVs, so you supply a map from FOV to core/patient/response:

    fov_map.csv:  fov,sample,patient,response      (response = Responder / Non-responder)

  python analysis/spatial/adapters/cosmx.py OUT.csv.gz EXPR.csv.gz META.csv.gz fov_map.csv [--um-per-px 0.12]

For the HNSCC ICB deposit (Dryad doi:10.5061/dryad.n2z34tn7x) unzip CosMx_RNA_TMA_Data_Expr.zip and
read the README to build fov_map.csv from its core annotations. Check the pixel size in the
deposit's README: 0.12 µm/px is the usual CosMx value, but confirm it.
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from celltyping import type_cells  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("out"); ap.add_argument("expr"); ap.add_argument("meta"); ap.add_argument("fov_map")
ap.add_argument("--um-per-px", type=float, default=0.12)
a = ap.parse_args()

expr = pd.read_csv(a.expr)
meta = pd.read_csv(a.meta)
fmap = pd.read_csv(a.fov_map)
key = ["fov", "cell_ID"]
expr = expr[expr["cell_ID"] != 0]                                   # cell_ID 0 = transcripts outside cells
genes = [c for c in expr.columns if c not in key and not re.match(r"^(NegPrb|Negative|SystemControl|FalseCode)", c)]
df = expr[key].merge(meta, on=key, how="inner").merge(fmap, on="fov", how="inner")
expr = expr.set_index(key).loc[pd.MultiIndex.from_frame(df[key]), genes]
labels, rep = type_cells(sp.csr_matrix(expr.to_numpy()), genes)
print(rep["counts"], file=sys.stderr)
xcol = next(c for c in meta.columns if c.lower().startswith("centerx_global"))
ycol = next(c for c in meta.columns if c.lower().startswith("centery_global"))
out = pd.DataFrame({"cell_id": df["fov"].astype(str) + "_" + df["cell_ID"].astype(str), "sample": df["sample"],
                    "x": df[xcol] * a.um_per_px, "y": df[ycol] * a.um_per_px, "cell_type": labels,
                    "patient": df["patient"], "response": df["response"]})
out.to_csv(a.out, index=False)
print(f"wrote {len(out):,} cells from {out['sample'].nunique()} samples to {a.out}", file=sys.stderr)
