#!/usr/bin/env python3
"""SYNTHETIC tissue cores in the standard cell-table format, for the spatial tab's demo.

12 invented 1-mm "TMA cores" (6 responders, 6 non-responders). By construction, responder cores
put CD8 T cells inside tumor nests and next to myeloid hubs; non-responder cores keep CD8 cells
out in the stroma. The point is to exercise the real pipeline (export_spatial.py) end to end and
show what a CD8–myeloid network or an immune-excluded pattern looks like on the page. It is not data.

  python analysis/spatial/make_synthetic_tissue.py OUT.csv.gz
"""
import sys

import numpy as np
import pandas as pd

rng = np.random.default_rng(2024)
STROMA = {"Fibroblast": 0.34, "Endothelial": 0.10, "CD4 T": 0.14, "B cell": 0.08, "Treg": 0.04, "Myeloid": 0.12}


def core(i: int, responder: bool) -> pd.DataFrame:
    R, spacing = 500.0, 16.0
    g = np.arange(-R, R, spacing)
    xx, yy = np.meshgrid(g, g)
    xy = np.c_[xx.ravel(), yy.ravel()] + rng.normal(0, 4, (xx.size, 2))
    xy = xy[np.hypot(xy[:, 0], xy[:, 1]) < R]
    n = len(xy)

    nests = rng.uniform(-330, 330, (rng.integers(4, 7), 2))
    nest_r = rng.uniform(90, 160, len(nests))
    d_nest = np.min(np.hypot(*(xy[:, None, :] - nests[None]).transpose(2, 0, 1)) - nest_r[None], axis=1)
    in_nest = d_nest < 0

    types = np.empty(n, dtype=object)
    types[in_nest] = np.where(rng.random(in_nest.sum()) < 0.9, "Tumor", "Myeloid")
    st = list(STROMA)
    p = np.array(list(STROMA.values())); p /= p.sum()
    types[~in_nest] = rng.choice(st, size=(~in_nest).sum(), p=p)

    # myeloid hubs: at nest margins in responders, deep in stroma in non-responders
    hubs = []
    for c, r in zip(nests, nest_r):
        ang = rng.uniform(0, 2 * np.pi)
        dist = r * (0.9 if responder else 1.0) + (0 if responder else rng.uniform(120, 200))
        hubs.append(c + dist * np.array([np.cos(ang), np.sin(ang)]))
    hubs = np.array(hubs)
    d_hub = np.min(np.hypot(*(xy[:, None, :] - hubs[None]).transpose(2, 0, 1)), axis=1)
    near_hub = (d_hub < 45) & (rng.random(n) < 0.55)
    types[near_hub] = "Myeloid"

    # CD8 T cells: placed by a weight map that differs by response (this IS the planted effect)
    n_cd8 = int(n * (rng.uniform(0.10, 0.15) if responder else rng.uniform(0.04, 0.07)))
    if responder:
        w = 4.0 * np.exp(-d_hub / 50.0) + 0.3 * in_nest + 0.03
    else:
        w = np.where(in_nest, 0.0, np.clip(d_nest / 150.0, 0, 1)) + 0.02 * (~in_nest)
    w = np.where(types == "Myeloid", 0.0, w)          # don't overwrite the hubs we just placed
    w = w / w.sum()
    pick = rng.choice(n, size=n_cd8, replace=False, p=w)
    types[pick] = "CD8 T"

    return pd.DataFrame({"cell_id": [f"C{i:02d}_{k}" for k in range(n)], "sample": f"core{i:02d}",
                         "patient": f"P{i:02d}", "x": xy[:, 0] + 500, "y": xy[:, 1] + 500, "cell_type": types,
                         "response": "Responder" if responder else "Non-responder"})


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic_cells.csv.gz"
    df = pd.concat([core(i, i <= 6) for i in range(1, 13)], ignore_index=True)
    df.to_csv(out, index=False)
    print(f"wrote {len(df):,} synthetic cells in 12 cores to {out}")
