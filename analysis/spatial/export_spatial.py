#!/usr/bin/env python3
"""Turn a standard per-cell table into the spatial tab's JSON.

Every spatial platform (Xenium, CosMx, IMC, mIF, ...) is first reduced by an adapter in
analysis/spatial/adapters/ to ONE table, one row per cell:

    cell_id, sample, x, y, cell_type[, patient][, response]
    x, y in micrometres; response = Responder / Non-responder (blank if none)

This script then computes, per sample:
  * cell-type composition and density per mm² (area = convex hull of the cells)
  * neighborhood enrichment: for every pair of cell types, how many pairs of cells lie
    within --radius µm of each other, compared with --perms random relabelings of the same
    cells (z = (observed − mean_perm) / sd_perm; positive = found together more than chance)
  * distance from each --focal cell (default CD8 T) to the nearest --target cell (default Tumor)
and writes site/data/spatial/<dataset-id>/{samples.json, nhood.json, cells/<sample>.json},
registering the dataset in site/data/spatial/index.json.

Only numpy, scipy and pandas are required. This is the same statistic as squidpy's
`gr.nhood_enrichment`, written out so you can read exactly what it does.

  python analysis/spatial/export_spatial.py cells.csv.gz --dataset-id xenium_melanoma \\
      --label "Xenium 5K melanoma" --platform "Xenium Prime 5K" [--radius 30] [--perms 200]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree

ROOT = Path(__file__).resolve().parents[2]


def neighborhood_z(xy: np.ndarray, labels: np.ndarray, k: int, radius: float, perms: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """Observed type-pair contact counts within `radius`, and z-scores against label permutations."""
    pairs = cKDTree(xy).query_pairs(radius, output_type="ndarray")
    if len(pairs) == 0:
        return np.zeros((k, k)), np.full((k, k), np.nan)

    def counts(lab):
        a, b = lab[pairs[:, 0]], lab[pairs[:, 1]]
        m = np.bincount(a * k + b, minlength=k * k).reshape(k, k)
        return m + m.T - np.diag(np.diag(m))           # undirected: count each pair once per side

    obs = counts(labels)
    null = np.stack([counts(rng.permutation(labels)) for _ in range(perms)])
    sd = null.std(axis=0)
    z = np.where(sd > 0, (obs - null.mean(axis=0)) / np.where(sd > 0, sd, 1), np.nan)
    return obs, z


def area_mm2(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return float("nan")
    return ConvexHull(xy).volume / 1e6                  # 2-D hull "volume" is area; µm² -> mm²


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("--dataset-id", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--platform", default="")
    ap.add_argument("--study", default="")
    ap.add_argument("--therapy", default="")
    ap.add_argument("--radius", type=float, default=30.0, help="contact radius in µm")
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--focal", default="CD8 T")
    ap.add_argument("--target", default="Tumor")
    ap.add_argument("--max-cells", type=int, default=25000, help="cells per sample kept for the map")
    ap.add_argument("--min-cells", type=int, default=200, help="skip samples smaller than this")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "site" / "data" / "spatial")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)

    df = pd.read_csv(a.table)
    need = {"cell_id", "sample", "x", "y", "cell_type"}
    if missing := need - set(df.columns):
        sys.exit(f"missing columns: {sorted(missing)} (have {list(df.columns)})")
    for c in ("patient", "response"):
        if c not in df:
            df[c] = None
    df = df.dropna(subset=["x", "y", "cell_type"])
    types = sorted(df["cell_type"].unique())
    k = len(types)
    tix = {t: i for i, t in enumerate(types)}
    out = a.out / a.dataset_id
    (out / "cells").mkdir(parents=True, exist_ok=True)

    samples, nhood = [], {}
    for sid, g in df.groupby("sample", sort=True):
        if len(g) < a.min_cells:
            print(f"  skip {sid}: {len(g)} cells", file=sys.stderr)
            continue
        xy = g[["x", "y"]].to_numpy(float)
        lab = g["cell_type"].map(tix).to_numpy()
        obs, z = neighborhood_z(xy, lab, k, a.radius, a.perms, rng)
        area = area_mm2(xy)
        n_by = np.bincount(lab, minlength=k)
        rec = {
            "sample": str(sid),
            "patient": None if pd.isna(g["patient"].iloc[0]) else str(g["patient"].iloc[0]),
            "response": None if pd.isna(g["response"].iloc[0]) or not str(g["response"].iloc[0]).strip() else str(g["response"].iloc[0]),
            "nCells": int(len(g)),
            "areaMm2": round(area, 4),
            "fraction": {t: round(n_by[i] / len(g), 4) for i, t in enumerate(types)},
            "density": {t: round(n_by[i] / area, 1) if area == area and area > 0 else None for i, t in enumerate(types)},
        }
        if a.focal in tix and a.target in tix and n_by[tix[a.focal]] and n_by[tix[a.target]]:
            d, _ = cKDTree(xy[lab == tix[a.target]]).query(xy[lab == tix[a.focal]])
            rec["focalToTarget"] = {"median": round(float(np.median(d)), 1),
                                    "within20": round(float(np.mean(d <= 20)), 4)}
        samples.append(rec)
        nhood[str(sid)] = [[None if np.isnan(v) else round(float(v), 2) for v in row] for row in z]

        keep = rng.choice(len(g), size=min(len(g), a.max_cells), replace=False)
        keep.sort()
        x0, y0 = xy[:, 0].min(), xy[:, 1].min()
        (out / "cells" / f"{slug(sid)}.json").write_text(json.dumps({
            "sample": str(sid), "nTotal": int(len(g)),
            "x": np.round(xy[keep, 0] - x0, 1).tolist(), "y": np.round(xy[keep, 1] - y0, 1).tolist(),
            "t": lab[keep].tolist()}, separators=(",", ":")))
        print(f"  {sid}: {len(g)} cells, {area:.3f} mm²", file=sys.stderr)

    meta = {"dataset": a.dataset_id, "label": a.label, "platform": a.platform, "study": a.study, "therapy": a.therapy,
            "synthetic": a.synthetic, "cellTypes": types, "radiusUm": a.radius, "perms": a.perms,
            "focal": a.focal, "target": a.target, "units": "µm", "source": Path(a.table).name}
    (out / "samples.json").write_text(json.dumps({"meta": meta, "samples": samples}, indent=1))
    (out / "nhood.json").write_text(json.dumps({"meta": meta, "types": types, "z": nhood}, separators=(",", ":")))

    idx_path = a.out / "index.json"
    index = json.loads(idx_path.read_text()) if idx_path.exists() else {"datasets": []}
    entry = {"id": a.dataset_id, "label": a.label, "platform": a.platform, "study": a.study, "therapy": a.therapy,
             "synthetic": a.synthetic, "path": f"data/spatial/{a.dataset_id}/",
             "samples": len(samples), "hasResponse": any(s["response"] for s in samples)}
    index["datasets"] = [d for d in index["datasets"] if d["id"] != a.dataset_id] + [entry]
    idx_path.write_text(json.dumps(index, indent=1))
    print(f"wrote {len(samples)} samples to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
