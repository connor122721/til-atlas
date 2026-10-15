"""Offline tests for the spatial adapters and exporter, on tiny files in the real instrument formats.
   python -m unittest discover -s analysis/spatial/tests      (needs numpy, scipy, pandas, h5py)
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SP = HERE.parent
GENES = ["PMEL", "MLANA", "SOX10", "CD8A", "CD8B", "CD3E", "CD68", "CD163", "LYZ", "COL1A1", "DCN", "ACTB"]
TRUTH = {"Tumor": ["PMEL", "MLANA", "SOX10"], "CD8 T": ["CD8A", "CD8B", "CD3E"],
         "Myeloid": ["CD68", "CD163", "LYZ"], "Fibroblast": ["COL1A1", "DCN"]}


def fake_counts(n, rng):
    kinds = rng.choice(list(TRUTH), size=n)
    X = rng.poisson(0.2, (n, len(GENES)))
    for i, k in enumerate(kinds):
        for g in TRUTH[k]:
            X[i, GENES.index(g)] += rng.poisson(8)
    return kinds, X


class TestSpatial(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import h5py
        import scipy.sparse as sp
        cls.tmp = Path(tempfile.mkdtemp()); rng = np.random.default_rng(0)
        # --- mock Xenium bundle (10x HDF5 layout: features x cells, CSC)
        b = cls.tmp / "bundle"; b.mkdir()
        cls.kinds, X = fake_counts(600, rng)
        M = sp.csc_matrix(X.T)
        bcs = [f"cell{i}" for i in range(600)]
        with h5py.File(b / "cell_feature_matrix.h5", "w") as f:
            m = f.create_group("matrix")
            m["data"], m["indices"], m["indptr"], m["shape"] = M.data, M.indices, M.indptr, np.array(M.shape)
            m["barcodes"] = np.array(bcs, dtype="S")
            feat = m.create_group("features")
            feat["name"] = np.array(GENES + ["NegControlProbe_1"], dtype="S")[: len(GENES)]
            feat["feature_type"] = np.array(["Gene Expression"] * len(GENES), dtype="S")
        pd.DataFrame({"cell_id": bcs, "x_centroid": rng.uniform(0, 800, 600), "y_centroid": rng.uniform(0, 800, 600)}) \
            .to_csv(b / "cells.csv.gz", index=False)
        # --- mock CosMx flat files, two FOVs = two cores
        kinds2, X2 = fake_counts(500, rng)
        fov = np.repeat([1, 2], 250); cid = np.r_[np.arange(1, 251), np.arange(1, 251)]
        e = pd.DataFrame(X2, columns=GENES); e.insert(0, "cell_ID", cid); e.insert(0, "fov", fov)
        e["NegPrb1"] = 0
        e.to_csv(cls.tmp / "expr.csv", index=False)
        pd.DataFrame({"fov": fov, "cell_ID": cid, "CenterX_global_px": rng.uniform(0, 6000, 500),
                      "CenterY_global_px": rng.uniform(0, 6000, 500)}).to_csv(cls.tmp / "meta.csv", index=False)
        pd.DataFrame({"fov": [1, 2], "sample": ["coreA", "coreB"], "patient": ["P1", "P2"],
                      "response": ["Responder", "Non-responder"]}).to_csv(cls.tmp / "fov.csv", index=False)
        py = sys.executable
        subprocess.run([py, SP / "adapters/xenium.py", cls.tmp / "xen.csv", f"{b}:s1"], check=True, capture_output=True)
        subprocess.run([py, SP / "adapters/cosmx.py", cls.tmp / "cos.csv", cls.tmp / "expr.csv", cls.tmp / "meta.csv",
                        cls.tmp / "fov.csv"], check=True, capture_output=True)
        subprocess.run([py, SP / "export_spatial.py", cls.tmp / "cos.csv", "--dataset-id", "t", "--label", "t",
                        "--out", cls.tmp / "site", "--perms", "50", "--min-cells", "50"], check=True, capture_output=True)

    def test_xenium_typing_recovers_truth(self):
        df = pd.read_csv(self.tmp / "xen.csv")
        self.assertEqual(len(df), 600)
        self.assertGreater((df["cell_type"].to_numpy() == self.kinds).mean(), 0.9)

    def test_cosmx_units_and_samples(self):
        df = pd.read_csv(self.tmp / "cos.csv")
        self.assertEqual(sorted(df["sample"].unique()), ["coreA", "coreB"])
        self.assertLessEqual(df["x"].max(), 6000 * 0.12 + 1)

    def test_export_schema(self):
        out = self.tmp / "site"
        s = json.loads((out / "t/samples.json").read_text())
        n = json.loads((out / "t/nhood.json").read_text())
        self.assertEqual(len(s["samples"]), 2)
        k = len(n["types"])
        self.assertEqual(len(n["z"]["coreA"]), k)
        self.assertTrue((out / "t/cells/coreA.json").exists())
        idx = json.loads((out / "index.json").read_text())
        self.assertTrue(idx["datasets"][0]["hasResponse"])


if __name__ == "__main__":
    unittest.main()
