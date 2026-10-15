"""python -m unittest discover -s analysis/bulk/tests   (numpy + pandas; no network)"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fetch_iatlas as fi  # noqa: E402
import geo_cohort as gc  # noqa: E402


class TestBulk(unittest.TestCase):
    def test_auc_matches_bruteforce(self):
        rng = np.random.default_rng(0)
        pos, neg = rng.normal(1, 1, 15), rng.normal(0, 1, 20)
        brute = np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
        auc, lo, hi, p = fi.auc_ci(pos, neg, rng, B=200)
        self.assertAlmostEqual(auc, round(brute, 3), places=3)
        self.assertTrue(lo <= auc <= hi and 0 < p < 1)

    def test_geo_cohort_end_to_end(self):
        tmp = Path(tempfile.mkdtemp())
        rng = np.random.default_rng(1)
        genes = ["CD8A", "CD8B", "GZMA", "PRF1", "MS4A1", "CD79A", "ACTB"]
        samples = [f"S{i}" for i in range(20)]
        resp = ["CR" if i < 8 else "PD" for i in range(20)]
        X = pd.DataFrame(rng.poisson(50, (len(genes), 20)).astype(float), index=genes, columns=samples)
        X.loc[["CD8A", "CD8B"], samples[:8]] *= 6          # responders: more CD8
        X.index.name = "gene"
        X.to_csv(tmp / "expr.tsv", sep="\t")
        pd.DataFrame({"sample": samples, "response": resp}).to_csv(tmp / "clin.csv", index=False)
        out = tmp / "bulk"; out.mkdir()
        (out / "index.json").write_text(json.dumps({"cohorts": []}))
        (out / "summary.json").write_text(json.dumps({"stats": {}}))
        gc.OUT = out
        gc.main(["--expr", str(tmp / "expr.tsv"), "--clinical", str(tmp / "clin.csv"), "--id", "t", "--label", "T",
                 "--cancer", "Melanoma", "--therapy", "TIL ACT"])
        st = json.loads((out / "summary.json").read_text())["stats"]["t"]["response"]
        self.assertGreater(st["cd8"]["auc"], 0.9)
        self.assertEqual(json.loads((out / "index.json").read_text())["cohorts"][0]["endpoints"]["response"], {"R": 8, "NR": 12})


if __name__ == "__main__":
    unittest.main()
