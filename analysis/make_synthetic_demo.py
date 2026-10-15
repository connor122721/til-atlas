#!/usr/bin/env python3
"""Write SYNTHETIC lab data so site/lab.html renders before you've run the real analysis.

The files use exactly the schema that analysis/01_score_signatures.R and
analysis/02_clonotypes.R export, so swapping in real data is a file copy.
Every file carries meta.synthetic = true and the page shows a banner while it does.
Nothing here is a finding: the responder/non-responder difference is put in by hand.

  python analysis/make_synthetic_demo.py site/data/lab
"""
import json
import math
import random
import sys
from pathlib import Path

R = random.Random(2020)
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "site/data/lab")
OUT.mkdir(parents=True, exist_ok=True)

# phenotype clusters: (name, umap center, spread, stem score mean, terminal score mean)
CLUSTERS = [
    ("Stem-like (TCF7+ DN)", (-6.0, 3.0), 1.4, 0.42, 0.08),
    ("Effector memory", (-1.5, 5.5), 1.6, 0.25, 0.18),
    ("Proliferating", (3.5, 6.0), 1.1, 0.15, 0.25),
    ("Effector (GZMB+)", (1.0, 0.5), 1.7, 0.12, 0.30),
    ("Terminal exhausted (CD39+ DP)", (5.5, -2.5), 1.5, 0.05, 0.48),
    ("Tissue-resident-like", (-3.0, -3.5), 1.3, 0.10, 0.32),
]
patients = [(f"P{i:02d}", "Responder" if i <= 7 else "Non-responder") for i in range(1, 15)]


def r3(x):
    return round(x, 3)


cells = {k: [] for k in ("x", "y", "stem", "term", "cluster", "patient")}
pat_rows = []
for pid, resp in patients:
    # put the effect in by hand: responders get more stem-like cells
    stem_w = R.uniform(0.18, 0.40) if resp == "Responder" else R.uniform(0.02, 0.14)
    w = [stem_w, 0.15, 0.08, 0.22, 0.25, 0.10]
    total = sum(w)
    w = [v / total for v in w]
    n = R.randint(320, 520)
    counts = [0] * len(CLUSTERS)
    for _ in range(n):
        k = R.choices(range(len(CLUSTERS)), w)[0]
        counts[k] += 1
        name, (cx, cy), sd, ms, mt = CLUSTERS[k]
        cells["x"].append(r3(R.gauss(cx, sd)))
        cells["y"].append(r3(R.gauss(cy, sd)))
        cells["stem"].append(r3(min(1, max(0, R.gauss(ms, 0.07)))))
        cells["term"].append(r3(min(1, max(0, R.gauss(mt, 0.07)))))
        cells["cluster"].append(k)
        cells["patient"].append(pid)
    pat_rows.append({"patient": pid, "response": resp, "nCells": n,
                     "stemFrac": r3(counts[0] / n), "termFrac": r3(counts[4] / n),
                     "persistence30": r3(min(0.99, max(0.01, R.gauss(0.55 if resp == "Responder" else 0.25, 0.15))))})

meta = {"synthetic": True, "source": "analysis/make_synthetic_demo.py (invented data)",
        "sample": "infusion product, CD8+ T cells", "clusters": [c[0] for c in CLUSTERS],
        "signatures": {"stem": "signatures/stem_like.txt", "term": "signatures/terminal.txt"}}
(OUT / "cells.json").write_text(json.dumps({"meta": meta, "patients": [p[0] for p in patients],
                                            "responses": dict(patients), **cells}, separators=(",", ":")))
(OUT / "patients.json").write_text(json.dumps({"meta": meta, "patients": pat_rows}, indent=1))

# clonotypes across compartments
AA = "ACDEFGHIKLMNPQRSTVWY"
COMPARTMENTS = ["tumor_pre", "product", "blood_d30", "tumor_post"]
clone_out = {}
for pid, resp in patients:
    clones = []
    for j in range(25):
        reactive = R.random() < (0.5 if resp == "Responder" else 0.2)
        stemmy = R.random() < (0.6 if resp == "Responder" else 0.2)
        pre = 10 ** R.uniform(-4, -1.3) * (2.5 if reactive else 0.6)
        prod = 10 ** R.uniform(-3.5, -1.2)
        blood = prod * (10 ** R.uniform(-0.3, 0.6) if stemmy else 10 ** R.uniform(-2.5, -0.5))
        post = (blood * 10 ** R.uniform(0, 0.8)) if (reactive and resp == "Responder") else blood * 10 ** R.uniform(-2, -0.3)
        cdr3 = "CAS" + "".join(R.choice(AA) for _ in range(R.randint(7, 11))) + "F"
        clones.append({"cdr3": cdr3, "reactive": reactive, "phenotype": "stem-like" if stemmy else "terminal",
                       "freq": dict(zip(COMPARTMENTS, [r3(math.log10(v)) if v > 1e-5 else None for v in (pre, prod, blood, post)]))})
    clone_out[pid] = {"response": resp, "clones": clones}
(OUT / "clones.json").write_text(json.dumps({"meta": {**meta, "compartments": COMPARTMENTS,
                                                      "freqScale": "log10 clone frequency"},
                                             "patients": clone_out}, separators=(",", ":")))
print(f"wrote synthetic lab data to {OUT}")
