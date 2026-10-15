#!/usr/bin/env python3
"""Wu et al., Nature 2020 (GEO GSE139555): single-cell RNA + TCR of T cells from tumor, normal adjacent
tissue (NAT) and blood of 14 patients (lung, endometrial, colorectal, renal).

Input: the processed AnnData that scirpy publishes on GitHub releases (raw counts, the authors' cluster
labels and clonotypes, UMAP). Downloaded automatically (563 MB) into data/raw/wu2020/.

Outputs (all from the real data):
  site/data/lab/wu2020/{cells,patients,clones,expansion}.json   Single-cell tab
  site/data/sim/repertoires.json                                 real tumor repertoires for the Clone simulator

Per CD8 T cell (authors' clusters 8.*): log-normalized expression, z-scored across CD8 cells, averaged over
the gene sets in analysis/signatures/ -> stem-like, terminal and tumor-reactivity scores. Each cell is then
called stem-like (S), effector (E) or terminal (X) for the simulator:
  X if terminal score >= its 75th percentile; S if stem-like score >= its 75th percentile (and not X); else E.
A clone is flagged "predicted tumor-reactive" when its tumor cells' mean reactivity score is in the top
quarter of expanded tumor clones (>= 3 cells). This is a transcriptional prediction, not a functional assay.

  python analysis/datasets/wu2020/prepare.py
Needs numpy, scipy, h5py. Peak memory ~3 GB.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis" / "lib"))
from h5ad_lite import obs_col, read_X, var_names  # noqa: E402

URL = "https://github.com/scverse/scirpy/releases/download/d0.1.0/wu2020.h5ad"
RAW = ROOT / "data" / "raw" / "wu2020" / "wu2020.h5ad"
LAB = ROOT / "site" / "data" / "lab"
SIM = ROOT / "site" / "data" / "sim"
SIG = ROOT / "analysis" / "signatures"
REACT = ["CXCL13", "ENTPD1", "ITGAE", "PDCD1", "TNFRSF9", "HAVCR2", "LAYN", "TOX"]
CANCER = {"Lung": "Lung", "Endo": "Endometrial", "Colon": "Colorectal", "Renal": "Kidney"}
COMP = ["Tumor", "NAT", "Blood"]
rng = np.random.default_rng(2020)


def read_set(p: Path) -> list[str]:
    return [l.strip() for l in p.read_text().splitlines() if l.strip() and not l.startswith("#")]


def r3(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 3)


def main():
    if not RAW.exists():
        RAW.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL}", file=sys.stderr)
        urllib.request.urlretrieve(URL, RAW)
    f = h5py.File(RAW, "r")
    patient = obs_col(f, "patient").astype(str)
    source = obs_col(f, "source").astype(str)
    cluster = obs_col(f, "cluster_orig").astype(str)
    clono = obs_col(f, "clonotype_orig")
    umap = f["obsm"]["X_umap_orig"][:]
    genes = var_names(f)

    cd8 = np.array([c.startswith("8.") for c in cluster])
    print(f"{cd8.sum():,} CD8 T cells of {len(cd8):,}", file=sys.stderr)
    X = read_X(f)[cd8]
    tot = np.asarray(X.sum(axis=1)).ravel()
    gi = {g: i for i, g in enumerate(genes)}

    sets = {"stem": read_set(SIG / "stem_like.txt"), "term": read_set(SIG / "terminal.txt"), "react": REACT}
    need = sorted({g for s in sets.values() for g in s if g in gi})
    E = X[:, [gi[g] for g in need]].toarray().astype(np.float32)
    E = np.log1p(E / tot[:, None] * 1e4)
    E = (E - E.mean(0)) / np.where(E.std(0) > 0, E.std(0), 1)
    col = {g: j for j, g in enumerate(need)}
    score = {k: E[:, [col[g] for g in s if g in col]].mean(1) for k, s in sets.items()}
    used = {k: [g for g in s if g in col] for k, s in sets.items()}
    del X, E

    p_cd8, s_cd8, c_cd8, k_cd8, u_cd8 = patient[cd8], source[cd8], cluster[cd8], clono[cd8], umap[cd8]
    q_term, q_stem = np.quantile(score["term"], 0.75), np.quantile(score["stem"], 0.75)
    state = np.where(score["term"] >= q_term, "X", np.where(score["stem"] >= q_stem, "S", "E"))
    clusters = sorted(set(c_cd8))
    patients = sorted(set(p_cd8), key=lambda p: (CANCER[p.rstrip("0123456789")], p))
    cancer_of = {p: CANCER[p.rstrip("0123456789")] for p in patients}

    # ---------------- cells.json (stratified subsample for the UMAP)
    keep = []
    for p in patients:
        ix = np.where(p_cd8 == p)[0]
        keep += list(rng.choice(ix, size=min(len(ix), 1800), replace=False))
    keep = np.sort(np.array(keep))
    out = LAB / "wu2020"
    out.mkdir(parents=True, exist_ok=True)
    meta = {"synthetic": False, "dataset": "wu2020", "study": "Wu et al. 2020", "therapy": "untreated (surgical resection)",
            "sample": "tumor, normal adjacent tissue and blood", "clusters": clusters, "groupBy": "cancer",
            "source": "GSE139555 via scirpy GitHub release d0.1.0 (wu2020.h5ad)", "signatures": used,
            "scoreNote": "mean z-scored log-normalized expression across CD8 T cells",
            "compartmentsOfCells": COMP}
    (out / "cells.json").write_text(json.dumps({
        "meta": meta, "patients": patients, "responses": {}, "groups": cancer_of,
        "x": np.round(u_cd8[keep, 0], 3).tolist(), "y": np.round(u_cd8[keep, 1], 3).tolist(),
        "stem": np.round(score["stem"][keep], 3).tolist(), "term": np.round(score["term"][keep], 3).tolist(),
        "cluster": [clusters.index(c) for c in c_cd8[keep]], "patient": p_cd8[keep].tolist(),
        "compartment": [COMP.index(s) for s in s_cd8[keep]]}, separators=(",", ":")))

    # ---------------- clones, expansion categories, simulator repertoires
    pat_rows, clone_out, expansion, reps = [], {}, [], []
    for p in patients:
        m = p_cd8 == p
        tumor = m & (s_cd8 == "Tumor")
        has_tcr = np.array([k is not None and str(k) not in ("nan", "None", "") for k in k_cd8])
        counts = {}
        for c in COMP:
            mm = m & (s_cd8 == c) & has_tcr
            ids, n = np.unique(k_cd8[mm].astype(str), return_counts=True)
            counts[c] = dict(zip(ids, n))
        totals = {c: sum(counts[c].values()) for c in COMP}
        tcl = counts["Tumor"]
        # expansion categories, weighted by tumor cells
        cat = {"Tumor singleton": 0, "Expanded in tumor only": 0, "Expanded in tumor and NAT": 0}
        in_blood = 0
        for cid, n in tcl.items():
            nn = counts["NAT"].get(cid, 0)
            k = "Tumor singleton" if n == 1 else ("Expanded in tumor and NAT" if nn >= 2 else "Expanded in tumor only")
            cat[k] += n
            in_blood += n if counts["Blood"].get(cid, 0) > 0 else 0
        tt = max(1, totals["Tumor"])
        expansion.append({"patient": p, "group": cancer_of[p], "tumorCells": int(totals["Tumor"]),
                          "categories": {k: r3(v / tt) for k, v in cat.items()},
                          "sharedWithBlood": r3(in_blood / tt) if totals["Blood"] else None,
                          "hasBlood": bool(totals["Blood"])})
        # per-clone phenotype and reactivity among tumor cells
        tix = np.where(tumor & has_tcr)[0]
        tk = k_cd8[tix].astype(str)
        per = {}
        for cid in np.unique(tk):
            ii = tix[tk == cid]
            st, n = np.unique(state[ii], return_counts=True)
            per[cid] = {"n": len(ii), "S": int(dict(zip(st, n)).get("S", 0)), "E": int(dict(zip(st, n)).get("E", 0)),
                        "X": int(dict(zip(st, n)).get("X", 0)), "react": float(score["react"][ii].mean()),
                        "stem": float(score["stem"][ii].mean())}
        expanded = [v for v in per.values() if v["n"] >= 3]
        thr = np.quantile([v["react"] for v in expanded], 0.75) if expanded else np.inf
        tN = max(1, totals["NAT"])
        for cid, v in per.items():
            v["reactive"] = bool(v["n"] >= 3 and v["react"] >= thr)
            # TCR-only proxy, independent of expression: expanded in tumor and >= 4x more frequent there than in NAT
            nn = counts["NAT"].get(cid, 0)
            v["enriched"] = bool(v["n"] >= 3 and (v["n"] / tt) >= 4 * ((nn + 0.5) / tN)) if totals["NAT"] else None
        top = sorted(per.items(), key=lambda kv: -kv[1]["n"])
        # patients.json row (unit = patient, tumor CD8 cells)
        tm = tumor
        pat_rows.append({"patient": p, "group": cancer_of[p], "response": None, "nCells": int(tm.sum()),
                         "stemFrac": r3((state[tm] == "S").mean()), "termFrac": r3((state[tm] == "X").mean()),
                         "persistence30": None,
                         "dualExpanded": r3(cat["Expanded in tumor and NAT"] / tt)})
        # clones.json (top 25 tumor clones)
        def phen(v):
            return "stem-like" if v["S"] >= max(v["E"], v["X"]) else ("terminal" if v["X"] >= v["E"] else "effector")
        clone_out[p] = {"response": None, "group": cancer_of[p], "clones": [{
            "cdr3": cid.split(".")[-1], "id": cid, "reactive": v["reactive"], "phenotype": phen(v),
            "freq": {c: (r3(np.log10(counts[c][cid] / totals[c])) if counts[c].get(cid) else None) for c in COMP if totals[c]}}
            for cid, v in top[:25]]}
        # simulator repertoire: top 80 tumor clones with S/E/X counts
        reps.append({"patient": p, "cancer": cancer_of[p], "tumorCD8": int(tm.sum()), "clonesTotal": len(per),
                     "hasBlood": bool(totals["Blood"]),
                     "clones": [{"id": cid, "S": v["S"], "E": v["E"], "X": v["X"], "reactive": v["reactive"], "enriched": v["enriched"],
                                 "inBlood": bool(counts["Blood"].get(cid, 0)), "inNAT": int(counts["NAT"].get(cid, 0))}
                                for cid, v in top[:80]]})
        print(f"  {p}: {int(tm.sum())} tumor CD8, {len(per)} tumor clones, "
              f"dual-expanded {cat['Expanded in tumor and NAT'] / tt:.0%}", file=sys.stderr)

    (out / "patients.json").write_text(json.dumps({"meta": meta, "patients": pat_rows}, indent=1))
    (out / "clones.json").write_text(json.dumps({"meta": {**meta, "compartments": COMP, "freqScale": "log10 clone frequency among CD8 T cells with a TCR in that compartment"},
                                                 "patients": clone_out}, separators=(",", ":")))
    (out / "expansion.json").write_text(json.dumps({"meta": meta, "patients": expansion}, indent=1))
    SIM.mkdir(parents=True, exist_ok=True)
    (SIM / "repertoires.json").write_text(json.dumps({
        "meta": {"source": "Wu et al. 2020 (GSE139555), tumor CD8 T cells", "stateRule": "X: terminal score >= 75th pct; S: stem-like >= 75th pct; else E",
                 "reactiveRule": "mean tumor-reactivity score in top quarter of tumor clones with >= 3 cells", "reactGenes": used["react"],
                 "enrichedRule": ">= 3 tumor cells and >= 4x higher frequency in tumor than in normal adjacent tissue (TCR only)"},
        "patients": reps}, separators=(",", ":")))

    idx_path = LAB / "index.json"
    index = json.loads(idx_path.read_text()) if idx_path.exists() else {"datasets": []}
    entry = {"id": "wu2020", "label": "Wu 2020 · tumor, NAT and blood T cells", "study": "Wu et al., Nature 2020", "therapy": "untreated",
             "sample": "tumor CD8 T cells", "synthetic": False, "path": "data/lab/wu2020/", "hasClones": True,
             "hasExpansion": True, "groupBy": "cancer"}
    index["datasets"] = [d for d in index["datasets"] if d["id"] != "wu2020" and not d.get("synthetic")] + [entry]
    idx_path.write_text(json.dumps(index, indent=1))
    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
