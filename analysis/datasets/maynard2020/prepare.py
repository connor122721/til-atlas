#!/usr/bin/env python3
"""Maynard et al., Cell 2020 (PRJNA591860): Smart-seq2 of lung cancer biopsies before targeted therapy
(treatment-naive, TN), at residual disease on therapy (RD) and at progression (PD), with B and T cell
receptors reconstructed by BraCeR/TraCeR. Processed AnnData from scirpy's GitHub release (919 MB).

Output: site/data/lab/maynard2020/receptors.json — per biopsy, BCR isotype mix and clonal expansion,
and TCR clonal expansion. Clones are defined within a patient by identical V gene + CDR3 amino acids
of the heavy (IGH) or beta (TRB) chain.

  python analysis/datasets/maynard2020/prepare.py
"""
from __future__ import annotations

import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis" / "lib"))
from h5ad_lite import obs_col  # noqa: E402

URL = "https://github.com/scverse/scirpy/releases/download/d0.1.0/maynard2020.h5ad"
RAW = ROOT / "data" / "raw" / "maynard2020" / "maynard2020.h5ad"
OUT = ROOT / "site" / "data" / "lab" / "maynard2020"
STAGE = {"naive": "Before treatment (TN)", "grouped_pr": "Residual disease on therapy (RD)", "grouped_pd": "Progression (PD)"}
ISO = {"IGHG": "IgG", "IGHA": "IgA", "IGHM": "IgM/IgD", "IGHD": "IgM/IgD"}


def iso(c):
    if c is None or str(c) in ("None", "nan"):
        return None
    return ISO.get(str(c)[:4], None)


def main():
    if not RAW.exists():
        RAW.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL}", file=sys.stderr)
        urllib.request.urlretrieve(URL, RAW)
    f = h5py.File(RAW, "r")
    g = lambda k: obs_col(f, k)
    loc, cgene, v, cdr3 = g("IR_VDJ_1_locus"), g("IR_VDJ_1_c_gene"), g("IR_VDJ_1_v_gene"), g("IR_VDJ_1_cdr3")
    stage, patient, sample = g("analysis"), g("patient_id"), g("sample_name")
    resp, site = g("best_response_status"), g("biopsy_site")
    vgene = np.array([None if x is None else str(x).split("*")[0] for x in v], dtype=object)

    rows = []
    for sid in sorted({s for s in sample if s is not None}):
        m = sample == sid
        st = stage[m][0]
        if st not in STAGE:
            continue
        rec = {"sample": sid, "patient": str(patient[m][0]).replace("_NAT", ""), "stage": st,
               "tissue": "normal adjacent" if str(patient[m][0]).endswith("_NAT") else str(site[m][0]),
               "response": None if resp[m][0] in (None, "nan") else str(resp[m][0]), "cells": int(m.sum())}
        for chain, key in (("IGH", "bcr"), ("TRB", "tcr")):
            mm = m & (loc == chain) & np.array([c is not None for c in cdr3])
            n = int(mm.sum())
            rec[key] = {"n": n}
            if n == 0:
                continue
            clones = Counter(zip(vgene[mm], cdr3[mm]))
            sizes = np.array(sorted(clones.values(), reverse=True))
            rec[key].update({"clones": int(len(sizes)), "expandedShare": round(float(sizes[sizes >= 2].sum() / n), 3),
                             "topClone": round(float(sizes[0] / n), 3)})
            if chain == "IGH":
                iso_counts = Counter(i for i in (iso(c) for c in cgene[mm]) if i)
                tot = sum(iso_counts.values()) or 1
                rec[key]["isotype"] = {k: round(iso_counts.get(k, 0) / tot, 3) for k in ("IgG", "IgA", "IgM/IgD")}
                rec[key]["isotypeN"] = dict(iso_counts)
        rows.append(rec)

    pooled = {}
    for st in STAGE:
        mm = (stage == st) & (loc == "IGH")
        c = Counter(i for i in (iso(x) for x in cgene[mm]) if i)
        pooled[st] = {"n": int(sum(c.values())), "isotype": dict(c)}
    vg = Counter(x for x in vgene[(loc == "IGH")] if x)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "receptors.json").write_text(json.dumps({
        "meta": {"dataset": "maynard2020", "study": "Maynard et al., Cell 2020", "therapy": "targeted therapy (EGFR/ALK/BRAF/ROS1 inhibitors)",
                 "source": "PRJNA591860 via scirpy GitHub release d0.1.0 (maynard2020.h5ad)", "stages": STAGE,
                 "cloneRule": "same V gene + CDR3 amino acids of IGH or TRB, within a biopsy",
                 "nBCR": int((loc == "IGH").sum()), "nTCR": int((loc == "TRB").sum())},
        "pooled": pooled, "topIGHV": vg.most_common(15), "samples": rows}, indent=1))
    print(f"wrote {len(rows)} biopsies; BCR {int((loc == 'IGH').sum())}, TCR {int((loc == 'TRB').sum())}", file=sys.stderr)


if __name__ == "__main__":
    main()
