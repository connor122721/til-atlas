#!/usr/bin/env python3
"""Score any bulk RNA-seq cohort you've downloaded (e.g. from GEO) the same way as the iAtlas cohorts,
and add it to the Bulk RNA tab.

Built for the TIL-therapy cohort of Lauss et al., Nat Commun 2017 (GEO GSE100797: pre-treatment melanoma
RNA-seq, 25 patients treated with TIL ACT; clinical response in the paper's Supplementary Data 4), but it
works for any matrix + clinical table.

Inputs
  --expr      genes x samples table (TSV/CSV, optionally .gz). First column = gene symbol. Counts, FPKM or TPM;
              values are log2(x+1)-transformed if they look linear.
  --clinical  CSV with columns: sample, response (Responder/Non-responder or CR/PR/SD/PD),
              optional benefit (Responder/Non-responder), os (months), os_event (1 = death, 0 = censored)
Example
  python analysis/bulk/geo_cohort.py --expr data/raw/GSE100797/GSE100797_ProcessedData.txt.gz \\
      --clinical data/processed/lauss2017_clinical.csv --id lauss2017 --label "Lauss 2017" \\
      --cancer Melanoma --therapy "TIL ACT" --citation "Lauss et al., Nat Commun 2017" --pmid 29170503
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from fetch_iatlas import OUT, auc_ci  # noqa: E402

RECIST = {"CR": True, "PR": True, "SD": False, "PD": False}


def to_resp(v):
    if pd.isna(v):
        return None
    s = str(v).strip()
    if s in RECIST:
        return RECIST[s]
    if s.lower() in ("responder", "r", "yes", "true", "1"):
        return True
    if s.lower() in ("non-responder", "nonresponder", "nr", "no", "false", "0"):
        return False
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--expr", "--clinical", "--id", "--label", "--cancer", "--therapy"):
        ap.add_argument(a, required=True)
    ap.add_argument("--citation", default="")
    ap.add_argument("--pmid", default=None)
    a = ap.parse_args(argv)
    sigs = json.loads((Path(__file__).with_name("signatures.json")).read_text())
    sep = "," if ".csv" in a.expr else "\t"
    ex = pd.read_csv(a.expr, sep=sep, index_col=0)
    ex.index = ex.index.astype(str)
    ex = ex[~ex.index.duplicated()].apply(pd.to_numeric, errors="coerce")
    cl = pd.read_csv(a.clinical, dtype={"sample": str})
    cl = cl[cl["sample"].isin(ex.columns)]
    if cl.empty:
        sys.exit(f"no clinical samples match expression columns (e.g. {list(ex.columns[:3])})")
    M = ex[cl["sample"].tolist()]
    if np.nanmax(M.to_numpy()) > 50:
        M = np.log2(M.clip(lower=0) + 1)
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
    scores, used = {}, {}
    for k, sig in sigs["computed"].items():
        g = [x for x in sig["genes"] if x in Z.index]
        used[k] = g
        scores[k] = Z.loc[g].mean(axis=0) if len(g) >= max(1, len(sig["genes"]) // 2) else pd.Series(np.nan, index=Z.columns)
    panel = {g: Z.loc[g] for g in sigs["panel"] if g in Z.index}
    cl["response_b"] = cl["response"].map(to_resp)
    cl["benefit_b"] = cl["benefit"].map(to_resp) if "benefit" in cl else None
    rng = np.random.default_rng(7)
    stats, counts = {}, {}
    for ep, col in (("response", "response_b"), ("benefit", "benefit_b")):
        if col not in cl or cl[col].isna().all():
            continue
        r = cl.set_index("sample")[col].dropna().astype(bool)
        if r.sum() < 3 or (~r).sum() < 3:
            continue
        counts[ep] = {"R": int(r.sum()), "NR": int((~r).sum())}
        stats[ep] = {}
        for k, v in scores.items():
            v = v.reindex(r.index).astype(float)
            ok = v.notna()
            if (ok & r).sum() >= 3 and (ok & ~r).sum() >= 3:
                au, lo, hi, p = auc_ci(v[ok & r].to_numpy(), v[ok & ~r].to_numpy(), rng)
                stats[ep][k] = {"auc": au, "lo": lo, "hi": hi, "p": p, "nR": int((ok & r).sum()), "nNR": int((ok & ~r).sum())}
    if not stats:
        sys.exit("need at least 3 responders and 3 non-responders")
    r3 = lambda x: None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), 3)
    lab = lambda v: None if v is None or pd.isna(v) else ("Responder" if v else "Non-responder")
    rows = [{"sample": s["sample"], "patient": s["sample"], "response": lab(s["response_b"]),
             "benefit": lab(s.get("benefit_b")), "recist": None if pd.isna(s["response"]) else str(s["response"]),
             "os": r3(s.get("os")) if "os" in s else None,
             "osEvent": bool(s["os_event"]) if "os_event" in s and pd.notna(s.get("os_event")) else None,
             "scores": {k: r3(v.get(s["sample"])) for k, v in scores.items()},
             "genes": {g: r3(z.get(s["sample"])) for g, z in panel.items()}} for _, s in cl.iterrows()]
    meta = {"id": a.id, "label": a.label, "cancer": a.cancer, "therapy": a.therapy, "name": a.label, "citation": a.citation,
            "pmid": a.pmid, "source": Path(a.expr).name, "n": len(rows), "endpoints": counts, "genesUsed": used}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{a.id}.json").write_text(json.dumps({"meta": meta, "stats": stats, "samples": rows}, separators=(",", ":")))
    idx = json.loads((OUT / "index.json").read_text())
    idx["cohorts"] = [c for c in idx["cohorts"] if c["id"] != a.id] + [{k: meta[k] for k in ("id", "label", "cancer", "therapy", "name", "citation", "pmid", "n", "endpoints")}]
    (OUT / "index.json").write_text(json.dumps(idx, indent=1))
    summ = json.loads((OUT / "summary.json").read_text())
    summ["stats"][a.id] = stats
    (OUT / "summary.json").write_text(json.dumps(summ, separators=(",", ":")))
    print(f"added {a.id}: {len(rows)} samples; " + "; ".join(f"{e} {c['R']} R / {c['NR']} NR" for e, c in counts.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
