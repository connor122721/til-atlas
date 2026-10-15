#!/usr/bin/env python3
"""Real bulk RNA-seq of immunotherapy cohorts -> the Bulk RNA tab.

Source: the CRI iAtlas-harmonized immune checkpoint cohorts that cBioPortal publishes in its
datahub repository on GitHub (ODbL licence). Each cohort has per-patient response, per-sample
expression (TPM) and iAtlas immune signature scores. Files are Git LFS objects, fetched from
media.githubusercontent.com, so this runs anywhere GitHub is reachable (laptop, cluster, Actions).

For every cohort it:
  1. keeps one pre-treatment sample per patient with a response label
  2. scores the gene sets in analysis/bulk/signatures.json (mean of within-cohort z-scored log2 TPM+1)
  3. pulls the matching iAtlas signature scores
  4. compares responders vs non-responders: AUC (Mann-Whitney), 95% bootstrap CI, two-sided p
and writes site/data/bulk/{index.json, summary.json, <cohort>.json}.

  python analysis/bulk/fetch_iatlas.py                  # all cohorts
  python analysis/bulk/fetch_iatlas.py --only mel_iatlas_liu_2019
Needs numpy + pandas. Downloads ~540 MB once into data/raw/iatlas/ (cached).
"""
from __future__ import annotations

import argparse
import io
import json
import math
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "iatlas"
OUT = ROOT / "site" / "data" / "bulk"
MEDIA = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public"
RAWGH = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public"

COHORTS = [  # id, short label, cancer  (citation, PMID and treatment are read from the study's own files)
    ("mel_iatlas_liu_2019", "Liu 2019", "Melanoma"),
    ("mel_iatlas_gide_2019", "Gide 2019", "Melanoma"),
    ("mel_iatlas_riaz_nivolumab_2017", "Riaz 2017", "Melanoma"),
    ("mel_iatlas_hugo_ucla_2016", "Hugo 2016", "Melanoma"),
    ("blca_iatlas_imvigor210_2017", "IMvigor210", "Bladder"),
    ("rcc_iatlas_immotion150_2018", "IMmotion150", "Kidney"),
    ("ccrcc_iatlas_choueiri_2016", "Choueiri 2016", "Kidney"),
    ("gbm_iatlas_prins_2019", "Cloughesy 2019", "Glioblastoma"),
    ("brca_iatlas_anders_2022", "Anders 2022", "Breast (TNBC)"),
    ("paad_iatlas_prince_2022", "PRINCE 2022", "Pancreas"),
]
ENDPOINTS = {"response": ("RESPONDER", "Objective response (CR/PR vs SD/PD)"),
             "benefit": ("CLINICAL_BENEFIT", "Clinical benefit (as defined by iAtlas for the cohort)")}


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url.split('/public/')[1]}", file=sys.stderr)
    req = urllib.request.Request(url, headers={"User-Agent": "til-atlas"})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    return dest


def read_cbio(path: Path) -> pd.DataFrame:
    """cBioPortal tab files: '#'-prefixed description rows, then a header row."""
    with open(path, encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False,
                     dtype={"PATIENT_ID": str, "SAMPLE_ID": str, "ENTITY_STABLE_ID": str})
    df.columns = [str(c) for c in df.columns]          # numeric sample IDs as column names -> str
    return df


def files_of(cid: str) -> dict:
    listing = {}
    for name in ("data_mrna_seq_tpm.txt", "data_mrna_seq_expression.txt"):
        try:
            with urllib.request.urlopen(f"{RAWGH}/{cid}/{name}", timeout=60) as r:
                if r.status == 200:
                    listing["expr"] = name
                    break
        except Exception:
            continue
    return listing


def auc_ci(pos: np.ndarray, neg: np.ndarray, rng, B=1000):
    """AUC = P(score_responder > score_nonresponder); ties count 1/2. CI by stratified bootstrap."""
    def auc(p, n):
        allv = np.r_[p, n]
        ranks = pd.Series(allv).rank().to_numpy()
        return (ranks[: len(p)].sum() - len(p) * (len(p) + 1) / 2) / (len(p) * len(n))
    a = auc(pos, neg)
    bs = [auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg))) for _ in range(B)]
    n1, n2 = len(pos), len(neg)
    U = a * n1 * n2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (abs(U - n1 * n2 / 2) - 0.5) / sd
    p = math.erfc(z / math.sqrt(2))
    return round(a, 3), round(float(np.quantile(bs, 0.025)), 3), round(float(np.quantile(bs, 0.975)), 3), float(f"{min(1.0, p):.3g}")


def truthy(v) -> bool | None:
    s = str(v).strip().lower()
    if s in ("true", "yes", "1"):
        return True
    if s in ("false", "no", "0"):
        return False
    return None


def process(cid, label, cancer, sigs, rng):
    d = RAW / cid
    pat = read_cbio(fetch(f"{MEDIA}/{cid}/data_clinical_patient.txt", d / "data_clinical_patient.txt"))
    smp = read_cbio(fetch(f"{MEDIA}/{cid}/data_clinical_sample.txt", d / "data_clinical_sample.txt"))
    gs = read_cbio(fetch(f"{MEDIA}/{cid}/data_gene_signatures.txt", d / "data_gene_signatures.txt"))
    expr_name = files_of(cid).get("expr")
    if not expr_name:
        print(f"  {cid}: no expression file, skipped", file=sys.stderr)
        return None
    ex = read_cbio(fetch(f"{MEDIA}/{cid}/{expr_name}", d / expr_name))
    meta_txt = fetch(f"{RAWGH}/{cid}/meta_study.txt", d / "meta_study.txt").read_text()
    md = dict(l.split(": ", 1) for l in meta_txt.splitlines() if ": " in l)

    # ---- one pre-treatment sample per patient that has expression
    for ep, (col, _) in ENDPOINTS.items():
        pat[ep] = pat[col].map(truthy) if col in pat else None
    keep_cols = ["PATIENT_ID", *ENDPOINTS] + [c for c in ("OS_MONTHS", "OS_STATUS", "PFS_MONTHS", "ICI_RX") if c in pat]
    s = smp.merge(pat[keep_cols], on="PATIENT_ID", how="left")
    if "SAMPLE_TREATMENT" in s:
        pre = s["SAMPLE_TREATMENT"].astype(str).str.lower().str.startswith("pre")
        s = s[pre | s["SAMPLE_TREATMENT"].isna()]
    expr_samples = [c for c in ex.columns if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    s = s[s["SAMPLE_ID"].isin(expr_samples)].drop_duplicates("PATIENT_ID")
    s = s[s[list(ENDPOINTS)].notna().any(axis=1)]
    usable = {ep: s[ep].dropna() for ep in ENDPOINTS}
    usable = {ep: v for ep, v in usable.items() if len(v) >= 8 and v.nunique() == 2 and v.value_counts().min() >= 3}
    if not usable:
        print(f"  {cid}: no endpoint with >=3 patients per group, skipped", file=sys.stderr)
        return None
    therapy = " / ".join(f"{k} ({v})" for k, v in s["ICI_RX"].fillna("none").value_counts().items()) if "ICI_RX" in s else ""

    # ---- expression: genes x samples, log2(x+1) if linear, z-score per gene within cohort
    ex = ex.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    M = ex[s["SAMPLE_ID"].tolist()].apply(pd.to_numeric, errors="coerce")
    if np.nanmax(M.to_numpy()) > 50:
        M = np.log2(M.clip(lower=0) + 1)
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)

    scores, used = {}, {}
    for k, sig in sigs["computed"].items():
        g = [x for x in sig["genes"] if x in Z.index]
        used[k] = g
        scores[k] = Z.loc[g].mean(axis=0) if len(g) >= max(1, len(sig["genes"]) // 2) else pd.Series(np.nan, index=Z.columns)
    gsi = gs.set_index("ENTITY_STABLE_ID")
    for k in sigs["iatlas"]:
        scores[f"ia_{k}"] = pd.to_numeric(gsi.loc[k, s["SAMPLE_ID"]], errors="coerce") if k in gsi.index else pd.Series(np.nan, index=Z.columns)
    panel = {g: Z.loc[g] for g in sigs["panel"] if g in Z.index}

    stats = {}
    for ep in usable:
        resp = s.set_index("SAMPLE_ID")[ep].dropna().astype(bool)
        stats[ep] = {}
        for k, v in scores.items():
            v = v.reindex(resp.index).astype(float)
            ok = v.notna()
            pos, neg = v[ok & resp].to_numpy(), v[ok & ~resp].to_numpy()
            if len(pos) >= 3 and len(neg) >= 3:
                a, lo, hi, p = auc_ci(pos, neg, rng)
                stats[ep][k] = {"auc": a, "lo": lo, "hi": hi, "p": p, "nR": int(len(pos)), "nNR": int(len(neg))}
    lab = lambda v: None if v is None or pd.isna(v) else ("Responder" if v else "Non-responder")

    r3 = lambda x: None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), 3)
    rows = []
    for _, r in s.iterrows():
        sid = r["SAMPLE_ID"]
        rows.append({
            "sample": sid, "patient": r["PATIENT_ID"], "response": lab(r.get("response")), "benefit": lab(r.get("benefit")),
            "recist": None if pd.isna(r.get("RESPONSE")) else str(r.get("RESPONSE")),
            "os": r3(r.get("OS_MONTHS")) if "OS_MONTHS" in r else None,
            "osEvent": (str(r.get("OS_STATUS", "")).startswith("1")) if "OS_STATUS" in r and pd.notna(r.get("OS_STATUS")) else None,
            "scores": {k: r3(v.get(sid)) for k, v in scores.items()},
            "genes": {g: r3(z.get(sid)) for g, z in panel.items()},
        })
    counts = {ep: {"R": int(usable[ep].astype(bool).sum()), "NR": int((~usable[ep].astype(bool)).sum())} for ep in usable}
    meta = {"id": cid, "label": label, "cancer": cancer, "therapy": therapy, "name": md.get("name", ""),
            "citation": md.get("citation", ""), "pmid": md.get("pmid", "").split(",")[0] or None,
            "source": f"cBioPortal datahub / CRI iAtlas ({cid})", "license": "ODbL 1.0",
            "expressionFile": expr_name, "n": len(rows), "endpoints": counts, "genesUsed": used}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{cid}.json").write_text(json.dumps({"meta": meta, "stats": stats, "samples": rows}, separators=(",", ":")))
    print(f"  {cid}: {len(rows)} patients; " + "; ".join(f"{e}: {c['R']} R / {c['NR']} NR" for e, c in counts.items()), file=sys.stderr)
    return meta, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args(argv)
    sigs = json.loads((Path(__file__).with_name("signatures.json")).read_text())
    rng = np.random.default_rng(7)
    index, summary = [], {}
    prev = json.loads((OUT / "index.json").read_text()) if a.only and (OUT / "index.json").exists() else None
    for c in COHORTS:
        if a.only and c[0] not in a.only:
            continue
        res = process(*c, sigs, rng)
        if res:
            meta, stats = res
            index.append({k: meta[k] for k in ("id", "label", "cancer", "therapy", "name", "citation", "pmid", "n", "endpoints")})
            summary[meta["id"]] = stats
    if prev:  # merge partial runs
        keep = [d for d in prev["cohorts"] if d["id"] not in {i["id"] for i in index}]
        index = keep + index
        summary = {**json.loads((OUT / "summary.json").read_text())["stats"], **summary}
    labels = {k: v["label"] for k, v in sigs["computed"].items()} | {f"ia_{k}": v for k, v in sigs["iatlas"].items()}
    sources = {k: v["source"] for k, v in sigs["computed"].items()} | {f"ia_{k}": "precomputed by CRI iAtlas" for k in sigs["iatlas"]}
    order = [c[0] for c in COHORTS]
    index.sort(key=lambda d: order.index(d["id"]))
    (OUT / "index.json").write_text(json.dumps({"cohorts": index, "signatures": labels, "signatureSources": sources,
                                                "endpoints": {k: v[1] for k, v in ENDPOINTS.items()},
                                                "panel": sigs["panel"], "license": "ODbL 1.0 (cBioPortal datahub, CRI iAtlas)"}, indent=1))
    (OUT / "summary.json").write_text(json.dumps({"stats": summary}, separators=(",", ":")))
    print(f"wrote {len(index)} cohorts to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
