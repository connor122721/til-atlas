#!/usr/bin/env python3
"""SYNTHETIC E-utilities responses (esearch + esummary JSON shapes) for offline tests and the demo.
Accessions start with GSE9000 and titles with [FIXTURE]; nothing here is a real dataset.

  python scraper/tests/make_geo_fixture.py OUTDIR
"""
import json
import sys
from pathlib import Path

ITEMS = {
    "til_sc": [
        ("Single-cell RNA and TCR sequencing of TIL infusion products from responders and non-responders",
         "Expression profiling by high throughput sequencing", 24, "2026/09/30"),
        ("Spatial transcriptomics (Visium) of melanoma metastases before and after adoptive cell transfer",
         "Expression profiling by high throughput sequencing", 12, "2026/09/12"),
        ("TCR repertoire of expanded tumor-infiltrating lymphocytes in NSCLC",
         "Other", 40, "2026/08/02"),
    ],
    "mel_spatial": [
        ("Xenium profiling of melanoma under anti-PD-1 therapy",
         "Expression profiling by high throughput sequencing", 18, "2026/10/01"),
        ("Spatial transcriptomics (Visium) of melanoma metastases before and after adoptive cell transfer",
         "Expression profiling by high throughput sequencing", 12, "2026/09/12"),
    ],
    "cd8_stem": [
        ("Single-cell atlas of TCF7+ stem-like CD8 T cells across solid tumors",
         "Expression profiling by high throughput sequencing", 60, "2026/07/21"),
    ],
}


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    acc_of = {}
    for qid, rows in ITEMS.items():
        uids, result = [], {}
        for title, gtype, n, pdat in rows:
            acc = acc_of.setdefault(title, f"GSE9000{len(acc_of) + 1:03d}")
            uid = "2" + acc[3:]
            uids.append(uid)
            result[uid] = {"uid": uid, "accession": acc, "title": f"[FIXTURE] {title}",
                           "summary": f"Synthetic test record. {title}.", "gdstype": gtype,
                           "n_samples": n, "pdat": pdat, "pubmedids": []}
        result["uids"] = uids
        (out / f"{qid}_esearch.json").write_text(json.dumps({"esearchresult": {"count": str(len(uids)), "idlist": uids}}))
        (out / f"{qid}_esummary.json").write_text(json.dumps({"result": result}))
    print(f"wrote GEO fixture to {out}")


if __name__ == "__main__":
    main()
