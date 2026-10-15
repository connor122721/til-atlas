#!/usr/bin/env python3
"""Weekly GEO watch: find new public GEO series that might belong in the dataset catalog.

For each query in config.json["geo"]["queries"]:
  1. esearch  db=gds, human series (gse[ETYP]) published in the last `lookback_days`
  2. esummary for the returned UIDs -> accession, title, summary, sample count, date, PubMed IDs
Items accumulate in site/data/geo_watch.json (newest first, with the date first seen), and
anything already listed in site/data/catalog.json is flagged so you only review what's new.

NCBI asks scripted users to identify themselves: set NCBI_EMAIL (and optionally NCBI_API_KEY)
as repository secrets; they are passed as the `email` / `api_key` parameters.

  python scraper/fetch_geo.py                    # live
  python scraper/fetch_geo.py --fixture DIR      # offline: DIR/<query_id>_esearch.json + _esummary.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_trials import get_json  # noqa: E402  (shared retry/backoff)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "site" / "data"

MODALITY = {
    "spatial": r"spatial|visium|xenium|cosmx|geomx|merfish|imaging mass cytometry|codex|slide-seq|stereo-seq",
    "single-cell": r"single[- ]cell|scrna|10x|smart-seq|cite-seq",
    "TCR": r"\btcr\b|t[- ]cell receptor|clonotyp|repertoire",
    "bulk RNA": r"bulk rna|rna-seq",
}
THERAPY = {
    "TIL-ACT": r"tumou?r[- ]infiltrating lymphocyte|\btil\b|adoptive (cell|t[- ]cell)",
    "Checkpoint": r"pd-1|pd-l1|ctla-4|checkpoint|nivolumab|pembrolizumab|ipilimumab",
    "CAR-T": r"\bcar[- ]t\b|chimeric antigen receptor",
}


def params(extra: dict) -> str:
    p = {"db": "gds", "retmode": "json", "tool": "til-atlas", **extra}
    if os.environ.get("NCBI_EMAIL"):
        p["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        p["api_key"] = os.environ["NCBI_API_KEY"]
    return urllib.parse.urlencode(p)


def search(q: dict, lookback: int, ua: str) -> list[str]:
    since = (dt.date.today() - dt.timedelta(days=lookback)).strftime("%Y/%m/%d")
    term = f'({q["term"]}) AND gse[ETYP] AND "Homo sapiens"[ORGN] AND ("{since}"[PDAT] : "3000"[PDAT])'
    data = get_json(f"{EUTILS}/esearch.fcgi?{params({'term': term, 'retmax': 200})}", ua)
    return data.get("esearchresult", {}).get("idlist", [])


def summarize(uids: list[str], ua: str) -> list[dict]:
    out = []
    for i in range(0, len(uids), 100):
        data = get_json(f"{EUTILS}/esummary.fcgi?{params({'id': ','.join(uids[i:i + 100])})}", ua)
        res = data.get("result", {})
        out += [res[u] for u in res.get("uids", [])]
        time.sleep(0.4)  # stay under 3 requests/second without an API key
    return out


def normalize(doc: dict, query_id: str) -> dict:
    acc = doc.get("accession") or f"GSE{doc.get('gse', '')}"
    text = f"{doc.get('title', '')} {doc.get('summary', '')} {doc.get('gdstype', '')}".lower()
    return {
        "accession": acc,
        "title": doc.get("title", "").strip(),
        "summary": re.sub(r"\s+", " ", doc.get("summary", "")).strip()[:600],
        "types": doc.get("gdstype", ""),
        "samples": doc.get("n_samples"),
        "published": (doc.get("pdat") or "").replace("/", "-"),
        "pubmed": [str(p) for p in doc.get("pubmedids", [])],
        "modalities": [k for k, rx in MODALITY.items() if re.search(rx, text)],
        "therapies": [k for k, rx in THERAPY.items() if re.search(rx, text)],
        "queries": [query_id],
        "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
    }


def load_fixture(fx: Path, qid: str, kind: str) -> dict:
    return json.loads((fx / f"{qid}_{kind}.json").read_text())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", type=Path)
    ap.add_argument("--out", type=Path, default=DATA)
    ap.add_argument("--today")
    ap.add_argument("--catalog", type=Path, default=DATA / "catalog.json", help="accessions already curated")
    args = ap.parse_args(argv)
    cfg = json.loads((HERE / "config.json").read_text())
    geo, ua = cfg["geo"], cfg["user_agent"]
    today = args.today or dt.date.today().isoformat()

    found: dict[str, dict] = {}
    for q in geo["queries"]:
        if args.fixture:
            uids = load_fixture(args.fixture, q["id"], "esearch")["esearchresult"]["idlist"]
            res = load_fixture(args.fixture, q["id"], "esummary")["result"]
            docs = [res[u] for u in res.get("uids", [])]
        else:
            uids = search(q, geo["lookback_days"], ua)
            docs = summarize(uids, ua) if uids else []
            time.sleep(0.4)
        print(f"  [{q['id']}] {len(uids)} series", file=sys.stderr)
        for d in docs:
            item = normalize(d, q["id"])
            if item["accession"] in found:
                found[item["accession"]]["queries"] = sorted(set(found[item["accession"]]["queries"]) | {q["id"]})
            else:
                found[item["accession"]] = item

    catalog = json.loads(args.catalog.read_text()) if args.catalog.exists() else {"datasets": []}
    known = {m.group(0) for ds in catalog["datasets"] for a in ds["accessions"] for m in re.finditer(r"GSE\d+", a["label"])}

    path = args.out / "geo_watch.json"
    prev = {i["accession"]: i for i in json.loads(path.read_text())["items"]} if path.exists() else {}
    for acc, item in found.items():
        item["firstSeen"] = prev.get(acc, {}).get("firstSeen", today)
        prev[acc] = {**prev.get(acc, {}), **item}
    for item in prev.values():
        item["inCatalog"] = item["accession"] in known
    items = sorted(prev.values(), key=lambda i: (i["firstSeen"], i["published"]), reverse=True)[: geo["max_items"]]
    new = [i for i in items if i["firstSeen"] == today]
    meta = {"generated": today, "fixture": bool(args.fixture), "lookbackDays": geo["lookback_days"],
            "queries": [{"id": q["id"], "label": q["label"]} for q in geo["queries"]], "newThisRun": len(new)}
    args.out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"meta": meta, "items": items}, indent=1))
    print(f"GEO watch: {len(found)} matched, {len(new)} new, {len(items)} kept", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
