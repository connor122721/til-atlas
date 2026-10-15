#!/usr/bin/env python3
"""Weekly ClinicalTrials.gov scraper for the TIL Atlas site.

What it does, in order:
  1. For each topic in config.json, page through /api/v2/studies (cursor = nextPageToken).
  2. Flatten each study's nested JSON into one small record (see `normalize`).
  3. Tag every record with indication / combination / engineering labels (tags.json).
  4. Diff against the previous snapshot -> new trials, status changes, removals.
  5. Write site/data/trials.json, append to history.json and changelog.json.

Standard library only, so the GitHub Action needs no `pip install`.

Usage:
  python scraper/fetch_trials.py                 # live run, writes site/data/
  python scraper/fetch_trials.py --fixture DIR   # offline run from saved API pages (tests, demos)
  python scraper/fetch_trials.py --dry-run       # fetch + summarize, write nothing
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://clinicaltrials.gov/api/v2/studies"
ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
DATA = ROOT / "site" / "data"
MAX_CHANGELOG = 104  # two years of weekly entries


# --------------------------------------------------------------------------- fetch
def get_json(url: str, user_agent: str, retries: int = 5) -> dict:
    """GET with exponential backoff on 429 / 5xx (the API rate-limits bursts)."""
    req = urllib.request.Request(url, headers={"User-Agent": user_agent, "Accept": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                wait = 2 ** attempt * 3
                print(f"  HTTP {e.code}; retrying in {wait}s", file=sys.stderr)
                time.sleep(wait)
                continue
            raise
        except urllib.error.URLError:
            if attempt < retries - 1:
                time.sleep(2 ** attempt * 3)
                continue
            raise
    raise RuntimeError("unreachable")


def fetch_topic(topic: dict, page_size: int, user_agent: str):
    """Yield raw study dicts for one topic, following the pagination cursor."""
    params = dict(topic["params"], pageSize=page_size, format="json", countTotal="true")
    page = 0
    while True:
        url = f"{API}?{urllib.parse.urlencode(params)}"
        data = get_json(url, user_agent)
        if page == 0:
            print(f"  [{topic['id']}] totalCount={data.get('totalCount')}", file=sys.stderr)
        yield from data.get("studies", [])
        token = data.get("nextPageToken")
        if not token:
            break
        params["pageToken"] = token
        page += 1
        time.sleep(0.5)  # be polite


def fetch_fixture(fixture_dir: Path, topic_id: str):
    """Offline mode: read <fixture_dir>/<topic_id>_page*.json in order."""
    for p in sorted(fixture_dir.glob(f"{topic_id}_page*.json")):
        yield from json.loads(p.read_text())["studies"]


# ----------------------------------------------------------------------- normalize
def _g(d: dict, *path, default=None):
    """Safe nested get: _g(study, 'protocolSection', 'statusModule', 'overallStatus')."""
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def normalize(study: dict) -> dict:
    """Flatten one API v2 study into the record the website reads.

    The v2 JSON nests everything under protocolSection.<module>. Field paths are
    documented at clinicaltrials.gov/data-api/about-api/study-data-structure.
    """
    ps = study.get("protocolSection", {})
    ident, status = ps.get("identificationModule", {}), ps.get("statusModule", {})
    design, sponsor = ps.get("designModule", {}), ps.get("sponsorCollaboratorsModule", {})

    locations = []
    for loc in _g(ps, "contactsLocationsModule", "locations", default=[]) or []:
        gp = loc.get("geoPoint") or {}
        locations.append({
            "facility": loc.get("facility"),
            "city": loc.get("city"),
            "country": loc.get("country"),
            "lat": round(gp["lat"], 3) if "lat" in gp else None,
            "lon": round(gp["lon"], 3) if "lon" in gp else None,
        })

    return {
        "id": ident.get("nctId"),
        "title": ident.get("briefTitle"),
        "acronym": ident.get("acronym"),
        "status": status.get("overallStatus"),
        "phases": design.get("phases") or [],
        "type": design.get("studyType"),
        "enrollment": _g(design, "enrollmentInfo", "count"),
        "start": _g(status, "startDateStruct", "date"),
        "primaryCompletion": _g(status, "primaryCompletionDateStruct", "date"),
        "firstPosted": _g(status, "studyFirstPostDateStruct", "date"),
        "lastUpdate": _g(status, "lastUpdatePostDateStruct", "date"),
        "sponsor": _g(sponsor, "leadSponsor", "name"),
        "sponsorClass": _g(sponsor, "leadSponsor", "class"),
        "conditions": _g(ps, "conditionsModule", "conditions", default=[]) or [],
        "keywords": _g(ps, "conditionsModule", "keywords", default=[]) or [],
        "interventions": [i.get("name") for i in _g(ps, "armsInterventionsModule", "interventions", default=[]) or []],
        "hasResults": bool(study.get("hasResults")),
        "locations": locations,
    }


def compile_tags(tag_cfg: dict) -> dict:
    return {
        group: {label: re.compile(rx, re.I) for label, rx in rules.items()}
        for group, rules in tag_cfg.items() if not group.startswith("_")
    }


def tag(rec: dict, rules: dict) -> dict:
    text = " | ".join([rec["title"] or "", *rec["conditions"], *rec["interventions"], *rec["keywords"]])
    rec["tags"] = {g: [lab for lab, rx in rs.items() if rx.search(text)] for g, rs in rules.items()}
    if not rec["tags"].get("indication"):
        rec["tags"]["indication"] = ["Other / basket"]
    return rec


# ---------------------------------------------------------------------------- diff
def diff(old: dict[str, dict], new: dict[str, dict]) -> dict:
    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    status_changes = [
        {"id": k, "from": old[k]["status"], "to": new[k]["status"]}
        for k in sorted(set(old) & set(new))
        if old[k]["status"] != new[k]["status"]
    ]
    return {"new": added, "removed": removed, "statusChanges": status_changes}


def count_by(records, key) -> dict:
    out: dict[str, int] = {}
    for r in records:
        out[r[key] or "UNKNOWN"] = out.get(r[key] or "UNKNOWN", 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


# ---------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixture", type=Path, help="read saved API pages instead of the network")
    ap.add_argument("--out", type=Path, default=DATA, help="output directory (default site/data)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--today", help="override run date (YYYY-MM-DD), for tests")
    args = ap.parse_args(argv)

    cfg = json.loads((HERE / "config.json").read_text())
    rules = compile_tags(json.loads((HERE / "tags.json").read_text()))
    today = args.today or dt.date.today().isoformat()

    trials: dict[str, dict] = {}
    for topic in cfg["topics"]:
        print(f"Fetching topic {topic['id']} …", file=sys.stderr)
        src = (fetch_fixture(args.fixture, topic["id"]) if args.fixture
               else fetch_topic(topic, cfg["page_size"], cfg["user_agent"]))
        for s in src:
            rec = tag(normalize(s), rules)
            if not rec["id"]:
                continue
            # A trial can match several topics; keep one record, remember every topic.
            prev = trials.get(rec["id"])
            rec["topics"] = sorted(set((prev or {}).get("topics", [])) | {topic["id"]})
            trials[rec["id"]] = rec

    records = sorted(trials.values(), key=lambda r: (r["firstPosted"] or ""), reverse=True)
    print(f"Collected {len(records)} unique trials", file=sys.stderr)
    if args.dry_run:
        print(json.dumps(count_by(records, "status"), indent=2))
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    trials_path, hist_path, log_path = (args.out / n for n in ("trials.json", "history.json", "changelog.json"))

    old = {}
    if trials_path.exists():
        old = {r["id"]: r for r in json.loads(trials_path.read_text()).get("trials", [])}
    change = diff(old, trials) if old else {"new": [], "removed": [], "statusChanges": [], "initial": True}

    meta = {
        "generated": today,
        "source": "ClinicalTrials.gov API v2",
        "topics": [{"id": t["id"], "label": t["label"], "params": t["params"]} for t in cfg["topics"]],
        "fixture": bool(args.fixture),
        "count": len(records),
    }
    trials_path.write_text(json.dumps({"meta": meta, "trials": records}, separators=(",", ":")))

    history = json.loads(hist_path.read_text()) if hist_path.exists() else []
    history = [h for h in history if h["date"] != today]  # re-runs on the same day replace
    history.append({"date": today, "total": len(records), "byStatus": count_by(records, "status")})
    hist_path.write_text(json.dumps(history, indent=1))

    log = json.loads(log_path.read_text()) if log_path.exists() else []
    log = [e for e in log if e["date"] != today]
    log.insert(0, {"date": today, **change})
    log_path.write_text(json.dumps(log[:MAX_CHANGELOG], indent=1))

    print(f"new={len(change['new'])} statusChanges={len(change['statusChanges'])} "
          f"removed={len(change['removed'])}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
