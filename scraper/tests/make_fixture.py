#!/usr/bin/env python3
"""Generate SYNTHETIC ClinicalTrials.gov API v2 pages for offline tests and the UI demo.

Everything here is invented: NCT IDs start with NCT9 (real IDs don't, yet), titles
start with "[FIXTURE]", sponsors are fictional. The JSON *shape* mirrors the real
v2 response (protocolSection.<module>...), so the scraper code path is identical.

  python scraper/tests/make_fixture.py OUTDIR [--week 2]

--week 2 applies a few changes (new trial, status flips, one removal) so the diff
logic has something to detect.
"""
import json
import random
import sys
from pathlib import Path

CITIES = [  # (facility, city, country, lat, lon) - fictional facility names, real coordinates
    ("Example Cancer Center A", "Bethesda", "United States", 38.984, -77.094),
    ("Example Cancer Center B", "Tampa", "United States", 27.951, -82.457),
    ("Example Cancer Center C", "Houston", "United States", 29.760, -95.370),
    ("Example Cancer Center D", "Louisville", "United States", 38.253, -85.759),
    ("Example Cancer Center E", "Lausanne", "Switzerland", 46.520, 6.633),
    ("Example Cancer Center F", "Copenhagen", "Denmark", 55.676, 12.568),
    ("Example Cancer Center G", "Amsterdam", "Netherlands", 52.370, 4.895),
    ("Example Cancer Center H", "Shanghai", "China", 31.230, 121.474),
    ("Example Cancer Center I", "Toronto", "Canada", 43.653, -79.383),
    ("Example Cancer Center J", "Tel Aviv", "Israel", 32.085, 34.782),
    ("Example Cancer Center K", "Manchester", "United Kingdom", 53.481, -2.243),
    ("Example Cancer Center L", "Seattle", "United States", 47.606, -122.332),
]
CONDS = ["Melanoma", "Cervical Cancer", "Non-small Cell Lung Cancer", "Head and Neck Squamous Cell Carcinoma",
         "Ovarian Cancer", "Colorectal Cancer", "Breast Cancer", "Sarcoma", "Solid Tumor"]
INTERVENTIONS = [
    ["Tumor Infiltrating Lymphocytes", "Cyclophosphamide", "Fludarabine", "Aldesleukin"],
    ["Tumor-infiltrating lymphocytes", "Pembrolizumab", "Cyclophosphamide", "Fludarabine"],
    ["Gene-edited TIL (PD-1 knockout)", "Fludarabine"],
    ["Membrane-bound IL-15 armored TIL", "Cyclophosphamide"],
    ["Neoantigen-selected TIL", "Aldesleukin", "Nivolumab"],
    ["TIL", "Ipilimumab", "Nivolumab", "Aldesleukin"],
]
STATUSES = ["RECRUITING"] * 5 + ["ACTIVE_NOT_RECRUITING"] * 2 + ["COMPLETED"] * 3 + \
           ["NOT_YET_RECRUITING", "TERMINATED", "WITHDRAWN", "UNKNOWN"]
PHASES = [["PHASE1"], ["PHASE1", "PHASE2"], ["PHASE2"], ["PHASE2"], ["PHASE3"]]
SPONSOR = [("Fictional University Hospital", "OTHER"), ("Example Bio Inc.", "INDUSTRY"),
           ("National Example Institute", "NIH"), ("Demo Therapeutics", "INDUSTRY")]


def study(i: int, rng: random.Random, topic: str) -> dict:
    year = rng.randint(2010, 2026)
    cond = rng.choice(CONDS)
    interv = rng.choice(INTERVENTIONS) if topic == "til" else ["Neoantigen TCR-T cells", "Aldesleukin"]
    locs = rng.sample(CITIES, rng.randint(1, 4))
    return {
        "protocolSection": {
            "identificationModule": {"nctId": f"NCT9{i:07d}",
                                     "briefTitle": f"[FIXTURE] {interv[0]} for {cond} (demo study {i})"},
            "statusModule": {
                "overallStatus": rng.choice(STATUSES),
                "startDateStruct": {"date": f"{year}-{rng.randint(1, 12):02d}"},
                "primaryCompletionDateStruct": {"date": f"{min(year + rng.randint(2, 6), 2032)}-06"},
                "studyFirstPostDateStruct": {"date": f"{year}-{rng.randint(1, 12):02d}-{rng.randint(1, 28):02d}"},
                "lastUpdatePostDateStruct": {"date": f"2026-{rng.randint(1, 9):02d}-{rng.randint(1, 28):02d}"},
            },
            "sponsorCollaboratorsModule": {"leadSponsor": dict(zip(("name", "class"), rng.choice(SPONSOR)))},
            "conditionsModule": {"conditions": [cond], "keywords": ["adoptive cell therapy"]},
            "designModule": {"studyType": "INTERVENTIONAL", "phases": rng.choice(PHASES),
                             "enrollmentInfo": {"count": rng.choice([12, 20, 30, 48, 60, 100, 180])}},
            "armsInterventionsModule": {"interventions": [{"type": "BIOLOGICAL", "name": n} for n in interv]},
            "contactsLocationsModule": {"locations": [
                {"facility": f, "city": c, "country": co, "geoPoint": {"lat": la, "lon": lo}}
                for f, c, co, la, lo in locs]},
        },
        "hasResults": rng.random() < 0.15,
    }


def main():
    out = Path(sys.argv[1])
    week = int(sys.argv[sys.argv.index("--week") + 1]) if "--week" in sys.argv else 1
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(1328)  # fixed seed -> identical fixture every time
    til = [study(i, rng, "til") for i in range(1, 49)]
    tcrt = [study(i, rng, "tcrt") for i in range(49, 57)]
    tcrt.append(til[3])  # one trial matching both topics -> tests de-duplication
    if week == 2:
        til[0]["protocolSection"]["statusModule"]["overallStatus"] = "COMPLETED"
        til[5]["protocolSection"]["statusModule"]["overallStatus"] = "TERMINATED"
        del til[10]
        til.append(study(99, random.Random(99), "til"))
    # split TIL into two pages so the pagination path is exercised
    (out / "til_page1.json").write_text(json.dumps({"studies": til[:30], "nextPageToken": "fixture-token"}))
    (out / "til_page2.json").write_text(json.dumps({"studies": til[30:]}))
    (out / "tcrt_page1.json").write_text(json.dumps({"studies": tcrt}))
    print(f"wrote week-{week} fixture to {out}")


if __name__ == "__main__":
    main()
