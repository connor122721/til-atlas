"""Offline tests: python -m unittest discover -s scraper/tests"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fetch_trials as ft  # noqa: E402


class TestScraper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        for wk in (1, 2):
            subprocess.run([sys.executable, HERE / "make_fixture.py", cls.tmp / f"wk{wk}", "--week", str(wk)],
                           check=True, capture_output=True)
        cls.out = cls.tmp / "out"
        ft.main(["--fixture", str(cls.tmp / "wk1"), "--out", str(cls.out), "--today", "2026-10-07"])
        cls.wk1 = json.loads((cls.out / "trials.json").read_text())
        ft.main(["--fixture", str(cls.tmp / "wk2"), "--out", str(cls.out), "--today", "2026-10-14"])
        cls.wk2 = json.loads((cls.out / "trials.json").read_text())
        cls.log = json.loads((cls.out / "changelog.json").read_text())
        cls.hist = json.loads((cls.out / "history.json").read_text())

    def test_pagination_and_dedup(self):
        # 48 TIL + 8 TCR-T, one TIL trial also in TCR-T -> 56 unique
        self.assertEqual(self.wk1["meta"]["count"], 56)
        both = [t for t in self.wk1["trials"] if t["topics"] == ["tcrt", "til"]]
        self.assertEqual(len(both), 1)

    def test_normalize_fields(self):
        t = self.wk1["trials"][0]
        for k in ("id", "title", "status", "phases", "sponsor", "locations", "tags", "firstPosted"):
            self.assertIn(k, t)
        self.assertTrue(t["id"].startswith("NCT9"))

    def test_tagging(self):
        rules = ft.compile_tags(json.loads((HERE.parent / "tags.json").read_text()))
        rec = ft.tag({"title": "TIL plus pembrolizumab in metastatic melanoma", "conditions": [],
                      "interventions": ["Fludarabine"], "keywords": []}, rules)
        self.assertIn("Melanoma", rec["tags"]["indication"])
        self.assertIn("anti-PD-1", rec["tags"]["combination"])
        self.assertIn("Lymphodepletion", rec["tags"]["combination"])

    def test_diff(self):
        latest = self.log[0]
        self.assertEqual(latest["date"], "2026-10-14")
        self.assertEqual(len(latest["new"]), 1)
        self.assertEqual(len(latest["removed"]), 1)
        self.assertGreaterEqual(len(latest["statusChanges"]), 1)
        self.assertTrue(self.log[1].get("initial"))

    def test_history_appends(self):
        self.assertEqual([h["date"] for h in self.hist], ["2026-10-07", "2026-10-14"])


if __name__ == "__main__":
    unittest.main()
