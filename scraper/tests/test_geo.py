"""Offline tests for the GEO watcher: python -m unittest discover -s scraper/tests"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fetch_geo as fg  # noqa: E402


class TestGeoWatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        subprocess.run([sys.executable, HERE / "make_geo_fixture.py", cls.tmp / "fx"], check=True, capture_output=True)
        cls.out = cls.tmp / "out"; cls.out.mkdir()
        # a catalog that already lists one fixture accession
        (cls.out / "catalog.json").write_text(json.dumps({"datasets": [
            {"accessions": [{"label": "GEO GSE9000003 (TCR)", "url": "x"}]}]}))
        fg.main(["--fixture", str(cls.tmp / "fx"), "--out", str(cls.out), "--today", "2026-10-12", "--catalog", str(cls.out / "catalog.json")])
        fg.main(["--fixture", str(cls.tmp / "fx"), "--out", str(cls.out), "--today", "2026-10-19", "--catalog", str(cls.out / "catalog.json")])
        cls.watch = json.loads((cls.out / "geo_watch.json").read_text())

    def test_dedup_across_queries(self):
        accs = [i["accession"] for i in self.watch["items"]]
        self.assertEqual(len(accs), len(set(accs)))
        visium = next(i for i in self.watch["items"] if "Visium" in i["title"])
        self.assertEqual(visium["queries"], ["mel_spatial", "til_sc"])

    def test_first_seen_is_kept(self):
        self.assertTrue(all(i["firstSeen"] == "2026-10-12" for i in self.watch["items"]))
        self.assertEqual(self.watch["meta"]["newThisRun"], 0)

    def test_catalog_flag(self):
        flagged = [i["accession"] for i in self.watch["items"] if i["inCatalog"]]
        self.assertEqual(flagged, ["GSE9000003"])

    def test_tagging(self):
        visium = next(i for i in self.watch["items"] if "Visium" in i["title"])
        self.assertIn("spatial", visium["modalities"])
        self.assertIn("TIL-ACT", visium["therapies"])


if __name__ == "__main__":
    unittest.main()
