"""The build script's output feeds the verification harness, end to end.

A fake beacon is loaded with exactly the sites of the built member VCF, then the
harness is run on the build directory. This is the check that the two halves
agree on file names, coordinates and what "present" and "absent" mean.
"""

import io
import json
import os
import shutil
import subprocess
import tempfile
import unittest
import urllib.parse
from contextlib import redirect_stdout

import beacon_check
import test_build_population_dataset as toy


@unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
class TestBuildThenCheck(unittest.TestCase):

    def setUp(self):
        helper = toy.TestBuild()
        helper.build()
        self.addCleanup(shutil.rmtree, helper.tmp, ignore_errors=True)
        self.out = helper.out
        with open(os.path.join(self.out, "build_summary.json")) as fh:
            self.dataset = json.load(fh)["dataset_id"]
        rows = subprocess.run(
            ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
             os.path.join(self.out, "members.vcf.gz")],
            capture_output=True, text=True, check=True).stdout.splitlines()
        self.held = {(self.dataset, c.removeprefix("chr"), int(p), r, a)
                     for c, p, r, a in (row.split("\t") for row in rows)}

    def beacon(self, held, everywhere=False):
        def fetch(url):
            q = dict(urllib.parse.parse_qsl(urllib.parse.urlsplit(url).query))
            key = (q["datasetIds"], q["referenceName"], int(q["start"]) + 1,
                   q["referenceBases"], q["alternateBases"])
            return 200, {}, {"responseSummary": {
                "exists": everywhere or key in held}}
        return fetch

    def check(self, fetch):
        out = os.path.join(self.out, "report.json")
        with redirect_stdout(io.StringIO()):
            code = beacon_check.main(
                ["--base-url", "http://b", "--build-dir", self.out,
                 "--per-individual", "10", "--out", out],
                fetch=fetch, sleep=lambda s: None)
        with open(out) as fh:
            return code, json.load(fh)

    def test_a_beacon_loaded_from_the_member_vcf_passes_every_gate(self):
        code, report = self.check(self.beacon(self.held))
        self.assertEqual(code, 0, report)
        self.assertEqual(report["individuals"]["member_yes_rate"], 1.0)

    def test_the_whole_panel_beacon_fails_on_the_built_lists(self):
        code, report = self.check(self.beacon(self.held, everywhere=True))
        self.assertEqual(code, 1)
        self.assertFalse(report["passed"])

    def test_a_beacon_that_dropped_a_member_site_fails(self):
        held = {k for k in self.held if k[2] != 300}
        code, report = self.check(self.beacon(held))
        self.assertEqual(code, 1, report)


if __name__ == "__main__":
    unittest.main()
