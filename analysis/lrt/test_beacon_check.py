"""Tests for the harness that checks a loaded beacon against genotype truth.

The beacon is replaced by a small fake, so the harness is tested on defects it
must catch: one that answers yes everywhere (the defect found on
afrigend-all-dev), one that is off by one base, and one that ignores which
dataset was asked about.
"""

import io
import json
import os
import tempfile
import unittest
import urllib.parse
from contextlib import redirect_stdout

import beacon_check

PRESENT = [("chr1", 100, "A", "G"), ("chr1", 300, "A", "G"), ("chr1", 700, "G", "A")]
ABSENT = [("chr1", 200, "C", "T"), ("chr1", 600, "T", "C")]


class FakeBeacon:
    """Answers from a fixed set of (dataset, chrom, POS, ref, alt)."""

    def __init__(
        self,
        held,
        *,
        yes_everywhere=False,
        off_by_one=False,
        ignore_dataset=False,
        rate_limits=0,
    ):
        self.held = set(held)
        self.yes_everywhere = yes_everywhere
        self.off_by_one = off_by_one
        self.ignore_dataset = ignore_dataset
        self.rate_limits = rate_limits
        self.requests = []

    def __call__(self, url):
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlsplit(url).query))
        self.requests.append(q)
        if self.rate_limits > 0:
            self.rate_limits -= 1
            return 429, {"Retry-After": "7"}, None
        start = int(q["start"])
        pos = start + 1 + (1 if self.off_by_one else 0)
        key = (
            q["datasetIds"],
            q["referenceName"],
            pos,
            q["referenceBases"],
            q["alternateBases"],
        )
        if self.ignore_dataset:
            hit = any(k[1:] == key[1:] for k in self.held)
        else:
            hit = key in self.held
        return 200, {}, {"responseSummary": {"exists": hit or self.yes_everywhere}}


def _held(dataset="YRI_M3_S1", chrom="1"):
    return [(dataset, chrom, p, r, a) for _, p, r, a in PRESENT]


def _client(fake, **kw):
    return beacon_check.BeaconClient(
        "http://beacon.test", fetch=fake, sleep=lambda s: None, **kw
    )


class TestQuery(unittest.TestCase):
    def test_start_is_the_zero_based_position(self):
        fake = FakeBeacon(_held())
        _client(fake).exists("YRI_M3_S1", "chr1", 100, "A", "G")
        self.assertEqual(fake.requests[0]["start"], "99")

    def test_chromosome_loses_its_chr_prefix(self):
        fake = FakeBeacon(_held())
        _client(fake).exists("YRI_M3_S1", "chr1", 100, "A", "G")
        self.assertEqual(fake.requests[0]["referenceName"], "1")

    def test_asks_for_the_dataset_and_assembly(self):
        fake = FakeBeacon(_held())
        _client(fake).exists("YRI_M3_S1", "chr1", 100, "A", "G")
        self.assertEqual(fake.requests[0]["datasetIds"], "YRI_M3_S1")
        self.assertEqual(fake.requests[0]["assemblyId"], "GRCh38")

    def test_reads_yes_and_no(self):
        client = _client(FakeBeacon(_held()))
        self.assertTrue(client.exists("YRI_M3_S1", "chr1", 100, "A", "G"))
        self.assertFalse(client.exists("YRI_M3_S1", "chr1", 200, "C", "T"))

    def test_an_answer_that_is_not_a_boolean_is_an_error_not_a_yes_or_no(self):
        # "true" the string is truthy and "" is falsy: neither is an answer.
        for value in ("true", "", None, 1):
            def fetch(url, value=value):
                return 200, {}, {"responseSummary": {"exists": value}}

            with self.assertRaises(beacon_check.BeaconError, msg=repr(value)):
                _client(fetch).exists("D", "chr1", 100, "A", "G")

    def test_a_malformed_answer_is_an_error_not_a_no(self):
        def fetch(url):
            return 200, {}, {"responseSummary": {}}

        with self.assertRaises(beacon_check.BeaconError):
            _client(fetch).exists("D", "chr1", 100, "A", "G")

    def test_a_server_error_is_raised_not_retried(self):
        calls = []

        def fetch(url):
            calls.append(url)
            return 500, {}, None

        with self.assertRaises(beacon_check.BeaconError):
            _client(fetch).exists("D", "chr1", 100, "A", "G")
        self.assertEqual(len(calls), 1)


class TestRateLimit(unittest.TestCase):
    def test_waits_the_time_the_server_asks_for_then_retries(self):
        sleeps = []
        fake = FakeBeacon(_held(), rate_limits=2)
        client = beacon_check.BeaconClient("http://b", fetch=fake, sleep=sleeps.append)
        self.assertTrue(client.exists("YRI_M3_S1", "chr1", 100, "A", "G"))
        self.assertEqual(sleeps, [7, 7])

    def test_gives_up_after_a_bounded_number_of_retries(self):
        # A tight retry loop against a limiter has taken a beacon offline
        # before. Bounded, and every wait honours Retry-After.
        sleeps = []
        fake = FakeBeacon(_held(), rate_limits=10_000)
        client = beacon_check.BeaconClient(
            "http://b", fetch=fake, sleep=sleeps.append, max_retries=3
        )
        with self.assertRaises(beacon_check.RateLimited):
            client.exists("YRI_M3_S1", "chr1", 100, "A", "G")
        self.assertEqual(len(fake.requests), 4)
        self.assertTrue(all(s >= 7 for s in sleeps))
        # No pointless wait after the final refusal: 3 retries, 3 waits.
        self.assertEqual(len(sleeps), 3)

    def test_a_missing_retry_after_still_waits(self):
        sleeps = []

        def fetch(url):
            return (
                (429, {}, None)
                if len(sleeps) < 1
                else (200, {}, {"responseSummary": {"exists": True}})
            )

        client = beacon_check.BeaconClient("http://b", fetch=fetch, sleep=sleeps.append)
        client.exists("D", "chr1", 1, "A", "G")
        self.assertGreater(sleeps[0], 0)


class TestKnownSites(unittest.TestCase):
    def test_a_correct_beacon_passes(self):
        result = beacon_check.check_known_sites(
            _client(FakeBeacon(_held())), "YRI_M3_S1", PRESENT, ABSENT
        )
        self.assertTrue(result["passed"], result)
        self.assertEqual(result["present_missed"], [])
        self.assertEqual(result["absent_wrongly_present"], [])

    def test_a_beacon_that_says_yes_everywhere_fails_on_the_absent_sites(self):
        # The defect found on afrigend-all-dev, as a regression guard: every
        # site of the parent panel answered yes whether or not a member
        # carried it.
        result = beacon_check.check_known_sites(
            _client(FakeBeacon(_held(), yes_everywhere=True)),
            "YRI_M3_S1",
            PRESENT,
            ABSENT,
        )
        self.assertFalse(result["passed"])
        self.assertEqual(len(result["absent_wrongly_present"]), len(ABSENT))
        self.assertEqual(result["present_missed"], [])

    def test_a_beacon_missing_a_member_site_fails_on_the_present_sites(self):
        held = _held()[:-1]
        result = beacon_check.check_known_sites(
            _client(FakeBeacon(held)), "YRI_M3_S1", PRESENT, ABSENT
        )
        self.assertFalse(result["passed"])
        self.assertEqual(result["present_missed"], [PRESENT[-1]])

    def test_an_empty_present_list_is_refused(self):
        with self.assertRaises(ValueError):
            beacon_check.check_known_sites(
                _client(FakeBeacon(_held())), "YRI_M3_S1", [], ABSENT
            )

    def test_an_empty_absent_list_is_refused(self):
        # With nothing absent the negative control tests nothing, and "pass"
        # would be a statement about a test that never ran.
        with self.assertRaises(ValueError):
            beacon_check.check_known_sites(
                _client(FakeBeacon(_held())), "YRI_M3_S1", PRESENT, []
            )


class TestPositiveControl(unittest.TestCase):
    def test_is_the_first_query_sent(self):
        fake = FakeBeacon(_held())
        beacon_check.check_known_sites(_client(fake), "YRI_M3_S1", PRESENT, ABSENT)
        first = fake.requests[0]
        self.assertEqual((first["start"], first["alternateBases"]), ("99", "G"))

    def test_a_failing_control_stops_everything_after_one_query(self):
        fake = FakeBeacon([])
        with self.assertRaises(beacon_check.ControlFailed):
            beacon_check.check_known_sites(_client(fake), "YRI_M3_S1", PRESENT, ABSENT)
        self.assertEqual(len(fake.requests), 1)

    def test_an_off_by_one_beacon_is_caught_by_the_control(self):
        # VCF POS used as the 0-based start: every real variant answers no,
        # which would read as "this population is safer".
        fake = FakeBeacon(_held(), off_by_one=True)
        with self.assertRaises(beacon_check.ControlFailed):
            beacon_check.check_known_sites(_client(fake), "YRI_M3_S1", PRESENT, ABSENT)

    def test_a_beacon_that_ignores_the_dataset_is_caught(self):
        # Queries scoped to a dataset that does not hold the site must say no.
        held = _held("OTHER_DATASET")
        fake = FakeBeacon(held, ignore_dataset=True)
        with self.assertRaises(beacon_check.ScopingFailed):
            beacon_check.check_dataset_scoping(_client(fake), "YRI_M3_S1", PRESENT[0])


class TestIndividuals(unittest.TestCase):
    def het(self):
        return {
            "m1": [PRESENT[0], PRESENT[1]],
            "m2": [PRESENT[2]],
            "c1": [ABSENT[0], PRESENT[0]],
            "c2": [ABSENT[1]],
        }

    def test_every_member_site_answers_yes_and_controls_do_not(self):
        result = beacon_check.check_individuals(
            _client(FakeBeacon(_held())),
            "YRI_M3_S1",
            self.het(),
            members=["m1", "m2"],
            controls=["c1", "c2"],
            per_individual=10,
            seed=1,
        )
        self.assertEqual(result["member_yes_rate"], 1.0)
        self.assertEqual(result["member_misses"], [])
        self.assertAlmostEqual(result["control_yes_rate"], 1 / 3)
        self.assertTrue(result["passed"])

    def test_one_member_miss_fails_the_gate(self):
        # A member's own het site must be in the beacon; a miss is a load or
        # coordinate bug, never biology.
        held = [k for k in _held() if k[2] != 300]
        result = beacon_check.check_individuals(
            _client(FakeBeacon(held)),
            "YRI_M3_S1",
            self.het(),
            members=["m1", "m2"],
            controls=["c1", "c2"],
            per_individual=10,
            seed=1,
        )
        self.assertFalse(result["passed"])
        self.assertEqual(result["member_misses"], [("m1", PRESENT[1])])

    def test_controls_all_saying_yes_fails_the_gate(self):
        # Controls that always answer yes mean the beacon cannot tell members
        # from non-members: the whole-panel defect, seen from the other side.
        result = beacon_check.check_individuals(
            _client(FakeBeacon(_held(), yes_everywhere=True)),
            "YRI_M3_S1",
            self.het(),
            members=["m1", "m2"],
            controls=["c1", "c2"],
            per_individual=10,
            seed=1,
        )
        self.assertEqual(result["control_yes_rate"], 1.0)
        self.assertFalse(result["passed"])

    def test_caps_the_sites_per_individual_reproducibly(self):
        het = {
            "m1": [("chr1", p, "A", "G") for p in range(1, 41)],
            "c1": [("chr1", 5000, "A", "G")],
        }
        held = [("YRI_M3_S1", "1", p, "A", "G") for p in range(1, 41)]
        a, b = FakeBeacon(held), FakeBeacon(held)
        for fake in (a, b):
            beacon_check.check_individuals(
                _client(fake),
                "YRI_M3_S1",
                het,
                members=["m1"],
                controls=["c1"],
                per_individual=5,
                seed=3,
            )
        self.assertEqual(len(a.requests), 6)
        self.assertEqual(a.requests, b.requests)

    def test_a_sample_without_a_het_table_is_an_error(self):
        with self.assertRaises(KeyError):
            beacon_check.check_individuals(
                _client(FakeBeacon(_held())),
                "YRI_M3_S1",
                {"m1": []},
                members=["m1", "m9"],
                controls=[],
                per_individual=5,
                seed=1,
            )

    def test_a_sample_with_no_het_sites_contributes_nothing_not_a_rate(self):
        result = beacon_check.check_individuals(
            _client(FakeBeacon(_held())),
            "YRI_M3_S1",
            {"m1": [PRESENT[0]], "m2": [], "c1": [ABSENT[0]]},
            members=["m1", "m2"],
            controls=["c1"],
            per_individual=5,
            seed=1,
        )
        self.assertEqual(result["member_queries"], 1)


class TestFilesAndCli(unittest.TestCase):
    def build_dir(self, tmp):
        os.makedirs(os.path.join(tmp, "het"))

        def write(name, rows):
            with open(os.path.join(tmp, name), "w") as fh:
                fh.write("".join("\t".join(map(str, r)) + "\n" for r in rows))

        write("present_sites.tsv", PRESENT)
        write("absent_sites.tsv", ABSENT)
        write("het/m1.tsv", [PRESENT[0]])
        write("het/c1.tsv", [ABSENT[0]])
        with open(os.path.join(tmp, "members.txt"), "w") as fh:
            fh.write("m1\n")
        with open(os.path.join(tmp, "controls.txt"), "w") as fh:
            fh.write("c1\n")
        with open(os.path.join(tmp, "build_summary.json"), "w") as fh:
            json.dump({"dataset_id": "YRI_M3_S1"}, fh)

    def test_reads_a_sites_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_dir(tmp)
            sites = beacon_check.read_sites(os.path.join(tmp, "present_sites.tsv"))
        self.assertEqual(sites, PRESENT)

    def test_a_malformed_sites_line_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.tsv")
            with open(path, "w") as fh:
                fh.write("chr1\tnot-a-number\tA\tG\n")
            with self.assertRaises(ValueError):
                beacon_check.read_sites(path)

    def run_cli(self, fake, tmp):
        out = os.path.join(tmp, "report.json")
        with redirect_stdout(io.StringIO()):
            code = beacon_check.main(
                [
                    "--base-url",
                    "http://beacon.test",
                    "--build-dir",
                    tmp,
                    "--per-individual",
                    "5",
                    "--seed",
                    "1",
                    "--out",
                    out,
                ],
                fetch=fake,
                sleep=lambda s: None,
            )
        with open(out) as fh:
            return code, json.load(fh)

    def test_cli_passes_a_correct_beacon_and_writes_a_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_dir(tmp)
            code, report = self.run_cli(FakeBeacon(_held()), tmp)
        self.assertEqual(code, 0)
        self.assertTrue(report["passed"])
        self.assertEqual(report["dataset_id"], "YRI_M3_S1")
        self.assertIn("queries_sent", report)

    def test_cli_fails_the_beacon_that_says_yes_everywhere(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_dir(tmp)
            code, report = self.run_cli(FakeBeacon(_held(), yes_everywhere=True), tmp)
        self.assertEqual(code, 1)
        self.assertFalse(report["passed"])

    def test_cli_reports_a_failed_control_instead_of_a_traceback(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_dir(tmp)
            code, report = self.run_cli(FakeBeacon([]), tmp)
        self.assertEqual(code, 1)
        self.assertIn("positive control", report["error"])


if __name__ == "__main__":
    unittest.main()
