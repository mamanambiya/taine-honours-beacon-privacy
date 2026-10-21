"""Tests for the seeded member/control split of a population's samples."""

import json
import os
import tempfile
import unittest

import make_split


def _ids(n, prefix="HG"):
    return [f"{prefix}{i:05d}" for i in range(n)]


class TestSplitSamples(unittest.TestCase):
    def test_exact_sizes(self):
        members, controls = make_split.split_samples(_ids(108), 40, 40, "YRI", 1)
        self.assertEqual((len(members), len(controls)), (40, 40))

    def test_members_and_controls_never_overlap(self):
        # A control that is also a member would answer "yes" for the wrong
        # reason, and every power figure built on it would be inflated.
        for seed in range(20):
            members, controls = make_split.split_samples(_ids(85), 40, 40, "MSL", seed)
            self.assertEqual(set(members) & set(controls), set())

    def test_every_chosen_sample_comes_from_the_population(self):
        ids = _ids(99)
        members, controls = make_split.split_samples(ids, 40, 40, "CEU", 1)
        self.assertTrue(set(members + controls) <= set(ids))

    def test_same_inputs_give_the_same_split(self):
        first = make_split.split_samples(_ids(108), 40, 40, "YRI", 1)
        second = make_split.split_samples(_ids(108), 40, 40, "YRI", 1)
        self.assertEqual(first, second)

    def test_input_order_does_not_change_the_split(self):
        # bcftools lists samples in file order, which can differ between two
        # cuts of the same panel. The split must not depend on it.
        ids = _ids(108)
        forward = make_split.split_samples(ids, 40, 40, "YRI", 1)
        backward = make_split.split_samples(list(reversed(ids)), 40, 40, "YRI", 1)
        self.assertEqual(forward, backward)

    def test_a_different_seed_gives_a_different_split(self):
        a = make_split.split_samples(_ids(108), 40, 40, "YRI", 1)
        b = make_split.split_samples(_ids(108), 40, 40, "YRI", 2)
        self.assertNotEqual(a[0], b[0])

    def test_the_same_seed_gives_different_splits_in_different_populations(self):
        # The seed is shared across the ten populations; mixing in the name is
        # what stops YRI and ESN from being split by the same permutation.
        ids = _ids(108)
        a = make_split.split_samples(ids, 40, 40, "YRI", 1)
        b = make_split.split_samples(ids, 40, 40, "ESN", 1)
        self.assertNotEqual(a[0], b[0])

    def test_outputs_are_sorted(self):
        members, controls = make_split.split_samples(_ids(108), 40, 40, "YRI", 1)
        self.assertEqual(members, sorted(members))
        self.assertEqual(controls, sorted(controls))

    def test_the_smallest_population_has_exactly_enough(self):
        members, controls = make_split.split_samples(_ids(80), 40, 40, "MSL", 1)
        self.assertEqual(set(members) | set(controls), set(_ids(80)))


class TestRefusals(unittest.TestCase):
    def test_too_few_samples(self):
        # 79 for an 80-sample design: silently shrinking the controls would
        # make this population incomparable with the others.
        with self.assertRaises(ValueError):
            make_split.split_samples(_ids(79), 40, 40, "MSL", 1)

    def test_empty_list(self):
        with self.assertRaises(ValueError):
            make_split.split_samples([], 40, 40, "YRI", 1)

    def test_duplicate_ids(self):
        ids = _ids(100) + [_ids(100)[0]]
        with self.assertRaises(ValueError):
            make_split.split_samples(ids, 40, 40, "YRI", 1)

    def test_non_positive_sizes(self):
        for members, controls in [(0, 40), (40, 0), (-1, 40)]:
            with self.assertRaises(ValueError):
                make_split.split_samples(_ids(100), members, controls, "YRI", 1)


class TestFilesAndProvenance(unittest.TestCase):
    def write_samples(self, tmp, ids, extra=""):
        path = os.path.join(tmp, "YRI_samples.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(ids) + "\n" + extra)
        return path

    def test_read_ignores_blank_lines_and_whitespace(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_samples(tmp, ["A", " B "], extra="\n\nC\n")
            self.assertEqual(make_split.read_samples(path), ["A", "B", "C"])

    def test_cli_writes_both_lists_and_a_provenance_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = self.write_samples(tmp, _ids(108))
            out = os.path.join(tmp, "out")
            make_split.main(
                [
                    "--population",
                    "YRI",
                    "--samples",
                    src,
                    "--members",
                    "40",
                    "--controls",
                    "40",
                    "--seed",
                    "1",
                    "--out-dir",
                    out,
                ]
            )
            members = make_split.read_samples(os.path.join(out, "YRI_members.txt"))
            controls = make_split.read_samples(os.path.join(out, "YRI_controls.txt"))
            with open(os.path.join(out, "YRI_split.json")) as fh:
                record = json.load(fh)
        self.assertEqual((len(members), len(controls)), (40, 40))
        self.assertEqual(record["population"], "YRI")
        self.assertEqual(record["seed"], 1)
        self.assertEqual(record["n_input_samples"], 108)
        self.assertEqual(len(record["input_sha256"]), 64)

    def test_the_input_hash_is_independent_of_listing_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self.write_samples(tmp, _ids(100))
            out_a = os.path.join(tmp, "a")
            make_split.main(
                [
                    "--population",
                    "YRI",
                    "--samples",
                    a,
                    "--members",
                    "40",
                    "--controls",
                    "40",
                    "--seed",
                    "1",
                    "--out-dir",
                    out_a,
                ]
            )
            b = self.write_samples(tmp, list(reversed(_ids(100))))
            out_b = os.path.join(tmp, "b")
            make_split.main(
                [
                    "--population",
                    "YRI",
                    "--samples",
                    b,
                    "--members",
                    "40",
                    "--controls",
                    "40",
                    "--seed",
                    "1",
                    "--out-dir",
                    out_b,
                ]
            )
            with open(os.path.join(out_a, "YRI_split.json")) as fh:
                ra = json.load(fh)
            with open(os.path.join(out_b, "YRI_split.json")) as fh:
                rb = json.load(fh)
        self.assertEqual(ra["input_sha256"], rb["input_sha256"])

    def test_cli_refuses_to_overwrite_an_existing_split(self):
        # A silent re-split after a dataset was built would orphan the control
        # list that dataset was verified against.
        with tempfile.TemporaryDirectory() as tmp:
            src = self.write_samples(tmp, _ids(100))
            out = os.path.join(tmp, "out")
            args = [
                "--population",
                "YRI",
                "--samples",
                src,
                "--members",
                "40",
                "--controls",
                "40",
                "--seed",
                "1",
                "--out-dir",
                out,
            ]
            make_split.main(args)
            with self.assertRaises(SystemExit):
                make_split.main(args)


if __name__ == "__main__":
    unittest.main()
