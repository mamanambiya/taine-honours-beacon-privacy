"""Tests for D_N taken straight from a panel's allele-count histogram."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from math import comb

import numpy as np

import empirical_dn
import lrt
import sfs_fit


def _histogram_from_frequencies(f, panel_individuals, seed):
    """Allele-count histogram of the sites a panel of N people would list."""
    rng = np.random.default_rng(seed)
    k = rng.binomial(2 * panel_individuals, f)
    return np.bincount(k[k >= 1], minlength=2 * panel_individuals + 1)


def _true_het_d_n(f, n_beacon):
    """E[f(1-f)(1-f)^(2n)] / E[f(1-f)] over known frequencies: the target."""
    weight = f * (1.0 - f)
    return float((weight * (1.0 - f) ** (2 * n_beacon)).sum() / weight.sum())


class TestHandComputedCase(unittest.TestCase):
    def test_matches_exact_fractions_on_a_six_haplotype_panel(self):
        # Panel of 3 people = M 6 haplotypes. Two sites carry 1 copy, one site
        # carries 3. A beacon of one person draws 2 haplotypes from the 4 that
        # are neither of the target's.
        m = 6
        hist = [0, 2, 0, 1, 0, 0, 0]
        weight = lambda k: Fraction(2 * k * (m - k), m * (m - 1))
        none_carry = lambda k: Fraction(comb(m - 1 - k, 2), comb(m - 2, 2))
        expected = (2 * weight(1) * none_carry(1) + 1 * weight(3) * none_carry(3)) / (
            2 * weight(1) + 1 * weight(3)
        )
        self.assertEqual(none_carry(1), 1)
        self.assertEqual(none_carry(3), Fraction(1, 6))
        self.assertAlmostEqual(
            empirical_dn.empirical_d_n(hist, 3, 1), float(expected), places=12
        )


class TestEdges(unittest.TestCase):
    def test_all_singletons_means_nobody_else_can_carry_it(self):
        # The target is the only person with the allele, so the beacon never
        # holds a copy, whatever its size.
        hist = np.zeros(2 * 50 + 1, dtype=int)
        hist[1] = 1000
        for n_beacon in (1, 10, 40):
            self.assertAlmostEqual(empirical_dn.empirical_d_n(hist, 50, n_beacon), 1.0)

    def test_a_site_missing_in_only_one_haplotype_is_always_carried(self):
        # k = M - 1: after the target's two haplotypes at most one reference
        # copy is left, so a beacon of two or more haplotypes must hold the
        # allele.
        panel = 50
        hist = np.zeros(2 * panel + 1, dtype=int)
        hist[2 * panel - 1] = 500
        self.assertEqual(empirical_dn.empirical_d_n(hist, panel, 1), 0.0)

    def test_sites_no_one_is_heterozygous_at_carry_no_weight(self):
        # Counts 0 and M can never be a heterozygous query site; the rest of the
        # histogram alone decides the answer.
        panel = 30
        base = np.zeros(2 * panel + 1, dtype=int)
        base[5] = 200
        base[20] = 100
        padded = base.copy()
        padded[0] = 10**6
        padded[2 * panel] = 10**6
        self.assertAlmostEqual(
            empirical_dn.empirical_d_n(padded, panel, 4),
            empirical_dn.empirical_d_n(base, panel, 4),
            places=12,
        )

    def test_a_larger_beacon_is_less_likely_to_miss_the_allele(self):
        rng = np.random.default_rng(1)
        hist = _histogram_from_frequencies(rng.beta(0.5, 2.0, 50_000), 500, 2)
        values = [empirical_dn.empirical_d_n(hist, 500, n) for n in (5, 20, 65, 200)]
        self.assertEqual(values, sorted(values, reverse=True))
        self.assertTrue(all(0.0 < v < 1.0 for v in values))

    def test_a_bigger_beacon_shrinks_the_separation_the_attack_relies_on(self):
        rng = np.random.default_rng(3)
        hist = _histogram_from_frequencies(rng.beta(0.5, 2.0, 50_000), 500, 4)
        d_n = empirical_dn.empirical_d_n(hist, 500, 20)
        d_n_minus_1 = empirical_dn.empirical_d_n(hist, 500, 19)
        self.assertLess(d_n, d_n_minus_1)


class TestAgreesWithTheParametricModelWhenItIsRight(unittest.TestCase):
    def test_recovers_the_het_conditioned_beta_d_n(self):
        a, b, panel = 0.5, 2.0, 1000
        f = np.random.default_rng(5).beta(a, b, 400_000)
        hist = _histogram_from_frequencies(f, panel, 6)
        for n_beacon in (20, 65, 200):
            parametric = lrt.d_n(a, b, n_beacon, heterozygous=True)
            empirical = empirical_dn.empirical_d_n(hist, panel, n_beacon)
            self.assertAlmostEqual(
                empirical, parametric, delta=0.01 * parametric + 0.003
            )

    def test_unaffected_by_how_many_empty_sites_the_panel_would_have_had(self):
        # A zero-truncated panel cannot know the K = 0 mass, and the answer must
        # not need it: het-conditioning already excludes those sites.
        a, b, panel = 0.5, 2.0, 400
        f = np.random.default_rng(7).beta(a, b, 100_000)
        hist = _histogram_from_frequencies(f, panel, 8)
        with_empties = hist.copy()
        with_empties[0] = 10**9
        self.assertEqual(
            empirical_dn.empirical_d_n(hist, panel, 20),
            empirical_dn.empirical_d_n(with_empties, panel, 20),
        )


class TestWhereASingleBetaFails(unittest.TestCase):
    def test_follows_a_two_population_spectrum_a_beta_cannot_describe(self):
        # Mostly very rare variants plus a block near fixation (reference-rare
        # sites). The V6HC-S_AFR fit hit its bound on exactly this shape.
        rng = np.random.default_rng(11)
        rare = rng.beta(0.05, 5.0, 300_000)
        common = rng.beta(8.0, 0.4, 60_000)
        f = np.concatenate([rare, common])
        panel = 800
        hist = _histogram_from_frequencies(f, panel, 12)
        n_beacon = 65
        truth = _true_het_d_n(f, n_beacon)

        empirical = empirical_dn.empirical_d_n(hist, panel, n_beacon)
        self.assertAlmostEqual(empirical, truth, delta=0.01)

        fit = sfs_fit.fit_beta_binomial(hist, panel)
        beta = lrt.d_n(fit["a"], fit["b"], n_beacon, heterozygous=True)
        self.assertGreater(abs(beta - truth), 3 * abs(empirical - truth))


class TestQueryBudget(unittest.TestCase):
    def test_budget_uses_d_n_and_d_n_minus_1_of_the_histogram(self):
        rng = np.random.default_rng(13)
        hist = _histogram_from_frequencies(rng.beta(0.5, 2.0, 80_000), 600, 14)
        d_n = empirical_dn.empirical_d_n(hist, 600, 65)
        d_n_minus_1 = empirical_dn.empirical_d_n(hist, 600, 64)
        expected = lrt.required_queries(d_n, d_n_minus_1, 1e-6, 0.05, 0.95)
        self.assertEqual(
            empirical_dn.required_queries_empirical(hist, 600, 65, 1e-6, 0.05, 0.95),
            expected,
        )


class TestRefusals(unittest.TestCase):
    def setUp(self):
        self.hist = np.zeros(2 * 20 + 1, dtype=int)
        self.hist[1:10] = 5

    def test_a_beacon_as_large_as_the_panel_is_refused(self):
        # Needs the target's two haplotypes plus the beacon's 2n from the panel.
        with self.assertRaisesRegex(ValueError, "beacon"):
            empirical_dn.empirical_d_n(self.hist, 20, 20)

    def test_the_largest_workable_beacon_is_one_below_the_panel(self):
        self.assertTrue(0.0 <= empirical_dn.empirical_d_n(self.hist, 20, 19) <= 1.0)

    def test_a_beacon_of_zero_is_refused(self):
        with self.assertRaisesRegex(ValueError, "positive"):
            empirical_dn.empirical_d_n(self.hist, 20, 0)

    def test_the_wrong_histogram_length_is_refused(self):
        with self.assertRaisesRegex(ValueError, "bins"):
            empirical_dn.empirical_d_n(self.hist[:-1], 20, 5)

    def test_negative_entries_are_refused(self):
        bad = self.hist.copy()
        bad[3] = -1
        with self.assertRaisesRegex(ValueError, "negative"):
            empirical_dn.empirical_d_n(bad, 20, 5)

    def test_a_histogram_with_no_heterozygous_sites_is_refused(self):
        only_ends = np.zeros(2 * 20 + 1, dtype=int)
        only_ends[0] = 100
        only_ends[40] = 100
        with self.assertRaisesRegex(ValueError, "heterozygous"):
            empirical_dn.empirical_d_n(only_ends, 20, 5)

    def test_a_panel_too_small_to_hold_a_target_and_a_beacon_is_refused(self):
        with self.assertRaisesRegex(ValueError, "panel"):
            empirical_dn.empirical_d_n(np.array([0, 1, 1]), 1, 1)


class TestCommandLine(unittest.TestCase):
    def test_sfs_fit_reports_the_empirical_figures_at_the_asked_beacon_sizes(self):
        rng = np.random.default_rng(17)
        panel = 300
        hist = _histogram_from_frequencies(rng.beta(0.5, 2.0, 60_000), panel, 18)
        tmp = tempfile.mkdtemp()
        try:
            table = os.path.join(tmp, "h.tsv")
            with open(table, "w") as fh:
                for k, c in enumerate(hist):
                    if c:
                        fh.write(f"{k}\t{c}\n")
            out = os.path.join(tmp, "fit.json")
            here = os.path.dirname(os.path.abspath(__file__))
            subprocess.run(
                [
                    sys.executable,
                    os.path.join(here, "sfs_fit.py"),
                    table,
                    "--n-individuals",
                    str(panel),
                    "--beacon-sizes",
                    "20,65",
                    "--out",
                    out,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with open(out) as fh:
                record = json.load(fh)
        finally:
            shutil.rmtree(tmp)

        rows = record["empirical_d_n"]
        self.assertEqual([r["n_beacon"] for r in rows], [20, 65])
        for row in rows:
            self.assertLess(row["d_n"], row["d_n_minus_1"])
            self.assertGreater(row["queries"], 0)
            self.assertIn("beta_fit_queries", row)
        self.assertEqual(rows[0]["d_n"], empirical_dn.empirical_d_n(hist, panel, 20))

    def _run(self, tmp, panel, extra):
        rng = np.random.default_rng(19)
        hist = _histogram_from_frequencies(rng.beta(0.5, 2.0, 20_000), panel, 20)
        table = os.path.join(tmp, "h.tsv")
        with open(table, "w") as fh:
            for k, c in enumerate(hist):
                if c:
                    fh.write(f"{k}\t{c}\n")
        out = os.path.join(tmp, "fit.json")
        here = os.path.dirname(os.path.abspath(__file__))
        return subprocess.run(
            [sys.executable, os.path.join(here, "sfs_fit.py"), table,
             "--n-individuals", str(panel), "--out", out, *extra],
            capture_output=True, text=True,
        ), out

    def test_a_panel_too_small_for_the_default_sizes_still_fits(self):
        tmp = tempfile.mkdtemp()
        try:
            done, out = self._run(tmp, 30, [])
            self.assertEqual(done.returncode, 0, done.stderr)
            with open(out) as fh:
                rows = json.load(fh)["empirical_d_n"]
        finally:
            shutil.rmtree(tmp)
        # 20 fits in a panel of 30; 65 does not and is dropped, not an error.
        self.assertEqual([r["n_beacon"] for r in rows], [20])

    def test_a_size_asked_for_that_does_not_fit_is_an_error(self):
        tmp = tempfile.mkdtemp()
        try:
            done, _ = self._run(tmp, 30, ["--beacon-sizes", "40"])
        finally:
            shutil.rmtree(tmp)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("does not fit", done.stderr)


if __name__ == "__main__":
    unittest.main()
