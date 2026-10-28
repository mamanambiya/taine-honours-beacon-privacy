"""Tests for the attack run against a loaded beacon's real answers.

The beacon is the set of sites its members carry, so a query answers yes exactly
when the queried allele is in that set. The tests pin the three things that can
silently go wrong: the site key (an off-by-one or allele mix-up turns every
answer into no), the all-yes probability (the only quantity the attack's
rejection rule depends on once the mismatch rate is negligible), and the budget
search (the number the study reports).
"""

import math
import os
import tempfile
import unittest

import real_beacon_attack as rba


class TestSiteKey(unittest.TestCase):
    def test_chr_prefix_and_bare_number_give_the_same_key(self):
        self.assertEqual(
            rba.site_key("chr20", 68049, "C", "A"), rba.site_key("20", 68049, "C", "A")
        )

    def test_position_ref_and_alt_all_change_the_key(self):
        base = rba.site_key("20", 100, "A", "G")
        self.assertNotEqual(base, rba.site_key("20", 101, "A", "G"))
        self.assertNotEqual(base, rba.site_key("20", 100, "C", "G"))
        self.assertNotEqual(base, rba.site_key("20", 100, "A", "T"))
        self.assertNotEqual(base, rba.site_key("21", 100, "A", "G"))

    def test_a_ref_alt_swap_is_a_different_site(self):
        self.assertNotEqual(
            rba.site_key("20", 100, "A", "G"), rba.site_key("20", 100, "G", "A")
        )

    def test_a_base_that_is_not_acgt_is_refused_not_mapped(self):
        with self.assertRaises(ValueError):
            rba.site_key("20", 100, "AT", "G")
        with self.assertRaises(ValueError):
            rba.site_key("20", 100, "A", "N")


class TestAllYesProbability(unittest.TestCase):
    def test_one_query_is_the_fraction_of_present_sites(self):
        self.assertAlmostEqual(rba.p_all_yes(10, 2, 1), 0.8)

    def test_two_queries_are_drawn_without_replacement(self):
        self.assertAlmostEqual(rba.p_all_yes(10, 2, 2), (8 / 10) * (7 / 9))

    def test_no_absent_site_means_certain_yes_at_any_budget(self):
        self.assertEqual(rba.p_all_yes(50, 0, 50), 1.0)

    def test_more_queries_than_present_sites_cannot_all_be_yes(self):
        self.assertEqual(rba.p_all_yes(10, 2, 9), 0.0)

    def test_more_queries_than_heterozygous_sites_is_refused(self):
        with self.assertRaises(ValueError):
            rba.p_all_yes(10, 2, 11)

    def test_zero_queries_is_certain(self):
        self.assertEqual(rba.p_all_yes(10, 2, 0), 1.0)

    def test_matches_the_binomial_when_the_pool_is_large(self):
        h, m, n = 150000, 4200, 95
        self.assertAlmostEqual(
            rba.p_all_yes(h, m, n), (1 - m / h) ** n, places=3
        )


class TestFalsePositiveCurve(unittest.TestCase):
    def test_it_is_the_mean_over_outsiders(self):
        outsiders = [(10, 2), (10, 0)]
        self.assertAlmostEqual(rba.fpr(outsiders, 1), (0.8 + 1.0) / 2)

    def test_it_never_increases_with_the_budget(self):
        outsiders = [(100, 5), (100, 20), (100, 1)]
        values = [rba.fpr(outsiders, n) for n in range(0, 40)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_an_outsider_with_no_absent_site_keeps_the_rate_up(self):
        # A person every query of whom answers yes is indistinguishable from a
        # member: the rate cannot fall below their share of the outsiders.
        outsiders = [(100, 0), (100, 50)]
        self.assertGreaterEqual(rba.fpr(outsiders, 100), 0.5)


class TestBudget(unittest.TestCase):
    def test_smallest_budget_at_or_below_alpha(self):
        outsiders = [(1000, 100)] * 5  # p(yes) about 0.9 per query
        n = rba.budget_for_fpr(outsiders, alpha=0.05, max_queries=500)
        self.assertLessEqual(rba.fpr(outsiders, n), 0.05)
        self.assertGreater(rba.fpr(outsiders, n - 1), 0.05)
        self.assertTrue(25 <= n <= 35, n)  # ln(0.05) / ln(0.9) = 28.4

    def test_unreachable_target_returns_none_not_the_cap(self):
        outsiders = [(100, 0)]
        self.assertIsNone(rba.budget_for_fpr(outsiders, alpha=0.05, max_queries=100))

    def test_alpha_outside_the_open_unit_interval_is_refused(self):
        for bad in (0.0, 1.0, -0.1, 1.5):
            with self.assertRaises(ValueError):
                rba.budget_for_fpr([(10, 1)], alpha=bad, max_queries=10)

    def test_an_empty_outsider_set_is_refused(self):
        with self.assertRaises(ValueError):
            rba.budget_for_fpr([], alpha=0.05, max_queries=10)


class TestBootstrap(unittest.TestCase):
    def test_the_interval_brackets_the_point_estimate(self):
        outsiders = [(1000, m) for m in range(40, 160, 4)]
        point = rba.budget_for_fpr(outsiders, 0.05, 500)
        low, high = rba.bootstrap_budget(outsiders, 0.05, 500, draws=200, seed=1)
        self.assertLessEqual(low, point)
        self.assertGreaterEqual(high, point)

    def test_the_same_seed_gives_the_same_interval(self):
        outsiders = [(1000, m) for m in range(40, 160, 4)]
        a = rba.bootstrap_budget(outsiders, 0.05, 500, draws=100, seed=7)
        b = rba.bootstrap_budget(outsiders, 0.05, 500, draws=100, seed=7)
        self.assertEqual(a, b)


class TestCounting(unittest.TestCase):
    def test_a_member_whose_het_sites_are_all_in_the_beacon_has_no_absent_site(self):
        sites = {rba.site_key("chr1", 100, "A", "G"), rba.site_key("chr1", 300, "C", "T")}
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "NA1.tsv")
            with open(path, "w") as f:
                f.write("chr1\t100\tA\tG\nchr1\t300\tC\tT\n")
            self.assertEqual(rba.count_person(path, sites), (2, 0))

    def test_an_absent_site_is_counted_once_per_line(self):
        sites = {rba.site_key("chr1", 100, "A", "G")}
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "NA2.tsv")
            with open(path, "w") as f:
                f.write("chr1\t100\tA\tG\nchr1\t200\tC\tT\nchr1\t201\tC\tT\n")
            self.assertEqual(rba.count_person(path, sites), (3, 2))

    def test_an_empty_table_is_refused_not_scored_as_certain(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "NA3.tsv")
            open(path, "w").close()
            with self.assertRaises(ValueError):
                rba.count_person(path, set())


class TestLoadBeaconSites(unittest.TestCase):
    def _write(self, rows):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "sites.csv")
        with open(path, "w") as f:
            f.write("reference_name,start,reference_bases,alternate_bases\n")
            for r in rows:
                f.write(",".join(str(x) for x in r) + "\n")
        return path

    def test_the_beacons_zero_based_start_becomes_the_one_based_position(self):
        # Beacon start 60180 is VCF POS 60181; a het table row carries POS.
        sites = rba.load_beacon_sites(self._write([("20", 60180, "A", "G")]))
        self.assertIn(rba.site_key("chr20", 60181, "A", "G"), sites)
        self.assertNotIn(rba.site_key("chr20", 60180, "A", "G"), sites)

    def test_columns_are_found_by_name_not_position(self):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "s.csv")
        with open(path, "w") as f:
            f.write("alternate_bases,reference_bases,start,reference_name\nG,A,99,1\n")
        self.assertIn(rba.site_key("1", 100, "A", "G"), rba.load_beacon_sites(path))

    def test_indels_are_left_out_rather_than_mis_keyed(self):
        sites = rba.load_beacon_sites(
            self._write([("1", 10, "A", "AT"), ("1", 20, "C", "T")])
        )
        self.assertEqual(len(sites), 1)


if __name__ == "__main__":
    unittest.main()
