"""Tests for fitting the Beta site-frequency spectrum to a panel's allele counts."""

import contextlib
import io
import json
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

import cohort
import lrt
import sfs_fit


def _panel_counts(a, b, n_individuals, n_variants, seed, min_count=1):
    """Alternate-allele counts of the variants a real panel would list.

    A panel VCF holds only polymorphic sites, so draws with fewer than
    ``min_count`` copies are never seen. That truncation is the whole reason
    ``fit_beta_binomial`` is not a plain beta-binomial fit.
    """
    rng = np.random.default_rng(seed)
    f = rng.beta(a, b, size=n_variants)
    k = rng.binomial(2 * n_individuals, f)
    return k[k >= min_count]


class TestBetaBinomialPmf(unittest.TestCase):
    def test_uniform_spectrum_gives_a_flat_count_distribution(self):
        # Beta(1, 1) over f makes every count 0..n equally likely: 1 / (n + 1).
        n_individuals = 5
        k = np.arange(0, 2 * n_individuals + 1)
        pmf = np.exp(sfs_fit.log_pmf(k, n_individuals, 1.0, 1.0))
        np.testing.assert_allclose(pmf, 1.0 / (2 * n_individuals + 1))

    def test_sums_to_one(self):
        n_individuals = 40
        k = np.arange(0, 2 * n_individuals + 1)
        for a, b in [(0.5, 2.0), (0.2, 20.0), (3.0, 1.5)]:
            self.assertAlmostEqual(
                float(np.exp(sfs_fit.log_pmf(k, n_individuals, a, b)).sum()),
                1.0,
                places=10,
            )

    def test_zero_count_is_exactly_d_n(self):
        # P(K = 0) is the chance nobody carries the allele, which is D_N. If the
        # fit and the attack model disagree here, every fitted parameter is
        # being fed into the wrong quantity.
        for a, b, n in [(0.5, 2.0, 1895), (0.2, 20.0, 40)]:
            self.assertAlmostEqual(
                float(np.exp(sfs_fit.log_pmf(np.array([0]), n, a, b))[0]),
                lrt.d_n(a, b, n),
                places=12,
            )

    def test_does_not_underflow_at_the_live_panel_size(self):
        k = np.arange(0, 2 * 1895 + 1)
        log_pmf = sfs_fit.log_pmf(k, 1895, 0.5, 2.0)
        self.assertTrue(np.isfinite(log_pmf).all())


class TestCountHistogram(unittest.TestCase):
    def test_bins_counts_by_number_of_alternate_copies(self):
        hist = sfs_fit.count_histogram([1, 1, 2, 4, 4, 4], n_individuals=2)
        np.testing.assert_array_equal(hist, [0, 2, 1, 0, 3])

    def test_length_is_two_n_plus_one(self):
        self.assertEqual(len(sfs_fit.count_histogram([1], n_individuals=10)), 21)

    def test_rejects_a_count_above_two_n(self):
        # More alternate copies than haplotypes means AN was not 2N: missing
        # genotypes, or the wrong N. A silent clip would hide that.
        with self.assertRaises(ValueError):
            sfs_fit.count_histogram([5], n_individuals=2)

    def test_rejects_a_negative_count(self):
        with self.assertRaises(ValueError):
            sfs_fit.count_histogram([-1], n_individuals=2)


class TestFitRecoversTheTruth(unittest.TestCase):
    def assertRecovers(
        self, a, b, n_individuals, n_variants, seed, min_count=1, tolerance=0.10
    ):
        counts = _panel_counts(
            a, b, n_individuals, n_variants, seed, min_count=min_count
        )
        hist = sfs_fit.count_histogram(counts, n_individuals)
        fit = sfs_fit.fit_beta_binomial(hist, n_individuals, min_count=min_count)
        self.assertTrue(fit["converged"], fit)
        self.assertAlmostEqual(fit["a"] / a, 1.0, delta=tolerance, msg=fit)
        self.assertAlmostEqual(fit["b"] / b, 1.0, delta=tolerance, msg=fit)
        return fit

    def test_moderately_rare_spectrum(self):
        self.assertRecovers(0.5, 2.0, n_individuals=500, n_variants=400_000, seed=11)

    def test_strongly_rare_spectrum(self):
        self.assertRecovers(0.2, 20.0, n_individuals=500, n_variants=400_000, seed=12)

    def test_at_the_live_panel_size(self):
        self.assertRecovers(0.5, 2.0, n_individuals=1895, n_variants=400_000, seed=13)

    def test_a_panel_filtered_to_minor_allele_count_two_or_more(self):
        # Many reference panels drop singletons. Fitting such a panel as though
        # singletons had been seen and missed would be badly biased; the
        # truncation point has to be a parameter.
        self.assertRecovers(
            0.5, 2.0, n_individuals=500, n_variants=400_000, seed=14, min_count=2
        )

    def test_reports_standard_errors_that_cover_the_truth(self):
        a, b = 0.5, 2.0
        fit = self.assertRecovers(a, b, n_individuals=500, n_variants=400_000, seed=15)
        self.assertGreater(fit["se_a"], 0.0)
        self.assertGreater(fit["se_b"], 0.0)
        self.assertLess(abs(fit["a"] - a), 4 * fit["se_a"])
        self.assertLess(abs(fit["b"] - b), 4 * fit["se_b"])

    def test_agrees_with_a_genotype_matrix_built_by_cohort(self):
        # The fit sees only allele counts; the attack sees genotypes. They must
        # be the same spectrum, or the fitted parameters do not describe the
        # cohorts the power curves were run on.
        a, b, n = 0.5, 2.0, 60
        genotypes = cohort.simulate_cohort(n, 60_000, a, b, seed=16)
        counts = cohort.alternate_allele_counts(genotypes)
        counts = counts[counts >= 1]
        fit = sfs_fit.fit_beta_binomial(sfs_fit.count_histogram(counts, n), n)
        self.assertAlmostEqual(fit["a"] / a, 1.0, delta=0.15, msg=fit)
        self.assertAlmostEqual(fit["b"] / b, 1.0, delta=0.15, msg=fit)


class TestTruncationMatters(unittest.TestCase):
    def test_ignoring_the_missing_zero_class_biases_the_fit(self):
        # Known-negative control: the same data fitted as if no sites had been
        # dropped (min_count=0) must come out wrong. If it did not, the
        # recovery tests above would not be telling truncation apart from
        # luck.
        a, b, n = 0.5, 2.0, 500
        counts = _panel_counts(a, b, n, 400_000, seed=11)
        hist = sfs_fit.count_histogram(counts, n)
        right = sfs_fit.fit_beta_binomial(hist, n, min_count=1)
        wrong = sfs_fit.fit_beta_binomial(hist, n, min_count=0)
        self.assertLess(abs(right["a"] / a - 1.0), 0.10)
        self.assertGreater(abs(wrong["a"] / a - 1.0), 0.15, wrong)


class TestRefusesWhatItCannotFit(unittest.TestCase):
    def test_empty_histogram(self):
        with self.assertRaises(ValueError):
            sfs_fit.fit_beta_binomial(np.zeros(21, dtype=int), 10)

    def test_histogram_of_the_wrong_length(self):
        with self.assertRaises(ValueError):
            sfs_fit.fit_beta_binomial(np.ones(15, dtype=int), 10)

    def test_counts_below_the_truncation_point(self):
        # A singleton in a histogram declared to start at 2 means the panel
        # was not filtered the way the caller says.
        hist = np.zeros(21, dtype=int)
        hist[1] = 5
        hist[3] = 5
        with self.assertRaises(ValueError):
            sfs_fit.fit_beta_binomial(hist, 10, min_count=2)

    def test_negative_histogram_entries(self):
        hist = np.ones(21, dtype=int)
        hist[4] = -1
        with self.assertRaises(ValueError):
            sfs_fit.fit_beta_binomial(hist, 10)

    def test_a_single_occupied_bin_is_flagged_not_trusted(self):
        # Every variant a singleton: the likelihood is flat along a ridge, so
        # any (a, b) it returns is the optimiser's starting point, not data.
        hist = np.zeros(21, dtype=int)
        hist[1] = 1000
        with self.assertRaises(ValueError):
            sfs_fit.fit_beta_binomial(hist, 10)

    def test_a_fit_pinned_to_the_search_bound_is_reported_as_such(self):
        # All-common variants push b toward zero and a toward infinity; the
        # result is a bound, not an estimate, and says so.
        n = 50
        rng = np.random.default_rng(3)
        k = rng.binomial(2 * n, rng.beta(3000.0, 30.0, size=50_000))
        k = k[k >= 1]
        fit = sfs_fit.fit_beta_binomial(sfs_fit.count_histogram(k, n), n)
        self.assertTrue(fit["at_boundary"], fit)


class TestDiagnostics(unittest.TestCase):
    def test_a_well_specified_panel_matches_its_fit(self):
        a, b, n = 0.5, 2.0, 500
        hist = sfs_fit.count_histogram(_panel_counts(a, b, n, 400_000, 21), n)
        fit = sfs_fit.fit_beta_binomial(hist, n)
        diag = sfs_fit.diagnostics(hist, n, fit["a"], fit["b"])
        self.assertAlmostEqual(
            diag["singleton_fraction_observed"],
            diag["singleton_fraction_fitted"],
            delta=0.01,
        )
        self.assertAlmostEqual(
            diag["rare_fraction_observed"], diag["rare_fraction_fitted"], delta=0.01
        )

    def test_flags_a_panel_that_a_single_beta_cannot_describe(self):
        # Real panels carry a second population of variants where the
        # REFERENCE allele is the rare one, so the alternate allele is near
        # fixation. One Beta cannot hold both humps. The fit still returns
        # numbers; the diagnostics are what says they should not be trusted.
        n = 500
        rng = np.random.default_rng(22)
        f = np.concatenate(
            [rng.beta(0.5, 2.0, size=300_000), rng.beta(30.0, 1.0, size=100_000)]
        )
        k = rng.binomial(2 * n, f)
        hist = sfs_fit.count_histogram(k[k >= 1], n)
        fit = sfs_fit.fit_beta_binomial(hist, n)
        diag = sfs_fit.diagnostics(hist, n, fit["a"], fit["b"])
        self.assertGreater(
            abs(
                diag["high_frequency_fraction_observed"]
                - diag["high_frequency_fraction_fitted"]
            ),
            0.02,
            diag,
        )

    def test_observed_fractions_are_counted_on_a_histogram_read_by_hand(self):
        # 200 haplotypes. Alternate allele frequency, by count: 1 -> 0.5%,
        # 100 -> 50%, 150 -> 75%, 199 -> 99.5%. Minor allele frequency is the
        # smaller of the two alleles, so 1 and 199 are both rare.
        n = 100
        hist = np.zeros(2 * n + 1, dtype=int)
        hist[1], hist[100], hist[150], hist[199] = 4, 3, 2, 1
        diag = sfs_fit.diagnostics(hist, n, 0.5, 2.0)
        self.assertAlmostEqual(diag["singleton_fraction_observed"], 0.4)
        self.assertAlmostEqual(diag["rare_fraction_observed"], 0.5)
        # Strictly above 50%: the count-100 sites sit exactly on it.
        self.assertAlmostEqual(diag["high_frequency_fraction_observed"], 0.3)

    def test_reports_a_binned_comparison_for_plotting(self):
        a, b, n = 0.5, 2.0, 200
        hist = sfs_fit.count_histogram(_panel_counts(a, b, n, 100_000, 23), n)
        diag = sfs_fit.diagnostics(hist, n, a, b)
        self.assertGreater(len(diag["bins"]), 3)
        for row in diag["bins"]:
            self.assertEqual(set(row), {"k_low", "k_high", "observed", "expected"})
        self.assertAlmostEqual(
            sum(r["observed"] for r in diag["bins"]), int(hist.sum()), places=6
        )
        self.assertAlmostEqual(
            sum(r["expected"] for r in diag["bins"]),
            int(hist.sum()),
            delta=1e-6 * hist.sum(),
        )


class TestWhatTheFitMeansForTheAttack(unittest.TestCase):
    def test_reproduces_the_documented_fixture_budget(self):
        # beacon-constraints.md: Beta(0.5, 2.0), N = 1895, delta = 1e-6 needs
        # 108,898 queries for 95% power at alpha 0.05. The fit's route to that
        # number must land on the same figure, or the "how far did the fixture
        # move" comparison is comparing two different calculations.
        self.assertEqual(
            sfs_fit.required_queries_for_fit(0.5, 2.0, 1895, delta=1e-6), 108_898
        )

    def test_budget_moves_with_the_fitted_spectrum(self):
        base = sfs_fit.required_queries_for_fit(0.5, 2.0, 1895, delta=1e-6)
        rarer = sfs_fit.required_queries_for_fit(0.2, 20.0, 1895, delta=1e-6)
        self.assertLess(rarer, base)

    def test_compare_to_fixture_reports_the_ratio(self):
        out = sfs_fit.compare_to_fixture(0.5, 2.0, 1895, delta=1e-6)
        self.assertEqual(out["fitted_queries"], out["fixture_queries"])
        self.assertAlmostEqual(out["ratio_fitted_to_fixture"], 1.0)


class TestHistogramFile(unittest.TestCase):
    def test_reads_a_two_column_table_and_ignores_comments(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hist.tsv")
            with open(path, "w") as fh:
                fh.write("# AC\tvariants\n1\t7\n2\t3\n\n4\t1\n")
            hist = sfs_fit.read_histogram(path, n_individuals=2)
        np.testing.assert_array_equal(hist, [0, 7, 3, 0, 1])

    def test_sums_repeated_counts_so_chromosomes_can_be_concatenated(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hist.tsv")
            with open(path, "w") as fh:
                fh.write("1\t5\n1\t6\n2\t1\n")
            hist = sfs_fit.read_histogram(path, n_individuals=2)
        np.testing.assert_array_equal(hist, [0, 11, 1, 0, 0])

    def test_rejects_a_malformed_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hist.tsv")
            with open(path, "w") as fh:
                fh.write("1\t5\nnot a number\n")
            with self.assertRaises(ValueError):
                sfs_fit.read_histogram(path, n_individuals=2)

    def test_command_line_writes_a_json_record_with_provenance(self):
        a, b, n = 0.5, 2.0, 100
        hist = sfs_fit.count_histogram(_panel_counts(a, b, n, 100_000, 31), n)
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "hist.tsv")
            out = os.path.join(tmp, "fit.json")
            with open(src, "w") as fh:
                for k, c in enumerate(hist):
                    if c:
                        fh.write(f"{k}\t{c}\n")
            with contextlib.redirect_stdout(io.StringIO()):
                sfs_fit.main(
                    [
                        src,
                        "--n-individuals",
                        str(n),
                        "--label",
                        "toy",
                        "--delta",
                        "1e-6",
                        "--out",
                        out,
                    ]
                )
            with open(out) as fh:
                record = json.load(fh)
        self.assertEqual(record["label"], "toy")
        self.assertEqual(record["n_individuals"], n)
        self.assertEqual(record["source"], src)
        self.assertAlmostEqual(record["fit"]["a"] / a, 1.0, delta=0.15)
        self.assertIn("diagnostics", record)
        self.assertIn("fixture_comparison", record)


_TOY_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5
chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\t0|0
chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1\t0|0\t0|0\t0|0
chr1\t300\t.\tA\tG,T\t.\tPASS\t.\tGT\t1|2\t1|0\t0|0\t0|0\t0|0
chr1\t400\t.\tAC\tA\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\t0|0
chr1\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0|.\t0|1\t0|0\t0|0\t0|0
chr1\t600\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
"""


@unittest.skipUnless(shutil.which("bcftools"), "bcftools not installed")
class TestExtractionScript(unittest.TestCase):
    """The step that runs on the cluster, checked on a VCF small enough to read."""

    script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "sfs_histogram.sh"
    )

    def run_script(self, n_individuals=5, **env):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = os.path.join(tmp, "toy.vcf")
            out = os.path.join(tmp, "hist.tsv")
            with open(vcf, "w") as fh:
                fh.write(_TOY_VCF)
            result = subprocess.run(
                [self.script, vcf, str(n_individuals), out],
                capture_output=True,
                text=True,
                env={**os.environ, **env},
            )
            table = ""
            if os.path.exists(out):
                with open(out) as fh:
                    table = fh.read()
        return result, table

    def test_counts_every_alt_allele_of_every_snp_site(self):
        result, table = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        # site 100: AC 1.  site 200: AC 3.  site 300 (multi-allelic): G has 2,
        # T has 1.  The indel is excluded, the site with a missing genotype
        # (AN = 9) is skipped, and the all-reference site has AC = 0.
        self.assertEqual(table, "1\t2\n2\t1\n3\t1\n")

    def test_reports_what_it_skipped(self):
        result, _ = self.run_script()
        self.assertIn("kept 4 ALT alleles", result.stderr)
        self.assertIn("skipped 1 at sites with AN != 10", result.stderr)
        self.assertIn("1 with AC = 0", result.stderr)

    def test_the_table_feeds_straight_into_the_fit_reader(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = os.path.join(tmp, "toy.vcf")
            out = os.path.join(tmp, "hist.tsv")
            with open(vcf, "w") as fh:
                fh.write(_TOY_VCF)
            subprocess.run(
                [self.script, vcf, "5", out], check=True, capture_output=True
            )
            hist = sfs_fit.read_histogram(out, n_individuals=5)
        np.testing.assert_array_equal(hist[:5], [0, 2, 1, 1, 0])

    def test_indels_come_back_when_asked_for(self):
        _, table = self.run_script(VTYPES="any")
        self.assertEqual(table, "1\t3\n2\t1\n3\t1\n")

    def test_wrong_n_skips_every_site_rather_than_guessing(self):
        # With N = 4 no site has AN = 8, so nothing is kept. Silence here would
        # be a fit on an empty table.
        result, table = self.run_script(n_individuals=4)
        self.assertEqual(table, "")
        self.assertIn("kept 0 ALT alleles", result.stderr)

    def test_refuses_a_non_numeric_n(self):
        result, _ = self.run_script(n_individuals="five")
        self.assertEqual(result.returncode, 2)


if __name__ == "__main__":
    unittest.main()
