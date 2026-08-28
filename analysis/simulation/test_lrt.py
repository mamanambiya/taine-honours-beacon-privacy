"""Tests for the Shringarpure & Bustamante membership-inference LRT.

Run:  python3 -m unittest discover -s analysis/simulation -t analysis/simulation

The implementation predates this suite, so "watch it fail first" is replaced by
MUTATION CHECKS: tests that assert on a NUMBER which changes when the model is
broken. That matters specifically here, because the defect this module already
suffered once — using D_N where D_{N-1} belongs — does not crash and does not
produce an obviously wrong curve. It produces a plausible curve with the wrong
required-n, which is the single number the study reports.

numpy is used only to cross-check a closed form against sampling. The module
under test is stdlib-only.
"""
import math
import unittest

import numpy as np

import lrt


class TestDN(unittest.TestCase):
    """D_N = E_f[(1 - f)^(2N)] for f ~ Beta(a, b)."""

    A, B = 0.5, 2.0

    def test_closed_form_matches_monte_carlo_integration(self):
        # The point of the closed form B(a, b+2N)/B(a, b) is that it replaces
        # an integral. If it disagrees with sampling, it is not that integral.
        n = 50
        rng = np.random.default_rng(20260821)
        f = rng.beta(self.A, self.B, size=2_000_000)
        draws = (1.0 - f) ** (2 * n)
        sampled = float(np.mean(draws))
        # Tolerance from the sample's own standard error, not a decimal place:
        # 2M draws cannot resolve 5e-5, so places=4 would flake on the seed.
        standard_error = float(np.std(draws, ddof=1)) / np.sqrt(draws.size)
        self.assertLess(abs(lrt.d_n(self.A, self.B, n) - sampled),
                        5.0 * standard_error)

    def test_decreases_as_the_beacon_grows(self):
        values = [lrt.d_n(self.A, self.B, n) for n in (10, 100, 1000)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_d_n_minus_1_exceeds_d_n(self):
        # Fewer genomes -> better chance nobody carries it. This ordering is
        # what makes the two hypotheses separate at all.
        self.assertGreater(lrt.d_n(self.A, self.B, 99),
                           lrt.d_n(self.A, self.B, 100))


class TestAscertainment(unittest.TestCase):
    """Queries come only from heterozygous sites, so the SFS is conditioned."""

    A, B, N = 0.5, 2.0, 50

    def test_conditioning_shifts_both_beta_parameters_by_one(self):
        # P(het) is proportional to 2f(1-f), so the ascertained density goes as
        # f^(a-1)(1-f)^(b-1) * f(1-f) = f^a (1-f)^b, i.e. Beta(a+1, b+1).
        self.assertAlmostEqual(
            lrt.d_n(self.A, self.B, self.N, heterozygous=True),
            lrt.d_n(self.A + 1.0, self.B + 1.0, self.N),
            places=12)

    def test_matches_monte_carlo_over_the_heterozygote_weighted_sfs(self):
        # Independent check that does NOT reuse the a+1/b+1 identity: sample
        # from Beta(a, b), resample proportional to the heterozygote weight
        # 2f(1-f), and integrate over what survives.
        rng = np.random.default_rng(11071987)
        f = rng.beta(self.A, self.B, size=4_000_000)
        weight = 2.0 * f * (1.0 - f)
        keep = rng.random(f.size) < (weight / weight.max())
        ascertained = f[keep]
        self.assertGreater(ascertained.size, 100_000)
        self.assertAlmostEqual(
            lrt.d_n(self.A, self.B, self.N, heterozygous=True),
            float(np.mean((1.0 - ascertained) ** (2 * self.N))),
            places=3)

    def test_conditioning_lowers_d_n(self):
        # Heterozygous sites are enriched for intermediate frequencies, so
        # "nobody carries it" becomes less likely. Using the unconditioned SFS
        # therefore UNDERSTATES the attack.
        self.assertLess(lrt.d_n(self.A, self.B, self.N, heterozygous=True),
                        lrt.d_n(self.A, self.B, self.N))


class TestHypotheses(unittest.TestCase):

    DELTA = 1e-3

    def setUp(self):
        self.d_n = lrt.d_n(0.5, 2.0, 100)
        self.d_n_minus_1 = lrt.d_n(0.5, 2.0, 99)

    def test_h1_is_driven_by_the_other_n_minus_1_genomes(self):
        # MUTATION SENTINEL. Passing D_N here is the exact defect this module
        # shipped with once. It does not raise; it silently narrows the gap.
        correct = lrt.theta_h1(self.d_n_minus_1, self.DELTA)
        defective = lrt.theta_h1(self.d_n, self.DELTA)
        self.assertNotEqual(correct, defective)
        self.assertAlmostEqual(correct, 1.0 - self.DELTA * self.d_n_minus_1,
                               places=15)

    def test_h1_predicts_more_yes_responses_than_h0(self):
        self.assertGreater(lrt.theta_h1(self.d_n_minus_1, self.DELTA),
                           lrt.theta_h0(self.d_n))


class TestStatistic(unittest.TestCase):

    DELTA = 1e-3

    def setUp(self):
        self.d_n = lrt.d_n(0.5, 2.0, 100)
        self.d_n_minus_1 = lrt.d_n(0.5, 2.0, 99)

    def stat(self, responses):
        return lrt.lrt_statistic(responses, self.d_n, self.d_n_minus_1,
                                 self.DELTA)

    def test_a_member_scores_lower_than_a_nonmember(self):
        # TRAP 1, THE SIGN. L = log P(H0) - log P(H1), so mostly-yes (a member)
        # must score LOWER. An implementation that infers membership on a HIGH
        # value yields an inverted power curve that still looks like a curve.
        self.assertLess(self.stat([1] * 100), self.stat([1] * 80 + [0] * 20))

    def test_only_the_yes_count_matters_not_their_order(self):
        # The yes-count is a sufficient statistic; position must not weight.
        self.assertAlmostEqual(self.stat([1, 1, 0, 0, 1]),
                               self.stat([0, 1, 1, 1, 0]), places=15)

    def test_more_yes_responses_lower_the_statistic_monotonically(self):
        values = [self.stat([1] * k + [0] * (100 - k))
                  for k in range(0, 101, 10)]
        self.assertEqual(values, sorted(values, reverse=True))

    def test_zero_delta_is_refused(self):
        # delta = 0 makes a single "no" infinitely strong evidence, so L blows
        # up rather than being merely large.
        with self.assertRaises(ValueError):
            lrt.lrt_statistic([1, 0], self.d_n, self.d_n_minus_1, 0.0)

    def test_delta_of_one_or_more_is_refused(self):
        with self.assertRaises(ValueError):
            lrt.lrt_statistic([1, 0], self.d_n, self.d_n_minus_1, 1.0)


class TestRequiredQueries(unittest.TestCase):

    A, B, DELTA = 0.5, 2.0, 1e-3

    def n_for(self, beacon_size, power=0.95):
        return lrt.required_queries(lrt.d_n(self.A, self.B, beacon_size),
                                    lrt.d_n(self.A, self.B, beacon_size - 1),
                                    self.DELTA, alpha=0.05, power=power)

    def test_a_larger_beacon_needs_more_queries(self):
        # Each member is a smaller share of a bigger beacon, so it hides better.
        self.assertLess(self.n_for(50), self.n_for(500))

    def test_higher_power_needs_more_queries(self):
        self.assertLess(self.n_for(100, power=0.80),
                        self.n_for(100, power=0.95))

    def test_confusing_d_n_for_d_n_minus_1_changes_the_reported_number(self):
        # THE MUTATION CHECK THAT MATTERS. Reverting to D_N in the H1 slot does
        # not crash and still returns a plausible query count, so the guard has
        # to assert the NUMBER moves. If this ever passes trivially, the two
        # D values have collapsed and the test is no longer protecting anything.
        # Measured 2026-08-21: the correction is worth 0% at delta<=0.05 and
        # bites only when the error rate is high and the beacon small, because
        # the two D values differ by under 1% and delta scales that difference.
        # So the probe runs in the regime where it CAN discriminate; asserting
        # at delta=1e-6 would pass vacuously and protect nothing.
        small, high_error = 10, 0.6
        d_n = lrt.d_n(self.A, self.B, small)
        d_n_minus_1 = lrt.d_n(self.A, self.B, small - 1)
        self.assertNotEqual(d_n, d_n_minus_1)

        correct = lrt.required_queries(d_n, d_n_minus_1, high_error)
        defective = lrt.required_queries(d_n, d_n, high_error)
        self.assertNotEqual(correct, defective)
        self.assertGreater(correct, defective)

    def test_raises_when_the_hypotheses_do_not_separate(self):
        # delta near 1 means a member looks like a non-member; there is no
        # answer, and the function must say so rather than return a number.
        d_n = lrt.d_n(self.A, self.B, 100)
        with self.assertRaises(ValueError):
            lrt.required_queries(d_n, d_n / 0.5, 0.99)


class TestAttackEndToEnd(unittest.TestCase):
    """Milestone 1: does the attack actually achieve the power it predicts?"""

    A, B, N = 0.5, 2.0, 50
    DELTA, ALPHA, POWER = 1e-3, 0.05, 0.95
    REPLICATES = 3000

    def setUp(self):
        self.d_n = lrt.d_n(self.A, self.B, self.N, heterozygous=True)
        self.d_n_minus_1 = lrt.d_n(self.A, self.B, self.N - 1,
                                   heterozygous=True)
        self.n_queries = lrt.required_queries(
            self.d_n, self.d_n_minus_1, self.DELTA,
            alpha=self.ALPHA, power=self.POWER)
        self.threshold = lrt.rejection_threshold(
            self.n_queries, self.d_n, self.d_n_minus_1, self.DELTA,
            alpha=self.ALPHA)

    def rejection_rate(self, in_beacon):
        rejections = 0
        for seed in range(self.REPLICATES):
            responses = lrt.simulate_responses(
                self.n_queries, self.d_n, self.d_n_minus_1, self.DELTA,
                in_beacon=in_beacon, seed=seed)
            statistic = lrt.lrt_statistic(responses, self.d_n,
                                          self.d_n_minus_1, self.DELTA)
            rejections += statistic < self.threshold
        return rejections / self.REPLICATES

    def test_required_query_count_is_usable(self):
        self.assertGreater(self.n_queries, 0)
        self.assertLess(self.n_queries, 100_000)

    def test_empirical_power_reaches_the_requested_power(self):
        # The closed form is a normal approximation to a Binomial, so it is
        # allowed to be a little conservative — but not to fall short.
        self.assertGreaterEqual(self.rejection_rate(in_beacon=True),
                                self.POWER - 0.05)

    def test_empirical_false_positive_rate_respects_alpha(self):
        self.assertLessEqual(self.rejection_rate(in_beacon=False),
                             self.ALPHA + 0.03)


if __name__ == "__main__":
    unittest.main()


class TestWhichCorrectionActuallyMatters(unittest.TestCase):
    """Characterisation: pin the RELATIVE size of the two model corrections.

    Both were found on 2026-08-21 by checking docs/methods-lrt.md against the
    paper. They are wildly asymmetric in effect, and recording that here stops
    the project from re-litigating it from memory.
    """

    A, B, DELTA = 0.5, 2.0, 1e-6
    BEACON = 1895                     # the live AfriGen-D panel

    def n_queries(self, heterozygous):
        return lrt.required_queries(
            lrt.d_n(self.A, self.B, self.BEACON, heterozygous=heterozygous),
            lrt.d_n(self.A, self.B, self.BEACON - 1, heterozygous=heterozygous),
            self.DELTA)

    def test_d_n_minus_1_correction_is_negligible_at_realistic_error_rates(self):
        # Correct model, but it does NOT move this study's headline number.
        d_n = lrt.d_n(self.A, self.B, self.BEACON)
        d_n_minus_1 = lrt.d_n(self.A, self.B, self.BEACON - 1)
        correct = lrt.required_queries(d_n, d_n_minus_1, self.DELTA)
        defective = lrt.required_queries(d_n, d_n, self.DELTA)
        self.assertLess(abs(correct - defective) / correct, 0.01)

    def test_ascertainment_correction_changes_the_answer_by_orders_of_magnitude(self):
        # THE one that matters. Heterozygous sites are enriched for
        # intermediate frequencies, so "nobody carries it" gets much rarer and
        # each query is far less informative.
        self.assertGreater(self.n_queries(heterozygous=True),
                           100 * self.n_queries(heterozygous=False))

    def test_corrected_attack_needs_more_queries_than_a_day_of_rate_limit(self):
        # Ties the model to docs/beacon-constraints.md: production allows 50
        # queries/hour/IP, so the corrected attack is not merely slow, it is
        # infeasible against the live beacon by a wide margin.
        queries_per_day = 50 * 24
        self.assertGreater(self.n_queries(heterozygous=True), queries_per_day)


class TestPaperEquationA35(unittest.TestCase):
    """The paper's literal A35 and our plus-form are the same expression.

    Verified 2026-08-21 against the Europe PMC full-text JATS, which serves
    MathML and so preserves the radicals that PMC's HTML-to-text rendering
    drops:

        curl -sL https://www.ebi.ac.uk/europepmc/webservices/rest/\\
            PMC4667107/fullTextXML

    A35 as printed:

        n = [ ( z_alpha*sqrt(D_N (1 - D_N))
                - z_{1-beta}*sqrt(delta*D_{N-1} (1 - delta*D_{N-1})) )
              / ( delta*D_{N-1} - D_N ) ]^2

    Note the MINUS and the negative denominator. ``required_queries`` uses the
    conventional plus-form with a positive denominator. THEY ARE IDENTICAL:
    ``z_alpha`` there is the LOWER-tail quantile (negative), so both the
    numerator and the denominator flip sign together, and the square would
    remove any residual sign regardless.

    This test exists so that a reader who checks the code against A35, sees
    the sign difference and concludes the code is wrong does not "fix" it.
    """

    A, B = 0.5, 2.0

    @staticmethod
    def _a35_verbatim(d_n, d_n_minus_1, delta, alpha, power):
        """A35 with the paper's own signs, transcribed rather than rearranged."""
        z_alpha = lrt._z(alpha)              # LOWER tail: negative
        z_power = lrt._z(power)
        numerator = (z_alpha * math.sqrt(d_n * (1.0 - d_n))
                     - z_power * math.sqrt(delta * d_n_minus_1
                                           * (1.0 - delta * d_n_minus_1)))
        denominator = delta * d_n_minus_1 - d_n
        return (numerator / denominator) ** 2

    def test_matches_required_queries_across_the_parameter_space(self):
        for beacon in (85, 100, 1000, 1895):
            for delta in (1e-6, 1e-3, 0.05):
                for alpha, power in ((0.05, 0.95), (0.01, 0.80)):
                    with self.subTest(N=beacon, delta=delta, alpha=alpha):
                        d_n = lrt.d_n(self.A, self.B, beacon)
                        d_n_1 = lrt.d_n(self.A, self.B, beacon - 1)
                        paper = self._a35_verbatim(d_n, d_n_1, delta,
                                                   alpha, power)
                        ours = lrt.required_queries(d_n, d_n_1, delta,
                                                    alpha=alpha, power=power)
                        # required_queries applies ceil(); compare on the
                        # same footing rather than against the raw real value.
                        self.assertEqual(max(1, math.ceil(paper)), ours)
