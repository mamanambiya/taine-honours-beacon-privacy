"""Tests for the genotype-level beacon and the attack against it."""
import unittest

import numpy as np

import cohort


class TestSimulateCohort(unittest.TestCase):

    def test_shape_is_individuals_by_variants(self):
        g = cohort.simulate_cohort(20, 500, 0.5, 2.0, seed=1)
        self.assertEqual(g.shape, (20, 500))

    def test_genotypes_are_allele_counts_zero_one_or_two(self):
        g = cohort.simulate_cohort(20, 500, 0.5, 2.0, seed=1)
        self.assertTrue(set(np.unique(g)).issubset({0, 1, 2}))

    def test_is_reproducible_from_the_seed(self):
        first = cohort.simulate_cohort(10, 100, 0.5, 2.0, seed=7)
        second = cohort.simulate_cohort(10, 100, 0.5, 2.0, seed=7)
        np.testing.assert_array_equal(first, second)

    def test_allele_frequencies_track_the_requested_beta(self):
        # A cohort whose SFS does not match the Beta the LRT assumes would make
        # every downstream power number meaningless.
        a, b = 0.5, 2.0
        g = cohort.simulate_cohort(400, 4000, a, b, seed=3)
        freq = cohort.alternate_allele_counts(g) / (2 * g.shape[0])
        self.assertAlmostEqual(float(freq.mean()), a / (a + b), delta=0.03)


class TestBeaconAnswers(unittest.TestCase):

    def test_yes_exactly_when_someone_carries_the_alternate_allele(self):
        g = np.array([[0, 1, 0],
                      [0, 0, 0],
                      [0, 0, 2]])
        np.testing.assert_array_equal(cohort.beacon_answers(g),
                                      np.array([False, True, True]))

    def test_censoring_suppresses_singletons_only(self):
        # One alternate COPY in the whole cohort is a singleton. Two copies —
        # whether one homozygote or two heterozygotes — is not.
        g = np.array([[0, 1, 1, 2],
                      [0, 0, 1, 0]])
        np.testing.assert_array_equal(
            cohort.beacon_answers(g, censor_singletons=True),
            np.array([False, False, True, True]))

    def test_censoring_never_creates_a_yes(self):
        g = cohort.simulate_cohort(30, 800, 0.5, 2.0, seed=11)
        plain = cohort.beacon_answers(g)
        censored = cohort.beacon_answers(g, censor_singletons=True)
        self.assertTrue(bool(np.all(plain | ~censored)))
        self.assertLess(int(censored.sum()), int(plain.sum()))


class TestHeterozygousPositions(unittest.TestCase):

    def test_selects_only_positions_with_exactly_one_alternate_copy(self):
        g = np.array([[0, 1, 2, 1]])
        np.testing.assert_array_equal(
            cohort.heterozygous_positions(g, 0), np.array([1, 3]))


class TestAttack(unittest.TestCase):
    """The end-to-end claim: a member scores lower than a non-member."""

    A, B, DELTA = 0.5, 2.0, 1e-3
    N, VARIANTS = 40, 6000

    def setUp(self):
        self.pool = cohort.simulate_cohort(self.N + 20, self.VARIANTS,
                                           self.A, self.B, seed=99)
        self.beacon = self.pool[:self.N]          # in the beacon
        self.outsiders = self.pool[self.N:]       # not in it

    def test_members_score_lower_than_nonmembers(self):
        members = [cohort.attack(self.beacon, self.beacon[i], 200,
                                 self.DELTA, seed=i, a=self.A, b=self.B)
                   for i in range(15)]
        outsiders = [cohort.attack(self.beacon, self.outsiders[i], 200,
                                   self.DELTA, seed=100 + i, a=self.A, b=self.B)
                     for i in range(15)]
        self.assertLess(float(np.mean(members)), float(np.mean(outsiders)))

    def test_separation_widens_with_more_queries(self):
        def gap(n_queries):
            m = np.mean([cohort.attack(self.beacon, self.beacon[i], n_queries,
                                       self.DELTA, seed=i, a=self.A, b=self.B) for i in range(10)])
            o = np.mean([cohort.attack(self.beacon, self.outsiders[i],
                                       n_queries, self.DELTA, seed=50 + i, a=self.A, b=self.B)
                         for i in range(10)])
            return float(o - m)
        self.assertGreater(gap(400), gap(50))


class TestSingletonCensoring(unittest.TestCase):
    """Milestone 2: what does the censoring defence actually cost the attacker?

    Measured as POWER AT A CONTROLLED FALSE-POSITIVE RATE, not as the gap
    between member and non-member scores. The attacker's LRT is derived for an
    UNCENSORED beacon, so under censoring it is mis-specified and its raw gap
    moves in either direction — measured 2026-08-21, censoring *widened* the
    gap by 33% under a common-variant SFS while doing nothing for real power.
    A mean gap is not a power measure.
    """

    DELTA, ALPHA = 1e-3, 0.05
    N_BEACON, N_OUT, VARIANTS = 40, 200, 20000
    RARE = (0.2, 20.0)          # mean frequency ~0.01; singletons common
    COMMON = (0.5, 2.0)         # mean frequency 0.2; singletons scarce

    def cohorts(self, sfs):
        pool = cohort.simulate_cohort(self.N_BEACON + self.N_OUT,
                                      self.VARIANTS, *sfs, seed=2026)
        return pool[:self.N_BEACON], pool[self.N_BEACON:]

    def power(self, sfs, n_queries, censor):
        beacon, outsiders = self.cohorts(sfs)
        power, _ = cohort.empirical_power(
            beacon, outsiders, n_queries, self.DELTA, *sfs,
            alpha=self.ALPHA, censor_singletons=censor)
        return power

    def test_the_attack_is_powerless_against_a_common_variant_spectrum(self):
        # With mean frequency 0.2 nearly every position answers yes whoever
        # asks, so the responses carry no membership information at all. Attack
        # power is a property of RARE-variant load — which is the mechanism the
        # study's ancestry hypothesis runs through.
        self.assertLess(self.power(self.COMMON, 300, censor=False), 0.10)

    def test_the_attack_succeeds_against_a_rare_variant_spectrum(self):
        self.assertGreaterEqual(self.power(self.RARE, 160, censor=False), 0.95)

    def test_censoring_raises_the_query_budget_the_attacker_needs(self):
        # Uncensored reaches full power by ~20 queries; censored needs ~160.
        # Compared at 40, where the margin is wide enough not to be noise.
        self.assertGreaterEqual(self.power(self.RARE, 40, censor=False), 0.95)
        self.assertLess(self.power(self.RARE, 40, censor=True), 0.75)

    def test_censoring_delays_the_attack_rather_than_preventing_it(self):
        # A cost, not a cure — which is the argument against treating a single
        # MAF threshold as a solution.
        self.assertGreaterEqual(self.power(self.RARE, 160, censor=True), 0.95)


if __name__ == "__main__":
    unittest.main()


class TestWhatTheAdversaryActuallyNeeds(unittest.TestCase):
    """Why two mutations to `attack` survive the suite: the behaviour is real.

    Both were checked by injecting the defect and observing that nothing
    changed. Rather than tightening the tests until they fail on a difference
    that does not exist, the invariance itself is asserted here.
    """

    A, B = 0.2, 20.0
    N_BEACON, N_OUT, VARIANTS = 40, 200, 20000

    def setUp(self):
        pool = cohort.simulate_cohort(self.N_BEACON + self.N_OUT,
                                      self.VARIANTS, self.A, self.B, seed=2026)
        self.beacon = pool[:self.N_BEACON]
        self.outsiders = pool[self.N_BEACON:]

    def power_assuming(self, a, b, delta=1e-3, n_queries=40):
        power, _ = cohort.empirical_power(self.beacon, self.outsiders,
                                          n_queries, delta, a, b, alpha=0.05)
        return power

    def test_empirical_power_does_not_depend_on_the_assumed_sfs(self):
        # L is a monotone linear function of the yes-count, and the SFS only
        # sets its coefficients. A threshold calibrated on observed
        # non-members therefore absorbs any D_N whatsoever.
        #
        # CONSEQUENCE FOR THE STUDY: the 885x ascertainment effect is a fact
        # about the PAPER'S PREDICTED BUDGET, not about what an adversary
        # achieves. A real attacker calibrates empirically and never needs D_N.
        # Operational risk must be measured, not read off required_queries.
        truth = self.power_assuming(self.A, self.B)
        self.assertEqual(truth, self.power_assuming(0.5, 2.0))
        self.assertEqual(truth, self.power_assuming(5.0, 5.0))

    def test_sequencing_error_is_inert_at_realistic_rates(self):
        # Same shape as the D_{N-1} result: structurally necessary, numerically
        # irrelevant until delta is large.
        self.assertEqual(self.power_assuming(self.A, self.B, delta=1e-6),
                         self.power_assuming(self.A, self.B, delta=1e-3))

    def test_a_high_error_rate_does_degrade_the_attack(self):
        # Controls the test above: if this did not drop, the delta pathway
        # would be unexercised and the inertness result would be vacuous.
        self.assertLess(self.power_assuming(self.A, self.B, delta=0.5),
                        self.power_assuming(self.A, self.B, delta=1e-3))
