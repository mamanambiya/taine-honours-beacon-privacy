import json, os, sys, tempfile, unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import binomial_delta as bd
import real_beacon_attack as rba


class BinomialDelta(unittest.TestCase):
    def test_delta_zero_reproduces_reject_on_any_no(self):
        outs = [(1000, 30), (900, 20), (1100, 45)]
        for n in (1, 5, 40):
            self.assertAlmostEqual(bd.fpr(outs, n, 0.0, 0.95), rba.fpr(outs, n), places=9)
        self.assertEqual(bd.budget(outs, 0.05, 0.0, 0.95, 500), rba.budget_for_fpr(outs, 0.05, 500))

    def test_cutoff_grows_with_delta_and_n(self):
        self.assertEqual(bd.accept_cutoff(100, 0.0, 0.95), 0)
        self.assertGreater(bd.accept_cutoff(100, 0.01, 0.95), 0)
        self.assertGreaterEqual(bd.accept_cutoff(300, 0.01, 0.95), bd.accept_cutoff(100, 0.01, 0.95))

    def test_hypergeometric_sums_to_one_and_edges(self):
        self.assertAlmostEqual(bd.p_at_most(10, 20, 5, 10), 1.0, places=9)
        self.assertAlmostEqual(bd.p_at_most(0, 20, 0, 10), 1.0)  # no absent sites: never a no
        self.assertAlmostEqual(bd.p_at_most(0, 10, 10, 1), 0.0)  # every site absent
        with self.assertRaises(ValueError):
            bd.p_at_most(0, 5, 1, 6)

    def test_bad_arguments(self):
        with self.assertRaises(ValueError):
            bd.accept_cutoff(10, 1.0, 0.95)
        with self.assertRaises(ValueError):
            bd.accept_cutoff(10, 0.01, 1.0)
        with self.assertRaises(ValueError):
            bd.fpr([], 5, 0.0, 0.95)

    def test_more_error_allowed_never_needs_fewer_queries(self):
        outs = [(2000, 60)] * 5
        b0 = bd.budget(outs, 0.05, 0.0, 0.95, 2000)
        b1 = bd.budget(outs, 0.05, 0.01, 0.95, 2000)
        self.assertGreaterEqual(b1, b0)

    def test_power_is_a_step_at_the_budget_when_delta_is_zero(self):
        outs = [(1000, 30), (900, 20), (1100, 45)]
        b = bd.budget(outs, 0.05, 0.0, 0.95, 500)
        self.assertEqual(bd.power_at(outs, b - 1, 0.0, 0.05), 0.0)
        self.assertEqual(bd.power_at(outs, b, 0.0, 0.05), 1.0)

    def test_power_with_delta_is_between_0_and_1_and_rises_with_n(self):
        outs = [(2000, 60)] * 5
        values = [bd.power_at(outs, n, 0.01, 0.05) for n in (60, 100, 200, 400)]
        self.assertTrue(all(0.0 <= v <= 1.0 for v in values))
        self.assertLessEqual(values[0], values[-1])
        self.assertGreater(values[-1], 0.0)


if __name__ == "__main__":
    unittest.main()
