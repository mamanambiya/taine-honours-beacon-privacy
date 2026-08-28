"""Tests for the frequency-spectrum contrast."""
import unittest

import contrast


class TestRareFraction(unittest.TestCase):

    def test_more_skewed_spectra_carry_more_rare_variants(self):
        # The axis has to be monotone or the sweep is unreadable.
        loads = [contrast.rare_fraction(a, b, 40)
                 for a, b in ((2.0, 2.0), (0.5, 2.0), (0.2, 20.0))]
        self.assertEqual(loads, sorted(loads))

    def test_a_symmetric_spectrum_has_almost_no_rare_variants(self):
        self.assertLess(contrast.rare_fraction(2.0, 2.0, 40), 0.01)


class TestContrast(unittest.TestCase):

    SHAPES = [("common", 0.5, 2.0), ("rare", 0.2, 20.0)]

    def test_reports_one_row_per_shape_with_a_valid_power(self):
        rows = contrast.contrast(self.SHAPES, n_variants=4000)
        self.assertEqual([r["label"] for r in rows], ["common", "rare"])
        for r in rows:
            self.assertGreaterEqual(r["power"], 0.0)
            self.assertLessEqual(r["power"], 1.0)

    def test_power_moves_with_the_spectrum(self):
        # THE STUDY'S PREMISE. If power did not move with the spectrum, no
        # ancestry effect could exist downstream, because the spectrum is the
        # only route by which a cohort enters the model.
        rows = contrast.contrast(self.SHAPES, n_variants=4000)
        by_label = {r["label"]: r["power"] for r in rows}
        self.assertGreater(by_label["rare"], by_label["common"])


class TestBudget(unittest.TestCase):

    def test_a_rare_spectrum_is_attackable_within_the_search_range(self):
        self.assertIsNotNone(contrast.budget_for_power(0.2, 20.0,
                                                       n_variants=8000))

    def test_a_common_spectrum_is_not(self):
        # Returning None is a result, not a failure: the attack cannot exploit
        # a spectrum with almost no rare variants at any budget searched.
        self.assertIsNone(contrast.budget_for_power(2.0, 2.0,
                                                    n_variants=8000))


if __name__ == "__main__":
    unittest.main()
