"""Does attack power depend on the shape of the frequency spectrum?

That is the study's question stripped of ancestry labels. The spectrum is the
only route by which a cohort enters the model, so if power does NOT move with
it, no ancestry difference can exist downstream — and the whole hypothesis is
dead regardless of which cohorts are used.

This runs entirely on simulated cohorts. It says nothing about African versus
European data until a real spectrum is fitted; what it establishes is whether
the effect is detectable AT ALL, and how large a difference in rare-variant
load is needed to produce one.
"""

import numpy as np

import cohort


def rare_fraction(a, b, n_individuals, threshold=0.01, seed=0, draws=200_000):
    """Share of variants below ``threshold`` minor allele frequency.

    A readable axis for the sweep: Beta parameters are hard to reason about,
    "12% of variants are rare" is not.
    """
    rng = np.random.default_rng(seed)
    f = rng.beta(a, b, size=draws)
    return float(np.mean(np.minimum(f, 1.0 - f) < threshold))


def contrast(shapes, n_beacon=40, n_outsiders=200, n_variants=20000,
             n_queries=40, delta=1e-3, alpha=0.05, seed=2026):
    """Empirical attack power for each ``(label, a, b)`` spectrum shape.

    Power is measured at a threshold calibrated on the observed non-member
    distribution, not from the closed form — that is what makes it a statement
    about an adversary rather than about the model.
    """
    rows = []
    for label, a, b in shapes:
        pool = cohort.simulate_cohort(n_beacon + n_outsiders, n_variants,
                                      a, b, seed=seed)
        beacon, outsiders = pool[:n_beacon], pool[n_beacon:]
        power, achieved_fpr = cohort.empirical_power(
            beacon, outsiders, n_queries, delta, a, b, alpha=alpha)
        rows.append({
            "label": label,
            "a": a,
            "b": b,
            "rare_fraction": rare_fraction(a, b, n_beacon),
            "power": power,
            "false_positive_rate": achieved_fpr,
        })
    return rows


def budget_for_power(a, b, target=0.95, n_beacon=40, n_outsiders=200,
                     n_variants=20000, delta=1e-3, alpha=0.05, seed=2026,
                     budgets=(5, 10, 20, 40, 80, 160, 320, 640)):
    """Smallest query budget in ``budgets`` reaching ``target`` power.

    Returns None when no budget in the range gets there — which is itself the
    result for a spectrum the attack cannot exploit.
    """
    pool = cohort.simulate_cohort(n_beacon + n_outsiders, n_variants, a, b,
                                  seed=seed)
    beacon, outsiders = pool[:n_beacon], pool[n_beacon:]
    for budget in budgets:
        power, _ = cohort.empirical_power(beacon, outsiders, budget, delta,
                                          a, b, alpha=alpha)
        if power >= target:
            return budget
    return None
