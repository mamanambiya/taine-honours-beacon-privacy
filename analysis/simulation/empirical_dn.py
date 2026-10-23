"""``D_N`` taken straight from a panel's allele-count histogram, with no Beta.

``lrt.d_n`` is an expectation over a fitted ``Beta(a, b)``. On the real
V6HC-S_AFR panel the fit lands on its search bound (``at_boundary = true``) and
the fitted counts are far from the observed ones, so a single Beta is the wrong
shape for that panel. This module skips the model: it treats the panel's own
sites as the population and asks, for a beacon of ``n`` people, how often none
of them carries the allele at a site where the attacked target is
heterozygous. It is the replacement ``docs/fitting-the-sfs.md`` ("Limits")
points to.

The calculation. The panel has ``M = 2N`` haplotypes and ``h_k`` sites with
``k`` copies of the alternate allele. A site is a query site when the target is
heterozygous there (docs/methods-lrt.md, trap 2). Under random pairing of the
panel's haplotypes::

    P(target heterozygous | k copies)  =  2 k (M - k) / (M (M - 1))

Given a heterozygous target, two haplotypes are fixed (one alternate, one
reference), leaving ``M - 2`` with ``k - 1`` alternates. A beacon of ``n`` other
people draws ``2n`` of them without replacement, so::

    P(none carries | k)  =  C(M - 1 - k, 2n) / C(M - 2, 2n)

and ``D_N`` is the het-weighted mean of that over the histogram::

    D_N = sum_k h_k w_k p_k / sum_k h_k w_k

**Why the K = 0 mass is not needed.** A panel lists only the sites it saw, so
``K = 0`` never appears and the Beta fit has to model the missing mass. Here the
weight ``w_k`` is zero at ``k = 0`` and at ``k = M``, so those bins drop out and
nothing unobserved enters. For the same reason the histogram is enough: counts
cannot tell a heterozygote from a homozygote, and the HWE weight above is the
only place that matters.

**Why drawing from the panel is unbiased.** The panel's haplotypes are draws
from the population, so ``2n`` of them drawn without replacement are again iid
draws. Given a site's population frequency ``f`` the expected ``p_k`` is
exactly ``(1 - f)^(2n)``, which is the quantity ``lrt.d_n`` integrates. The
estimate therefore agrees with the Beta answer when the Beta is right
(``TestAgreesWithTheParametricModelWhenItIsRight``) and keeps working when it
is not (``TestWhereASingleBetaFails``).

``D_{N-1}`` is the same function at ``n - 1``: under H1 the target is one of the
beacon's ``n`` people and the others number ``n - 1`` (``lrt.theta_h1``).

The beacon cannot be as large as the panel: the target's two haplotypes and the
beacon's ``2n`` must come from the panel, so ``n <= N - 1``.

numpy only.
"""

import math

import numpy as np

import lrt

__all__ = ["empirical_d_n", "required_queries_empirical"]

_lgamma = np.vectorize(math.lgamma, otypes=[float])


def empirical_d_n(histogram, panel_individuals, n_beacon):
    """``P(no beacon member carries the allele | target heterozygous there)``.

    ``histogram[k]`` is the number of variants with ``k`` alternate copies, over
    ``2 * panel_individuals`` haplotypes (the ``read_histogram`` layout).
    """
    hist = np.asarray(histogram, dtype=float)
    if panel_individuals < 2:
        raise ValueError(
            f"panel too small: {panel_individuals} individuals cannot hold a "
            "target and a beacon"
        )
    if n_beacon < 1:
        raise ValueError(f"n_beacon must be positive, got {n_beacon}")
    if n_beacon > panel_individuals - 1:
        raise ValueError(
            f"beacon of {n_beacon} does not fit in a panel of "
            f"{panel_individuals}: the target and the beacon must be different "
            f"people, so the beacon is at most {panel_individuals - 1}"
        )
    m = 2 * panel_individuals
    if hist.ndim != 1 or len(hist) != m + 1:
        raise ValueError(
            f"histogram must have 2N + 1 = {m + 1} bins, got shape {hist.shape}"
        )
    if (hist < 0).any():
        raise ValueError("histogram has negative entries")

    k = np.arange(m + 1, dtype=float)
    weight = hist * (2.0 * k * (m - k) / (m * (m - 1.0)))
    total = weight.sum()
    if total <= 0:
        raise ValueError(
            "no site can be heterozygous in a target: every variant sits at "
            "count 0 or at the full allele number"
        )

    draw = 2 * n_beacon
    free = m - 1.0 - k  # reference haplotypes left after the target
    fits = free >= draw
    # log C(free, draw) - log C(M - 2, draw); the draw! terms cancel.
    log_p = np.full(m + 1, -np.inf)
    log_p[fits] = (
        _lgamma(free[fits] + 1.0)
        - _lgamma(free[fits] - draw + 1.0)
        - (math.lgamma(m - 1.0) - math.lgamma(m - 1.0 - draw))
    )
    p = np.exp(log_p)
    return float((weight * p).sum() / total)


def required_queries_empirical(
    histogram, panel_individuals, n_beacon, delta, alpha=0.05, power=0.95
):
    """Queries for ``power`` at ``alpha`` using the histogram's own ``D_N``."""
    d_n = empirical_d_n(histogram, panel_individuals, n_beacon)
    d_n_minus_1 = empirical_d_n(histogram, panel_individuals, n_beacon - 1)
    return lrt.required_queries(d_n, d_n_minus_1, delta, alpha, power)
