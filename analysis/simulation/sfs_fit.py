"""Fit the Beta site-frequency spectrum the attack model assumes to a real panel.

``lrt.d_n`` is an expectation over ``f ~ Beta(a, b)``, and until now ``(a, b)``
has been a fixture, ``(0.5, 2.0)``. This module replaces the fixture with a fit
to the alternate-allele counts of a reference panel. It is Milestone 1b.

The model. Each site has a population frequency ``f ~ Beta(a, b)``; the panel's
``2N`` haplotypes then carry ``K ~ Binomial(2N, f)`` copies, so ``K`` is
beta-binomial. The fit works from the histogram of ``K`` over the sites the
panel lists, which is all that needs to leave the cluster.

**The panel lists only the sites it saw.** A VCF holds polymorphic sites, so
``K = 0`` never appears, and a panel filtered to a minimum allele count loses
more. The likelihood is therefore the beta-binomial conditioned on
``K >= min_count`` — the zero-truncated form when ``min_count = 1``. Fitting the
untruncated form to such data is badly biased (``TestTruncationMatters``). The
mass removed at ``K = 0`` is exactly ``lrt.d_n(a, b, N)``, so the truncation and
the attack model share one quantity.

numpy only. No scipy: the optimiser is a small Nelder-Mead, because the problem
is two parameters and smooth, and that keeps the package installable anywhere.
"""

import argparse
import json
import math
import sys

import numpy as np

import lrt

__all__ = [
    "log_pmf",
    "count_histogram",
    "read_histogram",
    "fit_beta_binomial",
    "diagnostics",
    "required_queries_for_fit",
    "compare_to_fixture",
    "FIXTURE",
]

# The spectrum every absolute figure used before this module was written.
FIXTURE = (0.5, 2.0)

# Search is in log-parameters. exp(+-7) spans 0.0009 .. 1096, far wider than any
# real spectrum; a result touching either end is a bound, not an estimate.
_LOG_BOUND = 7.0
_BOUNDARY_TOLERANCE = 0.05

_lgamma = np.vectorize(math.lgamma, otypes=[float])


def _log_beta(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def log_pmf(k, n_individuals, a, b):
    """``log P(K = k)`` for ``f ~ Beta(a, b)``, ``K | f ~ Binomial(2N, f)``.

    Closed form ``C(2N, k) * B(k + a, 2N - k + b) / B(a, b)``, in log space
    because the beta function underflows long before ``2N = 3790``.
    """
    k = np.asarray(k, dtype=float)
    n = 2.0 * n_individuals
    log_choose = math.lgamma(n + 1.0) - _lgamma(k + 1.0) - _lgamma(n - k + 1.0)
    return (
        log_choose
        + _lgamma(k + a)
        + _lgamma(n - k + b)
        - math.lgamma(n + a + b)
        - _log_beta(a, b)
    )


def count_histogram(counts, n_individuals):
    """Number of variants at each alternate-allele count ``0 .. 2N``.

    A count above ``2N`` means the allele number was not ``2N`` — missing
    genotypes, or the wrong ``N`` — and is refused rather than clipped.
    """
    counts = np.asarray(counts)
    if counts.size and counts.min() < 0:
        raise ValueError("allele counts cannot be negative")
    if counts.size and counts.max() > 2 * n_individuals:
        raise ValueError(
            f"allele count {int(counts.max())} exceeds 2N = "
            f"{2 * n_individuals}: wrong N, or sites with missing genotypes"
        )
    return np.bincount(counts.astype(int), minlength=2 * n_individuals + 1)


def read_histogram(path, n_individuals):
    """Read ``count<TAB>variants`` lines; repeated counts are summed.

    Summing is what lets per-chromosome tables be concatenated. ``#`` lines and
    blanks are ignored.
    """
    hist = np.zeros(2 * n_individuals + 1, dtype=np.int64)
    with open(path) as fh:
        for number, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            try:
                k, c = int(fields[0]), int(fields[1])
            except (ValueError, IndexError):
                raise ValueError(
                    f"{path}:{number}: expected 'count<TAB>variants', got {line!r}"
                )
            if not 0 <= k <= 2 * n_individuals:
                raise ValueError(
                    f"{path}:{number}: count {k} outside 0..{2 * n_individuals}"
                )
            if c < 0:
                raise ValueError(f"{path}:{number}: negative variant count")
            hist[k] += c
    return hist


def _validated(histogram, n_individuals, min_count):
    hist = np.asarray(histogram)
    if hist.ndim != 1 or len(hist) != 2 * n_individuals + 1:
        raise ValueError(
            f"histogram must have 2N + 1 = {2 * n_individuals + 1} bins, "
            f"got shape {hist.shape}"
        )
    if (hist < 0).any():
        raise ValueError("histogram has negative entries")
    if min_count > 0 and hist[:min_count].any():
        raise ValueError(
            f"histogram has variants below min_count = {min_count}: the panel "
            "was not filtered the way the caller says"
        )
    if hist.sum() == 0:
        raise ValueError("histogram is empty")
    if np.count_nonzero(hist) < 2:
        raise ValueError(
            "only one occupied bin: the likelihood is flat along a ridge, so "
            "no spectrum can be fitted"
        )
    return hist.astype(float)


def _log_normaliser(n_individuals, a, b, min_count):
    """``log P(K >= min_count)``; zero when nothing was truncated."""
    if min_count <= 0:
        return 0.0
    below = np.exp(log_pmf(np.arange(min_count), n_individuals, a, b)).sum()
    return math.log1p(-min(below, 1.0 - 1e-300))


def _neg_log_likelihood(x, hist, ks, n_individuals, min_count):
    clipped = np.clip(x, -_LOG_BOUND, _LOG_BOUND)
    # Outside the bounds the objective is the bound's value plus a steep
    # penalty, so the simplex is pushed back rather than left to wander.
    penalty = 1e4 * float(((x - clipped) ** 2).sum())
    a, b = math.exp(clipped[0]), math.exp(clipped[1])
    logp = log_pmf(ks, n_individuals, a, b)
    norm = _log_normaliser(n_individuals, a, b, min_count)
    return -(float((hist[ks.astype(int)] * logp).sum()) - hist.sum() * norm) + penalty


def _nelder_mead(fun, start, step=0.5, tolerance=1e-10, max_iterations=2000):
    """Minimal Nelder-Mead. Returns ``(x, value, converged)``."""
    n = len(start)
    simplex = [np.array(start, dtype=float)]
    for i in range(n):
        vertex = np.array(start, dtype=float)
        vertex[i] += step
        simplex.append(vertex)
    values = [fun(v) for v in simplex]

    for _ in range(max_iterations):
        order = np.argsort(values)
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        if (
            abs(values[-1] - values[0]) <= tolerance * (abs(values[0]) + 1.0)
            and max(np.abs(v - simplex[0]).max() for v in simplex[1:]) <= 1e-7
        ):
            return simplex[0], values[0], True

        centroid = np.mean(simplex[:-1], axis=0)
        reflected = centroid + (centroid - simplex[-1])
        f_reflected = fun(reflected)
        if f_reflected < values[0]:
            expanded = centroid + 2.0 * (centroid - simplex[-1])
            f_expanded = fun(expanded)
            if f_expanded < f_reflected:
                simplex[-1], values[-1] = expanded, f_expanded
            else:
                simplex[-1], values[-1] = reflected, f_reflected
        elif f_reflected < values[-2]:
            simplex[-1], values[-1] = reflected, f_reflected
        else:
            if f_reflected < values[-1]:
                contracted = centroid + 0.5 * (reflected - centroid)
            else:
                contracted = centroid + 0.5 * (simplex[-1] - centroid)
            f_contracted = fun(contracted)
            if f_contracted < min(f_reflected, values[-1]):
                simplex[-1], values[-1] = contracted, f_contracted
            else:
                for i in range(1, n + 1):
                    simplex[i] = simplex[0] + 0.5 * (simplex[i] - simplex[0])
                    values[i] = fun(simplex[i])
    order = int(np.argmin(values))
    return simplex[order], values[order], False


def _standard_errors(fun, x, step=1e-3):
    """Delta-method SEs of ``(a, b)`` from a finite-difference Hessian.

    ``nan`` when the Hessian is not positive definite, which is what a flat
    likelihood looks like.
    """
    hessian = np.zeros((2, 2))
    f0 = fun(x)
    for i in range(2):
        for j in range(i, 2):
            ei = np.zeros(2)
            ej = np.zeros(2)
            ei[i] = step
            ej[j] = step
            if i == j:
                value = (fun(x + ei) - 2.0 * f0 + fun(x - ei)) / step**2
            else:
                value = (
                    fun(x + ei + ej)
                    - fun(x + ei - ej)
                    - fun(x - ei + ej)
                    + fun(x - ei - ej)
                ) / (4 * step**2)
            hessian[i, j] = hessian[j, i] = value
    try:
        covariance = np.linalg.inv(hessian)
        variances = np.diag(covariance)
        if (variances <= 0).any():
            raise np.linalg.LinAlgError
    except np.linalg.LinAlgError:
        return float("nan"), float("nan")
    return (
        math.exp(x[0]) * math.sqrt(variances[0]),
        math.exp(x[1]) * math.sqrt(variances[1]),
    )


def fit_beta_binomial(histogram, n_individuals, min_count=1):
    """Maximum-likelihood ``(a, b)`` from a histogram of allele counts.

    ``min_count`` is the smallest count the panel can contain: ``1`` for a
    panel of polymorphic sites, ``2`` if singletons were filtered out, ``0`` to
    treat the histogram as complete (only correct if it genuinely is).

    Returns a dict. ``converged`` is the optimiser's own verdict;
    ``at_boundary`` says an estimate sits on the edge of the search range and is
    therefore a bound, not a measurement. Check both before quoting ``a`` or
    ``b``, and read ``diagnostics`` — a single Beta fits real spectra only
    approximately, and the fit will return numbers whether or not it should.
    """
    hist = _validated(histogram, n_individuals, min_count)
    ks = np.flatnonzero(hist)
    ks = ks[ks >= min_count]

    def objective(x):
        return _neg_log_likelihood(x, hist, ks, n_individuals, min_count)

    best = None
    for start in (
        [math.log(0.5), math.log(2.0)],
        [math.log(0.1), math.log(10.0)],
        [math.log(1.0), math.log(1.0)],
    ):
        x, value, converged = _nelder_mead(objective, start)
        if best is None or value < best[1]:
            best = (x, value, converged)
    x, value, converged = best

    x = np.clip(x, -_LOG_BOUND, _LOG_BOUND)
    at_boundary = bool((np.abs(np.abs(x) - _LOG_BOUND) < _BOUNDARY_TOLERANCE).any())
    se_a, se_b = _standard_errors(objective, x)
    return {
        "a": math.exp(x[0]),
        "b": math.exp(x[1]),
        "se_a": se_a,
        "se_b": se_b,
        "loglik": -value,
        "converged": bool(converged),
        "at_boundary": at_boundary,
        "n_variants": int(hist.sum()),
        "n_individuals": int(n_individuals),
        "min_count": int(min_count),
    }


def _expected(hist, n_individuals, a, b, min_count):
    ks = np.arange(len(hist))
    logp = log_pmf(ks, n_individuals, a, b) - _log_normaliser(
        n_individuals, a, b, min_count
    )
    expected = np.exp(logp) * float(hist[min_count:].sum())
    expected[:min_count] = 0.0
    return expected


def _bin_edges(low, high):
    """Roughly geometric inclusive bins covering ``low .. high``."""
    edges, k = [], max(low, 0)
    while k <= high:
        edges.append(k)
        k = max(k + 1, int(round(k * math.sqrt(2.0))))
    return edges


def diagnostics(histogram, n_individuals, a, b, min_count=1, rare_threshold=0.01):
    """How well ``Beta(a, b)`` describes the observed histogram.

    The fit is a two-parameter summary of a spectrum that is not exactly Beta.
    What matters for the attack is whether it gets the **rare end** right —
    singletons drive ``D_N`` — and whether the **near-fixed alternate alleles**
    (sites where the reference allele is the rare one) are absent from the
    model. Compare the observed and fitted fractions below before using the
    parameters.

    ``singleton_*`` is ``None`` when singletons were filtered out of the panel.
    """
    hist = np.asarray(histogram, dtype=float)
    n_hap = 2 * n_individuals
    total = float(hist[min_count:].sum())
    expected = _expected(hist, n_individuals, a, b, min_count)
    ks = np.arange(len(hist))
    minor = np.minimum(ks, n_hap - ks) / n_hap
    high = ks / n_hap > 0.5

    def share(weights, mask):
        return float(weights[mask].sum() / weights.sum())

    observed = hist.copy()
    observed[:min_count] = 0.0

    bins = []
    edges = _bin_edges(max(min_count, 0), n_hap)
    for i, low in enumerate(edges):
        upper = (edges[i + 1] - 1) if i + 1 < len(edges) else n_hap
        bins.append(
            {
                "k_low": int(low),
                "k_high": int(upper),
                "observed": float(observed[low : upper + 1].sum()),
                "expected": float(expected[low : upper + 1].sum()),
            }
        )

    singletons_kept = min_count <= 1
    return {
        "n_variants": int(total),
        "singleton_fraction_observed": (
            float(observed[1] / total) if singletons_kept else None
        ),
        "singleton_fraction_fitted": (
            float(expected[1] / total) if singletons_kept else None
        ),
        "rare_fraction_observed": share(observed, minor < rare_threshold),
        "rare_fraction_fitted": share(expected, minor < rare_threshold),
        "high_frequency_fraction_observed": share(observed, high),
        "high_frequency_fraction_fitted": share(expected, high),
        "bins": bins,
    }


def required_queries_for_fit(a, b, n_individuals, delta, alpha=0.05, power=0.95):
    """Queries for ``power`` at ``alpha`` under ``Beta(a, b)``, ascertained.

    The same route the documented 108,898 figure took: query sites are
    heterozygous in the target, so ``D_N`` is taken under the ``(a+1, b+1)``
    spectrum (docs/methods-lrt.md, trap 2).
    """
    d_n = lrt.d_n(a, b, n_individuals, heterozygous=True)
    d_n_minus_1 = lrt.d_n(a, b, n_individuals - 1, heterozygous=True)
    return lrt.required_queries(d_n, d_n_minus_1, delta, alpha, power)


def compare_to_fixture(
    a, b, n_individuals, delta, alpha=0.05, power=0.95, fixture=FIXTURE
):
    """How far the fitted spectrum moves the budget from the fixture's."""
    fitted = required_queries_for_fit(a, b, n_individuals, delta, alpha, power)
    reference = required_queries_for_fit(
        fixture[0], fixture[1], n_individuals, delta, alpha, power
    )
    return {
        "fixture": {"a": fixture[0], "b": fixture[1]},
        "fixture_queries": reference,
        "fitted_queries": fitted,
        "ratio_fitted_to_fixture": fitted / reference,
        "delta": delta,
        "alpha": alpha,
        "power": power,
        "n_individuals": n_individuals,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Fit Beta(a, b) to an allele-count histogram and cost the "
        "attack under it."
    )
    parser.add_argument(
        "histogram", help="two columns: alternate-allele count, variants"
    )
    parser.add_argument(
        "--n-individuals",
        type=int,
        required=True,
        help="N: individuals in the panel (2N haplotypes)",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="smallest count the panel can contain "
        "(1 = polymorphic sites; 2 = singletons removed)",
    )
    parser.add_argument(
        "--label", default="", help="what this panel is, recorded in the output"
    )
    parser.add_argument("--delta", type=float, default=1e-6)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--power", type=float, default=0.95)
    parser.add_argument("--out", help="write the JSON record here")
    args = parser.parse_args(argv)

    hist = read_histogram(args.histogram, args.n_individuals)
    fit = fit_beta_binomial(hist, args.n_individuals, args.min_count)
    record = {
        "label": args.label,
        "source": args.histogram,
        "n_individuals": args.n_individuals,
        "min_count": args.min_count,
        "fit": fit,
        "diagnostics": diagnostics(
            hist, args.n_individuals, fit["a"], fit["b"], args.min_count
        ),
        "fixture_comparison": compare_to_fixture(
            fit["a"], fit["b"], args.n_individuals, args.delta, args.alpha, args.power
        ),
    }

    diag = record["diagnostics"]
    print(
        f"{args.label or args.histogram}: a' = {fit['a']:.4g} "
        f"(se {fit['se_a']:.2g}), b' = {fit['b']:.4g} (se {fit['se_b']:.2g}) "
        f"from {fit['n_variants']:,} variants"
    )
    if not fit["converged"] or fit["at_boundary"]:
        print(
            f"WARNING: converged={fit['converged']} "
            f"at_boundary={fit['at_boundary']} — do not quote these "
            "parameters",
            file=sys.stderr,
        )
    print(
        "  high-frequency (alt > 50%) share: observed "
        f"{diag['high_frequency_fraction_observed']:.3f}, fitted "
        f"{diag['high_frequency_fraction_fitted']:.3f}"
    )
    comparison = record["fixture_comparison"]
    print(
        f"  queries for {args.power:.0%} power: "
        f"{comparison['fitted_queries']:,} fitted vs "
        f"{comparison['fixture_queries']:,} fixture "
        f"(x{comparison['ratio_fitted_to_fixture']:.3g})"
    )

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(record, fh, indent=2, default=float)
            fh.write("\n")
    return record


if __name__ == "__main__":
    main()
