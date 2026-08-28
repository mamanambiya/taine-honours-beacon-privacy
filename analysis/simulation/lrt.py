"""The Shringarpure & Bustamante (2015) membership-inference LRT.

Reference: doi:10.1016/j.ajhg.2015.09.010. See ../../docs/methods-lrt.md, which
reconstructs the statistic from the proposal because the equation objects did
not survive conversion out of Word.

**This is a reconstruction, not a transcription.** Every expression here must be
checked against the paper before any result derived from it is reported. The
tests pin internal consistency and the two traps; they cannot tell you the
model is the paper's model.

No numpy or scipy — the closed forms below are exact and `math.lgamma` is
enough. That keeps milestone 1 runnable anywhere, which is the point of doing it
first.
"""
import math
import random

__all__ = [
    "d_n",
    "theta_h0",
    "theta_h1",
    "lrt_statistic",
    "required_queries",
    "rejection_threshold",
    "simulate_responses",
]


def _log_beta(a, b):
    """log B(a, b), via lgamma so large 2N does not overflow."""
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def d_n(a, b, n_individuals, heterozygous=False):
    """P(none of the N genomes carries the alternate allele at a site).

    ``D_N = E_f[(1 - f)^(2N)]`` with ``f ~ Beta(a, b)``. This is an expectation
    over the site-frequency spectrum, NOT ``(1 - f_bar)^(2N)``: the SFS is the
    only route by which the cohort enters the test, so collapsing it to a mean
    would erase the very thing the study is measuring.

    For ``f ~ Beta(a, b)`` the expectation is exact::

        E[(1 - f)^m] = B(a, b + m) / B(a, b)

    ``heterozygous=True`` applies **trap 2**. Query sites are positions where
    the TARGET is heterozygous, and ``P(het | f) = 2f(1 - f)``, so the
    conditioned density is::

        Beta(a, b) * 2f(1 - f)  ∝  f^a (1 - f)^b  =  Beta(a + 1, b + 1)

    The factor 2 cancels in the normalisation. Het sites are enriched for
    intermediate frequencies, so the conditioned ``D_N`` is SMALLER — using the
    unconditioned SFS overstates ``D_N`` and therefore understates the attack's
    power, by an amount that depends on the SFS shape. That shape is exactly
    what differs between the African and European arms, so this is a confound,
    not a rounding error.
    """
    if a <= 0 or b <= 0:
        raise ValueError(f"Beta parameters must be positive, got a={a}, b={b}")
    if n_individuals <= 0:
        raise ValueError(f"n_individuals must be positive, got {n_individuals}")

    if heterozygous:
        a, b = a + 1.0, b + 1.0

    haplotypes = 2 * n_individuals
    return math.exp(_log_beta(a, b + haplotypes) - _log_beta(a, b))


def theta_h0(d_n_value):
    """P(beacon answers yes) when the target is NOT in it."""
    return 1.0 - d_n_value


def theta_h1(d_n_minus_1, delta):
    """P(beacon answers yes) when the target IS in it.

    Takes ``D_{N-1}``, NOT ``D_N``. Under H1 the target is one of the N and
    carries the allele at a query site by construction, so a "no" needs two
    things at once: the target's allele mis-called (``delta``) AND none of the
    **other N-1** individuals carrying it (``D_{N-1}``).

    Using ``D_N`` here — as an earlier revision of this module did — models the
    target as absent under both hypotheses and narrows the gap the test
    measures. See ../../docs/methods-lrt.md, "Verification status".
    """
    return 1.0 - delta * d_n_minus_1


def lrt_statistic(responses, d_n_value, d_n_minus_1, delta):
    """``L = log P(data | H0) - log P(data | H1)``.

    **Trap 1 — the sign.** SMALL ``L`` rejects ``H0``; small ``L`` means the
    target IS in the beacon. An implementation that infers membership on large
    ``L`` produces an inverted power curve that still looks like a power curve.

    Only ``sum(responses)`` is used — the yes-count is sufficient, which is what
    makes ``required_queries`` a closed form rather than a simulation.
    """
    if not 0.0 <= delta < 1.0:
        raise ValueError(f"delta must be in [0, 1), got {delta}")
    if delta == 0.0:
        raise ValueError(
            "delta = 0 makes L infinite for any 'no' response: under H1 a "
            "mismatch is then impossible, so a single no is infinitely strong "
            "evidence against membership. Use a small positive delta."
        )

    n = len(responses)
    yes = sum(responses)

    theta_0 = theta_h0(d_n_value)
    theta_1 = theta_h1(d_n_minus_1, delta)

    per_yes = math.log(theta_0 / theta_1)
    # log(D_N / (delta * D_{N-1})). This does NOT collapse to log(1/delta) —
    # that cancellation was the D_N-under-both-hypotheses bug showing through.
    per_no = math.log((1.0 - theta_0) / (1.0 - theta_1))

    return yes * per_yes + (n - yes) * per_no


def _z(p):
    """Standard-normal quantile. Acklam's rational approximation, |err| < 1.15e-9."""
    if not 0.0 < p < 1.0:
        raise ValueError(f"p must be in (0, 1), got {p}")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    p_low, p_high = 0.02425, 1 - 0.02425
    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > p_high:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def required_queries(d_n_value, d_n_minus_1, delta, alpha=0.05, power=0.95):
    """Queries needed for ``power`` at false-positive rate ``alpha``.

    ``L`` is a linear function of a Binomial count, so it is asymptotically
    Gaussian under both hypotheses. Rejecting for small ``L`` and equating the
    two expressions for the threshold gives::

        n = [ (z_{1-alpha} * sqrt(t0(1-t0)) + z_{power} * sqrt(t1(1-t1)))
              / (D_N - delta*D_{N-1}) ]^2

    The separation is ``D_N - delta*D_{N-1}``. An earlier revision wrote
    ``D_N * (1 - delta)``, which is what that expression collapses to only if
    both hypotheses share ``D_N`` — the bug, not a simplification.

    Note what cancels: the per-response log-ratios drop out entirely, so the
    answer depends only on the two yes-probabilities. Separation shrinks as the
    beacon grows, which is why a larger beacon hides its members better.

    VERIFIED against the paper's A35 (2026-08-21), transcribed from the Europe
    PMC full-text JATS, which serves MathML and so keeps the radicals that
    PMC's HTML-to-text rendering drops.

    A35 prints a MINUS between the two z terms over a NEGATIVE denominator
    ``(delta*D_{N-1} - D_N)``. This function uses the conventional plus-form
    over a positive denominator. **They are the same expression** — the
    paper's ``z_alpha`` is the lower-tail quantile, so numerator and
    denominator flip together and the square removes any residual sign. Do NOT
    "correct" this toward A35's literal signs; TestPaperEquationA35 asserts the
    two agree across the parameter space.
    """
    t0 = theta_h0(d_n_value)
    t1 = theta_h1(d_n_minus_1, delta)

    separation = d_n_value - delta * d_n_minus_1
    if separation <= 0:
        raise ValueError(
            f"hypotheses not separated: D_N={d_n_value:.6g} is not greater "
            f"than delta*D_(N-1)={delta * d_n_minus_1:.6g}"
        )

    numerator = (_z(1.0 - alpha) * math.sqrt(t0 * (1.0 - t0))
                 + _z(power) * math.sqrt(t1 * (1.0 - t1)))
    return max(1, math.ceil((numerator / separation) ** 2))


def rejection_threshold(n_queries, d_n_value, d_n_minus_1, delta, alpha=0.05):
    """The ``L`` below which ``H0`` is rejected, from the H0 sampling distribution.

    Reject for SMALL ``L`` (trap 1), so the critical region is the lower tail.
    """
    t0 = theta_h0(d_n_value)
    t1 = theta_h1(d_n_minus_1, delta)

    per_yes = math.log(t0 / t1)
    per_no = math.log((1.0 - t0) / (1.0 - t1))
    slope = per_yes - per_no          # change in L per additional yes

    mean = n_queries * (per_no + slope * t0)
    sd = abs(slope) * math.sqrt(n_queries * t0 * (1.0 - t0))
    return mean - _z(1.0 - alpha) * sd


def simulate_responses(n_queries, d_n_value, d_n_minus_1, delta, in_beacon,
                       seed=None):
    """Draw beacon answers for a target that is or is not a member.

    Membership is fixed by construction, which is what lets the sign of ``L`` be
    asserted rather than assumed.
    """
    theta = (theta_h1(d_n_minus_1, delta) if in_beacon
             else theta_h0(d_n_value))
    rng = random.Random(seed)
    return [1 if rng.random() < theta else 0 for _ in range(n_queries)]
