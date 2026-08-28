"""Build an experimental beacon from genotypes and run the attack against it.

Milestone 1 validated the LRT against its own closed form. This module closes
the loop: the beacon's answers come from an actual genotype matrix rather than
from ``D_N``, so the model and the data can disagree — which is the only way to
notice that they do.

numpy for the genotype matrix; the statistic itself stays in ``lrt``.
"""

import numpy as np

import lrt


def simulate_cohort(n_individuals, n_variants, a, b, seed):
    """Genotypes as alternate-allele counts, drawn under the Beta SFS.

    Each variant gets a frequency ``f ~ Beta(a, b)``, then each individual gets
    ``Binomial(2, f)`` copies — Hardy-Weinberg at that frequency.
    """
    rng = np.random.default_rng(seed)
    frequencies = rng.beta(a, b, size=n_variants)
    return rng.binomial(2, frequencies, size=(n_individuals, n_variants))


def alternate_allele_counts(genotypes):
    """Alternate copies observed at each variant, summed over individuals."""
    return genotypes.sum(axis=0)


def beacon_answers(genotypes, censor_singletons=False):
    """What the beacon replies at each position.

    ``censor_singletons`` is the defence Shringarpure & Bustamante tested: a
    beacon that always answers "no" where only ONE alternate copy exists in the
    whole cohort, since those are the positions that most nearly identify a
    single person.
    """
    counts = alternate_allele_counts(genotypes)
    return counts > (1 if censor_singletons else 0)


def heterozygous_positions(genotypes, individual):
    """Positions where ``individual`` carries exactly one alternate copy.

    The adversary queries only these, which is what makes the query set
    ascertained and the SFS conditioned — see docs/methods-lrt.md, Trap 2.
    """
    return np.flatnonzero(genotypes[individual] == 1)


def _rows_excluding_target(genotypes, target_genotypes):
    """The beacon minus the target, if the target is in it.

    Under H1 a "no" needs the target's own copy mis-called AND none of the
    OTHER N-1 carrying it, so the two groups must be kept apart.
    """
    matches = np.flatnonzero((genotypes == target_genotypes).all(axis=1))
    if matches.size == 0:
        return genotypes, False
    return np.delete(genotypes, matches[0], axis=0), True


def attack(genotypes, target_genotypes, n_queries, delta, seed, a, b,
           censor_singletons=False):
    """Query the beacon at the target's heterozygous sites, return ``L``.

    LOW ``L`` is evidence of membership (docs/methods-lrt.md, Trap 1).

    ``delta`` corrupts only the TARGET's contribution — that is what makes
    ``P(no | member) = delta * D_{N-1}`` rather than ``delta``. Corrupting the
    beacon's whole answer instead would model a different, easier attack.
    """
    rng = np.random.default_rng(seed)

    candidates = heterozygous_positions(target_genotypes.reshape(1, -1), 0)
    if candidates.size == 0:
        raise ValueError("target has no heterozygous sites to query")
    replace = candidates.size < n_queries
    queried = rng.choice(candidates, size=n_queries, replace=replace)

    others, target_in_beacon = _rows_excluding_target(genotypes,
                                                      target_genotypes)
    others_count = alternate_allele_counts(others)[queried]

    # The beacon counts alternate copies and answers against a threshold. A
    # sequencing error loses the TARGET's copies from that count, not the
    # beacon's whole answer — which is what makes P(no | member) equal
    # delta * D_{N-1} rather than delta.
    if target_in_beacon:
        intact = rng.random(n_queries) >= delta
        target_count = np.where(intact, target_genotypes[queried], 0)
    else:
        target_count = 0

    threshold = 1 if censor_singletons else 0
    responses = (others_count + target_count) > threshold

    n_beacon = genotypes.shape[0]
    d_n = lrt.d_n(a, b, n_beacon, heterozygous=True)
    d_n_minus_1 = lrt.d_n(a, b, n_beacon - 1, heterozygous=True)
    return lrt.lrt_statistic(responses.astype(int).tolist(), d_n,
                             d_n_minus_1, delta)


def empirical_power(beacon, outsiders, n_queries, delta, a, b, alpha=0.05,
                    censor_singletons=False):
    """Power at a false-positive rate fixed EMPIRICALLY, not from the model.

    The threshold comes from the observed non-member distribution rather than
    from the Gaussian approximation. That matters whenever the beacon does
    something the attacker's model does not know about — singleton censoring
    being exactly such a case, since the LRT is derived for an uncensored
    beacon and is therefore mis-specified under it.

    Under mis-specification the raw separation between member and non-member
    scores can move in EITHER direction and is not a measure of anything.
    Power at a controlled false-positive rate still is.

    Returns ``(power, achieved_false_positive_rate)``.
    """
    settings = dict(a=a, b=b, censor_singletons=censor_singletons)
    null = np.array([
        attack(beacon, outsiders[i], n_queries, delta, seed=5000 + i,
               **settings)
        for i in range(len(outsiders))])
    member = np.array([
        attack(beacon, beacon[i], n_queries, delta, seed=i, **settings)
        for i in range(len(beacon))])

    threshold = np.quantile(null, alpha)      # reject for SMALL L
    return float((member < threshold).mean()), float((null < threshold).mean())
