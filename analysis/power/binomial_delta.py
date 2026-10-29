"""Query budgets when the attacker's genome may disagree with the beacon.

``real_beacon_attack`` takes the mismatch rate ``delta`` as zero, so the test
is "reject on any no". Shringarpure and Bustamante instead let each queried
genotype be wrong with probability ``delta`` (sequencing or calling error) and
test the number of "no" answers. Here a person asks ``n`` distinct heterozygous
sites. A true member's number of nos is Binomial(n, delta); the rule is

    call "member" when nos <= c, with c the smallest value for which a member
    is accepted with probability at least ``power`` (default 0.95).

For an outsider with ``H`` heterozygous sites, ``m`` of them absent from the
beacon, the number of nos in ``n`` draws without replacement is hypergeometric,
so the false-positive rate is the mean over outsiders of P(nos <= c). With
``delta = 0`` then ``c = 0`` and this reproduces ``real_beacon_attack``.

Reads the ``outsider_counts`` already saved by ``real_beacon_attack``; sends no
queries. Standard library only.
"""

import json
import math
import os
import sys

__all__ = ["accept_cutoff", "p_at_most", "fpr", "budget", "power_at"]


def _log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _binom_cdf(c, n, delta):
    if delta == 0.0:
        return 1.0
    return sum(
        math.exp(_log_comb(n, j) + j * math.log(delta) + (n - j) * math.log1p(-delta))
        for j in range(0, min(c, n) + 1)
    )


def accept_cutoff(n, delta, power):
    """Smallest c with P(Binomial(n, delta) <= c) >= power."""
    if not 0.0 <= delta < 1.0:
        raise ValueError(f"delta must be in [0, 1), got {delta}")
    if not 0.0 < power < 1.0:
        raise ValueError(f"power must be in (0, 1), got {power}")
    c = 0
    while c < n and _binom_cdf(c, n, delta) < power:
        c += 1
    return c


def p_at_most(c, het, absent, n):
    """P(at most ``c`` of ``n`` distinct het sites are absent): hypergeometric."""
    if n > het:
        raise ValueError(f"{n} queries from {het} heterozygous sites is not possible")
    total = 0.0
    base = _log_comb(het, n)
    for j in range(0, min(c, n, absent) + 1):
        if n - j > het - absent:
            continue
        total += math.exp(_log_comb(absent, j) + _log_comb(het - absent, n - j) - base)
    return min(total, 1.0)


def fpr(outsiders, n, delta, power):
    if not outsiders:
        raise ValueError("no outsiders")
    c = accept_cutoff(n, delta, power)
    return sum(p_at_most(c, h, m, n) for h, m in outsiders) / len(outsiders)


def budget(outsiders, alpha, delta, power, max_queries):
    """Smallest n with false-positive rate <= alpha, else None."""
    limit = min(max_queries, min(h for h, _ in outsiders))
    for n in range(1, limit + 1):
        if fpr(outsiders, n, delta, power) <= alpha:
            return n
    return None


def power_at(outsiders, n, delta, alpha):
    """Empirical power at ``n`` queries with the false-positive rate held <= alpha.

    The cutoff is the largest ``c`` (call "member" when nos <= c) whose
    false-positive rate over ``outsiders`` is still <= ``alpha``; power is then
    P(Binomial(n, delta) <= c), the share of true members accepted. When even
    ``c = 0`` has a false-positive rate above ``alpha`` no test with that
    control exists at this ``n`` and the power is 0. With ``delta = 0`` the
    power is therefore 0 below the budget and 1 from it: a step.
    """
    if not outsiders:
        raise ValueError("no outsiders")
    c = -1
    while c < n:
        rate = sum(p_at_most(c + 1, h, m, n) for h, m in outsiders) / len(outsiders)
        if rate > alpha:
            break
        c += 1
    return 0.0 if c < 0 else _binom_cdf(c, n, delta)


def _main(argv):
    import argparse

    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--attack-dir", required=True, help="folder of attack_*.json")
    p.add_argument("--out", required=True)
    p.add_argument("--deltas", default="0,0.001,0.003,0.01")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--power", type=float, default=0.95)
    p.add_argument("--max-queries", type=int, default=2000)
    a = p.parse_args(argv)
    deltas = [float(x) for x in a.deltas.split(",")]
    rows = {}
    for name in sorted(os.listdir(a.attack_dir)):
        if not (name.startswith("attack_") and name.endswith(".json")):
            continue
        rep = json.load(open(os.path.join(a.attack_dir, name)))
        outs = [tuple(x) for x in rep["outsider_counts"]]
        rows[name[7:-5]] = {
            "outsiders": len(outs),
            "budget_by_delta": {
                str(d): budget(outs, a.alpha, d, a.power, a.max_queries)
                for d in deltas
            },
        }
    json.dump(
        {"alpha": a.alpha, "member_power": a.power, "rows": rows},
        open(a.out, "w"),
        indent=2,
    )
    print("delta ->", ", ".join(str(d) for d in deltas))
    for k, v in rows.items():
        print(k.ljust(18), [v["budget_by_delta"][str(d)] for d in deltas])


if __name__ == "__main__":
    _main(sys.argv[1:])
