"""The membership attack, run against a loaded beacon's real answers.

A boolean beacon answers yes to "does anyone here carry this allele?" exactly
when the allele is in the set of sites its members carry (``beacon_check.py``
verifies that on the loaded datasets). The attacker knows one genome and may
query only that person's heterozygous positions. With the mismatch rate
``delta`` negligible the likelihood-ratio test of Shringarpure and Bustamante
reduces to one rule: reject membership of anyone who gets a single "no". A
member never does, because every allele they carry is in the set. So power is
1 by construction, and the quantity that varies is the false-positive rate:
the share of non-members whose first ``n`` queried sites are all yes.

For one outsider with ``H`` heterozygous sites, ``m`` of which are absent from
the beacon, drawing ``n`` sites without replacement gives::

    P(all n yes) = prod_{i<n} (H - m - i) / (H - i)

so the false-positive rate at budget ``n`` is that, averaged over outsiders, and
no random queries need to be sent. ``budget_for_fpr`` returns the smallest ``n``
at or below the target rate; ``bootstrap_budget`` resamples the outsiders to put
an interval on it.

Standard library only, so it runs on the VM that holds the beacon.
"""

import json
import os
import random
import sys

__all__ = [
    "site_key",
    "p_all_yes",
    "fpr",
    "budget_for_fpr",
    "bootstrap_budget",
    "count_person",
]

_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def site_key(chrom, pos, ref, alt):
    """One integer per (chromosome, 1-based position, ref, alt)."""
    c = str(chrom)
    c = c[3:] if c.lower().startswith("chr") else c
    try:
        number = int(c)
    except ValueError:
        raise ValueError(f"not an autosome number: {chrom!r}")
    if ref not in _BASE or alt not in _BASE:
        raise ValueError(f"not a single A/C/G/T base: ref={ref!r} alt={alt!r}")
    return (((number << 29) | int(pos)) << 4) | (_BASE[ref] << 2) | _BASE[alt]


def p_all_yes(het_sites, absent_sites, n_queries):
    """P(every one of ``n_queries`` distinct het sites answers yes)."""
    if n_queries < 0 or n_queries > het_sites:
        raise ValueError(
            f"{n_queries} queries from {het_sites} heterozygous sites is not possible"
        )
    present = het_sites - absent_sites
    if n_queries > present:
        return 0.0
    p = 1.0
    for i in range(n_queries):
        p *= (present - i) / (het_sites - i)
    return p


def fpr(outsiders, n_queries):
    """Mean P(all yes) over ``outsiders``, a list of ``(het, absent)`` counts."""
    if not outsiders:
        raise ValueError("no outsiders")
    return sum(p_all_yes(h, m, n_queries) for h, m in outsiders) / len(outsiders)


def budget_for_fpr(outsiders, alpha, max_queries):
    """Smallest budget whose false-positive rate is at or below ``alpha``.

    ``None`` when no budget up to ``max_queries`` gets there: that is a result
    (some outsiders look like members at every budget), not an error.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    if not outsiders:
        raise ValueError("no outsiders")
    limit = min(max_queries, min(h for h, _ in outsiders))
    for n in range(1, limit + 1):
        if fpr(outsiders, n) <= alpha:
            return n
    return None


def bootstrap_budget(outsiders, alpha, max_queries, draws=1000, seed=1):
    """95% interval for the budget, resampling the outsiders with replacement.

    A resample with no reachable budget counts as ``max_queries``, so the upper
    end is a floor, not a guess, when that happens.
    """
    rng = random.Random(seed)
    values = []
    for _ in range(draws):
        sample = [outsiders[rng.randrange(len(outsiders))] for _ in outsiders]
        b = budget_for_fpr(sample, alpha, max_queries)
        values.append(max_queries if b is None else b)
    values.sort()
    return values[int(0.025 * draws)], values[min(draws - 1, int(0.975 * draws))]


def count_person(path, beacon_sites):
    """``(het sites, those absent from the beacon)`` for one person's table."""
    het = absent = 0
    with open(path) as handle:
        for line in handle:
            chrom, pos, ref, alt = line.rstrip("\n").split("\t")[:4]
            het += 1
            if site_key(chrom, pos, ref, alt) not in beacon_sites:
                absent += 1
    if het == 0:
        raise ValueError(f"{path}: no heterozygous sites, cannot score this person")
    return het, absent


def load_beacon_sites(path):
    """Read ``reference_name,start,reference_bases,alternate_bases`` rows.

    ``start`` is the beacon's 0-based coordinate; the key uses the 1-based
    position, the same as the het tables.
    """
    sites = set()
    with open(path) as handle:
        header = next(handle).rstrip("\n").split(",")
        want = ["reference_name", "start", "reference_bases", "alternate_bases"]
        cols = [header.index(c) for c in want]
        for line in handle:
            f = line.rstrip("\n").split(",")
            chrom, start, ref, alt = (f[c] for c in cols)
            if len(ref) == 1 and len(alt) == 1:
                sites.add(site_key(chrom, int(start) + 1, ref, alt))
    return sites


def _main(argv):
    import argparse

    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "--sites", required=True, help="beacon sites CSV (see load_beacon_sites)"
    )
    p.add_argument("--dataset", required=True)
    p.add_argument("--members", required=True, help="file with the beacon's members")
    p.add_argument("--outsiders", required=True, help="file with outsider sample ids")
    p.add_argument(
        "--het-dir",
        action="append",
        required=True,
        help="directory holding <sample>.tsv het tables (repeatable)",
    )
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--max-queries", type=int, default=2000)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    def table(sample):
        for d in a.het_dir:
            path = os.path.join(d, f"{sample}.tsv")
            if os.path.exists(path):
                return path
        raise FileNotFoundError(f"no het table for {sample} in {a.het_dir}")

    sites = load_beacon_sites(a.sites)
    ids = lambda f: [l.strip() for l in open(f) if l.strip()]
    members = [count_person(table(s), sites) for s in ids(a.members)]
    outs = [count_person(table(s), sites) for s in ids(a.outsiders)]
    members_with_a_no = sum(1 for _, m in members if m > 0)
    point = budget_for_fpr(outs, a.alpha, a.max_queries)
    low, high = bootstrap_budget(outs, a.alpha, a.max_queries)
    report = {
        "dataset_id": a.dataset,
        "beacon_sites": len(sites),
        "members": len(members),
        "members_with_any_absent_het_site": members_with_a_no,
        "member_power": 1.0 if members_with_a_no == 0 else None,
        "outsiders": len(outs),
        "alpha": a.alpha,
        "budget_for_fpr_alpha": point,
        "budget_95ci_bootstrap_over_outsiders": [low, high],
        "outsider_absent_fraction_mean": sum(m / h for h, m in outs) / len(outs),
        "fpr_by_budget": {
            str(n): fpr(outs, n)
            for n in (5, 10, 20, 40, 60, 80, 100, 150, 200, 300)
            if n <= min(h for h, _ in outs)
        },
        "outsider_counts": outs,
    }
    with open(a.out, "w") as handle:
        json.dump(report, handle, indent=2)
    print(
        a.dataset,
        "budget",
        point,
        "95% CI",
        (low, high),
        "member_power",
        report["member_power"],
        "outsiders",
        len(outs),
    )


if __name__ == "__main__":
    _main(sys.argv[1:])
