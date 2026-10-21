"""Seeded split of one population's samples into beacon members and held-out controls.

The attack needs two groups from the SAME population: people whose genotypes
are in the beacon (members) and people whose are not (controls). A beacon that
holds a whole population has no controls from that population, and using
another population's people as controls would confound ancestry with the thing
being measured.

Deterministic on purpose. The result depends on the SET of sample ids, the
population name and the seed, never on the order a VCF happens to list them, so
the same split can be regenerated and checked years later. Standard library only.
"""

import argparse
import hashlib
import json
import os
import random
import sys

__all__ = ["split_samples", "read_samples", "main"]


def split_samples(sample_ids, n_members, n_controls, population, seed):
    """Return ``(members, controls)``, each sorted, disjoint, exactly sized.

    Too few samples raises rather than shrinking a group: a population with
    fewer controls than the others is not comparable with them.
    """
    if n_members < 1 or n_controls < 1:
        raise ValueError(
            f"group sizes must be positive, got members={n_members}, "
            f"controls={n_controls}"
        )
    ids = list(sample_ids)
    if not ids:
        raise ValueError("no samples to split")
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate sample ids")
    if len(ids) < n_members + n_controls:
        raise ValueError(
            f"{population}: {len(ids)} samples cannot supply "
            f"{n_members} members and {n_controls} controls"
        )

    ordered = sorted(ids)
    # A string seed is hashed by Random, so this is stable across runs and
    # Python versions; mixing in the population keeps the ten splits independent.
    random.Random(f"{population}:{seed}").shuffle(ordered)
    members = sorted(ordered[:n_members])
    controls = sorted(ordered[n_members : n_members + n_controls])
    return members, controls


def read_samples(path):
    """Sample ids from a file, one per line; blank lines and padding ignored."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _sha256_of(ids):
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Split a population's samples into members and controls."
    )
    parser.add_argument("--population", required=True)
    parser.add_argument(
        "--samples", required=True, help="one sample id per line (bcftools query -l)"
    )
    parser.add_argument("--members", type=int, required=True)
    parser.add_argument("--controls", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args(argv)

    ids = read_samples(args.samples)
    members, controls = split_samples(
        ids, args.members, args.controls, args.population, args.seed
    )

    prefix = os.path.join(args.out_dir, args.population)
    targets = [
        f"{prefix}_members.txt",
        f"{prefix}_controls.txt",
        f"{prefix}_split.json",
    ]
    existing = [t for t in targets if os.path.exists(t)]
    if existing:
        sys.exit(
            f"refusing to overwrite an existing split: {existing}. "
            "A dataset built from it was verified against these lists."
        )

    os.makedirs(args.out_dir, exist_ok=True)
    for path, group in ((targets[0], members), (targets[1], controls)):
        with open(path, "w") as fh:
            fh.write("\n".join(group) + "\n")
    record = {
        "population": args.population,
        "seed": args.seed,
        "n_members": len(members),
        "n_controls": len(controls),
        "n_input_samples": len(ids),
        "input_sha256": _sha256_of(ids),
        "samples_file": args.samples,
    }
    with open(targets[2], "w") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")
    print(
        f"{args.population}: {len(members)} members, {len(controls)} "
        f"controls from {len(ids)} samples (seed {args.seed})"
    )
    return record


if __name__ == "__main__":
    main()
