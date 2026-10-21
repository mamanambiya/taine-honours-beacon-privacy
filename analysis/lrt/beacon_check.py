"""Check a loaded beacon dataset against genotype ground truth.

A beacon is only a usable oracle for the attack if its answers equal what the
genotypes say: "yes" exactly where a member carries the alternate allele. The
population beacons first loaded on afrigend-all-dev answered yes at sites no
member carried, so every power figure built on them would have been noise. This
harness is the gate that would have caught it, and it runs before any result is
quoted.

It uses the files ``build_population_dataset.sh`` writes: sites known present,
sites known absent, and each individual's heterozygous positions. Standard
library only.

Three controls run on every dataset, and each one exists because of a failure
that reads as a result:

* a positive control FIRST (a site a member carries must answer yes), so an
  off-by-one in the coordinates or a dataset that did not load stops the run
  before it can look like "this population is safer";
* a negative control (sites no member carries must answer no), which is the
  defect actually found;
* a scoping control (an unknown dataset id must answer no), so a beacon that
  ignores ``datasetIds`` cannot pass by answering over the whole panel.
"""

import argparse
import json
import os
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

__all__ = [
    "BeaconClient",
    "BeaconError",
    "RateLimited",
    "ControlFailed",
    "ScopingFailed",
    "read_sites",
    "check_known_sites",
    "check_dataset_scoping",
    "check_individuals",
    "main",
]

DEFAULT_WAIT_SECONDS = 60.0


class BeaconError(Exception):
    """The beacon answered something that is not a yes or a no."""


class RateLimited(BeaconError):
    """429 persisted through every permitted retry."""


class ControlFailed(Exception):
    """The positive control did not answer yes; nothing after it is trustworthy."""


class ScopingFailed(Exception):
    """A query for a dataset that does not exist answered yes."""


def _default_fetch(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            return (response.status, dict(response.headers), json.load(response))
    except urllib.error.HTTPError as error:
        return error.code, dict(error.headers), None


def _retry_after(headers):
    for key in ("Retry-After", "retry-after"):
        if key in headers:
            try:
                return max(float(headers[key]), 1.0)
            except (TypeError, ValueError):
                break
    return DEFAULT_WAIT_SECONDS


class BeaconClient:
    """One boolean question at a time, with a bounded, polite retry on 429.

    A retry loop that hammers a limiter has taken a beacon offline before, so
    the wait is always the server's own ``Retry-After`` and the number of
    retries is bounded.
    """

    def __init__(self, base_url, fetch=None, sleep=time.sleep, max_retries=5):
        self.base_url = base_url.rstrip("/")
        self.fetch = fetch or _default_fetch
        self.sleep = sleep
        self.max_retries = max_retries
        self.queries_sent = 0

    def exists(self, dataset_id, chrom, pos, ref, alt):
        """Does ``dataset_id`` hold this allele? ``pos`` is the VCF POS."""
        params = {
            "assemblyId": "GRCh38",
            "referenceName": str(chrom).removeprefix("chr"),
            "start": int(pos) - 1,  # the beacon is 0-based
            "referenceBases": ref,
            "alternateBases": alt,
            "datasetIds": dataset_id,
        }
        url = f"{self.base_url}/api/g_variants?{urllib.parse.urlencode(params)}"
        for attempt in range(self.max_retries + 1):
            self.queries_sent += 1
            status, headers, body = self.fetch(url)
            if status == 429:
                if attempt == self.max_retries:
                    raise RateLimited(
                        f"still rate limited after {self.max_retries} retries"
                    )
                self.sleep(_retry_after(headers or {}))
                continue
            if status != 200:
                raise BeaconError(f"HTTP {status} for {url}")
            try:
                answer = body["responseSummary"]["exists"]
            except (TypeError, KeyError):
                raise BeaconError(f"no responseSummary.exists in answer to {url}")
            if not isinstance(answer, bool):
                raise BeaconError(f"exists is {answer!r}, not a boolean")
            return answer
        raise RateLimited("unreachable")  # pragma: no cover


def read_sites(path):
    """``chrom<TAB>POS<TAB>REF<TAB>ALT`` lines as ``(chrom, int POS, ref, alt)``."""
    sites = []
    with open(path) as fh:
        for number, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            try:
                chrom, pos, ref, alt = fields[:4]
                sites.append((chrom, int(pos), ref, alt))
            except (ValueError, TypeError):
                raise ValueError(
                    f"{path}:{number}: expected "
                    f"'chrom<TAB>POS<TAB>REF<TAB>ALT', got {line!r}"
                )
    return sites


def check_known_sites(client, dataset_id, present, absent):
    """Sites a member carries must answer yes; sites nobody carries, no.

    Raises ``ControlFailed`` if the FIRST present site does not answer yes: that
    is the positive control, and it is checked before anything else is asked.
    """
    if not present:
        raise ValueError("no present sites: the positive control has nothing to test")
    if not absent:
        raise ValueError(
            "no absent sites: the negative control has nothing to "
            "test, so a pass would say nothing"
        )

    control = present[0]
    if not client.exists(dataset_id, *control):
        raise ControlFailed(
            f"positive control failed: {control} is carried by a member but "
            f"{dataset_id} answered no (coordinates, assembly or load)"
        )

    present_missed = [s for s in present[1:] if not client.exists(dataset_id, *s)]
    absent_wrongly_present = [s for s in absent if client.exists(dataset_id, *s)]
    return {
        "present_checked": len(present),
        "absent_checked": len(absent),
        "present_missed": present_missed,
        "absent_wrongly_present": absent_wrongly_present,
        "passed": not present_missed and not absent_wrongly_present,
    }


def check_dataset_scoping(client, dataset_id, site):
    """A dataset id that does not exist must answer no for a real site."""
    probe = f"{dataset_id}__NO_SUCH_DATASET"
    if client.exists(probe, *site):
        raise ScopingFailed(
            f"{probe} answered yes: the beacon ignores datasetIds, so every "
            "answer is about the whole panel and not about this population"
        )
    return {"probe": probe, "passed": True}


def _sample_sites(sites, k, seed, name):
    if len(sites) <= k:
        return list(sites)
    chosen = random.Random(f"{seed}:{name}").sample(range(len(sites)), k)
    return [sites[i] for i in sorted(chosen)]


def check_individuals(
    client, dataset_id, het_tables, members, controls, per_individual, seed
):
    """What the attack will ask: each person's own heterozygous positions.

    Every member site must answer yes (the member's own allele is in the
    dataset by construction), so a single miss is a load or coordinate bug, never
    biology. Controls must NOT all answer yes: if they do, the beacon cannot
    tell members from outsiders.
    """

    def run(names):
        queries = yes = 0
        misses = []
        for name in names:
            for site in _sample_sites(het_tables[name], per_individual, seed, name):
                queries += 1
                if client.exists(dataset_id, *site):
                    yes += 1
                else:
                    misses.append((name, site))
        return queries, yes, misses

    member_queries, member_yes, member_misses = run(members)
    control_queries, control_yes, _ = run(controls)

    member_rate = member_yes / member_queries if member_queries else None
    control_rate = control_yes / control_queries if control_queries else None
    passed = (
        member_queries > 0
        and control_queries > 0
        and not member_misses
        and control_rate < 1.0
    )
    return {
        "member_queries": member_queries,
        "member_yes_rate": member_rate,
        "member_misses": member_misses,
        "control_queries": control_queries,
        "control_yes_rate": control_rate,
        "passed": passed,
    }


def _read_list(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def main(argv=None, fetch=None, sleep=time.sleep):
    parser = argparse.ArgumentParser(
        description="Check a loaded beacon dataset against the build output."
    )
    parser.add_argument(
        "--base-url", required=True, help="for example http://127.0.0.1:8086"
    )
    parser.add_argument(
        "--build-dir",
        required=True,
        help="the OUTDIR build_population_dataset.sh wrote",
    )
    parser.add_argument(
        "--dataset", help="dataset id (default: from build_summary.json)"
    )
    parser.add_argument("--per-individual", type=int, default=50)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", help="write the JSON report here")
    args = parser.parse_args(argv)

    build = args.build_dir
    with open(os.path.join(build, "build_summary.json")) as fh:
        summary = json.load(fh)
    dataset_id = args.dataset or summary["dataset_id"]
    members = _read_list(os.path.join(build, "members.txt"))
    controls = _read_list(os.path.join(build, "controls.txt"))
    present = read_sites(os.path.join(build, "present_sites.tsv"))
    absent = read_sites(os.path.join(build, "absent_sites.tsv"))
    het = {
        name: read_sites(os.path.join(build, "het", f"{name}.tsv"))
        for name in members + controls
    }

    client = BeaconClient(args.base_url, fetch=fetch, sleep=sleep)
    report = {"dataset_id": dataset_id, "base_url": args.base_url, "passed": False}
    try:
        report["known_sites"] = check_known_sites(client, dataset_id, present, absent)
        report["scoping"] = check_dataset_scoping(client, dataset_id, present[0])
        report["individuals"] = check_individuals(
            client, dataset_id, het, members, controls, args.per_individual, args.seed
        )
        report["passed"] = (
            report["known_sites"]["passed"]
            and report["scoping"]["passed"]
            and report["individuals"]["passed"]
        )
    except (ControlFailed, ScopingFailed, RateLimited, BeaconError) as error:
        report["error"] = str(error)
    report["queries_sent"] = client.queries_sent

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=2, default=str)
            fh.write("\n")
    status = "PASS" if report["passed"] else "FAIL"
    print(
        f"{dataset_id}: {status} ({client.queries_sent} queries)"
        + (f" - {report['error']}" if "error" in report else "")
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
