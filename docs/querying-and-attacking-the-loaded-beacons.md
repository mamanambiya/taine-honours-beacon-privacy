# Querying and attacking the twelve loaded beacons

How to ask the loaded beacons a question, check that they answer correctly,
and run the membership attack on all twelve. Everything here was read or run
on 2026-10-28 against `afrigend-all-dev`; re-measure before relying on it.

The twelve datasets are `<POP>_M65_S1` for CEU ESN FIN GBR GWD IBS LWK MSL TSI
YRI and the pooled AFR and EUR. How they were built and loaded is in
[`rebuild-steps-for-taine.md`](rebuild-steps-for-taine.md); what the attack
measures is in [`../analysis/power/README.md`](../analysis/power/README.md);
the live beacon's limits are in [`beacon-constraints.md`](beacon-constraints.md).

## 1. Reach the beacon

The beacon API listens on `127.0.0.1:8086` on the VM only. From your Mac, use
the ssh alias `afrigend-all-dev` (host `192.168.101.167`, user `ubuntu`, through
the jumphost `mamana-jumphost-prod`, set in `~/.ssh/config`). Queries on this
beacon are throttled (see rebuild-steps D1 and D2 for the rate-limit settings),
and the API answers 429 when it throttles.

## 2. Ask one dataset one question

Send `datasetIds` with one dataset id. `start` is 0-based, so use the VCF
position minus 1, and `referenceName` has no `chr` prefix.

```bash
ssh afrigend-all-dev "curl -s 'http://127.0.0.1:8086/api/g_variants?assemblyId=GRCh38&referenceName=20&start=709680&referenceBases=G&alternateBases=C&datasetIds=CEU_M65_S1'" | python3 -m json.tool
```

Expect `responseSummary.exists: true` (chr20:709681 G>C is a site a CEU member
carries; it is the first line of `~/pop65_build/CEU/present_sites.tsv` on the VM).

**Read `responseSummary.exists`, not the `resultSets` list.** The response
always lists all twelve datasets, but only the dataset named in `datasetIds` is
evaluated; the rest show `exists: false`, `resultsCount: 0` because they were
not asked. Measured with the same variant:

| Query | `numTotalResults` | Datasets with `exists: true` |
| --- | --- | --- |
| no `datasetIds` | 12 | all 12 |
| `datasetIds=CEU_M65_S1` | 1 | CEU only |
| `datasetIds=ESN_M65_S1` | 1 | ESN only |

So the variant is in both CEU and ESN, and a dataset's `exists: false` in a
query that did not name it says nothing about that dataset. One query per
dataset. Without `datasetIds` the answer is over all twelve at once.

A site no member carries (`absent_sites.tsv` in each `~/pop65_build/<POP>/`)
should answer `exists: false`; run one as a negative control before trusting any
positive answer.

## 3. Check that each beacon answers correctly

`analysis/lrt/beacon_check.py` samples members' sites (must answer yes) and
sites no member carries (must answer no), one dataset per call, and fails if the
beacon ignores `datasetIds`. Run on the VM, from `~/pop65_build`:

```bash
cd ~/pop65_build && mkdir -p checks && for P in CEU ESN FIN GBR GWD IBS LWK MSL TSI YRI AFR EUR; do python3 beacon_check.py --base-url http://127.0.0.1:8086 --build-dir ~/pop65_build/$P --dataset ${P}_M65_S1 --out checks/beacon_check_$P.json > checks/beacon_check_$P.log 2>&1; echo "$P exit $? $(tail -1 checks/beacon_check_$P.log)"; done
```

Expect `exit 0` and a PASS line for each population (all twelve passed on
2026-10-27). This is the same loop as rebuild-steps Part E.

## 4. Run the attack on all twelve beacons

`analysis/power/real_beacon_attack.py` scores outsiders against the set of sites
each beacon holds, so it sends no queries for the scoring. The driver that ran
all twelve is `analysis/power/run_attack_vm.py`, copied here as it ran on the VM
(paths default to the VM's `~/pop65_build` and `~/pop65_build/attack`). Two
environment variables change that: `POP65_BUILD` (the build directory, read only)
and `POP65_ATTACK_OUT` (where the sites CSVs and the results are written, default
`<POP65_BUILD>/attack`). It imports `real_beacon_attack.py` from its own folder,
and needs a running `beacon-dev-mongodb` container.

What it does for each dataset:

1. Exports the beacon's sites to `attack/sites_<POP>.csv` with `mongoexport`
   (columns `reference_name,start,reference_bases,alternate_bases`; one row per
   variant), if the file is missing.
2. Takes the beacon's 65 members from `<POP>/members.txt`.
3. Scores two outsider sets, writing `attack/attack_<POP>_<set>.json`:
   - `controls`: the dataset's own 20 held-out controls;
   - `continent`: every member and control of every dataset on the same continent
     (AFR group ESN GWD LWK MSL YRI AFR; EUR group CEU FIN GBR IBS TSI EUR),
     minus this dataset's members. These include people who are members of other
     beacons, and the pooled sets overlap the populations' controls (see
     `data/README.md`).
4. Reads each person's heterozygous sites from `<POP>/het/<sample>.tsv`.

Run it on the VM:

```bash
cd ~/pop65_build/attack && python3 run_attack_vm.py > run.out 2>&1
```

Expect `ATTACK_ALL_DONE` at the end of `run.out` and 24 JSON files. Timing seen:
about 47 s for the 20 controls and 310 s for the 397 continent outsiders of ESN.
The last line of each dataset in `run.out` gives the query budget at 5 per cent
false positives.

### Check the scoring against the live API

`analysis/power/validate_api.py` samples 100 heterozygous sites from one outsider
and one member in each of two beacons (LWK and EUR), asks the live API with
`datasetIds`, and compares each answer with the set lookup. It sends 400 queries,
so the throttle must be raised first. On 2026-10-27 it reported 0 disagreements.

```bash
cd ~/pop65_build/attack && python3 validate_api.py
```

Expect the last line: `API vs set lookup: queries 400 disagreements 0`.

## 5. Where the results are

On the VM: `~/pop65_build/checks/` and `~/pop65_build/attack/`. On the Mac, in
`results/tables/` (git-ignored): `beacon_check_<POP>.json` for the ten
populations only (AFR and EUR were checked on the VM and are not copied here) and
`attack/`. The summary table and its limits are in
[`../analysis/power/README.md`](../analysis/power/README.md).
