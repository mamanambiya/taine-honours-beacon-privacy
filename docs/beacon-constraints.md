# Operational constraints on the live AfriGen-D Beacon

Everything below was measured against the running deployment or read from the
deployed source on **2026-08-21**. Re-measure before relying on any of it — the
beacon is under active development.

These constraints are not incidental. Two of them change the method as written
in the proposal, and one of them is a confound that could invert the study's
conclusion.

## 1. The rate limit makes the live-beacon attack infeasible as written

`beacon_project/settings_boolean.py:238-241` sets:

```python
BEACON_RATE_LIMITS = {
    'query':       config('RATELIMIT_QUERY_ENDPOINT', default='50/hour'),
    'variants':    config('RATELIMIT_QUERY_ENDPOINT', default='50/hour'),
    'individuals': config('RATELIMIT_QUERY_ENDPOINT', default='50/hour'),
}
```

**50 queries per hour, per IP.**

Shringarpure & Bustamante report membership inference at 95% power with ~250
SNPs, and their real-data evaluation queried whole genomes. The proposal's
design follows theirs: 20 genomes in the beacon and 20 absent, per cohort.

At 50 queries/hour, the paper's headline figure of ~250 queries for one target
genome takes 5 hours, and forty targets take over a week of continuous querying.

**That arithmetic understates the problem by about three orders of magnitude,
because ~250 is not the right query count for this beacon.** Measured
2026-08-21 against `analysis/simulation/lrt.py`, at the live panel size
`N = 1895`, `delta = 1e-6`, `alpha = 0.05`, 95% power, `f ~ Beta(0.5, 2.0)`:

| | Queries for one genome | At 50/hour |
| --- | --- | --- |
| Unconditioned SFS (wrong) | 123 | 2.5 hours |
| **Ascertainment-corrected** | **108,898** | **91 days** |

The correction is the `a' + 1, b' + 1` heterozygote weighting in
[`methods-lrt.md`](methods-lrt.md) — query sites are heterozygous in the target,
so their allele frequencies are enriched for intermediate values, `D_N` falls,
and the two hypotheses move much closer together. It is worth an **885x**
increase in the required number of queries.

So a single target genome is **roughly three months of uninterrupted querying**
at the production limit, not five hours. Forty targets is a decade. The live
attack is not slow; it is infeasible by a margin that no reasonable rate-limit
exemption closes.

Two things follow. First, Milestone 4 is spot validation only, and this is now
a measured statement rather than a cautious one. Second, **the rate limit is
doing more privacy work than anyone had quantified** — that is a finding for the
write-up, not merely a constraint on the method.

A35 is now **verified**, so this number is quotable. The radicals were never
missing from the paper — only from PMC's HTML-to-text rendering. Europe PMC
serves the JATS full text with MathML, where they are structurally present:

```bash
curl -sL "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4667107/fullTextXML"
```

See [`methods-lrt.md`](methods-lrt.md) for the transcription and for a sign in
A35 that looks wrong and is not.

**One limit survives A35, and it is about inputs rather than the formula.**
Every absolute figure above assumes a `Beta(a', b')` site-frequency spectrum
with `a' = 0.5, b' = 2.0` — a plausible shape, but **not fitted to
V6HC-S_AFR**. `D_N` is an expectation over that spectrum, and it is the only
route by which the cohort enters the model, so the absolute query count moves
with it.

The **885x ratio is robust** to this: it divides two numbers computed from the
same spectrum, so the assumption cancels. The absolute **108,898 is an
order-of-magnitude figure** until the SFS is fitted to the real panel. That
fitting is objective 5 of the proposal, so the study generates its own
correction — and the ordering matters, because objective 5 is currently
scheduled last while three earlier conclusions already depend on it.

The study cannot run against production as specified. The realistic options,
in the order they should be considered:

1. **Run the attack against a local instance.** The repo ships a full local
   stack (`compose/docker-compose.dev.yml` in `../afrigen-beacon-v2/`), which
   can be loaded with the same V6HC-S_AFR data and has the rate limit under
   your control. This is the correct home for the *power curve*, because a
   power curve needs thousands of queries and none of them are interesting to
   an operator.
2. **Use production only for validation at a handful of points** — enough to
   confirm the local instance answers identically, not enough to build a curve.
3. **Ask for a rate-limit exemption** for a named source IP and a fixed window.
   That is a decision for the supervisors and the AfriGen-D operators, not
   something to arrange ad hoc.

## 2. A naive attack loop can take the beacon offline for everyone

The rate limiter has a history of *latching*: a client that keeps polling after
being limited had its window repeatedly extended, locking it out for a full
hour rather than the intended window. That fault took the Beacon Network
aggregator offline in August 2026. It was fixed on the AfriGen-D side, but the
third-party ARDI beacon has not picked the fix up.

Two consequences:

- **Back off on 429. Do not retry-in-a-tight-loop.** An attack script that
  hammers through rate-limit responses is indistinguishable from the traffic
  pattern that caused a real outage.
- **The federation aggregator is the fragile one.** Querying
  `beacon-network-dev.afrigen-d.dev` fans out to every registered beacon and is
  bounded by the slowest one. Prefer querying member beacons directly.

## 3. GRCh38 only — GRCh37 is refused

Measured live:

```console
$ B=https://beacon.afrigen-d.org/api/g_variants
$ Q='assemblyId=GRCh37&referenceName=1&start=10000&end=20000'
$ curl -o /dev/null -w '%{http_code}\n' "$B?$Q"
501
```

The beacon answers `501 Not Implemented` for any build it does not hold. It
holds GRCh38 only (`beacon_api/capabilities.py`), and the refusal is
data-driven — load a GRCh37 dataset and it starts answering, with no code
change.

**This constrains the 1000G replication.** The Shringarpure & Bustamante study
used 1000G Phase 1, which is GRCh37. Coordinates cannot be mixed. Use the
GRCh38-aligned 1000G release (the NYGC high-coverage 30x call set) so that
both arms of the comparison are on one build, or lift over and document it.

The good news is it fails loudly, with a 501, rather than silently returning
"no".

## 4. Coordinates are 0-based half-open — this is the confound to watch

From `beacon_api/query_semantics.py`:

> Variants are stored in Beacon v2 **0-based half-open** coordinates
> (`start = POS - 1`), so a stored variant occupies `[start, start + 1)`.

A VCF gives you 1-based `POS`. Forwarding `POS` straight into `start` queries
one base to the right of the variant, and the beacon answers **no**.

This matters more here than in ordinary use. A systematic off-by-one produces a
uniform excess of "no" responses, which the LRT reads as *evidence the target is
not in the beacon*, which reads as *lower attack power*, which reads as *the
African cohort is less vulnerable*. **The bug and the hypothesis produce the
same result.** There is a regression suite for exactly this bug at
`../afrigen-beacon-v2/beacon_api/test_query_semantics.py`.

Guard against it with a positive control before every run: query a locus you
know is present and confirm `exists: true`. A negative result you have not
controlled is not a measurement.

## 5. The beacon publishes N, so the attacker model's assumption holds for free

`GET /api/datasets` returns, unauthenticated:

```json
{"id": "H3A_V6_AFR",
 "name": "H3A V6 African Reference Panel",
 "assemblyId": "GRCh38",
 "variantCount": 42331932,
 "sampleCount": 1895}
```

The attack model assumes the adversary knows the beacon size `N`, which is
usually argued as a reasonable assumption. Here it is not an assumption at all
— `N = 1895` is published. Worth stating explicitly in the write-up, because it
removes a caveat that the original paper had to carry.

Note the **size mismatch** this creates with the proposal's design: an
experimental 1000G beacon of size 500 against a live beacon of size 1,895
confounds ancestry with `N`, and attack power is strongly `N`-dependent through
`D_N`. Either match the sizes by subsampling the African panel, or model `N`
explicitly as a covariate. This is a real methodological decision, not a detail.

## 6. Query granularity is boolean — which is what the attack wants

Live response shape:

```json
{"meta": {"returnedGranularity": "boolean", ...},
 "responseSummary": {"exists": true, "numTotalResults": 1}}
```

`responseSummary.exists` is the yes/no oracle the LRT consumes. No granularity
negotiation, no authentication, no handover needed.

## 7. Every query against production is logged

Query logging writes to the `query_logs` collection and surfaces on the
AfriGen-D Grafana dashboard. It records the query parameters, the response, and
a truncated client address, with a TTL.

An attack run against production is therefore **visible and auditable**. That is
appropriate — but it means the run must be agreed with the operators and the
supervisors in advance, not discovered in the dashboard afterwards. Treat it as
a disclosed penetration test of a system your own group runs.

## 8. Allele frequencies: compute them from the panel, not from the beacon

The boolean response carries no allele frequency, and the AF field on the
stored variants is not populated for the full production collection.

For objective 5 (characterising the allele-frequency spectrum), compute the AFS
directly from the source VCFs on ILIFU rather than through the beacon:

```text
/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/
    V6HC-S_AFR/vcfs/V6HC-S_AFR_chr{1..22}_all.vcf.gz
```

This is the same data the beacon was loaded from, so the SFS you derive is the
SFS the beacon's `D_N` is governed by. Never SSH to the ILIFU login node — use a
compute node.

## Summary of what this changes in the proposal

| Proposal says | Reality | Consequence |
| ------------- | ----------------- | ---------------------------------- |
| Attack the live beacons | 50 queries/hour/IP | Power curves must come from a local instance; production for spot validation only |
| Replicate 1000G attack | Beacon is GRCh38-only, 501s on GRCh37 | Use the GRCh38 1000G call set, not Phase 1/3 as published |
| Matched sample size 500 | Live beacon is N=1,895 | Subsample, or model N as a covariate — otherwise ancestry and N are confounded |
| Query with VCF positions | Storage is 0-based half-open | Convert `POS-1`; a positive control per run is mandatory |
