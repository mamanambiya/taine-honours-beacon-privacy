# Work plan

| | |
| --- | --- |
| Project | Implementation, Evaluation and Extension of Beacon v2.0 for African data |
| Student | Taine Archbold, B.Sc. Bioinformatics Honours |
| Supervisors | Nicola Mulder, Ayton Meintjes |

## Aim

A Beacon answers one question: *do you hold this variant?* It returns yes or
no and nothing else, which is what lets institutions advertise data without
releasing records. The protocol is
[GA4GH Beacon v2](https://docs.genomebeacons.org/).

[Shringarpure and Bustamante
(2015)](https://doi.org/10.1016/j.ajhg.2015.09.010) showed an attacker holding
a target genome can determine whether that person is in a Beacon from enough
yes/no answers. Each answer leaks a little, and a likelihood-ratio test
accumulates it. Their numbers were measured on European data. The paper is
[open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC4667107/); read it
before week 2.

Measure whether the attack is easier, harder or the same against an African
cohort, and identify which property of the data accounts for any difference.

The mechanism to keep in mind: the attack works on variants that few people
carry, because a "yes" to a rare variant narrows the field sharply. African
genomes carry more rare variation, so the attack plausibly performs
differently — but nobody has measured it.

## Schedule

| Weeks | Task | Hand in |
| --- | --- | --- |
| 1–2 | Read in; run the beacon locally | — |
| 1 | Put questions A–C to us | Written answers |
| 2–4 | 1. Build the attack; check on simulated data | Power-vs-spectrum plot |
| 5–7 | 2. European arm | European power curve |
| 8–10 | 3. African arm | Matched curve and comparison |
| 11–12 | 4. Explain the difference | Figure: data features vs attack cost |
| 11–12 | 5. Check against live beacons | Validation table |
| 13–16 | Write up | Draft by week 14 |
| spare | 6. Real African panel | Only if time remains |

Start the 1000 Genomes download in week 1. Tasks 2 onwards depend on it.

## Setting up

The beacon itself is
[AfriGen-D/variant-checker-beacon](https://github.com/AfriGen-D/variant-checker-beacon)
— a Django service over MongoDB, run locally with Docker.

**Follow the tutorial:
[docs/developers/tutorial-zero-to-beacon.md](https://github.com/AfriGen-D/variant-checker-beacon/blob/main/docs/developers/tutorial-zero-to-beacon.md).**
It goes from an empty machine to a beacon answering real queries in about
fifteen minutes, and prints the output each step should produce so you can tell
when something has gone wrong. Every command on it was run end to end against
the repository, and the traps it documents — the gitignored environment file,
the production template that redirects every request into a void, the legacy
Mongo shell — each cost someone an hour to find.

Two things from it worth carrying in your head afterwards:

- The answer to a query is at `responseSummary.exists`. There is no top-level
  `exists` field, and older documents in that repository show clients parsing
  one; those report "not found" for everything.
- The beacon holds GRCh38 only, and refuses other builds with a 501 rather than
  answering "no". That refusal is doing you a favour — see task 2.

The [rest of the developer documentation](https://github.com/AfriGen-D/variant-checker-beacon/blob/main/docs/developers/README.md)
covers the architecture, the API contract and how the tests gate a merge.

## Your own repository

Your analysis lives in `taine-honours-beacon-privacy`:

```text
analysis/simulation/     the attack, the cohort simulator, their tests
analysis/lrt/            attack runs against real cohorts
analysis/power/          power curves and the arm comparison
data/                    cohort data — gitignored, record provenance in its README
results/figures/         plots
results/tables/          numeric output, as JSON or CSV
docs/                    method notes and constraints
```

Python 3 with [numpy](https://numpy.org/doc/stable/); no other dependency is
needed for the simulation. You will want
[bcftools](https://samtools.github.io/bcftools/bcftools.html) for the VCF work
in tasks 2 and 3. Tests use `unittest` from the standard library, matching the
beacon repository:

```bash
python3 -m unittest discover -s analysis/simulation -t analysis/simulation
```

Write the tests as you go, not afterwards. For this attack in particular, a
test that has never been observed failing is not evidence the code works — the
sign error described in task 1 passes every test that does not specifically
check for it.

## Questions for us — week 1

**A. Beacon size.** The live panel holds 1,895; the proposal specifies 500.
Beacon size alone changes attack difficulty, so unmatched arms measure size
rather than ancestry. Subsample both to one size, or model size in the
statistics?

**B. Live querying.** The beacon allows 50 queries/hour per machine.
Experiments run against your own copy; the live one is for spot checks. Arrange
with the AfriGen-D operators, or request an exemption?

**C. Ethics.** A re-identification attack against a service the group runs,
with logged and attributable queries. HREC application, or written operator
agreement?

## Task 1 — Build the attack, check it on simulated data

1. Implement the likelihood-ratio test in `analysis/simulation/lrt.py`.
   `docs/methods-lrt.md` sets it out; reconcile it against the paper as you go.
   You need `D_N` — the chance nobody in the beacon carries a given allele —
   as an expectation over the frequency spectrum, not a per-site constant. It
   has a closed form; compute it in log space with `math.lgamma`, because the
   beta function underflows once `N` reaches the thousands.
2. Simulate a cohort in `analysis/simulation/cohort.py` — a genotype matrix of
   individuals by variants, each entry 0, 1 or 2 alternate copies. Draw a
   frequency per variant from `numpy.random.Generator.beta`, then genotypes
   with `binomial(2, f)`. Take an explicit seed and record it with any output.
3. Run the attack. Query only the target's heterozygous positions — that
   ascertainment is part of the model, not a shortcut. The beacon answers yes
   if any member carries the allele. Feed the answers into the statistic and
   compare member against non-member scores.
4. Sweep the frequency spectrum in `analysis/simulation/contrast.py`, from
   common-dominated to rare-dominated, and measure how attack difficulty
   changes. Write the numbers to `results/tables/` rather than reading them off
   a terminal.

The *frequency spectrum* is the distribution of how common variants are across
a cohort — how much of the variation is rare versus shared. A `Beta(a, b)`
distribution is a convenient stand-in: small `a` and large `b` puts most of the
mass at low frequency, which is a rare-variant-rich cohort.

Step 4 is the precondition for the rest. The frequency spectrum is the only
route by which a cohort enters the model; if difficulty does not move with it,
no ancestry effect is possible.

Two requirements that carry into every later task:

- Report the **query count needed** for a given success rate, not success at a
  fixed count. Real cohorts saturate, and a fixed-budget comparison returns no
  difference even when costs differ by an order of magnitude.
- **Assert the sign of the statistic** against a case known by construction. A
  low value means the target is in the beacon. Reversed, the power curve
  inverts and still looks plausible.

Three weeks.

## Task 2 — European arm

1. Download the 1000 Genomes **GRCh38 NYGC 30x** call set — the
   [30x GRCh38 collection](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38),
   served from
   [the EBI FTP site](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/).
2. Extract the European samples into the genotype-matrix form from task 1.
3. Split into beacon members and controls, with a seeded split
   (`analysis/splits/make_split.py`). **Cut the members' VCF so it keeps only
   the sites that vary in the members**: `bcftools view -S members.txt -c1`,
   with allele counts recomputed for those people. A VCF cut to a subset of
   samples otherwise keeps every site of the parent panel, and a beacon built
   from it answers yes at sites no member carries, so the attack has nothing to
   measure. This was missed on the first load of the population beacons
   (see `docs/reload-population-beacons.md`). The controls stay out of the
   beacon.
4. Find the query budget identifying 95% of members at a 5% false-positive
   rate. That pair of numbers is a *power curve* point: queries on one axis,
   the fraction of members correctly identified on the other, with the
   false-accusation rate held fixed.
5. Repeat with singleton censoring enabled, to cost that defence. A censoring
   beacon refuses to answer for variants only one person carries — the
   positions that most nearly name someone. It is the defence the 2015 paper
   tested, so measuring its cost puts your numbers alongside theirs.

**Do not use Phase 1 or Phase 3 as published.** They are GRCh37; the beacon
accepts GRCh38 only. A build mismatch shifts coordinates, the beacon answers
"no" to variants it holds, and the result is indistinguishable from the effect
being measured.

Three weeks. Most of it data preparation.

## Task 3 — African arm

1. Apply the answer to question A before generating any curves.
2. Build the African beacon from the 1000 Genomes African samples, the same
   way as the European one: members only, sites that vary in the members only.
   Before any curve, run `analysis/lrt/beacon_check.py`. It asks the beacon
   for sites no member carries, which must answer no, and for each member's
   own heterozygous sites, which must all answer yes. A beacon that answers
   yes everywhere passes every positive check and fails this one.
3. Run the identical attack — same code, settings and false-positive rate.
4. Compare the query counts.

Hand in both curves on one figure and a statement of whether the attack is
easier against the African cohort, and by how much.

Three weeks.

## Task 4 — Explain the difference

Compute three properties for both cohorts and plot each against the measured
attack cost, to see which one tracks it:

- **Rare-variant load** — the share of variants below some frequency, say 1%.
  The attack's fuel.
- **Linkage disequilibrium** — how strongly nearby variants are inherited
  together. High linkage means answers are correlated, so each one adds less
  than it appears to.
- **F_ST** — how genetically differentiated the populations are. Sets how far
  apart the two spectra are in the first place.

Hand in the figure and a statement of the defence implication: whether a
minimum-frequency threshold suffices, or whether risk modelling must be
population-specific.

Two weeks.

## Task 5 — Check against live beacons

Spot checks only. The rate limit precludes an attack run.

1. Obtain the operators' agreement. Queries are logged.
2. Select positions with known answers.
3. Query your local beacon and the live one
   ([beacon.afrigen-d.org](https://beacon.afrigen-d.org/)); confirm agreement.
4. Repeat against [beacon.ardi.africa](https://beacon.ardi.africa/), a second
   African beacon run by a different group.

- Run a positive control first. Without it, "variant absent" and "script
  broken" are indistinguishable.
- Convert coordinates: VCF counts from 1, the beacon from 0. An off-by-one
  returns "no" throughout, which reads as the African cohort being safer.
- Back off on rate-limit responses. Retry loops have taken the service down.

Hand in a table: position, local answer, live answer, agreement. One week.

## Task 6 — Real African panel

Only if time remains. Measures the true frequency spectrum of the panel behind
the live beacon, making the absolute query counts specific to it.

Requires an ILIFU **compute-node** session — not the login node; see the
[ILIFU documentation](https://docs.ilifu.ac.za/). Ask us for the panel
location.

Not obtainable from the beacon: it returns no per-variant detail, so there is
no frequency to read. Nothing else depends on this task.

## Write-up

Draft to us by week 14.

- Literature review: record the search strategy. The novelty claim is a
  negative one and needs evidence of how thoroughly you looked.
- Methods: what was actually run, including every deviation from the proposal
  and its cause.
- Ethics statement covering the live beacon test.

## Contingencies

| If | Then |
| --- | --- |
| 1000G download too slow or large | Use one or two chromosomes; state it in the methods |
| Ethics approval delayed | Only task 5 needs it; continue with 2, 3 and 4 |
| No difference between arms | A result, not a failure — but first confirm the beacon sizes matched and the coordinates were correct |
| Time runs short | Drop task 6, then task 5. Tasks 1–4 are the dissertation |

## Meetings

Regular slot from week 1. Bring questions A–C to the first. Bring the European
curve before starting task 3.

## Reference

| File | Contents |
| --- | --- |
| `docs/methods-lrt.md` | The statistic, written out, with the traps flagged |
| `docs/beacon-constraints.md` | Rate limits, genome build, coordinates |
| `docs/background.md` | Reading list by purpose |
| `proposal/` | Your proposal |

External:

| Link | What it is |
| --- | --- |
| [doi:10.1016/j.ajhg.2015.09.010](https://doi.org/10.1016/j.ajhg.2015.09.010) | Shringarpure & Bustamante 2015 — the attack |
| [PMC4667107](https://pmc.ncbi.nlm.nih.gov/articles/PMC4667107/) | The same paper, open access |
| [docs.genomebeacons.org](https://docs.genomebeacons.org/) | GA4GH Beacon v2 specification |
| [30x GRCh38 collection](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) | The 1000 Genomes release to use |
| [EBI FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/) | Where that release is served from |
| [variant-checker-beacon](https://github.com/AfriGen-D/variant-checker-beacon) | The beacon implementation |
| [beacon.afrigen-d.org](https://beacon.afrigen-d.org/) | The live AfriGen-D beacon |
| [beacon.ardi.africa](https://beacon.ardi.africa/) | A second African beacon |
| [bcftools](https://samtools.github.io/bcftools/bcftools.html) | VCF manipulation |
| [numpy](https://numpy.org/doc/stable/) | Arrays and random sampling |
| [docs.ilifu.ac.za](https://docs.ilifu.ac.za/) | The cluster, for task 6 |
