# Plan: rebuild the 1000G population beacons so the attack can be run on them

Reading copy: https://claude.ai/artifact/4cXTpWSM432M8AJddyttQi (this file is
the source of truth).

Taine honours repo (`taine-honours-beacon-privacy`), session 04678b, 2026-10-21.
Nothing here has been run; the VM is untouched since the 09:00 audit (containers
stopped, as found).

## Context

The ten 1000G population beacons on `afrigend-all-dev` list every site in the
panel (6,240,864 records each), so "exists" means "in the panel", not "a member
carries it". 3 of 3 YRI sites with allele frequency 0 returned yes. The attack,
and any African-versus-European comparison, would measure nothing. Evidence and
the full audit: shared notes
`afrigend-all-dev-1000g-beacon-audit-2026-10-21.md`.

**Root cause in the process (read from afrigen-beacon-v2 @ c5283d2).** The rule
"a record exists iff a member carries the ALT allele" is enforced nowhere.
`vcf_to_beacon.py` reads no genotypes and writes a record for every ALT allele
passing QUAL/DP. `exists` is "any record matches". The ILIFU guide teaches
region slices only; the right recipe (`bcftools view -S samples.txt -c1`) sits
in `examples/README.md` and is not taught. Guide and tutorial checks are
positive-only. Our own `worklist-for-taine.md` task 2 has the same hole.

**Outcome wanted.** Per population, one beacon dataset of 40 members built only
from sites that vary in those 40, with 40 held-out controls whose genotypes stay
out of the beacon, verified against genotype ground truth before any attack run.

## Design (recommended defaults, all reversible)

| Choice | Default | Why |
| --- | --- | --- |
| Members / controls | 40 / 40 per population, seeded split | Matches the simulation's N=40; MSL (85) is the smallest, so every population can supply 80. Resolves the N confound inside the 1000G arms. |
| Replicates | 1 split now (`S1`); more seeds later | Keeps the first load small; add seeds if variance matters. |
| Dataset ids | `<POP>_M40_S1` | Old ids untouched, nothing to confuse. |
| Where | same Mongo (`beacon_db`) on `afrigend-all-dev`, new ids | New data ~2 GB (estimate from 11.4 GB / 62.4M docs); VM has 17 GB free. |
| Old datasets | kept until every gate passes, deleted only on a separate approval | Rollback stays trivial. |
| Scope | chr 20-22, SNPs, het sites 0/1 only | What Taine's outputs already cover. Stated as a limit. |
| Serving | API stays on loopback; Taine uses `ssh -L`; rate limit raised on this instance only | A raised limit must not be reachable from the public hostname. The dead 8080/8086 tunnel mismatch is left alone and reported. |
| Attack runs | answers computed from genotypes (`cohort.beacon_answers`), beacon used as the verified oracle on a sample | Avoids ~10^6 HTTP queries; the beacon's job is to prove it equals the genotype truth (this is the "local instance identical to truth" check in Milestone 4). |

## Phases

**Phase 0 - provenance and recipe (read-only, on ILIFU; you or Taine, MFA).** On
`/cbio/users/arctai001/vcf_samples/snps_<POP>.vcf.gz`: read the header
(`##source`, `##reference`, `##bcftools_viewCommand`), sample count, `AN` of a
few records. Resolves two open facts: is the source the IGSR GRCh38 re-mapping
(the example's file) or the NYGC 30x set the project chose (decision 3); and why
a YRI record carries `an=42` for 108 samples. **Gate G0: if the source is a
lift-over, stop and report before building anything.**

**Phase 1 - build on ILIFU (per population, 10 times; compute node, `sbatch`).**

1. Seeded split of `bcftools query -l` into `members.txt` (40) and
   `controls.txt` (40); the lists are public sample IDs and are committed to the
   repo.
2. `bcftools view -S members.txt -c1 -Oz` (sites that vary in the members; INFO
   AC/AN recomputed), index.
3. Run `sfs_histogram.sh` (PR #1) on that VCF with N=40. **Gate G1: "with AC =
   0" is 0 and its kept count equals `bcftools view -H | wc -l` of alleles.**
   The same table is the matched 1000G spectrum Milestone 1b asks for.
4. Beacon pipeline, `--skip-import`, with `dataset.yaml` (`id: <POP>_M40_S1`,
   description naming members, seed, source file, build) per the ILIFU guide
   section 4.4.
5. Het-position tables per member and per control (chr, POS, REF, ALT; `GT=0|1
   or 1|0`, SNP, biallelic) from the full-population VCF. About 800 small files
   in total.
6. Verification site lists: 100 sites present in the member VCF and 100 sites in
   the full VCF but absent from it (`bcftools isec -C`).
**Gate G2: per-population record counts differ, are far below 6,240,864, and AFR
counts exceed EUR.**

**Phase 2 - move and load (me, over ssh to `afrigend-all-dev`, after this plan
is approved).** Copy outputs via your Mac (guide Route A; `scp` to
`transfer.ilifu.ac.za`, MFA each time). Start the three stopped containers,
`mongoimport` each dataset, `register_dataset.py --samples 40`, flush Redis
(results are cached 5 minutes). Add `RATELIMIT_QUERY_ENDPOINT` to the instance
env and recreate only the API container, with the compose file backed up first.
**Gate G3: per-dataset `countDocuments` equals the Phase 1 count; no dataset
loaded twice (the guide's double-load warning); `/api/datasets` shows 10 new ids
with `sampleCount=40`.**

**Phase 3 - verify against truth (me, new harness in this repo).**

- G4: the 100 + 100 verification sites per population, scoped with `datasetIds`:
  all 100 present answer yes, all 100 absent answer no. Zero errors allowed.
- G5: every member's het sites answer yes (100%, any miss is a load or
  coordinate bug); controls' yes-rate is below 100% and recorded.
- G6: a positive control (known present, 0-based `start = POS-1`) opens every
  batch; a 429 is backed off, never retried in a loop.
**All of G0 to G6 must pass before any power number is quoted.**

**Phase 4 (not in this approval).** Run the M2/M3 power curves from genotype
answers, fitting `(a', b')` per population with `sfs_fit.py`; report
queries-to-95%-power, AFR versus EUR, matched at N=40.

## Files (Taine repo, one feature branch and PR, tests first)

- `analysis/splits/make_split.py`, with `test_make_split.py`: deterministic,
  disjoint, sizes exact, seed recorded; edge cases: population smaller than 80,
  duplicate IDs, empty list.
- `analysis/splits/<POP>_members.txt`, `<POP>_controls.txt` (public IDs).
- `analysis/lrt/build_population_dataset.sh`: Phase 1 steps 1-6 in one
  reproducible script, same `BCFTOOLS` override pattern as `sfs_histogram.sh`.
- `analysis/lrt/beacon_check.py` with a mock-HTTP test: Phase 3 gates, positive
  control first, 429 back-off, reads the het tables. Reuses `lrt` and `cohort`
  for expected rates.
- `docs/reload-population-beacons.md`: the runbook, plus the "what we missed"
  analysis. Plan copied to `docs/superpowers/plans/` and published as an
  artifact (your rule).
- Reused unchanged: `sfs_histogram.sh`, `sfs_fit.py`, `cohort.beacon_answers`,
  the beacon repo's `run_ilifu_pipeline.sh` and `register_dataset.py`.

## Risks and what I have not verified

- **Source provenance (G0)** and the `an=42` oddity are open; the whole design
  assumes the source VCFs are sound once subset.
- **The population VCFs' own recipe is not recorded anywhere I can read**; "cut
  without `-c1`" is inferred from the counts.
- **`bcftools +fill-tags` inside the tools container** is unverified (it
  appeared in the PR #1 runbook).
- **Timings and disk are estimates**; ILIFU array speed for these files is
  unmeasured.
- **chr 20-22 only** limits LD and rare-variant counts; the closed-form scaling
  and the V6HC fit cover the rest.
- **I touch `afrigend-all-dev` only in Phases 2 and 3**, and the ILIFU steps
  need you or Taine (MFA).

## Rollback (not part of the sequence)

New datasets: `deleteMany`/`deleteOne` by id (guide section 10.5), compose file
restored from its `.bak`, containers stopped. Old datasets stay untouched until
a separate approval. Nothing on ILIFU is shared.

## Upstream fixes (proposals only; not our repo, not edited here)

For the beacon repo owners: a guide section "a dataset from a subset of samples"
with `-S` and `-c1` and a negative-control step; a transformer warning when more
than a few percent of records have `AC=0` (it already reads INFO `AC`); a
dataset check comparing counts across datasets; a tutorial step "an in-panel
site your samples do not carry must return false".

## Verification of this work, end to end

G0 to G6 above, run in order, each with its output kept in `results/tables/`.
The reload is done when all pass on all 10 populations, and the old datasets are
deleted only after that, on your say.

## Refinements made while implementing (2026-10-21)

- **Biallelic SNPs only** in the member VCF (`-m2 -M2 -v snps`). The plan said
  "SNPs"; multi-allelic sites are also dropped because the beacon writes one
  record per ALT allele, so an ALT with no carrier could survive at a site
  another ALT varies at. The scope limit is stated in the build script and the
  runbook.
- **`build_population_dataset.sh` copies `members.txt` and `controls.txt` into
  its output**, so a build directory is self-contained for the check.
- **Gate G3's "Phase 1 count" is the beacon transform's `variants_processed`**,
  not the VCF site count, because the transform applies its own QUAL and DP
  filters.
- **Status:** the three code pieces and their tests are done and
  mutation-checked; Phase 0 (ILIFU, read-only) is the open gate; Phases 1 to 3
  have not been run.
