# Data

**Nothing in this directory is committed.** `.gitignore` excludes `data/**`
except this file. Record provenance here; keep the bytes out of git.

## Cohorts needed

| Arm | Source | Build | Notes |
| --- | --- | --- | --- |
| African (live beacon) | H3A V6 African Reference Panel, `V6HC-S_AFR` | GRCh38 | 1,895 samples, 42,331,932 variants — the dataset behind `beacon.afrigen-d.org` |
| African (experimental) | 1000 Genomes AFR | GRCh38 | For the size-matched replication |
| European (experimental) | 1000 Genomes EUR / CEU | GRCh38 | The arm Shringarpure & Bustamante measured |

## V6HC-S_AFR on ILIFU

```text
/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/
    V6HC-S_AFR/vcfs/V6HC-S_AFR_chr{1..22}_all.vcf.gz
```

Autosomes only, GRCh38. This is the same data the production beacon was loaded
from, so an allele-frequency spectrum computed here is the spectrum that governs
the beacon's `D_N`.

**Never SSH to the ILIFU login node.** Work on a compute node — see the
`hpc-clusters` guidance for how to get an interactive session.

## 1000 Genomes build

The beacon answers `501` for any build other than GRCh38, so both arms must be
on GRCh38. The published Phase 1 data that Shringarpure & Bustamante used is
GRCh37, and Phase 3 as originally released is GRCh37 too.

Use the GRCh38-aligned high-coverage (30x) call set, or lift over and document
the liftover — but do not mix builds between arms, because a coordinate
mismatch produces exactly the signal this study is trying to measure.

## A population cut from a larger panel

A VCF cut to a subset of samples keeps every site of the parent panel, including
the ones no sample in the subset carries. Load that into a beacon and it answers
yes at those sites: the dataset then says "this site is in the panel", not "a
member carries this allele", and attack power measures nothing.

Cut the members' VCF with `bcftools view -S members.txt -c1` (restricted to
biallelic SNPs, with `AN`, `AC` and `AF` recomputed for the members), and check
the result: no record with `AC = 0`, and one allele number, `2N`, at every site.
`analysis/lrt/build_population_dataset.sh` does this and refuses to continue
if it cannot show it.

## Recording what you used

For every dataset actually analysed, note here: the source URL or path, the
release/version string, the date retrieved, the sample count after any
filtering, and the checksum if one is published. A power comparison between
cohorts is only interpretable if both arms' provenance is pinned.

## Analysed datasets

### Loaded beacon datasets on afrigend-all-dev (2026-10-27)

| | |
| --- | --- |
| Datasets | `CEU FIN GBR IBS TSI ESN GWD LWK MSL YRI AFR EUR` as `<POP>_M65_S1`, each 65 members (controls: 20 held out, not loaded) |
| Source | NYGC 30x GRCh38 phased panel, chr 20 to 22 (`/cbio/dbs/1kg/b38/nygc_30x/1kGP_high_coverage_Illumina.chr{20,21,22}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`), 1,004 samples (the AFR and EUR lists from the earlier build), 3,713,694 records in the cut |
| Cleaning | biallelic SNPs, FILTER PASS or missing, no missing genotype in any member or control: 3,232,835 sites; allele number 130 |
| Pooled sets | `AFR` and `EUR` drawn only from samples that are not members of a population dataset (beacon issue #130) |
| Variants | 9,859,803 documents = sum of the 12 `member_sites` (CEU 610,835 to AFR 1,077,584) |
| Check | `analysis/lrt/beacon_check.py`: PASS on all 12, 4,451 queries each (known present 100/100 yes, known absent 100/100 no, unknown dataset id answers no) |
| Replaces | the 20-member `<POP>_M20_S1` datasets, deleted from the same beacon on 2026-10-27 (backup `~ubuntu/backup-m20-2026-10-27/` on the VM) |

Not checked: the earlier `afr_M20_S1` and `eur_M20_S1` builds used a
different pooled sample draw and were never loaded. The beacon is queried only
through `localhost:8086` from inside the VM (no public URL is needed or
exposed), and the query limits are raised to 1,000,000 per hour on that
loopback-only instance.

### D_N and required queries at the beacon's real N = 65 (2026-10-27)

`D_N` taken straight from each cohort's allele-count histogram
(`analysis/simulation/empirical_dn.py`), for a beacon of n people drawn from
the cohort, target heterozygous, `delta = 1e-6`, `alpha = 0.05`, 95% power.
The 1000G histograms are of the 501 AFR and 503 EUR samples the datasets were
drawn from, over the same three chromosomes. They are made by
`sfs_histogram.sh` on the cut source VCF with
`VIEW_ARGS="-f PASS,. -S <cohort sample list>"` (FILTER PASS or missing,
complete genotypes, SNPs): 1,720,042 AFR and 1,071,532 EUR ALT alleles with a
nonzero count, which with the monomorphic sites (1,512,793 and 2,161,303) make
the 3,232,835-site universe each dataset was cleaned to. The tables are not
committed (`data/` is git-ignored); the 2026-10-27 copies are in
`/scratch3/users/mamana/pop65_checks/hist_AFR.tsv` and `hist_EUR.tsv` on ILIFU,
and the sample lists are `analysis/splits/samples/{afr,eur}_samples.txt`.

| Cohort | n | `D_N` | Queries for 95% power |
| --- | --- | --- | --- |
| 1000G AFR (501) | 20 | 0.0753 | 34 |
| 1000G AFR (501) | 65 | 0.0280 | 95 |
| 1000G EUR (503) | 20 | 0.0490 | 53 |
| 1000G EUR (503) | 65 | 0.0211 | 126 |
| V6HC-S_AFR (1,895) | 20 | 0.0824 | 31 |
| V6HC-S_AFR (1,895) | 65 | 0.0338 | 78 |

At n = 65, the African cohort needs 95 queries and the European one 126: about
25 per cent fewer, the direction the study's hypothesis predicts. This is a
point estimate on one cohort pair over three chromosomes, without an interval,
and it is not yet a test of that hypothesis. The V6HC-S_AFR row (whole
genome, 1,895 people) reproduces the earlier figure and is the control for the
calculation. The model assumes Hardy-Weinberg pairing, independent sites and
unrelated people; whether the 501 and 503 are all unrelated was not checked.
The attack-simulation power (as opposed to this closed-form budget) at n = 65
has not been run on the loaded datasets.

### Relatedness and sample overlap in the loaded datasets (2026-10-27)

**Relatedness.** The `D_N` budgets above assume unrelated people, and the NYGC
metadata cannot show it (every `FamilyID` equals the sample id). So kinship was
estimated from genotypes: `plink 1.9 --genome` within each cohort on common
biallelic SNPs, complete genotypes, every 10th site then LD-pruned (15,621 AFR
and 6,793 EUR SNPs). Pairs with `PI_HAT` of 0.2 or more:

| Pair | `PI_HAT` | Dataset roles |
| --- | --- | --- |
| NA19331, NA19334 | 0.61 (Z1 0.51, Z2 0.35) | NA19334 is an LWK member; NA19331 was an LWK **control** (since replaced) and is a pooled AFR member |
| HG03484, HG03464 | 0.31 | both MSL members |
| HG00116, HG00120 | 0.27 | both GBR members |

That is one first-degree-like pair in 501 AFR samples and none in 503 EUR; 56
AFR and 544 EUR pairs sit at 0.1 to 0.2. The 0.1 to 0.2 band is not evidence of
close kin: the cohorts pool five sub-populations, which inflates `PI_HAT`, and
only a few thousand SNPs were used. Effect on the figures above: negligible
for `D_N` (three close pairs among about 1,000 people). It matters for one
control: NA19331 was a control for the LWK dataset but a close relative of an
LWK member, so a query for that person was not an independent non-member.
**Done:** NA19331 was replaced by NA19445 in `analysis/splits/LWK_controls.txt`
(a seeded draw among the 14 LWK samples outside the split with no `PI_HAT` of
0.1 or more against any LWK member or control; `LWK_split.json` records the
swap). Controls are never loaded, so the loaded LWK beacon is unchanged. The
cleaned universe is still 3,232,835 sites, so the other 84 het tables stay
valid; the new control has 151,376 heterozygous sites (NA19331 had 154,147).
The pooled AFR dataset still has NA19331 as a *member*; its relative NA19334 is
an LWK member, not an AFR one, so no related pair sits inside the pooled
beacon. Not done: excluding the two within-dataset pairs (HG03484 and HG03464,
both MSL members; HG00116 and HG00120, both GBR members) and rebuilding MSL, GBR
and the pooled sets (about 3 hours of builds plus reload); the effect on the
budget is judged negligible.

**Pooled datasets overlap the populations' controls.** The pooled `AFR` and
`EUR` members were drawn to be disjoint from every population *member* (so the
importer cannot collide on a sample id), but not from the populations'
*controls*, which are not loaded: 38 of the 65 AFR members and 41 of the 65
EUR members are controls of a population dataset. Each dataset's controls are
checked only against its own beacon, so nothing is wrong inside one beacon;
do not use a population's controls as non-members of the pooled beacon, or the
reverse, without removing the overlap first. Pooled controls are disjoint from
all population members (0), and 11 AFR and 12 EUR pooled controls are also
population controls.

### V6HC-S_AFR allele-count histogram (2026-10-23)

| | |
| --- | --- |
| Source | `/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/V6HC-S_AFR/vcfs/V6HC-S_AFR_chr{1..22}_all.vcf.gz`, GRCh38, autosomes |
| Release | `version_6` (the path's own label; no other release string was read) |
| Samples | 1,895 (`N`; allele number `2N = 3790`) |
| Counted | alternate alleles of biallelic and multi-allelic **SNPs** (`VTYPES=snps`), one entry per ALT allele, counts recomputed from genotypes with `bcftools +fill-tags`; no `-f PASS` filter |
| Result | 42,331,932 ALT alleles in total, 83,379 table rows across the 22 chromosomes (counts repeat per chromosome; `sfs_fit.py` sums them) |
| Tool | `analysis/simulation/sfs_histogram.sh` at `d5920f9`, `bcftools` from `beacon-tools.sif`, one Slurm array task per chromosome (array job 14452681, merge job 14452682) on ILIFU compute nodes |
| Table | `/cbio/users/mamana/sfs/v6hc_s_afr.tsv`, copied to the author's Mac; 667,077 bytes, `sha256 df1febf6173d4c771694636fee30754d4b25a70da6fda6bb1b6a8dac11bcc992` |
| Fit record | `results/tables/sfs_fit_v6hc_s_afr.json` (gitignored; the numbers are copied below; regenerate with the command in `analysis/simulation/docs/fitting-the-sfs.md`). Fitted locally on the author's Mac from the table above (its sha256 was re-read on 2026-10-27 and matches); not an ILIFU run, only the histogram step was |

**Checked 2026-10-23 (merge job log).** All 22 chromosomes report
`skipped 0 at sites with AN != 3790` and `0 with AC = 0`, so no site was dropped
for missing genotypes and no record without a carrier was counted. The 22 per-chromosome
`kept` counts (chr1 3,346,708 through chr21 560,238) sum to 42,331,932, the same
total the fit read from the merged table. A slice test on
`chr22:16,000,000-17,000,000` kept 9,383 ALT alleles.

**Index warning traced (read-only, 2026-10-23).** On all 22 chromosomes the
`.tbi` is older than its VCF by 35 seconds to about 3 minutes (for example
chr22: VCF 2026-03-23 15:30:43, index 15:30:08). The warning is therefore a
property of how the panel was written, not of one file. The histogram runs
read each chromosome end to end without `-r`, which does not need the index
(standard `bcftools` behaviour, not re-tested here), and every chromosome's
counts are consistent with the allele number, so the counts do not depend on
it. Region queries (`-r`) do use the index and are the thing to distrust.

**Still not checked.** The 42,331,932 total equals the variant count in the table
at the top of this file; whether the two were counted the same way (SNP ALT
alleles against all records) was not checked.

**Fit record, copied from the JSON (2026-10-23 14:50).** Zero-truncated
beta-binomial, `N = 1895`, `min_count = 1`, 42,331,932 variants:
`a' = 0.000912` (se 1.7e-7), `b' = 1.715` (se 6.0e-4), log-likelihood
-207,390,462.8. `converged = true` and `at_boundary = true`: the optimiser
stopped, but on the lower search bound for `a'`, so this is a boundary
solution, not an interior optimum. **The fit is poor** and must not be quoted
as a description of the panel:

| Count bin | Observed | Fitted |
| --- | --- | --- |
| singletons (k = 1) | 4,542,441 | 5,258,158 |
| doubletons (k = 2) | 6,774,109 | 2,630,980 |
| k = 3 | 3,395,856 | 1,754,455 |
| fraction with k <= 3 ("rare") | 0.688 | 0.522 |
| fraction singletons | 0.107 | 0.124 |
| top bin (k 2985 to 3790) | 639,156 | 252,131 |

The panel has far more doubletons and rare alleles than any single Beta
reproduces, and the fitted shape under-predicts the high-frequency tail. The
fit's query figure (3,810 for `delta = 1e-6`, 0.035 of the fixture's 108,898)
is therefore unreliable; use the histogram-based `D_N` table below instead.

**What the fit says.** The single-Beta fit sits on its search bound
(`a' = exp(-7)`, `at_boundary = true`), so its parameters and its query budget are
not to be quoted. `D_N` taken straight from the histogram
(`analysis/simulation/empirical_dn.py`, PR #9), target heterozygous, `delta = 1e-6`,
`alpha = 0.05`, 95% power:

| Beacon size | `D_N` | Queries |
| --- | --- | --- |
| 20 | 0.0824 | 31 |
| 65 | 0.0338 | 78 |
| 500 | 0.0067 | 404 |
| 1,000 | 0.0032 | 851 |

It assumes Hardy-Weinberg pairing, independent sites and unrelated panel samples.
