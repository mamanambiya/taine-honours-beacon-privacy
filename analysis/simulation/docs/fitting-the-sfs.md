# Fitting the site-frequency spectrum to a real panel

Milestone 1b. Replaces the `Beta(0.5, 2.0)` fixture behind every absolute query
count with a spectrum fitted to the panel itself.

[← back to the package README](../README.md)

## Why

`D_N`, the chance nobody in the beacon carries an allele, is an expectation over
the site-frequency spectrum. The spectrum is the only route by which a cohort
enters the attack model, so until it is fitted the absolute figures — the
108,898 queries per genome in `docs/beacon-constraints.md` — describe a fixture,
not the AfriGen-D panel. The 885x ascertainment ratio does not depend on it; the
absolute number does.

It cannot be read off the beacon. Production is Boolean-mode and returns no
per-variant frequency, so the counts come from the source VCFs on ILIFU.

## The two steps

Only a two-column table leaves the cluster. The fit itself runs anywhere with
numpy.

| Step | Where | Tool |
| --- | --- | --- |
| 1. Count alternate alleles, per VCF | ILIFU **compute node** | `sfs_histogram.sh` |
| 2. Fit `Beta(a', b')`, cost the attack | your machine | `sfs_fit.py` |

### Step 1 — the histogram (on a compute node)

`sfs_histogram.sh VCF N_INDIVIDUALS [OUT.tsv]` writes `count<TAB>variants`.
`bcftools` lives in the beacon tools container, so pass it through `BCFTOOLS`
(a shell function such as `beacon-sing` is not visible to a script). The
container, the `/cbio` bind and the compute-node request are exactly the ones in
[`ILIFU_DATA_LOADING_GUIDE.md`](../../../../afrigen-beacon-v2/docs/ILIFU_DATA_LOADING_GUIDE.md)
sections 4.1 and 4.2.1 (A1) of the beacon repository: follow those first.

```bash
BCFTOOLS="singularity exec --bind /cbio,/scratch3 /cbio/users/$USER/beacon-tools.sif bcftools"
export BCFTOOLS
PANEL=/cbio/dbs/refpanels/h3a_reference_panels/version_6/v6hc_s_african/V6HC-S_AFR/vcfs
mkdir -p /cbio/users/$USER/sfs
```

Try one small slice before the full run, and compare it with what you expect:

```bash
VIEW_ARGS="-r chr22:16000000-17000000" ./analysis/simulation/sfs_histogram.sh \
    $PANEL/V6HC-S_AFR_chr22_all.vcf.gz 1895 /cbio/users/$USER/sfs/slice.tsv
```

The last line it prints says how many ALT alleles it kept and how many sites it
skipped. **On the full panel, `skipped ... AN != 3790` should be 0 or tiny.** A
large number means missing genotypes or the wrong `N`, and the table would be a
biased sample of the panel.

Then every autosome. This is a stream over about 354 GB; run it as a batch job
or inside `tmux`, not interactively on the login node:

```bash
for chr in $(seq 1 22); do
    ./analysis/simulation/sfs_histogram.sh \
        $PANEL/V6HC-S_AFR_chr${chr}_all.vcf.gz 1895 \
        /cbio/users/$USER/sfs/v6hc_s_afr_chr${chr}.tsv
done
cat /cbio/users/$USER/sfs/v6hc_s_afr_chr*.tsv > /cbio/users/$USER/sfs/v6hc_s_afr.tsv
```

Copy `v6hc_s_afr.tsv` (a few kilobytes) back to your own machine.

**Not run on ILIFU.** The script is tested on a toy VCF and the 1000G example
file (`tests: TestExtractionScript`), and the loop above has not been run on the
cluster. Expect to adjust the loop to Slurm's conventions (a job array, one
chromosome per task, is the natural form).

**Defaults, and what to change.** `VTYPES=snps` keeps SNPs only; the 2015 attack
and its spectrum are about SNPs. Use `VTYPES=any` only if the query sites will
include indels. `VIEW_ARGS="-f PASS"` mirrors a quality filter. Whatever you
choose, **use the same site set the beacon holds**, because the spectrum that
governs `D_N` is the spectrum of the sites that can be queried; say what you
chose in `data/README.md`.

**The matched 1000G arms.** Subset to the samples first and give `N` as the
number of samples: `VIEW_ARGS="-S eur_samples.txt"` with `N` set to the line
count. Counts are recomputed after subsetting, so sites that are monomorphic in
the subset drop out, which is what a beacon built from that subset would hold.

### Step 2 — the fit

```bash
python3 analysis/simulation/sfs_fit.py v6hc_s_afr.tsv \
    --n-individuals 1895 --label "V6HC-S_AFR" --delta 1e-6 \
    --out results/tables/sfs_fit_v6hc_s_afr.json
```

It prints the fitted `a'`, `b'` with standard errors, and the budget for 95%
power at `alpha = 0.05` under the fit next to the fixture's. The JSON carries
the source path, `N`, `min_count`, the diagnostics and the comparison, which is
the provenance `data/README.md` asks for.

`--min-count` is the smallest allele count the panel can hold. `1` for a panel
of polymorphic sites (the default); `2` if singletons were filtered out. Read
the lowest occupied row of the table to check: if nothing sits at count 1 and
the panel is not supposed to be filtered, find out why before fitting.

## What the likelihood is

Each site has `f ~ Beta(a', b')`; the `2N` haplotypes then carry
`K ~ Binomial(2N, f)` copies, which is beta-binomial. A panel lists only the
sites it saw, so the likelihood is conditioned on `K >= min_count`. The mass
removed at `K = 0` is exactly `lrt.d_n(a', b', N)`, so the truncation and the
attack model are the same quantity.

**Ignoring the truncation is not a small error.** Fitted as if `K = 0` sites had
been observed and missed, a known `Beta(0.5, 2.0)` comes back more than 15% off
(`TestTruncationMatters`).

## Before quoting the numbers

The fit returns numbers whether or not they deserve trust. Check, in order:

1. **`converged` is true and `at_boundary` is false.** The search is over
   `exp(-7) .. exp(7)`; an estimate on the edge is a bound, not a measurement.
   The command prints a warning and the JSON records both.
2. **The diagnostics agree.** `diagnostics` compares observed and fitted
   fractions for singletons, for sites with minor allele frequency below 1%,
   and for sites where the alternate allele is above 50%. The singleton share
   is the one that matters most: singletons drive `D_N`.
3. **Sites where the reference allele is the rare one.** A real panel carries
   them, as alternate alleles near fixation. A single Beta has no room for that
   second population, so the fit absorbs it into the shape at the rare end. The
   high-frequency share in the diagnostics is the check, and
   `TestDiagnostics.test_flags_a_panel_that_a_single_beta_cannot_describe` shows
   it firing on a deliberately mis-specified panel.

### What the example file shows

Run end to end on the 9-individual 1000G example in the beacon repository
(4,231 variants), the pipeline finishes and **correctly refuses to be trusted**:
`a'` lands on the lower search bound (`at_boundary = true`), and the fitted rare
fraction is 0.009 against 0.022 observed. That file is a slice of a 2,504-sample
panel, so its sites were chosen for being polymorphic in the large panel, not in
these nine people; the zero-truncated model is the wrong one for it. It is a test
of the plumbing, not of the method, and a reminder that the truncation has to
match how the sites were chosen. For V6HC-S_AFR itself the sites are those
polymorphic in the panel, which is the case the model assumes.

## `D_N` without a Beta

`empirical_dn.empirical_d_n(histogram, N, n)` treats the panel's own sites as the
population: for a target heterozygous at a site with `k` of the panel's `M = 2N`
haplotypes, it asks how often none of a beacon of `n` other people carries the
allele. Under random pairing the weight of a `k`-copy site is
`2k(M-k) / (M(M-1))`, and the chance the beacon's `2n` haplotypes, drawn from the
`M - 2` that are not the target's, all miss the allele is
`C(M-1-k, 2n) / C(M-2, 2n)`. `D_N` is the weighted mean over the histogram.

- **No zero-bin needed.** The weight is zero at `k = 0` and `k = M`, so the
  `K = 0` mass the Beta has to model never enters.
- **Unbiased for the quantity `lrt.d_n` integrates.** Haplotypes drawn without
  replacement from iid draws are iid, so given a frequency `f` the expected
  miss probability is `(1 - f)^(2n)`. It reproduces the Beta answer when the
  Beta is right, and follows the data when it is not (`test_empirical_dn.py`).
- **Beacon size.** At most `N - 1`: the target and the beacon are different
  people from the same panel.
- **What it still assumes.** Hardy-Weinberg pairing (the histogram cannot check
  it), independence between sites, and unrelated panel samples. Relatedness or
  population structure inside the panel would make the panel's haplotypes
  correlated, so the miss probability would be too optimistic or pessimistic by
  an amount this module cannot see.

`sfs_fit.py --beacon-sizes 20,65,500,1000` prints and records the figures beside
the Beta's. Run on the full V6HC-S_AFR table (42,331,932 variants, `N = 1895`),
whose Beta fit sat on the bound (`at_boundary = true`):

| Beacon | `D_N` (histogram) | Queries, 95% power (histogram) | Queries (Beta fit, at the bound) |
| --- | --- | --- | --- |
| 20 | 0.0824 | 31 | 41 |
| 65 | 0.0338 | 78 | 131 |
| 500 | 0.0067 | 404 | 1,004 |
| 1,000 | 0.0032 | 851 | 2,009 |

`delta = 1e-6`, `alpha = 0.05`. The histogram column is the one to quote for this
panel; the Beta column is a bound. The "3,810 queries" the fit printed for a
beacon the size of the whole panel is not a usable figure and has no
histogram counterpart, since a beacon cannot be as large as its own panel.

## Limits

- **One Beta is an approximation.** A neutral spectrum falls as `1/f`, which a
  Beta reaches only as `a' -> 0`. If the fit sits near the lower bound on the
  real panel, that is the shape of the data, and the answer is `D_N` taken
  straight from the histogram ([below](#d_n-without-a-beta)), not a better
  optimiser. On V6HC-S_AFR the fit did sit on the bound.
- **The standard errors are from a finite-difference Hessian** and assume the
  model is right. They say how well a Beta is pinned down, not whether a Beta is
  the right shape.
- **Linkage is ignored.** Sites are treated as independent draws. That does not
  bias the point estimate of the spectrum, but it makes the standard errors
  optimistic.
- **Allele counts, not genotypes.** The histogram cannot tell a heterozygote
  from a homozygote, so nothing here checks Hardy-Weinberg.
