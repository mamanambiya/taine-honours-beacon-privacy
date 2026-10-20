# API reference

Every public function, what it returns, and how to choose parameters.

[← back to the package README](../README.md)

## `lrt.py` — the statistic

Standard library only. `math.lgamma` keeps `D_N` from underflowing at large
`N`, and there is no scipy dependency to install.

| Function | |
| --- | --- |
| `d_n(a, b, n_individuals, heterozygous=False)` | `D_N = E_f[(1-f)^(2N)]` for `f ~ Beta(a, b)`, exact as `B(a, b+2N) / B(a, b)`. `heterozygous=True` applies the ascertainment correction. |
| `theta_h0(d_n_value)` | `1 - D_N` — P(yes) for a non-member |
| `theta_h1(d_n_minus_1, delta)` | `1 - delta*D_(N-1)` — P(yes) for a member. Takes `D_(N-1)`, not `D_N`. |
| `lrt_statistic(responses, d_n_value, d_n_minus_1, delta)` | `L`. Uses only `sum(responses)`. |
| `required_queries(d_n_value, d_n_minus_1, delta, alpha=0.05, power=0.95)` | Closed-form query count. Verified against the paper's A35. |
| `rejection_threshold(n_queries, d_n_value, d_n_minus_1, delta, alpha=0.05)` | The `L` below which `H0` is rejected, from the Gaussian approximation. |
| `simulate_responses(n_queries, d_n, d_n_minus_1, delta, in_beacon, seed=None)` | Draw responses for a target whose membership is fixed by construction. |

`_z(p)` is Acklam's rational approximation to the normal quantile,
`|err| < 1.15e-9` — so the module needs no scipy for `norm.ppf`.

Two inputs are rejected rather than silently accepted:

- **`delta = 0`** — a single "no" would then be infinitely strong evidence
  against membership, and `L` diverges. Use a small positive value.
- **`D_N <= delta*D_(N-1)`** in `required_queries` — the hypotheses do not
  separate, so there is no answer. Raising beats returning a number.

## `cohort.py` — a beacon made of genotypes

| Function | |
| --- | --- |
| `simulate_cohort(n_individuals, n_variants, a, b, seed)` | Genotype matrix of alternate-allele counts. `f ~ Beta(a,b)` per variant, then `Binomial(2, f)` per individual — Hardy-Weinberg at that frequency. |
| `alternate_allele_counts(genotypes)` | Alternate copies per variant, summed over individuals. |
| `beacon_answers(genotypes, censor_singletons=False)` | What the beacon replies. Censoring suppresses sites carrying exactly one alternate copy. |
| `heterozygous_positions(genotypes, individual)` | The positions an adversary is allowed to query. |
| `attack(genotypes, target, n_queries, delta, seed, a, b, censor=False)` | Query the target's het sites, return `L`. |
| `empirical_power(beacon, outsiders, n_queries, delta, a, b, alpha=0.05, censor=False)` | `(power, achieved_fpr)` at a threshold calibrated on non-members. |

**`delta` corrupts only the target's own contribution**, never the beacon's
whole answer. That is what makes `P(no | member) = delta * D_(N-1)` — the
member's allele must be mis-called *and* none of the other `N-1` may carry it.
Corrupting the whole answer models a different, easier attack.

**`empirical_power` calibrates its threshold on the data, not the model.** This
matters whenever the beacon does something the attacker's model does not know
about. Singleton censoring is exactly that case: the LRT is derived for an
uncensored beacon, so under censoring it is mis-specified, the raw gap between
member and non-member scores can move in *either* direction, and only power at
a controlled false-positive rate still measures anything.

## `contrast.py` — does the spectrum matter at all?

| Function | |
| --- | --- |
| `rare_fraction(a, b, n_individuals, threshold=0.01, ...)` | Share of variants below 1% MAF. A readable axis — Beta parameters are not. |
| `contrast(shapes, ...)` | Empirical power for each `(label, a, b)` spectrum. |
| `budget_for_power(a, b, target=0.95, ...)` | Smallest query budget reaching `target` power, or `None` if none does. |

`None` is a result, not a failure: it means the attack cannot exploit that
spectrum at any budget searched.

## `sfs_fit.py` — fit the spectrum to a real panel

numpy only. Full walk-through, including the ILIFU step: [fitting-the-sfs.md](fitting-the-sfs.md).

| Function | |
| --- | --- |
| `log_pmf(k, n_individuals, a, b)` | `log P(K = k)` for `f ~ Beta(a, b)`, `K ~ Binomial(2N, f)`. `P(K = 0)` equals `lrt.d_n(a, b, N)`. |
| `count_histogram(counts, n_individuals)` | Variants per alternate-allele count `0 .. 2N`. Refuses a count above `2N`. |
| `read_histogram(path, n_individuals)` | Reads `count<TAB>variants`; repeated counts are summed, so per-chromosome tables concatenate. |
| `fit_beta_binomial(histogram, n_individuals, min_count=1)` | MLE `(a, b)` conditioned on `K >= min_count`. Returns `a, b, se_a, se_b, loglik, converged, at_boundary, n_variants`. Refuses an empty histogram, a wrong length, counts below `min_count`, and a single occupied bin. |
| `diagnostics(histogram, n_individuals, a, b, min_count=1)` | Observed vs fitted singleton, rare (MAF < 1%) and high-frequency (alt > 50%) shares, plus a binned table. |
| `required_queries_for_fit(a, b, n_individuals, delta, ...)` | Budget under the ascertained spectrum; `(0.5, 2.0)` at `N = 1895`, `delta = 1e-6` gives 108,898. |
| `compare_to_fixture(a, b, n_individuals, delta, ...)` | The fitted budget beside the fixture's, and their ratio. |

`python3 sfs_fit.py HIST.tsv --n-individuals N --out FIT.json` does the lot and
writes the provenance record. `sfs_histogram.sh` is the cluster-side step.

## `make_figures.py`

Six panels into `results/figures/sfs-explained.png`: what an SFS is, the
ascertainment shift, `D_N` against beacon size, predicted queries, measured
queries against rare-variant load, and the SFS-invariance of empirical power.
Every panel is computed from the modules above.

## What you get back

| Call | Returns |
| --- | --- |
| `simulate_cohort(...)` | `numpy.ndarray`, shape `(individuals, variants)`, dtype int, values `{0,1,2}` |
| `heterozygous_positions(...)` | `numpy.ndarray` of column indices |
| `beacon_answers(...)` | boolean `numpy.ndarray`, one per variant |
| `attack(...)` | one `float` — the statistic `L`. Low means member. |
| `empirical_power(...)` | `(power, achieved_fpr)`, both floats in `[0, 1]` |
| `d_n(...)` | `float` in `(0, 1)` |
| `required_queries(...)` | `int >= 1` |
| `contrast(...)` | `list[dict]` — `label, a, b, rare_fraction, power, false_positive_rate` |
| `budget_for_power(...)` | `int`, or `None` when no budget searched reaches the target |
| `make_figures.py` | `results/figures/sfs-explained.png` |

`contrast()` output is JSON-serialisable on purpose — write it to
`results/tables/` rather than reading numbers off a terminal.

## Choosing parameters

| | Meaning | Sensible range |
| --- | --- | --- |
| `a`, `b` | Beta shape of the SFS. Small `a`, large `b` = rare-variant-rich. | `Beta(0.2, 20)` is realistically rare; `Beta(5,5)` is a control that should fail |
| `n_individuals` | Beacon size `N`. Bigger hides members better. | 40 for fast iteration; 1,895 is the live panel |
| `n_variants` | Sites in the cohort. Caps how many het sites a target has. | 20,000 gives ~340 het sites at `Beta(0.2,20)` |
| `n_queries` | Query budget. **The thing to report.** | sweep it; do not fix it |
| `delta` | Mismatch rate. Inert below ~0.05; bites hard above. | `1e-6` (paper) to `1e-3` |
| `alpha` | False-positive rate. | 0.05 |
| `seed` | Always set it, always record it. | any int |

**Runtime is dominated by `n_outsiders x n_queries`**, because every outsider
must be scored to calibrate the threshold. The defaults run in seconds; raising
`n_variants` costs memory rather than time.
