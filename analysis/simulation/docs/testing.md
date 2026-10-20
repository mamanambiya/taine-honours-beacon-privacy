# Testing

85 tests, and why two deliberately injected defects are allowed to survive.

[← back to the package README](../README.md)

85 tests. `test_lrt.py` also carries `TestPaperEquationA35`, which asserts the
paper's literal A35 — printed with a minus over a negative denominator — agrees
with the conventional plus-form used here, across 24 parameter combinations.
They are the same expression; the test exists so nobody "corrects" the code
toward A35's literal signs.

**The implementation predates its tests, so "watch it fail first" is replaced by
mutation checks**: inject the defect, confirm the suite catches it, restore.
Five were run, and two of them **survived** — the `D_(N-1)` swap at realistic
`delta`, and dropping the ascertainment correction. Both survivals are genuine
invariances rather than coverage gaps, and both are now asserted explicitly in
`TestWhatTheAdversaryActuallyNeeds` so they read as documented behaviour.

A mutation that survives is a finding about the model. It is not a licence to
tighten the test until it fails.

## `test_sfs_fit.py`

38 tests, added with Milestone 1b. The fit is checked by **recovering a known
spectrum**: draw counts from a chosen `Beta(a, b)`, drop the sites a panel would
not list, fit, and require both parameters back within 10%, at `N = 500` and at
the live `N = 1895`. A known-negative control (`TestTruncationMatters`) fits the
same data as though nothing had been dropped and requires that to come out more
than 15% wrong, so the recovery tests are known to tell truncation apart from
luck.

Eight mutations were injected into `sfs_fit.py` and each was caught: truncation
removed, `b` swapped for `a` in the pmf, the below-`min_count` guard dropped,
ascertainment removed from the budget, the boundary flag disabled, `D_(N-1)`
replaced by `D_N`, the high-frequency comparison reversed, and minor allele
frequency replaced by alternate allele frequency. **The first run missed the last
two**: the diagnostics were only tested through a symmetric comparison, so
reversing their direction changed nothing. `test_observed_fractions_are_counted_on_a_histogram_read_by_hand`
now pins them against a histogram worked out by hand.

The shell step is tested on a six-site toy VCF holding a multi-allelic site, an
indel, a site with a missing genotype and an all-reference site, and is skipped
if `bcftools` is absent. That test caught a real defect: an empty `bash` array
under `set -u` aborts on the bash 3.2 that macOS ships, which broke `VTYPES=any`.
