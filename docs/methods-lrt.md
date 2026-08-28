# The Shringarpure & Bustamante likelihood-ratio test

The proposal states the LRT as an embedded equation object, which does not
survive conversion out of Word — `proposal/proposal-v3.txt` line 63 reads "The
LRT statistic is stated as:" followed by nothing, and line 65 reads "Where and
, DN being...". This file writes it out so the implementation has something to
be checked against.

**Verification status, 2026-08-21.** Checked against the open-access text at
[PMC4667107](https://pmc.ncbi.nlm.nih.gov/articles/PMC4667107/) by a second
reader, which found two errors in the previous revision of this file — both now
corrected and both flagged inline:

| | Status | Measured impact on required-`n` |
| --- | --- | --- |
| `H1` uses `D_{N-1}` | **corrected** — was `D_N` | **none** at `N=1895, delta=1e-6` (123 -> 123) |
| Ascertainment: `a = a'+1, b = b'+1` | **corrected** — trap was named, fix was not stated | **885x** (123 -> 108,898) |
| Direction: reject `H0` for small `L` | verified — *"rejected if Λ < t"* | — |
| `D_N` as an SFS expectation | verified | — |
| The required-`n` expression (paper's A35) | **verified** — transcribed from the JATS MathML, radicals intact | — |

The body below is still a reconstruction in places. Treat a line without a row
above as unchecked.

Notation is ASCII throughout (`x_i`, `mu_0`, `a'`) so that it copies cleanly
into code and renders identically in every viewer.

## Setup

The adversary holds a target genome as a VCF and can query a beacon of `N`
individuals (`2N` haplotypes).

A query is `q = {C, P, A}` — *"have you observed allele `A` as an alternate
allele at position `P` on chromosome `C`?"*

For a query set `Q = {q_1 ... q_n}` the beacon returns `R = {x_1 ... x_n}` with
each `x_i` in `{0, 1}`.

Two quantities carry the model:

- **`delta`** — the probability of a mismatch between the adversary's copy of
  the genome and the beacon's, from sequencing error. Small.
- **`D_N`** — the probability that *none* of the `N` genomes in the beacon
  carries the alternate allele at a given position.

`D_N` is an expectation over the site-frequency spectrum, not a per-site
constant. With alternate-allele frequency `f ~ Beta(a', b')`:

```text
D_N = E_f[ (1 - f)^(2N) ]
```

This is the only place the cohort enters the test, and it is the whole reason
the study can ask an ancestry question at all: change the SFS, and `D_N`
changes, and with it the power of the test.

## The two hypotheses

| | | `P(x_i = 1)` |
| --- | --- | --- |
| `H0` | target is **not** in the beacon | `theta_0 = 1 - D_N` |
| `H1` | target **is** in the beacon | `theta_1 = 1 - delta*D_{N-1}` |

**`H1` uses `D_{N-1}`, not `D_N`.** Under `H1` the target is one of the `N`, and
at a query site it carries the allele by construction. A "no" therefore needs
*two* things at once: the target's allele mis-called (`delta`), **and** none of
the **other `N-1`** individuals carrying it (`D_{N-1}`). Writing `delta*D_N`
here models the target as absent under both hypotheses.

**How much does this one matter? Measured, at this study's parameters: not at
all.** `N = 1895`, `delta = 1e-6`, `alpha = 0.05`, 95% power —

```text
required_queries with D_(N-1)  =  123
required_queries with D_N      =  123      difference: zero
```

`D_{N-1}/D_N = 1.000264` at `N = 1895`, and `delta` scales that difference, so
at `delta = 1e-6` it vanishes. A sweep puts the correction at 0% for
`delta <= 0.05` at every `N`, reaching -15.9% only at `N = 10, delta = 0.6`.

Keep the fix — it is the correct model and it matters for small beacons with
high error rates. But an earlier revision of this file claimed it "sets the
required number of queries, the number this study reports", and that does not
survive measurement. **The ascertainment correction below is the one that moves
the headline number, by 885x.** Both were found in the same review, which made
them look comparably serious; they are not, and only measuring told us so.

Under `H1` a "no" needs *both* a sequencing mismatch at the target's allele
*and* none of the other `N-1` individuals carrying it — so `H1` still predicts
more "yes" responses than `H0`. The whole test is a comparison of the observed
yes-rate against those two rates.

## The statistic

```text
B = log( D_N / (delta * D_{N-1}) )
C = log( delta*D_{N-1} * (1 - D_N) / ( D_N * (1 - delta*D_{N-1}) ) )

L = n*B + C * SUM_i x_i
```

This is the paper's form (Shringarpure & Bustamante 2015, eq. for `Λ`), verified
against the open-access text at PMC4667107 on 2026-08-21.

**An earlier revision of this file wrote `delta*D_N` in place of
`delta*D_{N-1}`,** which made the second logarithm collapse to `log(1/delta)`.
That cancellation looked like a simplification and was actually the error
showing through — it only exists when the two hypotheses share `D_N`. Code
implemented from that revision produces a plausible-looking curve and a wrong
required-`n`.

`L` is linear in `SUM_i x_i`, which makes the sufficient statistic explicit:
**the number of "yes" responses**. Everything else is a constant given `n`,
`N`, `delta` and the SFS. The binomial structure is what lets the required-`n`
calculation be
done in closed form rather than by simulation.

### Trap 1 — the sign

`L` as written is `log P(data | H0) - log P(data | H1)`. So **small `L` rejects
`H0`**, i.e. small `L` means *the target is in the beacon*. If your
implementation infers membership on large `L`, the power curve will come out
inverted and still look like a plausible curve. Assert the direction against a
synthetic case where the answer is known by construction, before running
anything real.

### Trap 2 — ascertainment bias in the query set

`Q` contains only positions where the target is **heterozygous**. The alternate
allele frequency at those positions is therefore *not* distributed
`Beta(a', b')` — heterozygous sites are enriched for intermediate frequencies.
The proposal flags this (line 62) and the paper handles it.

Using the unconditioned SFS to compute `D_N` for a query set drawn from
heterozygous sites will bias `D_N`, and it will bias it *differently* for
cohorts with different SFS shapes — which is precisely the comparison this study
is making. **This is the single most consequential modelling detail in the
project.**

### The fix, which an earlier revision of this file failed to state

The paper gives it directly: *"We can show that `f_i ~ beta(a, b)`, where
`a = a' + 1` and `b = b' + 1`"*. Conditioning on the target being heterozygous
weights the density by `P(het | f) = 2f(1-f)`:

```text
Beta(a', b') * 2f(1-f)   proportional to   f^a' * (1-f)^b'   =   Beta(a'+1, b'+1)
```

The factor 2 cancels in the normalisation. So the correction is exact and costs
nothing — compute `D_N` under `Beta(a'+1, b'+1)` rather than `Beta(a', b')`:

```text
D_N = B(a'+1, b'+1+2N) / B(a'+1, b'+1)
```

Het sites are enriched for intermediate frequencies, so the conditioned `D_N` is
**smaller**. Using the unconditioned SFS therefore *overstates* `D_N` and
*understates* the attack's power — and by a cohort-dependent amount.

A previous revision named this trap and stopped there, so a reader implementing
from the doc alone would have got it wrong while believing they had been warned.

## Power and the required number of queries

`L` is asymptotically Gaussian under both hypotheses, which gives a closed-form
relation between the number of queries `n`, the false-positive rate `alpha`, the
power `1 - beta`, and the beacon size `N`:

Transcribed from the paper's equation A35, via the Europe PMC JATS full text,
which serves MathML and therefore preserves the radicals that PMC's HTML-to-text
rendering drops:

```text
curl -sL "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4667107/fullTextXML"
```

Equation A34, immediately above it, states the relation the square gives:

```text
SQRT(n) * (delta*D_{N-1} - D_N)
    = z_alpha * SQRT( D_N*(1 - D_N) )
      - z_{1-beta} * SQRT( delta*D_{N-1}*(1 - delta*D_{N-1}) )
```

so

```text
n = [ ( z_alpha*SQRT(D_N (1-D_N)) - z_{1-beta}*SQRT(delta*D_{N-1} (1-delta*D_{N-1})) )
      / ( delta*D_{N-1} - D_N ) ]^2
```

### A sign that looks wrong and is not

A35 as printed has a **minus** between the two `z` terms and a **negative**
denominator (`delta*D_{N-1} - D_N`, and `D_N` is the larger). Implementations —
including `analysis/simulation/lrt.py` — use the standard **plus** with a
**positive** denominator. Both are correct: `z_alpha` as printed is the
*lower*-tail quantile, so it is negative, and the sign flips in numerator and
denominator together and cancels under the square.

Checked numerically rather than argued, at `N` = 85, 1000, 1895 and
`delta` = 1e-6, 1e-3 — the ratio of the two forms is `1.000000` in every case.

**Do not "fix" the plus-form toward the paper's literal signs.** It reads like a
discrepancy, and correcting it without flipping `z_alpha` to the lower tail at
the same time introduces a real error.

The proposal fixes `alpha = 0.05` a priori and reports power as a curve.

## Reference results to replicate

From Shringarpure & Bustamante (2015), as summarised in the proposal:

| Experiment | Setup | Result |
| ---------- | ------------------------- | ------------------------- |
| Synthetic validation | 1,000 diploid individuals, 500,000 SNPs | LRT framework validated against theoretical prediction |
| Real data (CEU, 1000G Phase 1) | 85 individuals, 65 in beacon; 20 in / 20 out queried | Membership inference at **95% power, 5% FPR, with ~250 SNPs** |
| Censoring test | Beacon always answers "no" to singleton queries | Measures the cost of that defence to attack power |

The synthetic validation is the right first milestone: it needs no beacon, no
network, and no data access, and it is where both traps above will surface.

## What the closed form predicts vs what an adversary achieves

**Measured 2026-08-21, and it qualifies every absolute figure in this file.**

`required_queries` is the paper's closed form. It depends on `D_N`, which is an
expectation over the assumed site-frequency spectrum. Change the SFS and the
predicted budget moves enormously — the ascertainment correction alone is worth
885x at `N = 1895`.

**Empirical attack power does not depend on the SFS at all.**

`L` is a monotone linear function of the yes-count, and the SFS only sets its
coefficients. An adversary who calibrates the rejection threshold on observed
non-members therefore absorbs any `D_N`. Measured against a simulated cohort
with the true spectrum, a wrong one, and an absurd one:

| Assumed SFS | Empirical power |
| --- | --- |
| `Beta(0.2, 20)` — the truth | 1.000 |
| `Beta(0.5, 2.0)` — wrong | 1.000 |
| `Beta(5, 5)` — absurd | 1.000 |

So the two numbers answer different questions:

- **`required_queries`** — how many queries the *paper's model* says are needed,
  given an assumed spectrum. Sensitive to the SFS; this is the 885x figure.
- **Empirical power** — what an adversary actually achieves. Insensitive to the
  SFS entirely.

**This matters most where the study is most likely to overclaim.** Quoting
108,898 queries as the operational cost of attacking the AfriGen-D beacon would
be reporting a property of our own modelling assumption as though it were a
property of the adversary. The rate-limit argument in
[`beacon-constraints.md`](beacon-constraints.md) rests on the predicted figure,
so it needs the empirical number before it can carry weight.

The same shape recurs throughout: `D_{N-1}` and `delta` are both structurally
necessary and numerically inert at realistic parameters, biting only at high
error rates. Structural necessity and numerical importance are independent, and
this project has now been caught assuming they travel together twice.

## What "attack power varies with ancestry" would look like

The comparison is the **number of queries required to reach 95% power at
`alpha = 0.05`**, matched on `N`, between an African-cohort beacon and a
European-cohort beacon.

- Fewer queries needed against the African cohort → African data is *more*
  exposed, and the published numbers understate the risk for AfriGen-D.
- More queries needed → the published numbers overstate it, and the case for
  population-aware risk modelling is made.

Either direction is a publishable result. The failure mode to avoid is a
difference that turns out to be driven by `N`, by an off-by-one in coordinates,
or by the ascertainment-bias handling — see
[`beacon-constraints.md`](beacon-constraints.md).
