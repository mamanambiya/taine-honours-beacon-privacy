# Getting started — step by step

From nothing to a reproduced finding. Nine steps, each with the output you
should see and what to do when you do not see it.

[← back to the package README](../README.md)

Run everything from the **repository root** unless a step says otherwise.

## Step 1 — Check Python

```bash
python3 --version
```

```text
Python 3.14.5
```

Anything 3.9 or newer works. `lrt.py` uses only the standard library.

## Step 2 — Check numpy

```bash
python3 -c "import numpy; print('numpy', numpy.__version__)"
```

```text
numpy 2.4.6
```

If that fails: `python3 -m pip install numpy`. Nothing else is needed — there
is no scipy dependency anywhere, deliberately, so this package runs on a
cluster login shell or a fresh laptop without a build toolchain.

## Step 3 — Run the tests before changing anything

```bash
python3 -m unittest discover -s analysis/simulation -t analysis/simulation
```

```text
Ran 47 tests in 2.6s

OK
```

**Do this first, every time.** A failure here means your environment differs
from the one the numbers below came from, and every later step would be
measuring that instead of the model.

## Step 4 — Compute your first `D_N`

`D_N` is the probability that **nobody** in a beacon of `N` people carries the
alternate allele at a site. Everything else is built on it.

```bash
cd analysis/simulation
python3 -c "
import lrt
print(round(lrt.d_n(0.5, 2.0, 40), 6))
print(round(lrt.d_n(0.5, 2.0, 40, heterozygous=True), 6))"
```

```text
0.147025
0.007657
```

The first is over all sites. The second applies the **ascertainment
correction** — you only query sites where your target is heterozygous, and
those are enriched for intermediate frequencies, so `D_N` falls. Here by a
factor of 19.

If both numbers are identical you passed `heterozygous` wrongly; the whole
885x result depends on this line behaving differently.

## Step 5 — Build a cohort

```bash
python3 -c "
import cohort
g = cohort.simulate_cohort(60, 20000, 0.2, 20.0, seed=1)
print(g.shape)
print(cohort.heterozygous_positions(g, 0).size)"
```

```text
(60, 20000)
339
```

60 individuals, 20,000 variants, values `0/1/2` alternate copies. Individual 0
is heterozygous at 339 of them — those are the only positions an adversary is
allowed to query.

**Change the seed and the numbers change.** Record whichever you use.

## Step 6 — Run one attack

```bash
python3 -c "
import cohort
pool = cohort.simulate_cohort(60, 20000, 0.2, 20.0, seed=1)
b, o = pool[:40], pool[40:]
print(round(cohort.attack(b, b[0], 80, 1e-3, seed=1, a=0.2, b=20.0), 2))
print(round(cohort.attack(b, o[0], 80, 1e-3, seed=2, a=0.2, b=20.0), 2))"
```

```text
-13.23
113.66
```

First is a **member**, second an **outsider**. The member scores far lower.

**If yours come out the other way round, stop.** `L` is
`log P(data|H0) - log P(data|H1)`, so low means member. An implementation that
infers membership on high `L` produces an inverted power curve that still
looks like a power curve. That is trap 1 — see [traps.md](traps.md).

## Step 7 — Turn scores into power

One score tells you nothing on its own. Power needs a threshold calibrated on
people you know are **not** in the beacon.

```bash
python3 -c "
import cohort
pool = cohort.simulate_cohort(240, 20000, 0.2, 20.0, seed=2026)
b, o = pool[:40], pool[40:]
for q in (10, 20, 40, 80):
    p, f = cohort.empirical_power(b, o, q, 1e-3, 0.2, 20.0, alpha=0.05)
    print(f'{q:3d} queries -> power {p:.3f}  fpr {f:.3f}')"
```

```text
 10 queries -> power 0.000  fpr 0.000
 20 queries -> power 1.000  fpr 0.040
 40 queries -> power 1.000  fpr 0.020
 80 queries -> power 1.000  fpr 0.020
```

Power turns on sharply between 10 and 20 queries. The achieved false-positive
rate sits at or below the 0.05 requested — conservative, because a 40-member
quantile is coarse.

## Step 8 — Measure a defence

Singleton censoring makes the beacon answer "no" wherever exactly one
alternate copy exists in the whole cohort.

```bash
python3 -c "
import cohort
pool = cohort.simulate_cohort(240, 20000, 0.2, 20.0, seed=2026)
b, o = pool[:40], pool[40:]
for q in (20, 40, 80, 160):
    plain = cohort.empirical_power(b, o, q, 1e-3, 0.2, 20.0)[0]
    cens  = cohort.empirical_power(b, o, q, 1e-3, 0.2, 20.0,
                                   censor_singletons=True)[0]
    print(f'{q:3d} queries -> plain {plain:.3f}  censored {cens:.3f}')"
```

```text
 20 queries -> plain 1.000  censored 0.400
 40 queries -> plain 1.000  censored 0.475
 80 queries -> plain 1.000  censored 0.975
160 queries -> plain 1.000  censored 1.000
```

Roughly **8x** the budget, then it stops mattering. A cost, not a cure — which
is the argument against treating a single frequency threshold as a solution.

## Step 9 — Draw the figures

```bash
cd ../..
python3 analysis/simulation/make_figures.py
```

```text
wrote .../results/figures/sfs-explained.png
```

Six panels: what a frequency spectrum is, the ascertainment shift, `D_N`
against beacon size, predicted queries, measured queries against rare-variant
load, and the invariance of empirical power to the assumed spectrum.

Every panel is computed from the modules you just ran.

## What to do next

| | |
| --- | --- |
| Understand what you just ran | [how-it-works.md](how-it-works.md) |
| Change a parameter and know what it does | [api.md](api.md) |
| Avoid the three plausible wrong answers | [traps.md](traps.md) |
| See what all of this established | [findings.md](findings.md) |

## When a number does not match

Work down this list in order.

1. **Did step 3 pass?** If not, nothing else is meaningful.
2. **Is the seed the same?** Every function taking `seed` is deterministic
   given it, and different given anything else.
3. **Is `heterozygous=True` where it belongs?** It changes `D_N` by more than
   an order of magnitude, and forgetting it is the largest single error
   available here.
4. **Is your member scoring lower than your outsider?** If not, see trap 1.
5. **Different numpy version?** The generator is stable across 2.x, but a
   major-version change can move the draws. Say which version produced a
   number you report.
