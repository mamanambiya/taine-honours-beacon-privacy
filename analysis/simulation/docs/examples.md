# Worked examples

Every block below was run to produce the output shown.

[← back to the package README](../README.md)

Real output, produced by running exactly what is shown.

## Predicted cost against the live panel

```python
import lrt

d  = lrt.d_n(0.5, 2.0, 1895, heterozygous=True)
d1 = lrt.d_n(0.5, 2.0, 1894, heterozygous=True)
print(lrt.required_queries(d, d1, delta=1e-6))
```

```text
D_N        = 2.489414e-05
D_(N-1)    = 2.491384e-05
queries    = 108898
```

108,898 queries at the production limit of 50/hour is **91 days for one
genome**. Drop `heterozygous=True` and it reads 123 — that difference is the
885x.

## One attack, one target

```python
import cohort

pool = cohort.simulate_cohort(60, 20000, 0.2, 20.0, seed=1)
beacon, outsiders = pool[:40], pool[40:]

member  = cohort.attack(beacon, beacon[0],    80, 1e-3, seed=1, a=0.2, b=20.0)
outside = cohort.attack(beacon, outsiders[0], 80, 1e-3, seed=2, a=0.2, b=20.0)
```

```text
genotype matrix : (60, 20000) (individuals x variants), values [0, 1, 2]
het sites for target 0 : 339
L member   =    -13.23   <- LOWER means in the beacon
L outsider =    113.66
```

The member scores far below the outsider. If yours come out the other way
round, you have hit trap 1.

## Power against query budget

```python
pool = cohort.simulate_cohort(240, 20000, 0.2, 20.0, seed=2026)
beacon, out = pool[:40], pool[40:]
for q in (10, 20, 40, 80):
    print(q, cohort.empirical_power(beacon, out, q, 1e-3, 0.2, 20.0, alpha=0.05))
```

```text
   10 queries -> power 0.000   fpr 0.040
   20 queries -> power 1.000   fpr 0.040
   40 queries -> power 1.000   fpr 0.020
   80 queries -> power 1.000   fpr 0.020
```

Note how sharply it turns on, and that the achieved FPR sits at or below the
0.05 asked for. With 40 members the quantile is coarse, so it is conservative.

## What singleton censoring costs

```python
for q in (20, 40, 80, 160):
    plain    = cohort.empirical_power(beacon, out, q, 1e-3, 0.2, 20.0)[0]
    censored = cohort.empirical_power(beacon, out, q, 1e-3, 0.2, 20.0,
                                      censor_singletons=True)[0]
```

```text
   20 queries -> plain 1.000   censored 0.400
   40 queries -> plain 1.000   censored 0.475
   80 queries -> plain 1.000   censored 0.975
  160 queries -> plain 1.000   censored 1.000
```

Roughly an 8x increase in the budget the attacker needs, then it stops
mattering. A cost, not a cure.
