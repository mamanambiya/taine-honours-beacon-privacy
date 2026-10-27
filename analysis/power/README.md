# Attack power on the loaded beacons

`real_beacon_attack.py` runs the membership attack against the answers of the
12 `<POP>_M65_S1` beacons loaded on `afrigend-all-dev`, not against a
simulated cohort. Tests: `python3 -m unittest test_real_beacon_attack` in this
folder (26 tests, mutation-checked).

## Method

A boolean beacon answers yes to "does anyone here carry this allele?" exactly
when the allele is in the set of sites its members carry. The attacker holds
one genome and may query only that person's heterozygous positions. With the
mismatch rate `delta` negligible, the likelihood-ratio test reduces to one
rule: reject membership of anyone who receives a single "no".

- **Members are always detected.** Every allele a member carries is in the
  beacon, so a member never receives a "no": power is 1 (65 of 65 members, in
  all 12 beacons).
- **What varies is the false-positive rate.** It is the share of non-members
  whose first `n` queried heterozygous sites all answer yes. For one outsider
  with `H` heterozygous sites, `m` of them absent from the beacon, drawing `n`
  without replacement gives `P(all yes) = prod_{i<n} (H - m - i) / (H - i)`.
  This is exact, so no query is sent; the rate at budget `n` is that value
  averaged over outsiders.
- **The reported number** is the smallest `n` whose false-positive rate is
  5 per cent or lower, with a 95 per cent interval from resampling the
  outsiders (1,000 draws, seed 1).
- **Check that the lookup is the beacon.** 400 live API queries (a member and
  an outsider in each of two beacons) agreed with the set lookup every time.

## Result (2026-10-27, beacon of 65, FPR 5 per cent)

Two outsider sets: the beacon's own 20 held-out controls (same sub-population,
wide intervals) and everyone with a heterozygous table on the same continent
who is not a member of the beacon (about 395 people, mostly other
sub-populations).

| Beacon | Member sites | 20 controls | Same-continent outsiders |
| --- | --- | --- | --- |
| ESN | 1,002,769 | 138 (133-144) | 92 (89-94) |
| GWD | 1,027,561 | 122 (116-128) | 82 (81-84) |
| LWK | 1,038,569 | 142 (134-152) | 82 (79-84) |
| MSL | 1,032,590 | 136 (129-144) | 97 (95-99) |
| YRI | 1,017,091 | 128 (122-134) | 95 (93-98) |
| **AFR (pooled)** | 1,077,584 | 116 (111-121) | 109 (107-111) |
| CEU | 610,835 | 147 (138-156) | 139 (136-143) |
| FIN | 565,814 | 248 (232-265) | 117 (113-122) |
| GBR | 584,895 | 187 (178-196) | 134 (131-136) |
| IBS | 664,903 | 123 (114-133) | 144 (142-147) |
| TSI | 609,724 | 160 (152-167) | 132 (129-135) |
| **EUR (pooled)** | 627,468 | 138 (125-150) | 150 (147-152) |

Queries needed, summarised:

| | African populations (5) | European populations (5) | Pooled AFR | Pooled EUR |
| --- | --- | --- | --- | --- |
| Same-continent outsiders | mean 89.6 (82-97) | mean 133.2 (117-144) | 109 | 150 |
| 20 controls | mean 133.2 (122-142) | mean 173.0 (123-248) | 116 | 138 |

- On same-continent outsiders every African beacon needs fewer queries than
  every European one, about two thirds as many (89.6 against 133.2; pooled 109
  against 150, ratio 0.73). An exact permutation of the ten population labels
  (252 splits) gives a two-sided p of 2/252 = 0.008.
- On the 20 controls the direction is the same (133 against 173) but the ranges
  overlap (IBS needs 123) and the permutation p is 18/252 = 0.071.
- The mean share of an outsider's heterozygous sites absent from the beacon
  (the empirical `D_N`) is 0.0284 for pooled AFR and 0.0211 for pooled EUR,
  against 0.0280 and 0.0211 from the closed form in `data/README.md`. The
  closed-form budgets there (95 and 126) are lower than the measured 109 and
  150 because the closed form assumes every outsider has the same `D_N`; real
  outsiders do not, and the few who resemble members for longest set the
  budget at 5 per cent.

## What this does and does not show

- It is a measurement on these 12 beacons, each drawn once (seed 1). The
  intervals come from resampling outsiders only; they do not include the
  variation from which 65 people were drawn into a beacon, which was not
  repeated.
- The outsiders are shared: the same people serve as outsiders for several
  beacons, so the ten population budgets are not independent observations, and
  the permutation p is conditional on this design.
- Same-continent outsiders are mostly from other sub-populations, which are
  easier to tell apart from the members than people from the beacon's own
  sub-population. The 20-control column is the same-sub-population test, and it
  is the weaker one (20 people per beacon; FIN's 248 is an outlier).
- Beacon size also differs: the African beacons list 1.0 million sites and the
  European ones 0.6 million, drawn from one panel, chromosomes 20 to 22 only.
  This design measures that ancestry-linked difference; it does not separate
  ancestry from the panel's sampling or from sub-population structure.
- Two within-beacon pairs of second-degree relatives remain (see
  `data/README.md`). The attack was not repeated with them removed.
- No query was sent for the scoring itself, and the answers assume
  `delta` is negligible. A larger `delta` or a censored beacon would change
  the rule and is not covered.
