# How the simulation works

The six steps, and the two layers they are split across.

[← back to the package README](../README.md)

Six steps. Each is one function, and every one is separately testable.

## 1. Draw a cohort

```text
for each of n_variants:      f ~ Beta(a, b)              its allele frequency
for each individual:         genotype ~ Binomial(2, f)   0, 1 or 2 alt copies
```

Result is an integer matrix, `individuals x variants`, values in `{0,1,2}`.
Hardy-Weinberg at each variant's own frequency. `simulate_cohort` takes an
explicit `seed`; record it with any output.

## 2. Split it

The first rows become **the beacon**. The rest are **outsiders** — people
genuinely not in it, which is what the false-positive rate is measured on.
Nothing distinguishes the two groups except which slice they landed in.

## 3. Choose the query set

For a target, take the positions where **that target is heterozygous**
(genotype `== 1`) and sample `n_queries` of them. This is not a modelling
convenience — it is what an adversary holding a VCF would do, and it is what
makes the sampled spectrum `Beta(a+1, b+1)` rather than `Beta(a, b)`.

## 4. Let the beacon answer

The beacon counts alternate copies at each queried position and answers yes
when the count clears a threshold — `> 0` normally, `> 1` when censoring
singletons.

The target's own contribution is handled separately, because it has to be:

```text
others_count = alt copies among everyone EXCEPT the target
target_count = the target's own copies, ZEROED with probability delta
answer       = (others_count + target_count) > threshold
```

That is what makes `P(no | member) = delta * D_(N-1)`: the member's allele must
be mis-called **and** none of the other `N-1` may carry it. Corrupting the
beacon's whole answer instead would model a different, easier attack.

## 5. Score the responses

Feed the yes/no vector to `lrt_statistic`. It uses only the count of yes
answers. **Low `L` means the target is in the beacon.**

## 6. Turn scores into power

Score every outsider, take the `alpha` quantile of that distribution as the
rejection threshold, then score every member and count how many fall below it.

**The threshold is calibrated on observed non-members, not taken from the
model.** That is deliberate: it keeps the measurement honest when the beacon
does something the attacker's model does not know about, and it is why the
result is a statement about an adversary rather than about the formula.

```text
  simulate_cohort ─┬─> beacon    ─┐
                   └─> outsiders ─┤
                                  ├─> attack() per person ─> L scores
   het positions of target  ──────┘                            │
                                                               ▼
                                        threshold = quantile(outsider L, alpha)
                                        power     = P(member L < threshold)
```

## The two layers, and why there are two

| | module | answers |
| --- | --- | --- |
| **predicted** | `lrt.py` | what the paper's closed form says the power should be |
| **measured** | `cohort.py` | what the attack actually achieves against a genotype matrix |

`lrt.py` computes `D_N` from an assumed spectrum. `cohort.py` builds an actual
cohort, lets the beacon answer from allele counts, and calibrates its rejection
threshold on observed non-members.

**Keeping them separate is the point.** A closed form can agree with itself
while both it and the model are wrong. Where the two layers disagree, the
disagreement is the finding — see *What this established* below.
