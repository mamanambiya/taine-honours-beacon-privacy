# The attack

What a Beacon membership-inference attack is, and which inputs it needs.

[← back to the package README](../README.md)

A beacon answers exactly one question — *"do you hold allele A at position P on
chromosome C?"* — with **yes** or **no**. Nothing else: no counts, no
frequencies, no genotypes.

An adversary who already holds a target genome queries the beacon at the
positions where that target is **heterozygous**, asking for the alternate
allele each time. A member of the beacon should answer yes almost every time; a
non-member answers yes only when somebody else happens to carry the same
allele. Enough queries separate the two.

| | | `P(answer = yes)` |
| --- | --- | --- |
| `H0` | target is **not** in the beacon | `1 - D_N` |
| `H1` | target **is** in the beacon | `1 - delta * D_(N-1)` |

`D_N` is the chance that *none* of the `N` genomes in the beacon carries the
allele at that site. `delta` is the mismatch rate between the adversary's copy
of the genome and the beacon's, from sequencing and variant-calling differences.

The log-likelihood ratio is

```text
L = n*B + C * SUM_i x_i        B = log( D_N / (delta*D_(N-1)) )
                               C = log( delta*D_(N-1)*(1-D_N) / (D_N*(1-delta*D_(N-1))) )
```

`L` depends on the responses only through `SUM_i x_i`, **the number of yes
answers**. That is what makes the required-query count a closed form rather
than a simulation.

## What comes from where

Only two of the five inputs come from the beacon.

| Input | Source |
| --- | --- |
| yes/no per query — `x_i` | **the beacon**, `responseSummary.exists` |
| `N`, the beacon size | **the beacon**, `/api/datasets` → `sampleCount` |
| which positions to query | the target's VCF — its heterozygous sites |
| which allele to ask for | the target's VCF — the alternate allele there |
| the site-frequency spectrum | a reference panel, or an assumed `Beta(a', b')` |

**The beacon serves no allele frequencies.** Verified against production: the
response carries `id`, `name`, `setType`, `exists`, `resultsCount` — no more.
The SFS has to come from elsewhere, and fitting it to the real panel is task 6
of the work plan.
