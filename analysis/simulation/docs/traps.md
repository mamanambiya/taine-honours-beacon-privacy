# The three traps

Each produces a plausible wrong answer rather than an error.

[← back to the package README](../README.md)

Each one produces a plausible-looking wrong answer rather than an error.

**1 — the sign of `L`.** `L = log P(data|H0) - log P(data|H1)`, so **small `L`
means the target IS in the beacon**. Infer membership on large `L` and the
power curve comes out inverted and still looks like a power curve. Assert the
direction against a case where membership is known by construction.

**2 — ascertainment.** Queries come only from het sites, and
`P(het|f) = 2f(1-f)`, so the sampled spectrum is `Beta(a'+1, b'+1)`, not
`Beta(a', b')`. Measured: het-conditioning lowers `D_N` by a factor of 0.006 to
0.59 depending on spectrum and `N`. Using the unconditioned SFS overstates
`D_N` and understates the attack's power — *by a cohort-dependent amount*,
which is a confound for this study specifically.

**3 — `D_(N-1)` under `H1`.** A member's "no" needs their own allele mis-called
**and** none of the *other* `N-1` carrying it. Writing `delta*D_N` models the
target as absent under both hypotheses. It also makes the second logarithm
collapse to `log(1/delta)` — a cancellation that looks like a simplification
and is the bug showing through.
