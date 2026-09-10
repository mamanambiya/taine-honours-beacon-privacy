# What the simulation established

Measured with this code. And what it does not establish.

[← back to the package README](../README.md)

Measured with this code; the commands are in the modules' own docstrings.

**The study's premise holds.** Attack cost moves across the entire range with
rare-variant load alone — 10 queries at 93% rare, 640 at 15%, and unattackable
below that. Since the spectrum is the only route by which a cohort enters the
model, a flat result here would have killed the hypothesis for any cohorts.

**Report the required budget, not power at a fixed budget.** Real cohorts sit
at the high-rare end where power saturates at 1.000. An African-versus-European
comparison run as "power at 40 queries" would report 1.000 against 1.000 and
detect nothing, while the underlying costs differ by an order of magnitude.

**Singleton censoring is a cost, not a cure.** It raises the budget roughly 8x
(20 queries to 160 at 95% power, rare spectrum, `N=40`) and then stops
mattering.

**The two model corrections are wildly asymmetric.** At `N=1895, delta=1e-6`:
the `D_(N-1)` fix changes the required count by **0%**, the ascertainment fix
by **885x**. Both were found in the same review, which made them look
comparably serious. Only measuring separated them.

**Empirical power is invariant to the assumed SFS.** Attackers assuming the
true `Beta(0.2,20)`, a wrong `Beta(0.5,2)` and an absurd `Beta(5,5)` all achieve
identical power, because `L` is monotone in the yes-count and a threshold
calibrated on observed non-members absorbs any `D_N`. **So the 885x describes
the paper's predicted budget, not what an adversary achieves.** Operational
risk has to be measured, never read off `required_queries`.

## What this does NOT establish

- **Nothing about real cohorts.** Every `Beta(a, b)` here is chosen to span a
  range, not fitted to 1000G or V6HC-S_AFR. Mapping the axis onto real
  populations needs a fit — task 6.
- **No absolute query count against the live beacon.** These are `N=40`
  simulated beacons; the live panel holds 1,895.
- **No ancestry claim.** The spectrum is a mechanism, not an ancestry label.
  That two populations differ in rare-variant load is well established; that
  the difference moves attack power measurably is what tasks 2 and 3 are for.
