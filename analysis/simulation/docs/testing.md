# Testing

47 tests, and why two deliberately injected defects are allowed to survive.

[← back to the package README](../README.md)

47 tests. `test_lrt.py` also carries `TestPaperEquationA35`, which asserts the
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
