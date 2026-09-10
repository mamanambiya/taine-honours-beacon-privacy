# `analysis/simulation`

The Shringarpure & Bustamante (2015) membership-inference attack against a
GA4GH Beacon, and a measurement of how its power depends on the cohort.

Reference: [doi:10.1016/j.ajhg.2015.09.010](https://doi.org/10.1016/j.ajhg.2015.09.010)
· [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC4667107/)

## Quick start

From the repository root:

```bash
python3 -m unittest discover -s analysis/simulation -t analysis/simulation
```

```text
Ran 47 tests in 2.8s

OK
```

```bash
python3 analysis/simulation/make_figures.py
```

```text
wrote .../results/figures/sfs-explained.png
```

Python 3 with numpy. `lrt.py` on its own is standard library only — no scipy
anywhere, so there is nothing to install beyond numpy.

## Documentation

| | |
| --- | --- |
| [the-attack.md](docs/the-attack.md) | What the attack is, the two hypotheses, and which inputs come from the beacon rather than from you |
| [how-it-works.md](docs/how-it-works.md) | The six steps of the simulation, and why it is split into a predicted layer and a measured one |
| [examples.md](docs/examples.md) | Worked examples with real output — one attack, a power curve, the cost of censoring |
| [api.md](docs/api.md) | Every public function, what it returns, and how to choose parameters |
| [traps.md](docs/traps.md) | Three ways to get a plausible wrong answer |
| [findings.md](docs/findings.md) | What this established in full, and what it does not |
| [testing.md](docs/testing.md) | 47 tests, mutation checks, and why two injected defects survive on purpose |

## The files

| | |
| --- | --- |
| `lrt.py` | The statistic. `d_n`, `theta_h0`, `theta_h1`, `lrt_statistic`, `required_queries`, `rejection_threshold`, `simulate_responses` |
| `cohort.py` | A beacon built from genotypes. `simulate_cohort`, `beacon_answers`, `heterozygous_positions`, `attack`, `empirical_power` |
| `contrast.py` | The spectrum sweep. `rare_fraction`, `contrast`, `budget_for_power` |
| `make_figures.py` | Six explanatory panels into `results/figures/` |
| `test_*.py` | 47 tests |

## The one-paragraph version

A beacon answers *"do you hold allele A at position P on chromosome C?"* with
yes or no. An adversary holding a target genome queries the positions where
that target is heterozygous. A member answers yes almost always; a non-member
answers yes only when somebody else happens to carry the same allele. A
likelihood-ratio test on the **count of yes answers** separates the two, and
how many queries that takes depends almost entirely on how much rare variation
the cohort carries — which is what makes this an ancestry question.

## Where this sits

The repository's top-level `docs/` describes **the study** — the method as
reconstructed from the paper, the live beacon's constraints, the findings.
This folder describes **the code**.

| | |
| --- | --- |
| `../../docs/methods-lrt.md` | The statistic written out, with each line's verification status |
| `../../docs/beacon-constraints.md` | Rate limits, genome build, coordinates |
| `../../worklist-for-taine.md` | The work plan these tasks belong to |
