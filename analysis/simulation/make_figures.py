"""Figures explaining the site-frequency spectrum and what it does to the attack.

Every panel is computed from lrt.py / cohort.py — the same code the study uses,
not illustrative curves drawn by hand.

    python3 analysis/simulation/make_figures.py
"""
import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import cohort
import lrt

OUT = pathlib.Path(__file__).resolve().parents[2] / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

GREY, INK = "#8a8a8a", "#222222"
COLS = ["#1b6ca8", "#c1462d", "#2e8b57", "#7d3c98"]


def beta_pdf(f, a, b):
    logp = ((a - 1) * np.log(f) + (b - 1) * np.log1p(-f)
            - (math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)))
    return np.exp(logp)


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=10, color=INK, pad=8)
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.18, linewidth=0.6)


fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
fig.suptitle("The site-frequency spectrum, and what it does to a Beacon re-identification attack",
             fontsize=13, color=INK, y=0.98)

# ---- 1. what an SFS IS -------------------------------------------------
ax = axes[0][0]
rng = np.random.default_rng(7)
for (a, b, lab), c in zip([(5.0, 5.0, "Beta(5, 5) — common-dominated"),
                           (0.5, 2.0, "Beta(0.5, 2) — intermediate"),
                           (0.2, 20.0, "Beta(0.2, 20) — rare-dominated")], COLS):
    draws = rng.beta(a, b, 200_000)
    ax.hist(draws, bins=np.linspace(0, 1, 51), histtype="step", lw=1.9,
            color=c, label=lab, density=True)
ax.set_ylim(0, 12)
ax.legend(fontsize=7.2, frameon=False)
style(ax, "1. What an SFS IS\nthe distribution of allele frequency across SITES",
      "alternate allele frequency  f", "share of variants")
ax.text(0.55, 10.6, "most variants RARE\n= left-heavy", fontsize=7, color=COLS[2])

# ---- 2. ascertainment shift -------------------------------------------
ax = axes[0][1]
f = np.linspace(0.002, 0.998, 900)
a0, b0 = 0.5, 2.0
ax.plot(f, beta_pdf(f, a0, b0), color=GREY, lw=1.9, label=f"all sites: Beta({a0}, {b0})")
ax.plot(f, beta_pdf(f, a0 + 1, b0 + 1), color=COLS[1], lw=1.9,
        label=f"het sites only: Beta({a0+1}, {b0+1})")
ax.fill_between(f, beta_pdf(f, a0, b0), beta_pdf(f, a0 + 1, b0 + 1),
                where=beta_pdf(f, a0, b0) > beta_pdf(f, a0 + 1, b0 + 1),
                color=GREY, alpha=0.12)
ax.set_yscale("log")
ax.legend(fontsize=7.2, frameon=False)
style(ax, "2. Ascertainment: you only query HET sites\nP(het|f) = 2f(1-f) reweights the spectrum",
      "alternate allele frequency  f", "density (log)")

# ---- 3. D_N vs beacon size --------------------------------------------
ax = axes[0][2]
Ns = np.unique(np.logspace(0.5, 3.6, 60).astype(int))
for (a, b, lab), c in zip([(5.0, 5.0, "common"), (0.5, 2.0, "intermediate"),
                           (0.2, 20.0, "rare")], COLS):
    ax.plot(Ns, [lrt.d_n(a, b, int(n)) for n in Ns], color=c, lw=1.9, label=lab)
    ax.plot(Ns, [lrt.d_n(a, b, int(n), heterozygous=True) for n in Ns],
            color=c, lw=1.4, ls="--", alpha=0.75)
ax.axvline(1895, color=INK, lw=0.9, ls=":")
ax.text(1895, ax.get_ylim()[1], " N=1895\n live panel", fontsize=6.8, va="top", color=INK)
ax.set_xscale("log"); ax.set_yscale("log")
ax.legend(fontsize=7.2, frameon=False, title="solid = all sites\ndashed = het-ascertained",
          title_fontsize=6.5)
style(ax, "3. $D_N$: chance NOBODY in the beacon\ncarries the allele", "beacon size  N", "$D_N$ (log)")

# ---- 4. required queries ----------------------------------------------
ax = axes[1][0]
for (a, b, lab), c in zip([(0.5, 2.0, "intermediate"), (0.2, 20.0, "rare")],
                          [COLS[1], COLS[2]]):
    plain, het = [], []
    for n in Ns:
        d, d1 = lrt.d_n(a, b, int(n)), lrt.d_n(a, b, int(n) - 1)
        dh, dh1 = (lrt.d_n(a, b, int(n), heterozygous=True),
                   lrt.d_n(a, b, int(n) - 1, heterozygous=True))
        plain.append(lrt.required_queries(d, d1, 1e-6))
        het.append(lrt.required_queries(dh, dh1, 1e-6))
    ax.plot(Ns, plain, color=c, lw=1.9, label=f"{lab} — all sites")
    ax.plot(Ns, het, color=c, lw=1.4, ls="--", label=f"{lab} — het-corrected")
ax.axvline(1895, color=INK, lw=0.9, ls=":")
ax.set_xscale("log"); ax.set_yscale("log")
ax.legend(fontsize=7.2, frameon=False)
style(ax, "4. Predicted queries for 95% power\nthe ascertainment gap is the 885x",
      "beacon size  N", "queries needed (log)")

# ---- 5. measured: power vs rare-variant load ---------------------------
ax = axes[1][1]
shapes = [("Beta(5,5)", 5.0, 5.0), ("Beta(1,3)", 1.0, 3.0), ("Beta(0.5,2)", 0.5, 2.0),
          ("Beta(0.5,10)", 0.5, 10.0), ("Beta(0.2,20)", 0.2, 20.0), ("Beta(0.1,40)", 0.1, 40.0)]
import contrast
rows = contrast.contrast([(l, a, b) for l, a, b in shapes])
rare = [r["rare_fraction"] for r in rows]
budget = [contrast.budget_for_power(r["a"], r["b"]) for r in rows]
plotted = [(x, y) for x, y in zip(rare, budget) if y is not None]
unatt = [x for x, y in zip(rare, budget) if y is None]
ax.plot([p[0] for p in plotted], [p[1] for p in plotted], "o-", color=COLS[0], lw=1.9, ms=6)
for x in unatt:
    ax.scatter([x], [1400], marker="x", color=COLS[1], s=60, zorder=5, clip_on=False)
ax.set_ylim(6, 1100)
ax.annotate("attack fails at\nany budget searched", xy=(max(unatt) if unatt else 0.03, 1100),
            xytext=(0.30, 620), fontsize=7, color=COLS[1],
            arrowprops=dict(arrowstyle="->", color=COLS[1], lw=0.9))
ax.set_yscale("log")
style(ax, "5. MEASURED: queries to reach 95% power\nsimulated cohorts, N=40",
      "fraction of variants rare (MAF < 1%)", "queries needed (log)")

# ---- 6. the invariance -------------------------------------------------
ax = axes[1][2]
A, B = 0.2, 20.0
pool = cohort.simulate_cohort(240, 20000, A, B, seed=2026)
beacon, out = pool[:40], pool[40:]
budgets = [5, 10, 20, 40, 80]
for (a, b, lab), c in zip([(A, B, "true Beta(0.2,20)"), (0.5, 2.0, "wrong Beta(0.5,2)"),
                           (5.0, 5.0, "absurd Beta(5,5)")], COLS):
    ys = [cohort.empirical_power(beacon, out, q, 1e-3, a, b, alpha=0.05)[0] for q in budgets]
    ax.plot(budgets, ys, "o-", color=c, lw=1.9, ms=5, alpha=0.85, label=lab)
ax.set_ylim(-0.05, 1.08)
ax.legend(fontsize=7.2, frameon=False, title="SFS the ATTACKER assumes", title_fontsize=6.5)
style(ax, "6. What an adversary ACHIEVES\nidentical — the assumed SFS cancels",
      "query budget", "empirical power")

fig.tight_layout(rect=[0, 0.01, 1, 0.955])
path = OUT / "sfs-explained.png"
fig.savefig(path, dpi=155)
print("wrote", path)
