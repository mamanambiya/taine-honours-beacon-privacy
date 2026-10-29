"""Figures from the attack results on the 12 loaded beacons.

Reads ``attack_<POP>_<set>.json`` (written by ``run_attack_vm.py``) from
``--attack-dir`` and writes two PNGs to ``--out-dir``:

* ``attack_budget.png``: queries needed for a 5 per cent false-positive rate,
  per beacon, with the bootstrap interval, for the beacon's own controls and
  for same-continent outsiders.
* ``attack_fpr_curves.png``: false-positive rate against the number of queried
  heterozygous sites, same-continent outsiders, one line per beacon.
* ``attack_absent_fraction.png``: for each beacon, the share of an outsider's
  heterozygous sites that are absent from the beacon (box plot over outsiders).
* ``attack_budget_vs_absent.png``: queries for 5 per cent against the mean of
  that share, with the single-rate approximation ``ln 0.05 / ln(1 - a)``.
* ``attack_budget_vs_size.png``: queries for 5 per cent against the number of
  sites the beacon lists.
* ``attack_budget_vs_delta.png`` (only with ``--binomial``): queries for 5 per
  cent when each answer may be wrong with probability delta, from
  ``binomial_delta.py``.
* ``attack_power_curves.png`` (only with ``--binomial``): empirical power against
  queries with the false-positive rate held at 5 per cent, for several delta.

Needs matplotlib. Usage::

    python3 make_real_figures.py --attack-dir results/tables/attack --out-dir docs/figures
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AFRICAN = ["ESN", "GWD", "LWK", "MSL", "YRI", "AFR"]
EUROPEAN = ["CEU", "FIN", "GBR", "IBS", "TSI", "EUR"]
ORDER = AFRICAN + EUROPEAN
COLOUR = {"AFR": "#c4621a", "EUR": "#2a6f9e"}


def colour(pop):
    return COLOUR["AFR"] if pop in AFRICAN else COLOUR["EUR"]


def load(attack_dir, pop, which):
    with open(os.path.join(attack_dir, f"attack_{pop}_{which}.json")) as handle:
        return json.load(handle)


def budget_figure(attack_dir, out_path):
    fig, ax = plt.subplots(figsize=(10, 4.6))
    width = 0.38
    for i, pop in enumerate(ORDER):
        for j, (which, hatch) in enumerate((("controls", "//"), ("continent", None))):
            d = load(attack_dir, pop, which)
            b = d["budget_for_fpr_alpha"]
            lo, hi = d["budget_95ci_bootstrap_over_outsiders"]
            x = i + (j - 0.5) * width
            ax.bar(x, b, width, color=colour(pop), hatch=hatch,
                   alpha=0.55 if which == "controls" else 1.0,
                   edgecolor="white", linewidth=0.6)
            ax.errorbar(x, b, yerr=[[b - lo], [hi - b]], color="black", capsize=2,
                        linewidth=0.8)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER)
    ax.set_ylabel("queries for a 5% false-positive rate")
    ax.set_title("Attack budget on the 12 loaded beacons (65 members each)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLOUR["AFR"]),
        plt.Rectangle((0, 0), 1, 1, color=COLOUR["EUR"]),
        plt.Rectangle((0, 0), 1, 1, facecolor="grey", alpha=0.55, hatch="//",
                      edgecolor="white"),
        plt.Rectangle((0, 0), 1, 1, color="grey"),
    ]
    ax.legend(handles, ["African beacons", "European beacons",
                        "own 20 controls (hatched)", "same-continent outsiders"],
              frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
OFFSETS = [(6, 4), (6, -13), (-28, 6)]  # keep neighbouring labels apart


def curve_figure(attack_dir, out_path):
    """Two panels (African, European), one colour and marker per beacon."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, group, title in ((axes[0], AFRICAN, "African beacons"),
                             (axes[1], EUROPEAN, "European beacons")):
        for i, pop in enumerate(group):
            d = load(attack_dir, pop, "continent")
            points = sorted((int(k), v) for k, v in d["fpr_by_budget"].items())
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            ax.plot(xs, ys, marker=MARKERS[i], markersize=5, linewidth=1.5,
                    color=PALETTE[i], label=f"{pop}  ({d['budget_for_fpr_alpha']})")
        ax.axhline(0.05, color="black", linewidth=0.8, linestyle=":")
        ax.set_yscale("log")
        ax.set_xlabel("queried heterozygous sites (n)")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9, title="beacon (queries at 5%)",
                  loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("false-positive rate (outsiders answering yes to all n)")
    axes[0].text(2, 0.057, "5%", ha="left", fontsize=8)
    fig.suptitle("False-positive rate against queries, same-continent outsiders")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _fractions(attack_dir, pop):
    d = load(attack_dir, pop, "continent")
    return d, [absent / het for het, absent in d["outsider_counts"]]


def absent_fraction_figure(attack_dir, out_path):
    fig, ax = plt.subplots(figsize=(10, 4.6))
    data = [_fractions(attack_dir, pop)[1] for pop in ORDER]
    box = ax.boxplot(data, tick_labels=ORDER, patch_artist=True, showfliers=False,
                     medianprops={"color": "black"})
    for patch, pop in zip(box["boxes"], ORDER):
        patch.set_facecolor(colour(pop))
        patch.set_alpha(0.8)
    ax.set_ylabel("share of an outsider's heterozygous sites\nabsent from the beacon")
    ax.set_title("Why the budgets differ: outsiders' sites missing from each beacon")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_figure(attack_dir, out_path, x_of, xlabel, title, curve=False):
    import math

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for pop in ORDER:
        d, _ = _fractions(attack_dir, pop)
        x = x_of(d)
        y = d["budget_for_fpr_alpha"]
        lo, hi = d["budget_95ci_bootstrap_over_outsiders"]
        ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt="o", color=colour(pop),
                    capsize=2, markersize=7)
        ax.annotate(pop, (x, y), textcoords="offset points",
                    xytext=OFFSETS[ORDER.index(pop) % 3], fontsize=9)
    if curve:
        xs = [i / 1000 for i in range(15, 50)]
        ax.plot(xs, [math.log(0.05) / math.log(1 - v) for v in xs], color="grey",
                linestyle="--", linewidth=1,
                label="single-rate approximation  ln 0.05 / ln(1 - a)")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("queries for a 5% false-positive rate")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def delta_figure(binomial_path, out_path):
    """Budget against delta, same-continent outsiders, one line per beacon."""
    with open(binomial_path) as handle:
        rows = json.load(handle)["rows"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, group, title in ((axes[0], AFRICAN, "African beacons"),
                             (axes[1], EUROPEAN, "European beacons")):
        for i, pop in enumerate(group):
            by_delta = rows[f"{pop}_continent"]["budget_by_delta"]
            deltas = sorted(by_delta, key=float)
            pts = [(k, by_delta[k]) for k in deltas if by_delta[k] is not None]
            ax.plot([deltas.index(k) for k, _ in pts], [v for _, v in pts],
                    marker=MARKERS[i], markersize=5, linewidth=1.5,
                    color=PALETTE[i], label=pop)
        ax.set_xticks(range(len(deltas)))
        ax.set_xticklabels(deltas)
        ax.set_yscale("log")
        ax.set_xlabel("mismatch rate delta")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("queries for a 5% false-positive rate")
    fig.suptitle("Budget when each answer may be wrong (member acceptance 95%); "
                 "a missing point is above 2000 queries")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def power_figure(attack_dir, out_path):
    """Power against queries at 5% false-positive rate: one panel per delta."""
    import binomial_delta as bd

    grid = list(range(20, 1201, 20))
    pops = ["YRI", "AFR", "CEU", "EUR"]
    deltas = [0.0, 0.001, 0.01]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharey=True)
    for ax, delta in zip(axes, deltas):
        for i, pop in enumerate(pops):
            outs = [tuple(x) for x in load(attack_dir, pop, "continent")["outsider_counts"]]
            ys = [bd.power_at(outs, n, delta, 0.05) for n in grid]
            ax.plot(grid, ys, color=PALETTE[i], linewidth=1.6, label=pop)
        ax.set_title(f"delta = {delta:g}")
        ax.set_xlabel("queried heterozygous sites (n)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("empirical power (false-positive rate held at 5%)")
    axes[0].legend(frameon=False, title="beacon")
    fig.suptitle("Power against queries, same-continent outsiders")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--attack-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--binomial", help="binomial_budgets.json from binomial_delta.py")
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    budget_figure(a.attack_dir, os.path.join(a.out_dir, "attack_budget.png"))
    curve_figure(a.attack_dir, os.path.join(a.out_dir, "attack_fpr_curves.png"))
    absent_fraction_figure(a.attack_dir,
                           os.path.join(a.out_dir, "attack_absent_fraction.png"))
    scatter_figure(a.attack_dir, os.path.join(a.out_dir, "attack_budget_vs_absent.png"),
                   lambda d: d["outsider_absent_fraction_mean"],
                   "mean share of an outsider's sites absent from the beacon",
                   "Budget against the share of missing sites", curve=True)
    scatter_figure(a.attack_dir, os.path.join(a.out_dir, "attack_budget_vs_size.png"),
                   lambda d: d["beacon_sites"] / 1e6,
                   "sites the beacon lists (millions)",
                   "Budget against beacon size")
    if a.binomial:
        delta_figure(a.binomial, os.path.join(a.out_dir, "attack_budget_vs_delta.png"))
        power_figure(a.attack_dir, os.path.join(a.out_dir, "attack_power_curves.png"))
    print("wrote", sorted(os.listdir(a.out_dir)))


if __name__ == "__main__":
    main()
