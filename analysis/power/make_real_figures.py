"""Figures from the attack results on the 12 loaded beacons.

Reads ``attack_<POP>_<set>.json`` (written by ``run_attack_vm.py``) from
``--attack-dir`` and writes two PNGs to ``--out-dir``:

* ``attack_budget.png``: queries needed for a 5 per cent false-positive rate,
  per beacon, with the bootstrap interval, for the beacon's own controls and
  for same-continent outsiders.
* ``attack_fpr_curves.png``: false-positive rate against the number of queried
  heterozygous sites, same-continent outsiders, one line per beacon.

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


def curve_figure(attack_dir, out_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    for pop in ORDER:
        d = load(attack_dir, pop, "continent")
        points = sorted((int(k), v) for k, v in d["fpr_by_budget"].items())
        ax.plot([p[0] for p in points], [p[1] for p in points], marker="o",
                markersize=3, linewidth=1.3, color=colour(pop),
                linestyle="-" if pop in AFRICAN else "--",
                label=pop)
    ax.axhline(0.05, color="black", linewidth=0.8)
    ax.text(300, 0.057, "5%", ha="right", fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("queried heterozygous sites (n)")
    ax.set_ylabel("false-positive rate (outsiders answering yes to all n)")
    ax.set_title("False-positive rate against queries, same-continent outsiders")
    ax.legend(frameon=False, ncol=3, fontsize=8, title="solid: African, dashed: European")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--attack-dir", required=True)
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    budget_figure(a.attack_dir, os.path.join(a.out_dir, "attack_budget.png"))
    curve_figure(a.attack_dir, os.path.join(a.out_dir, "attack_fpr_curves.png"))
    print("wrote", sorted(os.listdir(a.out_dir)))


if __name__ == "__main__":
    main()
