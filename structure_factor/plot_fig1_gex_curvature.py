#!/usr/bin/env python
"""Fig. 1(c,d) bottom — excess curvature (1/k_BT) d^2 g_ex/dx^2 = (Gamma-1)/(x_Fe x_B)
vs x_B, four (T,P) panels. Reads data/csv/fig1_gex_curvature.csv (CSV + matplotlib only)."""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "data" / "csv" / "fig1_gex_curvature.csv"
FIGDIR = ROOT / "figures"
PANELS = [(1500, 0), (1800, 0), (1500, 10), (1800, 10)]
cmap = plt.get_cmap("viridis"); norm = Normalize(vmin=1200, vmax=2600)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.titlesize": 12,
                     "axes.labelsize": 13, "xtick.labelsize": 10, "ytick.labelsize": 11,
                     "axes.linewidth": 1.0, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})


def load():
    d = defaultdict(lambda: {"xb": [], "y": []})
    with open(CSV) as fh:
        for row in csv.DictReader(fh):
            key = (int(float(row["T_K"])), int(float(row["P_GPa"])))
            d[key]["xb"].append(float(row["x_B"])); d[key]["y"].append(float(row["d2gex_over_kBT"]))
    return d


def main():
    d = load(); ylim = (-13.0, 21.0)
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.8), sharey=True, constrained_layout=True)
    for k, (ax, (T, P)) in enumerate(zip(axs, PANELS)):
        xb = np.array(d[(T, P)]["xb"]); y = np.array(d[(T, P)]["y"])
        order = np.argsort(xb); xb, y = xb[order], y[order]
        ax.axhspan(ylim[0], 0.0, color="navy", alpha=0.05)
        ax.axhline(0.0, color="k", lw=1.0, ls=(0, (6, 4)))
        ax.plot(xb, y, "o-", color=cmap(norm(T)), lw=1.8, ms=6, mec="white", mew=0.7)
        ax.set_xlabel(r"$x_B$"); ax.set_xlim(0.05, 0.95); ax.set_xticks(np.arange(0.1, 0.91, 0.2))
        ax.set_ylim(*ylim); ax.minorticks_on()
        ax.set_title(f"$T = {int(T)}$ K,  $P = {int(P)}$ GPa")
        ax.text(0.04, 0.95, f"({chr(97 + k)})", transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")
    axs[0].set_ylabel(r"$(1/k_B T)\,\partial^2 g_{ex}/\partial x^2$" "\n[S(k) route]")
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / "fig1_gex_curvature.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
