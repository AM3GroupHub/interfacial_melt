#!/usr/bin/env python
"""Fig. 3 — reciprocal S_cc(k) for x_B in {0.6,0.7,0.8,0.9}, 0 vs 10 GPa, 1500 K.
Reads data/csv/fig3_scc_k.csv (CSV + matplotlib only)."""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "data" / "csv" / "fig3_scc_k.csv"
FIGDIR = ROOT / "figures"
ORDER = [0.6, 0.7, 0.8, 0.9]
FAMILY = {0.6: ("#9ecae1", "#2171b5"), 0.7: ("#a1d99b", "#238b45"),
          0.8: ("#fdd0a2", "#e6550d"), 0.9: ("#bcbddc", "#6a51a3")}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 22,
                     "axes.labelsize": 26, "axes.titlesize": 26,
                     "xtick.labelsize": 22, "ytick.labelsize": 22, "legend.fontsize": 23,
                     "axes.linewidth": 1.4, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})


def load():
    """(x_B,P) -> dict(k=array, scc=array, s0=float)."""
    d = defaultdict(lambda: {"k": [], "scc": [], "s0": np.nan})
    with open(CSV) as fh:
        for row in csv.DictReader(fh):
            key = (float(row["x_B"]), int(float(row["P_GPa"])))
            d[key]["k"].append(float(row["k"])); d[key]["scc"].append(float(row["S_cc"]))
            if row["S_cc0_M2"]:
                d[key]["s0"] = float(row["S_cc0_M2"])
    for v in d.values():
        v["k"] = np.array(v["k"]); v["scc"] = np.array(v["scc"])
    return d


def thin(uk, k_lo=1.0, step=3):
    m = np.ones(len(uk), bool); hi = np.where(uk > k_lo)[0]
    m[hi] = (np.arange(len(hi)) % step == 0); return m


def main():
    d = load()
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True, sharey="row")
    for k, (ax, xt) in enumerate(zip(axs.flat, ORDER)):
        light, dark = FAMILY[xt]; ytop = 1.8 if k // 2 == 1 else 0.6
        r0, r10 = d[(xt, 0)], d[(xt, 10)]
        m0, m10 = thin(r0["k"]), thin(r10["k"])
        ax.scatter(r10["k"][m10], r10["scc"][m10], color=dark, marker="^", s=190,
                   edgecolor="none", alpha=0.9, zorder=3)
        ax.scatter(r0["k"][m0], r0["scc"][m0], facecolor="none", edgecolor=light,
                   marker="o", s=205, lw=2.2, alpha=0.9, zorder=4)

        def mark(val, color, marker, s, lw, z, slot):
            if not np.isfinite(val):
                return
            if val <= ytop:
                ax.scatter([0.0], [val], marker=marker, s=s, color=color, linewidths=lw, zorder=z)
            else:
                yc = ytop * (0.95 - 0.12 * slot)
                ax.scatter([0.0], [yc], marker=marker, s=s, color=color, linewidths=lw,
                           zorder=z, clip_on=False)
                ax.text(0.12, yc, f"{val:.2g}", color=color, fontsize=14, fontweight="bold", va="center")
        mark(r10["s0"], "#67000d", "+", 340, 3.2, 6, 1)
        if xt != 0.9:
            mark(r0["s0"], "#e31a1c", "x", 210, 3.0, 7, 0)
        ax.set_title(f"$x_B = {xt:.1f}$", fontsize=27, fontweight="bold", pad=14)
        ax.set_xlim(-0.08, 2.55); ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        if k // 2 == 1:
            ax.set_ylim(0, 1.8); ax.set_yticks([0.0, 0.6, 1.2, 1.8])
        else:
            ax.set_ylim(0, 0.6); ax.set_yticks([0.0, 0.2, 0.4, 0.6])
        if k == 2:
            h = [Line2D([], [], marker="x", ls="", color="#e31a1c", mew=3.0, ms=15, label=r"$S_{cc}(0)$, 0 GPa"),
                 Line2D([], [], marker="+", ls="", color="#67000d", mew=3.0, ms=17, label=r"$S_{cc}(0)$, 10 GPa")]
            ax.legend(handles=h, loc="upper right", fontsize=22, handletextpad=0.4, labelspacing=0.45)
        if k == 3:
            h = [Line2D([], [], marker="o", ls="", mfc="none", mec="0.25", mew=2.0, ms=17, label="0 GPa"),
                 Line2D([], [], marker="^", ls="", color="0.25", mec="none", ms=18, label="10 GPa")]
            ax.legend(handles=h, loc="upper right", fontsize=22, handletextpad=0.4, labelspacing=0.45)
    for ax in axs[-1]:
        ax.set_xlabel(r"$|k|$ [Å$^{-1}$]")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"$S_{cc}(|k|)$")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.10, wspace=0.11, hspace=0.30)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / "fig3_scc_k_T1500.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
