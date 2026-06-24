#!/usr/bin/env python
"""SM — per-composition reciprocal S_cc(k) with the M2 OZ-like fit (red dashed) and
S_cc(0) (red x), one figure per (T,P). Reads data/csv/sm_state_diagnostics.csv."""
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
CSV = ROOT / "data" / "csv" / "sm_state_diagnostics.csv"
FIGDIR = ROOT / "figures" / "state_diagnostics"
SHELL = "#3b78c2"; RED = "#e31a1c"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 19, "axes.labelsize": 23,
                     "axes.titlesize": 23, "xtick.labelsize": 18, "ytick.labelsize": 18,
                     "legend.fontsize": 19, "axes.linewidth": 1.3, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})


def load():
    """(T,P) -> {x_B -> dict(k, scc, a, b, s0)} sorted by x_B."""
    by = defaultdict(lambda: defaultdict(lambda: {"k": [], "scc": [], "a": np.nan, "b": np.nan, "s0": np.nan}))
    with open(CSV) as fh:
        for row in csv.DictReader(fh):
            tp = (int(float(row["T_K"])), int(float(row["P_GPa"]))); xb = float(row["x_B"])
            rec = by[tp][xb]; rec["k"].append(float(row["k"])); rec["scc"].append(float(row["S_cc"]))
            if row["a"]:
                rec["a"] = float(row["a"]); rec["b"] = float(row["b"]); rec["s0"] = float(row["S_cc0"])
    return by


def make_figure(T, P, comps, ncol=6):
    items = sorted(comps.items())                       # (x_B, rec)
    nrow = int(np.ceil((len(items) + 1) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 3.1 * nrow), squeeze=False)
    axes = list(axs.flat)
    for ax in axes:
        ax.set_visible(False)
    pool = [np.nanmax(np.array(r["scc"])[np.array(r["k"]) <= 2.5]) for _, r in items]
    top = (max([v for v in pool if np.isfinite(v)] or [0.6])) * 1.06
    for i, (xb, r) in enumerate(items):
        ax = axes[i]; ax.set_visible(True)
        k = np.array(r["k"]); scc = np.array(r["scc"])
        ax.scatter(k, scc, s=62, color=SHELL, edgecolor="none", alpha=0.5, zorder=3)
        x_fe = 1.0 - xb
        ax.axhline(x_fe * xb, color="0.72", lw=1.6, zorder=2)
        ax.set_ylim(0, top); ax.set_xlim(-0.06, 2.55); ax.set_xticks([0.0, 1.0, 2.0])
        if np.isfinite(r["a"]) and np.isfinite(r["b"]):
            kl = np.linspace(0.0, 0.8, 80)
            ax.plot(kl, r["a"] / (1.0 + r["b"] * kl ** 2), ls="--", color=RED, lw=2.4, zorder=4)
        if np.isfinite(r["s0"]) and 0.0 <= r["s0"] <= top:
            ax.scatter([0], [r["s0"]], marker="x", s=170, color=RED, linewidths=2.8, zorder=6, clip_on=False)
        ax.set_title(f"$x_B = {xb:.2f}$", fontweight="bold", pad=8)
        if i % ncol == 0:
            ax.set_ylabel(r"$S_{cc}(|k|)$")
        if i // ncol == nrow - 1:
            ax.set_xlabel(r"$|k|$ [Å$^{-1}$]")
    leg = axes[len(items)]; leg.set_visible(True); leg.axis("off")
    leg.legend(handles=[
        Line2D([], [], marker="o", ls="", color=SHELL, alpha=0.5, mec="none", ms=14, label=r"$S_{cc}(k)$"),
        Line2D([], [], marker="x", ls="", color=RED, mew=2.8, ms=16, label=r"$S_{cc}(0)$"),
        Line2D([], [], ls="--", color=RED, lw=2.4, label="fit"),
        Line2D([], [], ls="-", color="0.72", lw=1.6, label=r"ideal $x_{Fe}x_B$"),
    ], loc="center", frameon=False, labelspacing=0.9, handletextpad=0.6)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"sm_state_diagnostics_T{int(T)}_P{int(P)}.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig); print(f"[fig] {out}")


def main():
    by = load()
    for (T, P), comps in sorted(by.items()):
        make_figure(T, P, comps)


if __name__ == "__main__":
    main()
