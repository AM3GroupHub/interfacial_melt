#!/usr/bin/env python
"""Reproduce manuscript Fig. 4 with weighted OZ fits and 95% bootstrap bands."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

import plot_github_scc_last0p5_3x4 as fitlib


ROOT = Path(__file__).resolve().parent
CSV = ROOT / "data" / "fig3_scc_k.csv"
OUT_DIR = ROOT / "output"
OUT_STEM = OUT_DIR / "fig3_scc_k_T1500_with_fit_uncertainty"
ORDER = [0.6, 0.7, 0.8, 0.9]
FAMILY = {
    0.6: ("#9ecae1", "#2171b5"),
    0.7: ("#a1d99b", "#238b45"),
    0.8: ("#fdd0a2", "#e6550d"),
    0.9: ("#bcbddc", "#6a51a3"),
}
KMAX_FIT = 0.8
TEMPERATURE = 1500

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 22,
        "axes.labelsize": 26,
        "axes.titlesize": 26,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "legend.fontsize": 23,
        "axes.linewidth": 1.4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)


def load_upstream_csv():
    """Load the same exported CSV consumed by the upstream plotting script."""
    data = defaultdict(lambda: {"k": [], "scc": [], "s0": np.nan})
    with open(CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            key = (float(row["x_B"]), int(float(row["P_GPa"])))
            data[key]["k"].append(float(row["k"]))
            data[key]["scc"].append(float(row["S_cc"]))
            if row["S_cc0_M2"]:
                data[key]["s0"] = float(row["S_cc0_M2"])
    for record in data.values():
        record["k"] = np.asarray(record["k"], float)
        record["scc"] = np.asarray(record["scc"], float)
    return data


def selected_fits():
    """Use the same weighted residual bootstrap and seeds as our 3x4 figures."""
    records = fitlib.load_records()
    all_keys = [
        (pressure, temperature, composition)
        for pressure in fitlib.PRESSURES
        for composition in fitlib.COMPOSITIONS
        for temperature in fitlib.TEMPERATURES
    ]
    seeds = np.random.SeedSequence(fitlib.RNG_SEED).generate_state(len(all_keys))
    seed_by_key = {key: int(seed) for key, seed in zip(all_keys, seeds)}
    fits = {}
    for pressure in fitlib.PRESSURES:
        for composition in ORDER:
            key = (pressure, TEMPERATURE, composition)
            fits[(composition, pressure)] = fitlib.bootstrap_fit(
                records[key], np.random.default_rng(seed_by_key[key]), KMAX_FIT
            )
            print(f"[fit] P={pressure:2d} T={TEMPERATURE} x_B={composition:.1f}")
    return fits


def thin(uk, k_lo=1.0, step=3):
    mask = np.ones(len(uk), bool)
    high = np.where(uk > k_lo)[0]
    mask[high] = np.arange(len(high)) % step == 0
    return mask


def mark_intercept(ax, value, color, marker, size, linewidth, zorder, slot, ytop):
    if not np.isfinite(value):
        return
    if value <= ytop:
        ax.scatter(
            [0.0], [value], marker=marker, s=size, color=color,
            linewidths=linewidth, zorder=zorder,
        )
    else:
        yclipped = ytop * (0.95 - 0.12 * slot)
        ax.scatter(
            [0.0], [yclipped], marker=marker, s=size, color=color,
            linewidths=linewidth, zorder=zorder, clip_on=False,
        )
        ax.text(
            0.12, yclipped, f"{value:.2g}", color=color, fontsize=14,
            fontweight="bold", va="center",
        )


def add_fit(ax, result, color):
    ax.fill_between(
        result["kfit"], result["ylo"], result["yhi"],
        color=color, alpha=0.14, linewidth=0, zorder=1,
    )
    ax.plot(result["kfit"], result["yfit"], color=color, linewidth=2.1, zorder=2)


def main():
    data = load_upstream_csv()
    fits = selected_fits()
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True, sharey="row")

    for panel, (ax, composition) in enumerate(zip(axes.flat, ORDER)):
        light, dark = FAMILY[composition]
        ytop = 1.8 if panel // 2 == 1 else 0.6
        record0 = data[(composition, 0)]
        record10 = data[(composition, 10)]
        mask0, mask10 = thin(record0["k"]), thin(record10["k"])

        add_fit(ax, fits[(composition, 10)], dark)
        add_fit(ax, fits[(composition, 0)], light)
        ax.scatter(
            record10["k"][mask10], record10["scc"][mask10], color=dark,
            marker="^", s=190, edgecolor="none", alpha=0.9, zorder=3,
        )
        ax.scatter(
            record0["k"][mask0], record0["scc"][mask0], facecolor="none",
            edgecolor=light, marker="o", s=205, linewidth=2.2, alpha=0.9, zorder=4,
        )

        mark_intercept(ax, fits[(composition, 10)]["a"], "#67000d", "+", 340, 3.2, 6, 1, ytop)
        if composition != 0.9:
            mark_intercept(ax, fits[(composition, 0)]["a"], "#e31a1c", "x", 210, 3.0, 7, 0, ytop)

        ax.set_title(rf"$x_B = {composition:.1f}$", fontsize=27, fontweight="bold", pad=14)
        ax.set_xlim(-0.08, 2.55)
        ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        if panel // 2 == 1:
            ax.set_ylim(0, 1.8)
            ax.set_yticks([0.0, 0.6, 1.2, 1.8])
        else:
            ax.set_ylim(0, 0.6)
            ax.set_yticks([0.0, 0.2, 0.4, 0.6])

        if panel == 2:
            handles = [
                Line2D([], [], marker="x", linestyle="", color="#e31a1c", markeredgewidth=3.0,
                       markersize=15, label=r"$S_{cc}(0)$, 0 GPa"),
                Line2D([], [], marker="+", linestyle="", color="#67000d", markeredgewidth=3.0,
                       markersize=17, label=r"$S_{cc}(0)$, 10 GPa"),
            ]
            ax.legend(handles=handles, loc="upper right", fontsize=22, handletextpad=0.4, labelspacing=0.45)
        if panel == 3:
            handles = [
                Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="0.25",
                       markeredgewidth=2.0, markersize=17, label="0 GPa"),
                Line2D([], [], marker="^", linestyle="", color="0.25", markeredgecolor="none",
                       markersize=18, label="10 GPa"),
                Line2D([], [], color="0.25", linewidth=2.1, label="OZ fit"),
                Patch(facecolor="0.5", alpha=0.14, edgecolor="none", label="95% fit CI"),
            ]
            ax.legend(handles=handles, loc="upper right", fontsize=16.5, handletextpad=0.4, labelspacing=0.35)

    for ax in axes[-1]:
        ax.set_xlabel(r"$|k|$ [Å$^{-1}$]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$S_{cc}(|k|)$")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.10, wspace=0.11, hspace=0.30)
    panel_label_positions = [(0.025, 0.995), (0.555, 0.995), (0.025, 0.515), (0.555, 0.515)]
    for panel, (xpos, ypos) in enumerate(panel_label_positions):
        fig.text(
            xpos,
            ypos,
            f"({chr(ord('a') + panel)})",
            fontsize=plt.rcParams["axes.labelsize"],
            fontweight="normal",
            ha="left",
            va="center",
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {OUT_STEM.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
