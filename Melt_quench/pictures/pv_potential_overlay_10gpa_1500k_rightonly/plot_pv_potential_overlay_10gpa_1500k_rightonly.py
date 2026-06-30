import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data" / "csv"
TEMPERATURE = 1500.0

PV_SOURCE = DATA_DIR / "pv_composition_10gpa_1500k.csv"
ENERGY_SOURCE = DATA_DIR / "potential_energy_10gpa_1500k.csv"

POINT_COLORS = {
    "Potential": "#1f4e79",
    "Potential+PV": "#8b1e3f",
}


def read_energy_table(path, value_field):
    data = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if float(row["temperature_K"]) != TEMPERATURE:
                continue
            data[row["formula"]] = float(row[value_field])
    return data


def read_pv_table(path):
    data = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if float(row["temperature_K"]) != TEMPERATURE:
                continue
            data[row["formula"]] = float(row["PV_norm_eV_per_atom"])
    return data


def fit_values(x, y):
    coeff = np.polyfit(x, y, min(5, len(x) - 1))
    fit_x = np.linspace(0.1, 0.9, 300)
    fit_y = np.polyval(coeff, fit_x)
    return fit_x, fit_y


def main():
    pv = read_pv_table(PV_SOURCE)
    pot = read_energy_table(ENERGY_SOURCE, "PotEng_norm_eV_per_atom")

    rows = []
    with open(ENERGY_SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if float(row["temperature_K"]) != TEMPERATURE:
                continue
            rows.append((float(row["x_B"]), row["formula"]))
    rows.sort(key=lambda item: item[0])

    x = np.array([item[0] for item in rows], dtype=float)
    formulas = [item[1] for item in rows]
    pot_values = np.array([pot[formula] for formula in formulas], dtype=float)
    pot_plus_pv = pot_values + np.array([pv[formula] for formula in formulas], dtype=float)

    fig, ax = plt.subplots(figsize=(2.4, 1.6))

    ax.plot(
        x,
        pot_values,
        linestyle="none",
        marker="^",
        markersize=2.8,
        color=POINT_COLORS["Potential"],
        markeredgecolor="black",
        markeredgewidth=0.25,
        clip_on=False,
        zorder=3,
        label="Average potential",
    )
    ax.plot(
        x,
        pot_plus_pv,
        linestyle="none",
        marker="s",
        markersize=2.6,
        color=POINT_COLORS["Potential+PV"],
        markeredgecolor="black",
        markeredgewidth=0.25,
        clip_on=False,
        zorder=3,
        label="Potential + PV",
    )

    fit_x1, fit_y1 = fit_values(x, pot_values)
    fit_x2, fit_y2 = fit_values(x, pot_plus_pv)
    ax.plot(fit_x1, fit_y1, color=POINT_COLORS["Potential"], linestyle="--", linewidth=0.8, zorder=4, label="Potential fit")
    ax.plot(fit_x2, fit_y2, color=POINT_COLORS["Potential+PV"], linestyle="-", linewidth=0.8, zorder=4, label="Potential+PV fit")

    ax.axvline(0.8, color="gray", linestyle="--", linewidth=0.25, zorder=2)
    y_max = max(float(np.max(pot_values)), float(np.max(pot_plus_pv)))
    y_min = min(float(np.min(pot_values)), float(np.min(pot_plus_pv)))
    y_text = y_max - 0.03 * (y_max - y_min)
    ax.text(0.785, y_text - 0.2, r"FeB$_4$", fontsize=6, ha="right", va="center")

    ax.set_xlim(0.1, 0.9)
    ax.set_xticks(np.arange(0.1, 1.0, 0.1))
    ax.set_xlabel(r"$x$ in Fe$_{1-x}$B$_x$", fontsize=6)
    ax.set_ylabel("Average energy [eV/atom]", fontsize=6)
    ax.tick_params(labelsize=6, direction="in", width=0.25, length=2, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.grid(False)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")

    fig.tight_layout(pad=0.35)
    png_path = ROOT / "pv_potential_overlay_10GPa_1500K_4b_rightonly.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(ROOT / "pv_potential_overlay_10GPa_1500K_4b_rightonly.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "pv_potential_overlay_10GPa_1500K_4b_rightonly.eps", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
