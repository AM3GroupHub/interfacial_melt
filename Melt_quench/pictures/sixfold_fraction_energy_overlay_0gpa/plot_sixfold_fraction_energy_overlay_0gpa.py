import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SIXFOLD_SUMMARY = ROOT / "data" / "csv" / "sixfold_fraction_summary.csv"
ENERGY_0GPA = ROOT / "data" / "csv" / "normalized_energy_0GPa_liquid.csv"

TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]
COMPOSITIONS = [
    ("Fe160B240", "Fe0.4B0.6", 0.60),
    ("Fe120B280", "Fe0.3B0.7", 0.70),
    ("Fe135B405", "Fe0.25B0.75", 0.75),
    ("Fe72B288", "Fe0.2B0.8", 0.80),
    ("Fe81B459", "Fe0.15B0.85", 0.85),
    ("Fe54B486", "Fe0.1B0.9", 0.90),
]
PRESSURES = ["0 GPa", "10 GPa"]
BAR_COLORS = {"0 GPa": "#2b6cb0", "10 GPa": "#c53030"}
BAR_HATCH = {"0 GPa": "", "10 GPa": "//"}
LINE_STYLE = {
    "0 GPa": {
        "color": "#1f4e79",
        "linestyle": "--",
        "label": "0 GPa potential energy fit",
    },
}


def read_fraction_summary(path):
    data = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["composition"], float(row["temperature_K"]), row["pressure"])
            data[key] = {
                "mean": float(row["mean_fraction"]),
                "sem": float(row["sem"]),
            }
    return data


def read_csv_by_key(path, value_field):
    data = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["formula"], float(row["temperature_K"]))
            data[key] = float(row[value_field])
    return data


def read_csv_with_x(path, value_field):
    data = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["formula"], float(row["temperature_K"]))
            data[key] = (float(row["x_B"]), float(row[value_field]))
    return data


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    return np.polyfit(x, y, degree)


def subscript_label(label):
    return label.replace("Fe", r"Fe$_{").replace("B", r"}$B$_{") + r"}$"


def plot_temperature(temp, fraction_data, pot_0gpa, pot_0gpa_err, fit_source):
    x = np.array([x_b for _, _, x_b in COMPOSITIONS], dtype=float)
    width = 0.014
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax2 = ax.twinx()

    for idx, pressure in enumerate(PRESSURES):
        offset = (-0.5 + idx) * width
        means = []
        errs = []
        for composition, _, _ in COMPOSITIONS:
            stats = fraction_data[(composition, temp, pressure)]
            means.append(stats["mean"])
            errs.append(stats["sem"])
        ax.bar(
            x + offset,
            means,
            width=width,
            color=BAR_COLORS[pressure],
            edgecolor="black",
            linewidth=0.6,
            hatch=BAR_HATCH[pressure],
            alpha=0.85,
            label=pressure,
            zorder=3,
        )
        ax.errorbar(
            x + offset,
            means,
            yerr=errs,
            fmt="none",
            ecolor="black",
            elinewidth=0.8,
            capsize=3,
            capthick=0.8,
            zorder=4,
        )

    full_x = []
    full_y = []
    for (_, row_temp), (x_b, value) in fit_source.items():
        if row_temp == temp:
            full_x.append(x_b)
            full_y.append(value)
    pot0_y = [pot_0gpa[(composition, temp)] for composition, _, _ in COMPOSITIONS]
    pot0_err = [pot_0gpa_err[(composition, temp)] for composition, _, _ in COMPOSITIONS]
    coeff = fit_curve(np.array(full_x, dtype=float), np.array(full_y, dtype=float))
    fit_x = np.linspace(float(np.min(x)), float(np.max(x)), 300)
    fit_y = np.polyval(coeff, fit_x)

    ax2.plot(
        x,
        pot0_y,
        linestyle="none",
        marker="o",
        color=LINE_STYLE["0 GPa"]["color"],
        markeredgecolor="black",
        markeredgewidth=0.6,
        markersize=4.5,
        zorder=7,
    )
    ax2.errorbar(
        x,
        pot0_y,
        yerr=pot0_err,
        fmt="none",
        ecolor="black",
        elinewidth=0.8,
        capsize=3,
        capthick=0.8,
        zorder=8,
    )
    ax2.plot(
        fit_x,
        fit_y,
        color=LINE_STYLE["0 GPa"]["color"],
        linestyle=LINE_STYLE["0 GPa"]["linestyle"],
        linewidth=1.8,
        zorder=6,
        label=LINE_STYLE["0 GPa"]["label"],
    )

    right_min = float(np.min(fit_y)) - 0.05
    right_max = float(np.max(fit_y)) + 0.05
    ax.set_xticks(x)
    ax.set_xticklabels([subscript_label(label) for _, label, _ in COMPOSITIONS], rotation=25, ha="right", fontsize=16)
    ax.set_xlim(float(np.min(x)) - 0.03, float(np.max(x)) + 0.03)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"Fraction of six-coordinated B", fontsize=16)
    ax.set_title(f"{temp:.0f} K", fontsize=16)
    ax.tick_params(labelsize=16)
    ax.grid(True, axis="y", alpha=0.25, zorder=0)
    ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=16)

    ax2.tick_params(labelsize=16)
    ax2.legend(frameon=False, loc="upper right", fontsize=16)
    ax2.set_ylabel("Average potential energy (eV/atom)", fontsize=16)
    ax2.set_ylim(right_min, right_max)

    fig.tight_layout()
    png_path = ROOT / f"sixfold_fraction_energy_overlay_{temp:.0f}K.png"
    pdf_path = ROOT / f"sixfold_fraction_energy_overlay_{temp:.0f}K.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def main():
    fraction_data = read_fraction_summary(SIXFOLD_SUMMARY)
    pot_0gpa = read_csv_by_key(ENERGY_0GPA, "PotEng_norm_eV_per_atom")
    pot_0gpa_err = read_csv_by_key(ENERGY_0GPA, "PotEng_norm_err_eV_per_atom")
    fit_source = read_csv_with_x(ENERGY_0GPA, "PotEng_norm_eV_per_atom")
    outputs = []
    for temp in TEMPERATURES:
        outputs.extend(plot_temperature(temp, fraction_data, pot_0gpa, pot_0gpa_err, fit_source))
    for path in outputs:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
