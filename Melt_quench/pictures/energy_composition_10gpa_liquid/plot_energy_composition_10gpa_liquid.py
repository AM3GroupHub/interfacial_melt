import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "csv" / "normalized_energy_10GPa_liquid.csv"
PLOT_TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]
PRESSURE_LABEL = "10 GPa"
ENERGY_COLUMNS = [
    ("PotEng_norm_eV_per_atom", "PotEng_norm_err_eV_per_atom", "Potential energy"),
    ("Enthalpy_norm_eV_per_atom", "Enthalpy_norm_err_eV_per_atom", "Enthalpy"),
    ("Free_norm_eV_per_atom", "Free_norm_err_eV_per_atom", "Free energy"),
]


def read_rows():
    rows = []
    with open(SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature not in PLOT_TEMPERATURES:
                continue
            normalized = {key: float(value) if key not in {"formula"} else value for key, value in row.items()}
            rows.append(normalized)
    return rows


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    return xs, np.polyval(coeff, xs)


def plot_temperature(rows, temperature):
    temp_rows = [row for row in rows if row["temperature_K"] == temperature]
    temp_rows.sort(key=lambda row: row["x_B"])
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8), sharex=True)
    color = plt.cm.plasma(PLOT_TEMPERATURES.index(temperature) / max(len(PLOT_TEMPERATURES) - 1, 1))
    x = np.array([row["x_B"] for row in temp_rows], dtype=float)
    for ax, (y_key, err_key, title) in zip(axes, ENERGY_COLUMNS):
        y = np.array([row[y_key] for row in temp_rows], dtype=float)
        yerr = np.array([row[err_key] for row in temp_rows], dtype=float)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            color=color,
            ecolor="black",
            elinewidth=1.8,
            capthick=1.8,
            capsize=5,
            barsabove=True,
            zorder=5,
        )
        fit_x, fit_y = fit_curve(x, y)
        ax.plot(fit_x, fit_y, "--", color="black", linewidth=2.0, alpha=0.85)
        ax.set_title(title, fontsize=16)
        ax.set_xlabel(r"B fraction, $x_B = B / (Fe + B)$", fontsize=16)
        ax.set_ylabel("Normalized energy (eV/atom)", fontsize=16)
        ax.tick_params(labelsize=16)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f"{PRESSURE_LABEL}, {temperature:.0f} K", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    png_path = ROOT / f"energy_composition_10GPa_liquid_{temperature:.0f}K_deg5.png"
    pdf_path = ROOT / f"energy_composition_10GPa_liquid_{temperature:.0f}K_deg5.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def main():
    rows = read_rows()
    outputs = []
    for temperature in PLOT_TEMPERATURES:
        outputs.extend(plot_temperature(rows, temperature))
    for path in outputs:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
