import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "csv" / "pv_composition_10GPa_liquid.csv"
PRESSURE_LABEL = "10 GPa"
PLOT_TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]


def read_rows():
    rows = []
    with open(SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature not in PLOT_TEMPERATURES:
                continue
            rows.append(
                {
                    "temperature_K": temperature,
                    "x_B": float(row["x_B"]),
                    "PV_norm_eV_per_atom": float(row["PV_norm_eV_per_atom"]),
                    "PV_norm_err_eV_per_atom": float(row["PV_norm_err_eV_per_atom"]),
                }
            )
    return rows


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    return xs, np.polyval(coeff, xs)


def main():
    rows = read_rows()
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(PLOT_TEMPERATURES)))
    for color, temperature in zip(colors, PLOT_TEMPERATURES):
        temp_rows = [row for row in rows if row["temperature_K"] == temperature]
        x = np.array([row["x_B"] for row in temp_rows], dtype=float)
        y = np.array([row["PV_norm_eV_per_atom"] for row in temp_rows], dtype=float)
        yerr = np.array([row["PV_norm_err_eV_per_atom"] for row in temp_rows], dtype=float)
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
            label=f"{temperature:.0f} K",
        )
        fit_x, fit_y = fit_curve(x, y)
        ax.plot(fit_x, fit_y, "--", color=color, linewidth=1.8, alpha=0.85)
    ax.set_xlabel(r"B fraction, $x_B = B / (Fe + B)$", fontsize=16)
    ax.set_ylabel("Endpoint-normalized PV term (eV/atom)", fontsize=16)
    ax.set_title(f"{PRESSURE_LABEL}: endpoint-normalized PV term vs composition", fontsize=16)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=16, frameon=True)
    fig.tight_layout()
    png_path = ROOT / "pv_composition_10GPa_liquid.png"
    pdf_path = ROOT / "pv_composition_10GPa_liquid.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
