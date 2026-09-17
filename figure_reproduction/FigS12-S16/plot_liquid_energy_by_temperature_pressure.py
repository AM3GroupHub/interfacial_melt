import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parent
SOURCES = {
    "0 GPa": ROOT / "data" / "normalized_energy_0GPa_liquid.csv",
    "10 GPa": ROOT / "data" / "normalized_energy_10GPa_liquid.csv",
}
TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]
ENERGY_COLUMNS = [
    ("PotEng_norm_eV_per_atom", "PotEng_norm_err_eV_per_atom", "Potential energy"),
    ("Enthalpy_norm_eV_per_atom", "Enthalpy_norm_err_eV_per_atom", "Enthalpy"),
    ("Free_norm_eV_per_atom", "Free_norm_err_eV_per_atom", r"$H - TS_{vib}$"),
]
FONT = 14
ROW_COLORS = {"0 GPa": "#4b657d", "10 GPa": "#b84f41"}


def read_rows(path, pressure):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row["pressure"] = pressure
            row["temperature_K"] = float(row["temperature_K"])
            row["x_B"] = float(row["x_B"])
            for value_col, err_col, _ in ENERGY_COLUMNS:
                row[value_col] = float(row[value_col])
                row[err_col] = float(row[err_col])
            rows.append(row)
    return rows


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(float(np.min(x)), float(np.max(x)), 300)
    return xs, np.polyval(coeff, xs)


def plot_temperature(ax, rows, pressure, value_col, err_col, title=None):
    data = [row for row in rows if row["pressure"] == pressure]
    data.sort(key=lambda row: row["x_B"])
    x = np.array([row["x_B"] for row in data], dtype=float)
    y = np.array([row[value_col] for row in data], dtype=float)
    yerr = np.array([row[err_col] for row in data], dtype=float)
    color = ROW_COLORS[pressure]
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        color=color,
        ecolor="black",
        elinewidth=1.4,
        capsize=4,
        capthick=1.4,
        barsabove=True,
        markersize=4.2,
        zorder=4,
    )
    xs, ys = fit_curve(x, y)
    ax.plot(xs, ys, color="black", linestyle="--", linewidth=1.6, zorder=3)
    if title is not None:
        ax.set_title(title, fontsize=FONT)
    ax.tick_params(labelsize=FONT, direction="in", width=0.8, length=4, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.grid(True, alpha=0.25)
    ax.set_xlim(0.08, 0.92)
    ax.set_xticks(np.arange(0.2, 1.0, 0.2))


def main():
    all_rows = []
    for pressure, path in SOURCES.items():
        all_rows.extend(read_rows(path, pressure))

    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)

    for temperature in TEMPERATURES:
        temp_rows = [row for row in all_rows if row["temperature_K"] == temperature]
        fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.3333333333), sharex=True, sharey="col")

        for col, (value_col, err_col, title) in enumerate(ENERGY_COLUMNS):
            plot_temperature(axes[0, col], temp_rows, "0 GPa", value_col, err_col, title=title)
            plot_temperature(axes[1, col], temp_rows, "10 GPa", value_col, err_col, title=None)
            axes[0, col].yaxis.set_major_locator(MultipleLocator(0.05))
            axes[1, col].yaxis.set_major_locator(MultipleLocator(0.05))

        axes[0, 0].set_ylabel("0 GPa\nEnergy (eV/atom)", fontsize=FONT)
        axes[1, 0].set_ylabel("10 GPa\nEnergy (eV/atom)", fontsize=FONT)
        for ax in axes[1, :]:
            ax.set_xlabel(r"$x$ in Fe$_{1-x}$B$_x$", fontsize=FONT)

        fig.suptitle(f"{temperature:.0f} K", fontsize=FONT)
        fig.tight_layout(rect=[0, 0, 1, 0.95], pad=0.5, w_pad=0.45, h_pad=0.35)
        png = out_dir / f"energy_0GPa_10GPa_{temperature:.0f}K_combined.png"
        pdf = out_dir / f"energy_0GPa_10GPa_{temperature:.0f}K_combined.pdf"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {png}")


if __name__ == "__main__":
    main()
