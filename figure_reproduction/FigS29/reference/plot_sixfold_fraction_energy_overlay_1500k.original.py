import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data" / "csv"
SIXFOLD_SUMMARY = DATA_DIR / "sixfold_fraction_summary_1500k.csv"
POTENTIAL_DATA = DATA_DIR / "potential_energy_0gpa_1500k.csv"
TEMPERATURE = 1500.0
COMPOSITIONS = [
    ("Fe120B280", "Fe0.3B0.7", 0.70),
    ("Fe135B405", "Fe0.25B0.75", 0.75),
    ("Fe72B288", "Fe0.2B0.8", 0.80),
    ("Fe81B459", "Fe0.15B0.85", 0.85),
    ("Fe54B486", "Fe0.1B0.9", 0.90),
]
BAR_STYLE = {
    "0 GPa": {"color": "#4b657d", "label": "0 GPa fraction"},
    "10 GPa": {"color": "#b84f41", "label": "10 GPa fraction"},
}
LINE_COLOR = "#e69500"
LEFT_BG = "#66ccff"
RIGHT_BG = "#f8dce2"

AXIS_LABEL_FONT = 8
Y_TICK_LABEL_FONT = 8
X_TICK_LABEL_FONT = 7
LEGEND_FONT = 6


def apply_font_sizes():
    plt.rcParams.update(
        {
            "font.size": AXIS_LABEL_FONT,
            "axes.labelsize": AXIS_LABEL_FONT,
            "xtick.labelsize": X_TICK_LABEL_FONT,
            "ytick.labelsize": Y_TICK_LABEL_FONT,
            "legend.fontsize": LEGEND_FONT,
        }
    )


def read_fraction_summary():
    data = {}
    with open(SIXFOLD_SUMMARY, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["composition"], float(row["temperature_K"]), row["pressure"])
            data[key] = (float(row["mean_fraction"]), float(row["sem"]))
    return data


def read_potential_data():
    points = {}
    all_points = []
    with open(POTENTIAL_DATA, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature != TEMPERATURE:
                continue
            formula = row["formula"]
            x_b = float(row["x_B"])
            y = float(row["PotEng_norm_eV_per_atom"])
            err = float(row["PotEng_norm_err_eV_per_atom"])
            points[formula] = (x_b, y, err)
            all_points.append((x_b, y))
    return points, all_points


def fit_curve(points):
    x = np.array([item[0] for item in points], dtype=float)
    y = np.array([item[1] for item in points], dtype=float)
    coeff = np.polyfit(x, y, min(5, len(x) - 1))
    xs = np.linspace(0.60, 0.90, 300)
    ys = np.polyval(coeff, xs)
    return xs, ys


def map_xb_to_plot_x(values):
    values = np.asarray(values, dtype=float)
    return (values - 0.70) * ((len(COMPOSITIONS) - 1) / 0.20)


def add_horizontal_gradient(ax, x0, x1, y0, y1, left_color, right_color):
    n = 256
    left_rgba = np.array(to_rgba(left_color), dtype=float)
    right_rgba = np.array(to_rgba(right_color), dtype=float)
    ramp = np.linspace(0.0, 1.0, n)
    colors = left_rgba[None, :] * (1.0 - ramp[:, None]) + right_rgba[None, :] * ramp[:, None]
    image = colors[np.newaxis, :, :]
    ax.imshow(image, extent=[x0, x1, y0, y1], origin="lower", aspect="auto", zorder=0)


def main():
    apply_font_sizes()
    fraction = read_fraction_summary()
    potential, all_points = read_potential_data()

    x = np.arange(len(COMPOSITIONS), dtype=float)
    bar_width = 0.30
    fig, ax = plt.subplots(figsize=(3.25, 1.8))
    ax2 = ax.twinx()
    x_min, x_max = -0.45, len(COMPOSITIONS) - 0.55
    split_x = 1.5
    add_horizontal_gradient(ax, x_min, split_x, 0.0, 1.0, LEFT_BG, "white")
    add_horizontal_gradient(ax, split_x, x_max, 0.0, 1.0, "white", RIGHT_BG)

    for idx, pressure in enumerate(["0 GPa", "10 GPa"]):
        offset = (-0.5 + idx) * bar_width
        means = []
        for composition, _, _ in COMPOSITIONS:
            mean, _ = fraction[(composition, TEMPERATURE, pressure)]
            means.append(mean)
        ax.bar(
            x + offset,
            means,
            width=bar_width,
            color=BAR_STYLE[pressure]["color"],
            edgecolor="black",
            linewidth=0.25,
            alpha=0.95,
            label=BAR_STYLE[pressure]["label"],
            zorder=2,
        )

    # Match the original script: the fit uses all potential points,
    # while only the selected high-B compositions are drawn as markers.
    xs, ys = fit_curve(all_points)
    comp_x = np.arange(len(COMPOSITIONS), dtype=float)
    comp_y = [potential[composition][1] for composition, _, _ in COMPOSITIONS]
    fit_plot_x = map_xb_to_plot_x(xs)

    ax2.plot(fit_plot_x, ys, color="black", linestyle="--", linewidth=0.8, zorder=5)
    ax2.plot(
        comp_x,
        comp_y,
        linestyle="none",
        marker="^",
        color=LINE_COLOR,
        markerfacecolor=LINE_COLOR,
        markeredgecolor="black",
        markeredgewidth=0.25,
        markersize=5.6,
        zorder=6,
    )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, 1.0)
    ax2.set_ylim(min(comp_y) - 0.05, max(comp_y) + 0.05)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            r"Fe$_{0.3}$B$_{0.7}$",
            r"Fe$_{0.25}$B$_{0.75}$",
            r"Fe$_{0.2}$B$_{0.8}$",
            r"Fe$_{0.15}$B$_{0.85}$",
            r"Fe$_{0.1}$B$_{0.9}$",
        ],
        rotation=0,
        ha="center",
    )
    ax.set_ylabel("Fraction of B(CN$=6$)")
    ax2.set_ylabel("Potential energy [eV/atom]")
    ax.tick_params(axis="x", labelsize=X_TICK_LABEL_FONT, direction="in", width=0.25, length=2, top=True)
    ax.tick_params(axis="y", labelsize=Y_TICK_LABEL_FONT, direction="in", width=0.25, length=2, right=False)
    ax2.tick_params(axis="y", labelsize=Y_TICK_LABEL_FONT, direction="in", width=0.25, length=2, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.25)
    ax.axvline(split_x, color="gray", linestyle="--", linewidth=0.25, zorder=1.5)

    bar_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=BAR_STYLE[pressure]["color"], edgecolor="none", label=BAR_STYLE[pressure]["label"])
        for pressure in ["0 GPa", "10 GPa"]
    ]
    line_handle = plt.Line2D(
        [0],
        [0],
        color=LINE_COLOR,
        marker="^",
        markerfacecolor=LINE_COLOR,
        markeredgecolor="black",
        markeredgewidth=0.25,
        linestyle="none",
        markersize=3,
        label="0 GPa potential",
    )
    fit_handle = plt.Line2D([0], [0], color="black", linestyle="--", linewidth=0.8, label="0 GPa potential fit")
    ax.legend(handles=bar_handles + [line_handle, fit_handle], frameon=False, loc="upper right", ncol=1)

    fig.tight_layout(pad=0.4)
    png_path = ROOT / "sixfold_fraction_energy_overlay_1500K_3b.png"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(ROOT / "sixfold_fraction_energy_overlay_1500K_3b.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "sixfold_fraction_energy_overlay_1500K_3b.eps", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
