import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "csv" / "free_energy_0gpa_1500_1800.csv"
TEMPERATURES = [1500.0, 1800.0]
PANEL_COLORS = {
    1500.0: "#1f6b3b",
    1800.0: "#5b2a86",
}
LEFT_CURVE_COLOR = "#356fb3"
RIGHT_CURVE_COLOR = "#c94f66"
LEFT_FILL_COLOR = "#d9e8f6"
RIGHT_FILL_COLOR = "#f8dce2"
X_LIMITS = (0.1, 0.9)
GRAY_GUIDE_X = [1.0 / 3.0, 0.5, 0.25, 6.0 / 29.0]


def read_rows():
    rows = []
    with open(SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature not in TEMPERATURES:
                continue
            rows.append(
                {
                    "temperature_K": temperature,
                    "x_B": float(row["x_B"]),
                    "free_norm": float(row["Free_norm_eV_per_atom"]),
                }
            )
    return rows


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(X_LIMITS[0], X_LIMITS[1], 300)
    return coeff, xs, np.polyval(coeff, xs)


def closest_inflection_point(coeff, x_ref, x_min, x_max):
    second_derivative = np.polyder(coeff, 2)
    roots = np.roots(second_derivative)
    real_roots = [float(root.real) for root in roots if abs(root.imag) < 1e-8 and x_min <= root.real <= x_max]
    if not real_roots:
        real_roots = [float(root.real) for root in roots if abs(root.imag) < 1e-8]
    if not real_roots:
        return None
    return min(real_roots, key=lambda value: abs(value - x_ref))


def prepare_panel(rows, temperature):
    data = [row for row in rows if row["temperature_K"] == temperature]
    data.sort(key=lambda row: row["x_B"])
    x = np.array([row["x_B"] for row in data], dtype=float)
    y = np.array([row["free_norm"] for row in data], dtype=float)
    color = PANEL_COLORS[temperature]
    coeff, xs, ys = fit_curve(x, y)
    inflection_x = closest_inflection_point(coeff, 0.6, float(np.min(x)), float(np.max(x)))
    if inflection_x is None:
        inflection_x = 0.6
    return {
        "temperature": temperature,
        "x": x,
        "y": y,
        "color": color,
        "coeff": coeff,
        "xs": xs,
        "ys": ys,
        "inflection_x": inflection_x,
        "baseline": min(float(np.min(y)), float(np.min(ys))) - 0.01,
        "ymax": max(float(np.max(y)), float(np.max(ys))),
    }


def plot_panel(ax, panel, baseline):
    x = panel["x"]
    y = panel["y"]
    color = panel["color"]
    coeff = panel["coeff"]
    xs = panel["xs"]
    ys = panel["ys"]
    inflection_x = panel["inflection_x"]

    ax.plot(x, y, "o", color=color, markersize=2.8, zorder=5, clip_on=False)
    left_mask = xs <= inflection_x
    right_mask = xs >= inflection_x

    ax.fill_between(xs[left_mask], ys[left_mask], baseline, color=LEFT_FILL_COLOR, alpha=0.7, zorder=1)
    ax.fill_between(xs[right_mask], ys[right_mask], baseline, color=RIGHT_FILL_COLOR, alpha=0.7, zorder=1)
    ax.plot(xs[left_mask], ys[left_mask], "-", color=LEFT_CURVE_COLOR, linewidth=0.8, alpha=0.95, zorder=3)
    ax.plot(xs[right_mask], ys[right_mask], "-", color=RIGHT_CURVE_COLOR, linewidth=0.8, alpha=0.95, zorder=3)

    for x_guide in GRAY_GUIDE_X:
        y_guide = float(np.polyval(coeff, x_guide))
        ax.plot([x_guide, x_guide], [baseline, y_guide], color="gray", linestyle=":", linewidth=1, alpha=0.9, zorder=2)

    ax.set_xlabel(r"$x$ in Fe$_{1-x}$B$_x$", fontsize=6)
    ax.set_xlim(*X_LIMITS)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.tick_params(labelsize=5, direction="in", width=0.25, length=2, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.axvline(0.8, color="gray", linestyle="--", linewidth=1, alpha=0.9)
    ax.text(0.735, 0.155, r"FeB$_4$", fontsize=6, ha="center", va="center")

    legend_x0 = 0.08
    legend_x1 = 0.17
    text_x = 0.20
    y_rows = [0.93, 0.84, 0.75]
    ax.plot(
        [(legend_x0 + legend_x1) / 2.0],
        [y_rows[0]],
        marker="o",
        color=color,
        markersize=2.6,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(text_x, y_rows[0], f"0 GPa, {panel['temperature']:.0f} K", transform=ax.transAxes, fontsize=6, color=color, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[1], y_rows[1]], color=LEFT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[1], "Convex", transform=ax.transAxes, fontsize=6, color=LEFT_CURVE_COLOR, va="center")
    ax.plot(
        [legend_x0, legend_x1],
        [y_rows[2], y_rows[2]],
        color=RIGHT_CURVE_COLOR,
        linewidth=0.8,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        text_x,
        y_rows[2],
        "Concave",
        transform=ax.transAxes,
        fontsize=6,
        color=RIGHT_CURVE_COLOR,
        va="center",
    )


def main():
    rows = read_rows()
    fig, axes = plt.subplots(1, 2, figsize=(3.25, 1.6), sharey=True)
    panels = [prepare_panel(rows, temperature) for temperature in TEMPERATURES]
    global_ymin = min(panel["baseline"] for panel in panels)
    global_ymax = max(max(panel["ymax"] for panel in panels), 0.17)
    for ax, panel in zip(axes, panels):
        plot_panel(ax, panel, global_ymin)
    for ax in axes:
        ax.set_ylim(global_ymin, global_ymax)
    axes[0].set_ylabel(r"$H-TS_{vib}$ [eV/atom]", fontsize=6)
    fig.tight_layout(pad=0.45, w_pad=0.7)
    out_path = ROOT / "free_energy_0GPa_1500K_1800K.png"
    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    fig.savefig(ROOT / "free_energy_0GPa_1500K_1800K.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "free_energy_0GPa_1500K_1800K.eps", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
