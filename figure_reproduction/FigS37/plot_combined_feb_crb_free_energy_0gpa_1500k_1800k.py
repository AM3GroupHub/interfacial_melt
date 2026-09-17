import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FEB_SOURCE = ROOT / "data" / "normalized_energy_0GPa_liquid.csv"
CRB_SOURCE = ROOT / "data" / "normalized_energy_CrB_0GPa_1500_1800.csv"
TEMPERATURES = [1500.0, 1800.0]
KB_EV_PER_K = 8.617333262145e-5
PANEL_COLORS = {1500.0: "#1f6b3b", 1800.0: "#5b2a86"}
LEFT_CURVE_COLOR = "#356fb3"
RIGHT_CURVE_COLOR = "#c94f66"
LEFT_FILL_COLOR = "#d9e8f6"
RIGHT_FILL_COLOR = "#f8dce2"
X_LIMITS = (0.1, 0.9)
FE_GUIDE_X = [6.0 / 29.0, 1.0 / 3.0, 0.5, 0.25]
ENERGY_Y_TICK_STEP = 0.1
FIGSIZE = (3.25, 3.28)
ANNOTATION_FONTSIZE = 6
AXIS_LABEL_FONTSIZE = 6
TICK_FONTSIZE = 6


def configurational_free_energy(temperature, x):
    entropy_term = 0.0
    if x > 0.0:
        entropy_term += x * math.log(x)
    if x < 1.0:
        entropy_term += (1.0 - x) * math.log(1.0 - x)
    return KB_EV_PER_K * temperature * entropy_term


def read_normalized_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature not in TEMPERATURES:
                continue
            x_b = float(row["x_B"])
            free_vib_norm = float(row["Free_norm_eV_per_atom"])
            rows.append(
                {
                    "temperature_K": temperature,
                    "x_B": x_b,
                    "free_total_norm": free_vib_norm + configurational_free_energy(temperature, x_b),
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


def tick_bounds_including_zero(ymin, ymax, step):
    lower = min(ymin, 0.0)
    upper = max(ymax, 0.0)
    return math.floor(lower / step) * step, math.ceil(upper / step) * step


def prepare_panel(rows, temperature):
    data = [row for row in rows if row["temperature_K"] == temperature]
    data.sort(key=lambda row: row["x_B"])
    x = np.array([row["x_B"] for row in data], dtype=float)
    y = np.array([row["free_total_norm"] for row in data], dtype=float)
    coeff, xs, ys = fit_curve(x, y)
    inflection_x = closest_inflection_point(coeff, 0.6, float(np.min(x)), float(np.max(x)))
    if inflection_x is None:
        inflection_x = 0.6
    return {
        "temperature": temperature,
        "x": x,
        "y": y,
        "color": PANEL_COLORS[temperature],
        "coeff": coeff,
        "xs": xs,
        "ys": ys,
        "inflection_x": inflection_x,
        "baseline": min(float(np.min(y)), float(np.min(ys))) - 0.01,
        "ymax": max(float(np.max(y)), float(np.max(ys))),
    }


def plot_energy_panel(ax, panel, baseline, label_text, show_fe_guides, show_x_label, xlabel):
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

    if show_fe_guides:
        for x_guide in FE_GUIDE_X:
            y_guide = float(np.polyval(coeff, x_guide))
            ax.plot([x_guide, x_guide], [baseline, y_guide], color="gray", linestyle=":", linewidth=1, alpha=0.9, zorder=2)

    ax.set_xlim(*X_LIMITS)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.tick_params(labelsize=TICK_FONTSIZE, direction="in", width=0.25, length=2, top=True, right=True)
    if not show_x_label:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.axvline(0.8, color="gray", linestyle="--", linewidth=1, alpha=0.9)
    ax.text(0.735, 0.155, label_text, fontsize=ANNOTATION_FONTSIZE, ha="center", va="center")

    legend_x0, legend_x1, text_x = 0.08, 0.17, 0.20
    y_rows = [0.93, 0.84, 0.75]
    ax.plot([(legend_x0 + legend_x1) / 2.0], [y_rows[0]], marker="o", color=color, markersize=2.6, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[0], f"0 GPa, {panel['temperature']:.0f} K", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=color, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[1], y_rows[1]], color=LEFT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[1], "Convex", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=LEFT_CURVE_COLOR, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[2], y_rows[2]], color=RIGHT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[2], "Concave", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=RIGHT_CURVE_COLOR, va="center")


def row_limits(panels):
    global_ymin = min(panel["baseline"] for panel in panels)
    global_ymax = max(max(panel["ymax"] for panel in panels), 0.17)
    y_min, y_max = tick_bounds_including_zero(global_ymin, global_ymax, ENERGY_Y_TICK_STEP)
    yticks = np.arange(y_min, y_max + ENERGY_Y_TICK_STEP * 0.5, ENERGY_Y_TICK_STEP)
    return y_min, y_max, yticks


def apply_row_limits(axes, y_min, y_max, yticks):
    for ax in axes:
        ax.set_ylim(y_min, y_max)
        ax.set_yticks(yticks)
    axes[1].tick_params(labelleft=False)


def main():
    feb_rows = read_normalized_rows(FEB_SOURCE)
    crb_rows = read_normalized_rows(CRB_SOURCE)
    feb_panels = [prepare_panel(feb_rows, temperature) for temperature in TEMPERATURES]
    crb_panels = [prepare_panel(crb_rows, temperature) for temperature in TEMPERATURES]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=FIGSIZE,
        sharex="col",
        gridspec_kw={"height_ratios": [1, 1]},
    )

    feb_ymin, feb_ymax, feb_yticks = row_limits(feb_panels)
    crb_ymin, crb_ymax, crb_yticks = row_limits(crb_panels)

    for col, panel in enumerate(feb_panels):
        plot_energy_panel(axes[0, col], panel, feb_ymin, r"FeB$_4$", True, False, "")
    for col, panel in enumerate(crb_panels):
        plot_energy_panel(axes[1, col], panel, crb_ymin, r"CrB$_4$", False, True, r"$x$ in M$_{1-x}$B$_x$")

    apply_row_limits(axes[0, :], feb_ymin, feb_ymax, feb_yticks)
    apply_row_limits(axes[1, :], crb_ymin, crb_ymax, crb_yticks)
    axes[0, 0].set_ylabel(r"$H-T(S_{vib}+S_{mix})$ [eV/atom]", fontsize=AXIS_LABEL_FONTSIZE)
    axes[1, 0].set_ylabel(r"$H-T(S_{vib}+S_{mix})$ [eV/atom]", fontsize=AXIS_LABEL_FONTSIZE)

    fig.tight_layout(pad=0.35, w_pad=0.35, h_pad=0.2)
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    out_base = out_dir / "combined_FeB_CrB_free_energy_0GPa_1500K_1800K_with_config_entropy"
    fig.savefig(out_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_base.with_suffix('.png')}")
    print(f"Wrote {out_base.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
