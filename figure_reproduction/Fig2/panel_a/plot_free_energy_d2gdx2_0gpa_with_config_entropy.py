import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ENERGY_SOURCE = ROOT / "data" / "normalized_energy_0GPa_liquid.csv"
D2G_SOURCE = ROOT / "data" / "dgmix_recon_0GPa.csv"
OUTPUT_DIR = ROOT / "output"
TEMPERATURES = [1500.0, 1800.0]
KB_EV_PER_K = 8.617333262145e-5
PANEL_COLORS = {1500.0: "#1f6b3b", 1800.0: "#5b2a86"}
LEFT_CURVE_COLOR = "#356fb3"
RIGHT_CURVE_COLOR = "#c94f66"
LEFT_FILL_COLOR = "#d9e8f6"
RIGHT_FILL_COLOR = "#f8dce2"
X_LIMITS = (0.1, 0.9)
GRAY_GUIDE_X = [6.0 / 29.0, 1.0 / 3.0, 0.5, 0.25]
ENERGY_Y_TICK_STEP = 0.1
D2G_BOX_ASPECT = 1.0 / 2.4
D2G_FIT_DEGREE = 5
D2G_FIT_POINTS = 300
D2G_MARKERSIZE = 5


def configurational_free_energy(temperature, x):
    """Return k_B T [x ln x + (1-x) ln(1-x)] in eV/atom."""
    entropy_term = 0.0
    if x > 0.0:
        entropy_term += x * math.log(x)
    if x < 1.0:
        entropy_term += (1.0 - x) * math.log(1.0 - x)
    return KB_EV_PER_K * temperature * entropy_term


def read_energy_rows():
    rows = []
    with open(ENERGY_SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["temperature_K"])
            if temperature not in TEMPERATURES:
                continue
            x_b = float(row["x_B"])
            free_vib_norm = float(row["Free_norm_eV_per_atom"])
            free_config_norm = configurational_free_energy(temperature, x_b)
            rows.append({
                "temperature_K": temperature,
                "x_B": x_b,
                "free_vib_norm": free_vib_norm,
                "free_config_norm": free_config_norm,
                "free_total_norm": free_vib_norm + free_config_norm,
            })
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


def prepare_energy_panel(rows, temperature):
    data = [row for row in rows if row["temperature_K"] == temperature]
    data.sort(key=lambda row: row["x_B"])
    x = np.array([row["x_B"] for row in data], dtype=float)
    y = np.array([row["free_total_norm"] for row in data], dtype=float)
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


def plot_energy_panel(ax, panel, baseline):
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

    ax.set_xlim(*X_LIMITS)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.tick_params(labelsize=6, direction="in", width=0.25, length=2, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.axvline(0.8, color="gray", linestyle="--", linewidth=1, alpha=0.9)
    ax.text(0.735, 0.155, r"FeB$_4$", fontsize=6, ha="center", va="center")
    legend_x0, legend_x1, text_x = 0.08, 0.17, 0.20
    y_rows = [0.93, 0.84, 0.75]
    ax.plot([(legend_x0 + legend_x1) / 2.0], [y_rows[0]], marker="o", color=color, markersize=2.6, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[0], f"0 GPa, {panel['temperature']:.0f} K", transform=ax.transAxes, fontsize=6, color=color, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[1], y_rows[1]], color=LEFT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[1], "Convex", transform=ax.transAxes, fontsize=6, color=LEFT_CURVE_COLOR, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[2], y_rows[2]], color=RIGHT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[2], "Concave", transform=ax.transAxes, fontsize=6, color=RIGHT_CURVE_COLOR, va="center")


def read_d2g_rows():
    rows = []
    with open(D2G_SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            temperature = float(row["T_K"])
            if temperature not in TEMPERATURES:
                continue
            rows.append({
                "temperature_K": temperature,
                "x_B": float(row["x_B"]),
                "d2g_dx2": float(row["d2g_dx2"]),
                "is_slab": int(float(row["is_slab_4_6"])) == 1,
            })
    return rows


def plot_d2g_panel(ax, rows, temperature):
    data = [row for row in rows if row["temperature_K"] == temperature]
    data.sort(key=lambda row: row["x_B"])
    x = np.array([row["x_B"] for row in data], dtype=float)
    y = np.array([row["d2g_dx2"] for row in data], dtype=float)
    color = PANEL_COLORS[temperature]
    degree = min(D2G_FIT_DEGREE, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(float(np.min(x)), float(np.max(x)), D2G_FIT_POINTS)
    ys = np.polyval(coeff, xs)
    ax.plot(xs, ys, "-", color=color, linewidth=0.8, alpha=0.95, zorder=2)
    ax.plot(x, y, marker="o", linestyle="none", markerfacecolor="white", markeredgecolor=color, markeredgewidth=0.6, markersize=D2G_MARKERSIZE, zorder=4, clip_on=False)
    ax.axhline(0.0, color="black", linewidth=0.5, linestyle="--", zorder=1)
    ax.set_xlim(0.095, 0.905)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.tick_params(labelsize=6, direction="in", width=0.25, length=2, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    return float(np.min(y)), float(np.max(y))


def main():
    energy_rows = read_energy_rows()
    d2g_rows = read_d2g_rows()
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(3.25, 2.5),
        sharex="col",
        gridspec_kw={"height_ratios": [1.6, 0.9]},
    )
    panels = [prepare_energy_panel(energy_rows, temperature) for temperature in TEMPERATURES]
    global_ymin = min(panel["baseline"] for panel in panels)
    global_ymax = max(max(panel["ymax"] for panel in panels), 0.17)
    energy_ymin, energy_ymax = tick_bounds_including_zero(global_ymin, global_ymax, ENERGY_Y_TICK_STEP)
    energy_yticks = np.arange(energy_ymin, energy_ymax + ENERGY_Y_TICK_STEP * 0.5, ENERGY_Y_TICK_STEP)
    d2g_panel_data = []
    for col, panel in enumerate(panels):
        plot_energy_panel(axes[0, col], panel, energy_ymin)
        axes[0, col].set_ylim(energy_ymin, energy_ymax)
        axes[0, col].set_yticks(energy_yticks)
        d2g_panel_data.append(plot_d2g_panel(axes[1, col], d2g_rows, panel["temperature"]))
    d2g_min = min(item[0] for item in d2g_panel_data)
    d2g_max = max(item[1] for item in d2g_panel_data)
    d2g_ymin = d2g_min - 3
    d2g_ymax = d2g_max + 3
    for ax in axes[1, :]:
        ax.set_ylim(d2g_ymin, d2g_ymax)
        ax.axhspan(d2g_ymin, 0.0, color="#E8E8E8", alpha=1, zorder=0)
    axes[0, 1].tick_params(labelleft=False)
    axes[1, 1].tick_params(labelleft=False)
    for ax in axes[1, :]:
        ax.set_box_aspect(D2G_BOX_ASPECT)
    axes[0, 0].set_ylabel(r"$H-T(S_{vib}+S_{mix})$ [eV/atom]", fontsize=6)
    axes[1, 0].set_ylabel(r"$\partial^2g_{ex}/\partial x^2$ [$k_BT$]", fontsize=6)
    axes[1, 0].set_xlabel(r"$x$ in Fe$_{1-x}$B$_x$", fontsize=6)
    axes[1, 1].set_xlabel(r"$x$ in Fe$_{1-x}$B$_x$", fontsize=6)
    fig.tight_layout(pad=0.35, w_pad=0.35, h_pad=0.15)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / "free_energy_d2g_0GPa_1500K_1800K_with_config_entropy.pdf"
    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "free_energy_d2g_0GPa_1500K_1800K_with_config_entropy.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
