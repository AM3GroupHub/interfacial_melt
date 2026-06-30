import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "csv" / "energy_composition_b_1500k_raw.csv"

TEMPERATURE = 1500.0
REFERENCE_FORMULA = "B"
PRESSURE_LABEL = "0 GPa"
OUTPUT_STEM = "energy_composition_B_1500K"
FIT_DEGREE = 6
N_FIT_SAMPLES = 300

P0 = -8.043070824
P1 = -6.434663479166667
ENT0 = -7.850557452
ENT1 = -6.2431970833333335


def parse_formula(formula):
    counts = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def mole_ratio_of_reference(formula, reference_formula):
    counts = parse_formula(formula)
    reference_count = counts.get(reference_formula, 0)
    total_atoms = sum(counts.values())
    if total_atoms == 0:
        raise ValueError(f"Invalid formula: {formula}")
    return reference_count / total_atoms, total_atoms


def read_rows():
    rows = []
    with open(SOURCE, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if float(row["temperature_K"]) != TEMPERATURE:
                continue
            x_ref, divisor = mole_ratio_of_reference(row["formula"], REFERENCE_FORMULA)
            rows.append(
                {
                    "formula": row["formula"],
                    "x_ref": x_ref,
                    "divisor": divisor,
                    "pot_eng": float(row["Average_PotEng_eV"]),
                    "enthalpy": float(row["Average_Enthalpy_eV"]),
                }
            )
    if not rows:
        raise ValueError("No valid rows found for 1500 K")
    return rows


def normalize_rows(rows):
    normalized = []
    for row in rows:
        x_ref = row["x_ref"]
        divisor = row["divisor"]
        pot_norm = row["pot_eng"] / divisor - x_ref * P1 - (1.0 - x_ref) * P0
        ent_norm = row["enthalpy"] / divisor - x_ref * ENT1 - (1.0 - x_ref) * ENT0
        normalized.append(
            {
                "formula": row["formula"],
                "x_ref": x_ref,
                "pot_norm": pot_norm,
                "ent_norm": ent_norm,
            }
        )
    return normalized


def plot_segmented_curve_by_curvature(ax, poly, x_range, n_samples=N_FIT_SAMPLES):
    xs = np.linspace(x_range[0], x_range[1], n_samples)
    ys = poly(xs)
    first_derivative = np.polyder(poly, 1)
    second_derivative = np.polyder(poly, 2)
    y_prime = first_derivative(xs)
    y_double_prime = second_derivative(xs)
    curvatures = y_double_prime / ((1.0 + y_prime ** 2) ** 1.5)
    sign_changes = np.where(np.diff(np.sign(curvatures)))[0]

    segments = []
    start_idx = 0
    for change_idx in sign_changes:
        segments.append((start_idx, change_idx + 1))
        start_idx = change_idx + 1
    segments.append((start_idx, len(xs)))

    for start, end in segments:
        if start >= end:
            continue
        x_segment = xs[start:end]
        y_segment = ys[start:end]
        mid_idx = min((start + end) // 2, len(curvatures) - 1)
        kappa_mid = curvatures[mid_idx]
        color = "red" if kappa_mid < 0 else ("blue" if kappa_mid > 0 else "gray")
        ax.plot(x_segment, y_segment, color=color, linewidth=2.5, alpha=0.9)

    ax.legend(
        [Line2D([0], [0], color="blue", lw=2.5), Line2D([0], [0], color="red", lw=2.5)],
        ["κ > 0 (Convex)", "κ < 0 (Concave)"],
        fontsize=10,
        frameon=True,
    )


def plot_one_panel(ax, x, y, color_points, title, ylabel):
    ax.scatter(x, y, s=90, alpha=0.80, edgecolors="black", linewidth=1.0, color=color_points)
    ax.plot(x, y, "-", alpha=0.25, color=color_points)
    if len(x) >= 4:
        poly = np.poly1d(np.polyfit(x, y, FIT_DEGREE))
        plot_segmented_curve_by_curvature(ax, poly, (float(np.min(x)), float(np.max(x))))
    ax.set_xlabel(f"Mole ratio of {REFERENCE_FORMULA} (x)", fontsize=12, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)


def main():
    rows = normalize_rows(read_rows())
    rows.sort(key=lambda row: row["x_ref"])
    x = np.array([row["x_ref"] for row in rows], dtype=float)
    pot = np.array([row["pot_norm"] for row in rows], dtype=float)
    ent = np.array([row["ent_norm"] for row in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    plot_one_panel(
        axes[0],
        x,
        pot,
        color_points="tab:blue",
        title=f"Potential Energy vs x({REFERENCE_FORMULA}) at {int(TEMPERATURE)} K",
        ylabel="Normalized Potential Energy (eV)",
    )
    plot_one_panel(
        axes[1],
        x,
        ent,
        color_points="green",
        title=f"Enthalpy vs x({REFERENCE_FORMULA}) at {int(TEMPERATURE)} K and {PRESSURE_LABEL}",
        ylabel="Normalized Enthalpy (eV)",
    )

    fig.tight_layout()
    eps_path = ROOT / f"{OUTPUT_STEM}.eps"
    pdf_path = ROOT / f"{OUTPUT_STEM}.pdf"
    png_path = ROOT / f"{OUTPUT_STEM}.png"
    fig.savefig(eps_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {eps_path}")


if __name__ == "__main__":
    main()
