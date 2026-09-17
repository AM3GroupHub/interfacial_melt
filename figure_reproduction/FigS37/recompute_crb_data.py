import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE_ROOT = ROOT / "data" / "source_crb"
SYSTEM_LABEL = "Cr-B liquid"
PLOT_TEMPERATURES = [1500.0, 1800.0]
ENDPOINT_FILES = {
    "Cr": ("Cr250", SOURCE_ROOT / "Cr" / "Cr" / "NPT" / "Cr250_G.txt"),
    "B": ("B144", SOURCE_ROOT / "B" / "B_small" / "NPT" / "B144_G.txt"),
}
ENERGY_COLUMNS = [
    ("Average_PotEng", "Average_PotEng_err", "Potential energy"),
    ("Average_Enthalpy", "Average_Enthalpy_err", "Enthalpy"),
    ("Free_Energy_H_minus_TS", "Free_Energy_err", "Free energy"),
]
KB_EV_PER_K = 8.617333262145e-5
PANEL_COLORS = {1500.0: "#1f6b3b", 1800.0: "#5b2a86"}
LEFT_CURVE_COLOR = "#356fb3"
RIGHT_CURVE_COLOR = "#c94f66"
LEFT_FILL_COLOR = "#d9e8f6"
RIGHT_FILL_COLOR = "#f8dce2"
X_LIMITS = (0.1, 0.9)
GRAY_GUIDE_X = [6.0 / 29.0, 1.0 / 3.0, 0.5, 0.25]
ENERGY_Y_TICK_STEP = 0.1
FIGSIZE = (3.25, 2.05)
ANNOTATION_FONTSIZE = 6
AXIS_LABEL_FONTSIZE = 6
TICK_FONTSIZE = 6
LEGEND_FONTSIZE = 6


def parse_formula(formula):
    counts = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def clean_header(line):
    names = line.lstrip("# ").split()
    return [name.replace("(K)", "").replace("(eV)", "").replace("(eV/K)", "") for name in names]


def read_energy_table(path):
    header = None
    rows = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("# Temperature"):
                header = clean_header(line)
                continue
            if line.startswith("#") or header is None:
                continue
            parts = line.split()
            if len(parts) < len(header):
                continue
            row = {}
            for key, value in zip(header, parts):
                try:
                    row[key] = float(value)
                except ValueError:
                    row[key] = value
            rows[row["Temperature"]] = row
    return rows


def endpoint_per_atom(formula, path):
    atom_count = sum(parse_formula(formula).values())
    table = read_energy_table(path)
    output = {}
    for temperature, row in table.items():
        output[temperature] = {}
        for energy_key, err_key, _ in ENERGY_COLUMNS:
            output[temperature][energy_key] = row[energy_key] / atom_count
            output[temperature][err_key] = row[err_key] / atom_count
    return output


def output_key(energy_key):
    return energy_key.replace("Average_", "").replace("Free_Energy_H_minus_TS", "Free")


def normalize_record(formula, row, endpoints):
    counts = parse_formula(formula)
    n_cr = counts.get("Cr", 0)
    n_b = counts.get("B", 0)
    n_total = n_cr + n_b
    record = {
        "formula": formula,
        "temperature_K": row["Temperature"],
        "Cr_count": n_cr,
        "B_count": n_b,
        "x_B": n_b / n_total,
    }
    for energy_key, err_key, _ in ENERGY_COLUMNS:
        cr_e = endpoints["Cr"][row["Temperature"]][energy_key]
        b_e = endpoints["B"][row["Temperature"]][energy_key]
        cr_err = endpoints["Cr"][row["Temperature"]][err_key]
        b_err = endpoints["B"][row["Temperature"]][err_key]
        residual = row[energy_key] - n_cr * cr_e - n_b * b_e
        residual_err = math.sqrt(row[err_key] ** 2 + (n_cr * cr_err) ** 2 + (n_b * b_err) ** 2)
        short = output_key(energy_key)
        record[f"{short}_norm_eV_per_atom"] = residual / n_total
        record[f"{short}_norm_err_eV_per_atom"] = residual_err / n_total
        record[f"raw_{short}_eV"] = row[energy_key]
        record[f"raw_{short}_err_eV"] = row[err_key]
    return record


def collect_records():
    endpoints = {key: endpoint_per_atom(formula, path) for key, (formula, path) in ENDPOINT_FILES.items()}
    endpoint_formulas = {formula for formula, _ in ENDPOINT_FILES.values()}
    records = []
    for path in sorted(SOURCE_ROOT.glob("**/*_G.txt")):
        formula = path.stem[:-2] if path.stem.endswith("_G") else path.stem
        if formula in endpoint_formulas:
            continue
        table = read_energy_table(path)
        for temperature in PLOT_TEMPERATURES:
            if temperature in table:
                records.append(normalize_record(formula, table[temperature], endpoints))
    return sorted(records, key=lambda item: (item["temperature_K"], item["x_B"]))


def configurational_free_energy(temperature, x):
    entropy_term = 0.0
    if x > 0.0:
        entropy_term += x * math.log(x)
    if x < 1.0:
        entropy_term += (1.0 - x) * math.log(1.0 - x)
    return KB_EV_PER_K * temperature * entropy_term


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
    y = np.array(
        [row["Free_norm_eV_per_atom"] + configurational_free_energy(temperature, row["x_B"]) for row in data],
        dtype=float,
    )
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

    ax.set_xlim(*X_LIMITS)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.tick_params(labelsize=TICK_FONTSIZE, direction="in", width=0.25, length=2, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    ax.axvline(0.8, color="gray", linestyle="--", linewidth=1, alpha=0.9)
    ax.text(0.735, 0.155, r"CrB$_4$", fontsize=ANNOTATION_FONTSIZE, ha="center", va="center")
    legend_x0, legend_x1, text_x = 0.08, 0.17, 0.20
    y_rows = [0.93, 0.84, 0.75]
    ax.plot([(legend_x0 + legend_x1) / 2.0], [y_rows[0]], marker="o", color=color, markersize=2.6, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[0], f"0 GPa, {panel['temperature']:.0f} K", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=color, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[1], y_rows[1]], color=LEFT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[1], "Convex", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=LEFT_CURVE_COLOR, va="center")
    ax.plot([legend_x0, legend_x1], [y_rows[2], y_rows[2]], color=RIGHT_CURVE_COLOR, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(text_x, y_rows[2], "Concave", transform=ax.transAxes, fontsize=ANNOTATION_FONTSIZE, color=RIGHT_CURVE_COLOR, va="center")


def write_csv(records):
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "normalized_energy_CrB_0GPa_1500_1800.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    return out_path


def main():
    records = collect_records()
    if not records:
        raise ValueError("No Cr-B composition records found")
    csv_path = write_csv(records)

    print(f"Wrote {csv_path}")
    print(f"Records: {len(records)}")


if __name__ == "__main__":
    main()
