import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CRYSTAL = ROOT / "data" / "source" / "crystal"
LIQUID = ROOT / "data" / "source" / "liquid"
PRESSURE_LABEL = "10 GPa"
PLOT_TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]


def parse_formula(formula):
    counts = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def clean_header(line):
    return [name.replace("(K)", "").replace("(eV)", "") for name in line.lstrip("# ").split()]


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


def endpoint_pv_per_atom(formula, path):
    atom_count = sum(parse_formula(formula).values())
    table = read_energy_table(path)
    output = {}
    for temperature, row in table.items():
        output[temperature] = {
            "PV_eV_per_atom": row["Average_PV"] / atom_count,
            "PV_err_eV_per_atom": row["Average_PV_err"] / atom_count,
        }
    return output


def normalize_pv(row, counts, endpoints):
    n_fe = counts.get("Fe", 0)
    n_b = counts.get("B", 0)
    n_total = n_fe + n_b
    temperature = row["Temperature"]
    fe_pv = endpoints["Fe"][temperature]["PV_eV_per_atom"]
    b_pv = endpoints["B"][temperature]["PV_eV_per_atom"]
    fe_err = endpoints["Fe"][temperature]["PV_err_eV_per_atom"]
    b_err = endpoints["B"][temperature]["PV_err_eV_per_atom"]
    residual = row["Average_PV"] - n_fe * fe_pv - n_b * b_pv
    residual_err = math.sqrt(row["Average_PV_err"] ** 2 + (n_fe * fe_err) ** 2 + (n_b * b_err) ** 2)
    return residual / n_total, residual_err / n_total


def collect_records():
    endpoints = {
        "Fe": endpoint_pv_per_atom("Fe250", CRYSTAL / "Fe250.txt"),
        "B": endpoint_pv_per_atom("B144", CRYSTAL / "B144.txt"),
    }
    records = []
    for path in sorted(LIQUID.glob("*_G.txt")):
        formula = path.stem.removesuffix("_G")
        counts = parse_formula(formula)
        n_total = counts.get("Fe", 0) + counts.get("B", 0)
        x_b = counts.get("B", 0) / n_total
        table = read_energy_table(path)
        for temperature in PLOT_TEMPERATURES:
            if temperature not in table:
                continue
            row = table[temperature]
            pv_norm, pv_norm_err = normalize_pv(row, counts, endpoints)
            records.append(
                {
                    "formula": formula,
                    "temperature_K": temperature,
                    "Fe_count": counts.get("Fe", 0),
                    "B_count": counts.get("B", 0),
                    "x_B": x_b,
                    "Average_PV_eV": row["Average_PV"],
                    "Average_PV_err_eV": row["Average_PV_err"],
                    "PV_eV_per_atom": row["Average_PV"] / n_total,
                    "PV_err_eV_per_atom": row["Average_PV_err"] / n_total,
                    "PV_norm_eV_per_atom": pv_norm,
                    "PV_norm_err_eV_per_atom": pv_norm_err,
                }
            )
    return sorted(records, key=lambda item: (item["temperature_K"], item["x_B"]))


def fit_curve(x, y):
    degree = min(5, len(x) - 1)
    coeff = np.polyfit(x, y, degree)
    xs = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    return xs, np.polyval(coeff, xs)


def write_csv(records):
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "pv_composition_10GPa_liquid.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    return out_path


def plot(records):
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(PLOT_TEMPERATURES)))
    for color, temperature in zip(colors, PLOT_TEMPERATURES):
        rows = [row for row in records if row["temperature_K"] == temperature]
        x = np.array([row["x_B"] for row in rows], dtype=float)
        y = np.array([row["PV_norm_eV_per_atom"] for row in rows], dtype=float)
        yerr = np.array([row["PV_norm_err_eV_per_atom"] for row in rows], dtype=float)
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
        xs, ys = fit_curve(x, y)
        ax.plot(xs, ys, "--", color=color, linewidth=1.8, alpha=0.85)
    ax.set_xlabel(r"B fraction, $x_B = B / (Fe + B)$", fontsize=16)
    ax.set_ylabel("Endpoint-normalized PV term (eV/atom)", fontsize=16)
    ax.set_title(f"{PRESSURE_LABEL}: endpoint-normalized PV term vs composition", fontsize=16)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=16, frameon=True)
    fig.tight_layout()
    out_path = ROOT / "pv_composition_10GPa_liquid.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / "pv_composition_10GPa_liquid.pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    records = collect_records()
    csv_path = write_csv(records)
    print(f"Wrote {csv_path}")
    print(f"Records: {len(records)}")


if __name__ == "__main__":
    main()
