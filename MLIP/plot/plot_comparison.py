import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data" / "csv"

FONT_FAMILY = "Arial"
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 20
FONT_SIZE_TICK = 20
FONT_SIZE_LEGEND = 20
FONT_WEIGHT_TITLE = "bold"
FONT_WEIGHT_LABEL = "normal"
FONT_WEIGHT_TICK = "normal"
FONT_WEIGHT_LEGEND = "normal"

FIGURES = [
    (DATA_DIR / "energy_per_atom_comparison.csv", "Energy", "eV/atom", "Energy_plot"),
    (DATA_DIR / "forces_comparison.csv", "Forces", "eV/A", "Forces_plot"),
]


plt.rcParams.update(
    {
        "font.family": FONT_FAMILY,
        "font.size": FONT_SIZE_LABEL,
        "axes.titlesize": FONT_SIZE_TITLE,
        "axes.titleweight": FONT_WEIGHT_TITLE,
        "axes.labelsize": FONT_SIZE_LABEL,
        "axes.labelweight": FONT_WEIGHT_LABEL,
        "xtick.labelsize": FONT_SIZE_TICK,
        "ytick.labelsize": FONT_SIZE_TICK,
        "legend.fontsize": FONT_SIZE_LEGEND,
        "figure.titlesize": FONT_SIZE_TITLE,
    }
)


def read_numeric_pairs(csv_file):
    rows = []
    removed = 0
    with open(csv_file, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if len(fieldnames) < 2:
            raise ValueError(f"Expected at least two columns in {csv_file}")
        col_true, col_pred = fieldnames[0], fieldnames[1]
        for row in reader:
            try:
                true_val = float(row[col_true])
                pred_val = float(row[col_pred])
            except (TypeError, ValueError):
                removed += 1
                continue
            rows.append((true_val, pred_val))
    if not rows:
        raise ValueError(f"No valid numeric rows found in {csv_file}")
    return col_true, col_pred, np.array(rows, dtype=float), removed


def compute_metrics(true_vals, pred_vals):
    diff = pred_vals - true_vals
    mse = float(np.mean(diff ** 2))
    rmse = float(np.sqrt(mse))
    mae = float(np.mean(np.abs(diff)))
    return rmse, mae


def plot_comparison(csv_file, title, unit, output_stem):
    print(f"\nProcessing file: {csv_file.name} ...")
    col_true, col_pred, data, removed = read_numeric_pairs(csv_file)
    true_vals = data[:, 0]
    pred_vals = data[:, 1]

    print(f"  - Detected columns: {col_true}, {col_pred}")
    if removed:
        print(f"  - Removed {removed} invalid rows.")
    print("  - Data preview (first 5 rows):")
    for true_val, pred_val in data[:5]:
        print(f"    {true_val:.10f}, {pred_val:.10f}")

    rmse, mae = compute_metrics(true_vals, pred_vals)
    print("  - Metrics:")
    print(f"    RMSE: {rmse:.6f} {unit}")
    print(f"    MAE:  {mae:.6f} {unit}")

    fig, ax = plt.subplots(figsize=(6, 6), dpi=100)
    ax.scatter(
        true_vals,
        pred_vals,
        marker="x",
        s=70,
        alpha=1.0,
        linewidths=1.2,
        color="steelblue",
        label="Data",
        rasterized=True,
    )

    min_val = min(float(np.min(true_vals)), float(np.min(pred_vals)))
    max_val = max(float(np.max(true_vals)), float(np.max(pred_vals)))
    margin = (max_val - min_val) * 0.05
    if margin == 0.0:
        margin = 1.0
    limits = [min_val - margin, max_val + margin]
    ax.plot(limits, limits, "k--", alpha=1.0, zorder=0, label="Ideal")

    ax.set_xlabel(
        f"True {title} ({unit})",
        fontsize=FONT_SIZE_LABEL,
        fontweight=FONT_WEIGHT_LABEL,
        fontfamily=FONT_FAMILY,
    )
    ax.set_ylabel(
        f"Predicted {title} ({unit})",
        fontsize=FONT_SIZE_LABEL,
        fontweight=FONT_WEIGHT_LABEL,
        fontfamily=FONT_FAMILY,
    )
    ax.set_title(
        f"{title} Comparison\nRMSE = {rmse:.4f} {unit} | MAE = {mae:.4f} {unit}",
        fontsize=FONT_SIZE_TITLE,
        fontweight=FONT_WEIGHT_TITLE,
        fontfamily=FONT_FAMILY,
    )

    ax.tick_params(axis="both", labelsize=FONT_SIZE_TICK)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily(FONT_FAMILY)
        label.set_fontweight(FONT_WEIGHT_TICK)

    ax.legend(
        prop={
            "family": FONT_FAMILY,
            "size": FONT_SIZE_LEGEND,
            "weight": FONT_WEIGHT_LEGEND,
        }
    )
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()

    pdf_path = ROOT / f"{output_stem}.pdf"
    png_path = ROOT / f"{output_stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  - Saved: {pdf_path.name}")
    print(f"  - Saved: {png_path.name}")


def main():
    for csv_file, title, unit, output_stem in FIGURES:
        plot_comparison(csv_file, title, unit, output_stem)


if __name__ == "__main__":
    main()
