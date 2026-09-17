import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
PRESSURE_DIRS = {
    "0 GPa": ROOT / "data" / "coordination" / "0GPa",
    "10 GPa": ROOT / "data" / "coordination" / "10GPa",
}
TEMPERATURES = [1200.0, 1500.0, 1800.0]
COMPOSITIONS = [
    ("Fe160B240", "Fe0.4B0.6"),
    ("Fe120B280", "Fe0.3B0.7"),
    ("Fe135B405", "Fe0.25B0.75"),
    ("Fe72B288", "Fe0.2B0.8"),
    ("Fe81B459", "Fe0.15B0.85"),
    ("Fe54B486", "Fe0.1B0.9"),
]
PRESSURE_COLORS = {"0 GPa": "#2b6cb0", "10 GPa": "#c53030"}
PRESSURE_HATCH = {"0 GPa": "", "10 GPa": "//"}


def read_frame_fraction(coord_csv):
    with open(coord_csv, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    total = len(rows)
    six = sum(int(row["coordination_B"]) == 6 for row in rows)
    return six / total if total else np.nan


def collect_data():
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for pressure, base in PRESSURE_DIRS.items():
        for composition_dir in base.iterdir():
            if not composition_dir.is_dir():
                continue
            if composition_dir.name.endswith("_coor"):
                composition = composition_dir.name.removesuffix("_coor")
            else:
                composition = composition_dir.name
            for temp_dir in composition_dir.iterdir():
                if not temp_dir.is_dir() or not temp_dir.name.endswith("_dir"):
                    continue
                temperature = float(temp_dir.name.removesuffix("_dir"))
                if temperature not in TEMPERATURES:
                    continue
                frame_dirs = sorted([p for p in temp_dir.iterdir() if p.is_dir() and p.name.startswith("frame_")])
                for frame_dir in frame_dirs:
                    coord_csv = frame_dir / "b_coordination.csv"
                    if coord_csv.exists():
                        data[temperature][composition][pressure].append(read_frame_fraction(coord_csv))
    return data


def summarize(samples):
    arr = np.array(samples, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return np.nan, np.nan
    mean = float(np.mean(arr))
    err = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
    return mean, err


def subscript_label(label):
    return label.replace("Fe", r"Fe$_{").replace("B", r"}$B$_{") + r"}$"


def plot_temperature(temperature, data):
    x = np.arange(len(COMPOSITIONS))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11, 4.2))

    for idx, pressure in enumerate(PRESSURE_DIRS):
        offset = (-0.5 + idx) * width
        means = []
        errs = []
        for composition, _ in COMPOSITIONS:
            samples = data.get(temperature, {}).get(composition, {}).get(pressure, [])
            mean, err = summarize(samples)
            means.append(mean)
            errs.append(err)
        bars = ax.bar(
            x + offset,
            means,
            width=width,
            label=pressure,
            color=PRESSURE_COLORS[pressure],
            edgecolor="black",
            linewidth=0.6,
            hatch=PRESSURE_HATCH[pressure],
            alpha=0.85,
            zorder=3,
        )
        ax.errorbar(
            x + offset,
            means,
            yerr=errs,
            fmt="none",
            ecolor="black",
            elinewidth=0.8,
            capsize=3,
            capthick=0.8,
            zorder=4,
        )

    ax.set_xticks(x)
    ax.set_xticklabels([subscript_label(label) for _, label in COMPOSITIONS], rotation=25, ha="right", fontsize=16)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"Fraction of six-coordinated B", fontsize=16)
    ax.set_title(f"{temperature:.0f} K", fontsize=16)
    ax.tick_params(labelsize=16)
    ax.grid(True, axis="y", alpha=0.25, zorder=0)
    ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=16)
    fig.tight_layout()
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / f"sixfold_fraction_{temperature:.0f}K.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"sixfold_fraction_{temperature:.0f}K.pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def write_summary_csv(data):
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "sixfold_fraction_summary.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["temperature_K", "composition", "pressure", "mean_fraction", "sem", "n_frames"])
        for temperature in TEMPERATURES:
            for composition, _ in COMPOSITIONS:
                for pressure in PRESSURE_DIRS:
                    samples = data.get(temperature, {}).get(composition, {}).get(pressure, [])
                    mean, err = summarize(samples)
                    writer.writerow([temperature, composition, pressure, f"{mean:.8f}", f"{err:.8f}", len(samples)])
    return out_path


def main():
    data = collect_data()
    summary_csv = write_summary_csv(data)
    print(f"Wrote {summary_csv}")


if __name__ == "__main__":
    main()
