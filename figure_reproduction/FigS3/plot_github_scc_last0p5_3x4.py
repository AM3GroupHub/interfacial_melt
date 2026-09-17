"""Reproduce Figure S3 (1500 K) from the archived S_cc(k) data."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import numpy.core
import numpy.core.multiarray
import numpy.core.numeric
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data" / "scc_reciprocal"
OUT_DIR = ROOT / "output"
SAFE_DATA_PATH = ROOT / "data" / "last0p5_records_flat.csv"
SUMMARY_PATH = OUT_DIR / "last0p5_weighted_oz_bootstrap_summary.csv"
DEFAULT_KMAX_FIT = 0.8
KMAX_PLOT = 2.55
N_BOOTSTRAP = 500
RNG_SEED = 20260824
MAX_WORKERS = 4
TEMPERATURES = [1500, 1800]
PLOT_TEMPERATURES = [1500]  # Keep both fitting temperatures for the original seeds.
PRESSURES = [0, 10]
COMPOSITIONS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90]
PRESSURE_STYLES = {0: ("#4C88C7", "o"), 10: ("#d95f02", "^")}


def install_numpy_pickle_aliases():
    """Read NPZ files written by NumPy versions using numpy._core."""
    sys.modules.setdefault("numpy._core", numpy.core)
    sys.modules.setdefault("numpy._core.multiarray", numpy.core.multiarray)
    sys.modules.setdefault("numpy._core.numeric", numpy.core.numeric)


def oz(k, a, b):
    return a / (1.0 + b * k * k)


def parse_state(path):
    match = re.fullmatch(r"T(\d+)_P(\d+)_last0\.5", path.stem)
    if not match:
        raise ValueError(f"Unexpected last0.5 filename: {path.name}")
    return int(match.group(1)), int(match.group(2))


def source_paths():
    return sorted(DATA_DIR.glob("T*_P*_last0.5.npz"))


def extract_safe_data():
    paths = source_paths()
    if not paths:
        raise FileNotFoundError(
            "NPZ extraction is optional and no NPZ files are packaged. "
            "Use the supplied data/last0p5_records_flat.csv for reproduction."
        )
    install_numpy_pickle_aliases()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = ["P_GPa", "T_K", "x_B", "n_frames", "k", "S_cc", "weight"]
    with open(SAFE_DATA_PATH, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for path in paths:
            temperature, pressure = parse_state(path)
            with np.load(path, allow_pickle=True) as data:
                state_records = data["records"].tolist()
            for record in state_records:
                for k, scc, weight in zip(record["uk"], record["scc"], record["w"]):
                    writer.writerow({
                        "P_GPa": pressure,
                        "T_K": temperature,
                        "x_B": f"{float(record['x_B']):.8g}",
                        "n_frames": int(record["n_frames"]),
                        "k": f"{float(k):.10g}",
                        "S_cc": f"{float(scc):.10g}",
                        "weight": f"{float(weight):.10g}",
                    })
    print(f"[extract] {SAFE_DATA_PATH}")


def ensure_safe_data():
    # Use the archived table directly; no upstream checkout or NPZ cache is needed.
    if not SAFE_DATA_PATH.is_file():
        raise FileNotFoundError(f"Missing packaged fit data: {SAFE_DATA_PATH}")


def load_records():
    ensure_safe_data()
    records = {}
    with open(SAFE_DATA_PATH, newline="") as fh:
        for row in csv.DictReader(fh):
            pressure = int(row["P_GPa"])
            temperature = int(row["T_K"])
            if pressure not in PRESSURES or temperature not in TEMPERATURES:
                continue
            key = (pressure, temperature, round(float(row["x_B"]), 2))
            record = records.setdefault(key, {"uk": [], "scc": [], "w": [], "x_B": float(row["x_B"])})
            record["n_frames"] = int(row["n_frames"])
            record["uk"].append(float(row["k"]))
            record["scc"].append(float(row["S_cc"]))
            record["w"].append(float(row["weight"]))
    for record in records.values():
        record["uk"] = np.asarray(record["uk"], float)
        record["scc"] = np.asarray(record["scc"], float)
        record["w"] = np.asarray(record["w"], float)
    expected = len(PRESSURES) * len(TEMPERATURES) * len(COMPOSITIONS)
    if len(records) != expected:
        raise ValueError(f"Expected {expected} records, found {len(records)}")
    return records


def select_fit_data(record, kmax_fit):
    k = np.asarray(record["uk"], float)
    scc = np.asarray(record["scc"], float)
    weight = np.asarray(record.get("w", np.ones_like(k)), float)
    good = np.isfinite(k) & np.isfinite(scc) & (scc > 0) & (k < kmax_fit)
    return k[good], scc[good], weight[good]


def initial_guess(k, scc):
    a0 = float(scc[np.argmin(k)])
    with np.errstate(divide="ignore", invalid="ignore"):
        bg = np.nanmedian((a0 / scc - 1.0) / (k * k))
    b0 = bg if np.isfinite(bg) and bg > 0 else 1.0
    return [a0, b0]


def fit_weighted(k, scc, weight, p0=None):
    guess = initial_guess(k, scc) if p0 is None else p0
    popt, _ = curve_fit(
        oz,
        k,
        scc,
        p0=guess,
        sigma=1.0 / weight,
        bounds=([0.0, 0.0], [np.inf, np.inf]),
        maxfev=10000,
    )
    return np.asarray(popt, float)


def bootstrap_fit(record, rng, kmax_fit):
    k, scc, weight = select_fit_data(record, kmax_fit)
    popt = fit_weighted(k, scc, weight)
    fitted = oz(k, *popt)
    standardized = (scc - fitted) * weight
    standardized -= np.mean(standardized)
    samples = []
    failures = 0
    for _ in range(N_BOOTSTRAP):
        simulated = fitted + rng.choice(standardized, len(k), replace=True) / weight
        try:
            sample = fit_weighted(k, simulated, weight, p0=popt)
        except (RuntimeError, ValueError, FloatingPointError):
            failures += 1
            continue
        if np.all(np.isfinite(sample)):
            samples.append(sample)
        else:
            failures += 1
    if len(samples) < N_BOOTSTRAP // 2:
        raise RuntimeError(f"Too few bootstrap fits: {len(samples)}/{N_BOOTSTRAP}")

    samples = np.asarray(samples, float)
    kfit = np.linspace(0.0, kmax_fit, 180)
    curves = oz(kfit[:, None], samples[:, 0], samples[:, 1])
    ylo, yhi = np.percentile(curves, [2.5, 97.5], axis=1)
    a_ci = np.percentile(samples[:, 0], [2.5, 50.0, 97.5])
    b_ci = np.percentile(samples[:, 1], [2.5, 50.0, 97.5])
    residual = scc - fitted
    return {
        "k": k,
        "scc": scc,
        "weight": weight,
        "a": float(popt[0]),
        "b": float(popt[1]),
        "kfit": kfit,
        "yfit": oz(kfit, *popt),
        "ylo": ylo,
        "yhi": yhi,
        "a_ci": a_ci,
        "b_ci": b_ci,
        "rmse": float(np.sqrt(np.mean(residual * residual))),
        "success": len(samples),
        "failures": failures,
    }


def analyze_one(key, record, seed, kmax_fit):
    rng = np.random.default_rng(seed)
    return key, bootstrap_fit(record, rng, kmax_fit)


def analyze_records(records, kmax_fit):
    keys = []
    for pressure in PRESSURES:
        for composition in COMPOSITIONS:
            for temperature in TEMPERATURES:
                keys.append((pressure, temperature, composition))
    seeds = np.random.SeedSequence(RNG_SEED).generate_state(len(keys))
    analyses = {}
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [
            executor.submit(analyze_one, key, records[key], int(seed), kmax_fit)
            for key, seed in zip(keys, seeds)
        ]
        for future in as_completed(futures):
            key, result = future.result()
            analyses[key] = result
            pressure, temperature, composition = key
            print(f"[fit] P={pressure:2d} T={temperature} x_B={composition:.2f}", flush=True)
    return analyses


def thin_mask(k, threshold=1.0, step=3):
    mask = np.ones(len(k), dtype=bool)
    high = np.flatnonzero(k > threshold)
    mask[high] = np.arange(len(high)) % step == 0
    return mask


def fit_is_unstable(result):
    a_ratio = result["a_ci"][2] / max(result["a"], 1.0e-12)
    return a_ratio > 10.0 or result["b_ci"][2] > 1.0e4


def panel_limits(records, analyses, temperature, composition):
    values = [0.0]
    for pressure in PRESSURES:
        key = (pressure, temperature, composition)
        record = records[key]
        result = analyses[key]
        k = np.asarray(record["uk"], float)
        scc = np.asarray(record["scc"], float)
        visible = np.isfinite(scc) & (k <= KMAX_PLOT)
        values.extend(scc[visible].tolist())
        if not fit_is_unstable(result):
            values.extend(result["yhi"].tolist())
    ymin = min(0.0, float(np.nanmin(values)))
    ymax = max(0.2, float(np.nanmax(values)))
    padding = 0.08 * (ymax - ymin)
    return ymin - padding, ymax + padding


def plot_panel(ax, records, analyses, temperature, composition, kmax_fit):
    ideal = composition * (1.0 - composition)
    ylim = panel_limits(records, analyses, temperature, composition)
    ax.axhline(ideal, color="0.70", lw=1.0, ls=":", zorder=0)
    ax.axvspan(0.0, kmax_fit, color="0.75", alpha=0.10, lw=0, zorder=0)
    unstable_slot = 0
    for pressure in PRESSURES:
        color, marker = PRESSURE_STYLES[pressure]
        key = (pressure, temperature, composition)
        record = records[key]
        result = analyses[key]
        k = np.asarray(record["uk"], float)
        scc = np.asarray(record["scc"], float)
        mask = thin_mask(k)
        if pressure == 0:
            ax.scatter(
                k[mask], scc[mask], s=28, facecolor="none", edgecolor=color,
                marker=marker, linewidth=0.9, alpha=0.52,
            )
        else:
            ax.scatter(
                k[mask], scc[mask], s=30, color=color, marker=marker,
                edgecolor="none", alpha=0.42,
            )
        ax.fill_between(result["kfit"], result["ylo"], result["yhi"], color=color, alpha=0.14, lw=0)
        ax.plot(result["kfit"], result["yfit"], color=color, lw=1.45)
        if fit_is_unstable(result):
            ax.text(
                0.97,
                0.96 - 0.11 * unstable_slot,
                rf"{pressure} GPa: $S_{{cc}}(0)={result['a']:.2g}$",
                transform=ax.transAxes,
                color=color,
                fontsize=12,
                ha="right",
                va="top",
                zorder=7,
                bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.82, "pad": 1.5},
            )
            unstable_slot += 1
    ax.set_title(rf"$x_B={composition:.2f}$", fontsize=15, pad=4)
    ax.set_xlim(-0.04, KMAX_PLOT)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 1, 2])
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)


def add_legend_panel(ax, kmax_fit):
    ax.axis("off")
    handles = [
        Line2D(
            [], [], color=color, lw=1.8, marker=marker, ms=8,
            markerfacecolor="none" if pressure == 0 else color,
            label=f"{pressure} GPa",
        )
        for pressure, (color, marker) in PRESSURE_STYLES.items()
    ]
    handles.extend([
        Patch(facecolor="0.5", alpha=0.16, label="95% fit CI"),
        Line2D([], [], color="0.65", lw=1.2, ls=":", label=r"ideal $x_{Fe}x_B$"),
        Patch(facecolor="0.75", alpha=0.16, label=rf"fit window $k<{kmax_fit:g}$"),
        Patch(facecolor="white", edgecolor="0.35", label="off-scale value: top-right"),
    ])
    ax.legend(handles=handles, loc="center", frameon=False, fontsize=14, labelspacing=0.7)


def plot_temperature(records, analyses, temperature, kmax_fit, output_suffix):
    fig, axes = plt.subplots(3, 4, figsize=(14.5, 8.7), sharex=True)
    for ax, composition in zip(axes.flat, COMPOSITIONS):
        plot_panel(ax, records, analyses, temperature, composition, kmax_fit)
    add_legend_panel(axes.flat[-1], kmax_fit)
    for row in range(3):
        axes[row, 0].set_ylabel(r"$S_{cc}(|k|)$", fontsize=15)
    for col in range(4):
        if axes[2, col].axison:
            axes[2, col].set_xlabel(r"$|k|$ [$\mathrm{\AA}^{-1}$]", fontsize=15)
    fig.suptitle(f"Fe-B concentration structure factors, T = {temperature} K, last 50%", fontsize=19)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.08, top=0.92, wspace=0.25, hspace=0.28)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = OUT_DIR / f"scc_last0p5_T{temperature}_pressure_compare_3x4{output_suffix}.png"
    out_pdf = OUT_DIR / f"scc_last0p5_T{temperature}_pressure_compare_3x4{output_suffix}.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure] {out_png}")
    return out_png, out_pdf


def summary_rows(records, analyses, kmax_fit):
    for pressure in PRESSURES:
        for composition in COMPOSITIONS:
            for temperature in TEMPERATURES:
                key = (pressure, temperature, composition)
                record = records[key]
                result = analyses[key]
                yield {
                    "P_GPa": pressure,
                    "T_K": temperature,
                    "x_B": f"{composition:.2f}",
                    "n_frames": int(record["n_frames"]),
                    "n_fit_points": len(result["k"]),
                    "kmax_fit": kmax_fit,
                    "a_Scc0": f"{result['a']:.10g}",
                    "a_ci025": f"{result['a_ci'][0]:.10g}",
                    "a_ci50": f"{result['a_ci'][1]:.10g}",
                    "a_ci975": f"{result['a_ci'][2]:.10g}",
                    "b": f"{result['b']:.10g}",
                    "b_ci025": f"{result['b_ci'][0]:.10g}",
                    "b_ci50": f"{result['b_ci'][1]:.10g}",
                    "b_ci975": f"{result['b_ci'][2]:.10g}",
                    "rmse": f"{result['rmse']:.10g}",
                    "bootstrap_success": result["success"],
                    "bootstrap_failures": result["failures"],
                }


def write_summary(records, analyses, kmax_fit, summary_path):
    rows = list(summary_rows(records, analyses, kmax_fit))
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[summary] {summary_path}")


def main(kmax_fit=DEFAULT_KMAX_FIT, output_suffix=""):
    records = load_records()
    analyses = analyze_records(records, kmax_fit)
    for temperature in PLOT_TEMPERATURES:
        plot_temperature(records, analyses, temperature, kmax_fit, output_suffix)
    summary_path = (
        SUMMARY_PATH
        if not output_suffix
        else OUT_DIR / f"last0p5_weighted_oz_bootstrap_summary{output_suffix}.csv"
    )
    write_summary(records, analyses, kmax_fit, summary_path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kmax-fit", type=float, default=DEFAULT_KMAX_FIT)
    parser.add_argument(
        "--output-suffix",
        default="",
        help="Suffix added before .png/.pdf and to the summary CSV filename.",
    )
    args = parser.parse_args()
    if args.kmax_fit <= 0.0:
        parser.error("--kmax-fit must be positive")
    if args.output_suffix and not args.output_suffix.startswith("_"):
        args.output_suffix = "_" + args.output_suffix
    return args


if __name__ == "__main__":
    if sys.argv[1:] == ["--extract-safe"]:
        extract_safe_data()
    else:
        args = parse_args()
        main(args.kmax_fit, args.output_suffix)
