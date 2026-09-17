import csv
import re
from collections import defaultdict
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "output"
PRESSURE_OUT_DIR = ROOT / "output"
FIT_SUMMARY = ROOT / "output" / "fit_summary.csv"
KMAX_FIT = 1.2
K_PLOT_MAX = 3.05


def oz(k, a, b):
    return a / (1.0 + b * k * k)


def fit_oz(k, s):
    k = np.asarray(k, float)
    s = np.asarray(s, float)
    good = np.isfinite(k) & np.isfinite(s) & (s > 0) & (k <= KMAX_FIT)
    u, y = k[good], s[good]
    if len(u) < 3:
        return np.nan, np.nan
    a0 = float(y[np.argmin(u)])
    with np.errstate(divide="ignore", invalid="ignore"):
        bg = np.nanmedian((a0 / y - 1.0) / (u * u))
    b0 = bg if np.isfinite(bg) and bg > 0 else 1.0
    popt, _ = curve_fit(oz, u, y, p0=[a0, b0], bounds=([0.0, 0.0], [np.inf, np.inf]), maxfev=10000)
    return float(popt[0]), float(popt[1])


def parse_pressure(csv_path):
    parent_match = re.search(r"data-(\d+(?:\.\d+)?)GPa", csv_path.parent.parent.name)
    if parent_match:
        return float(parent_match.group(1))
    stem_match = re.search(r"_(\d+(?:\.\d+)?)GPa$", csv_path.stem)
    return float(stem_match.group(1)) if stem_match else np.nan


def parse_composition(stem):
    match = re.search(r"single_Fe_B_scc_(\d+(?:\.\d+)?)_", stem)
    return float(match.group(1)) if match else np.nan


def pressure_label(pressure):
    return f"{pressure:g}GPa"


def output_stem(case):
    return f"single_Fe_B_scc_{case['composition']:.1f}_{pressure_label(case['pressure'])}"


def load_case(csv_path):
    k, scc = [], []
    meta = {}
    scc0 = np.nan
    filename_composition = parse_composition(csv_path.stem)
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["route"] != "reciprocal":
                continue
            k.append(float(row["k"]))
            scc.append(float(row["S_cc"]))
            if row["S_cc0"]:
                scc0 = float(row["S_cc0"])
            if not meta:
                meta = {
                    "species_a": row["species_a"],
                    "species_b": row["species_b"],
                    "x_a": float(row["x_A"]),
                    "x_b": float(row["x_B_pair"]),
                    "n_frames": int(row["n_frames"]),
                    "frac_tag": row["frac_tag"],
                }
    order = np.argsort(k)
    return {
        "path": csv_path,
        "pressure": parse_pressure(csv_path),
        "composition": meta.get("x_b", filename_composition),
        "filename_composition": filename_composition,
        "scc0": scc0,
        "k": np.asarray(k)[order],
        "scc": np.asarray(scc)[order],
        **meta,
    }


def plot_case(case):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    a, b = fit_oz(case["k"], case["scc"])
    out_png = OUT_DIR / f"{output_stem(case)}_compare_fit.png"
    out_pdf = OUT_DIR / f"{output_stem(case)}_compare_fit.pdf"
    xprod = case["x_a"] * case["x_b"]

    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    ax.axhline(xprod, color="0.55", lw=1.3, ls="--", label=rf"$x_{{{case['species_a']}}}x_{{{case['species_b']}}}$")
    ax.scatter(case["k"], case["scc"], s=34, color="#4C88C7", alpha=0.75, edgecolor="none", label="reciprocal")
    if np.isfinite(a):
        kfit = np.linspace(0.0, KMAX_FIT, 200)
        ax.plot(kfit, oz(kfit, a, b), color="red", lw=1.8, ls="--", label="reciprocal OZ fit")
        ax.scatter([0.0], [a], marker="x", s=80, color="red", lw=2.0)
    ax.axvspan(0, KMAX_FIT, color="0.8", alpha=0.13, lw=0)
    ax.set_xlim(-0.05, K_PLOT_MAX)
    ax.set_ylim(0.0, max(0.6, float(np.nanmax(case["scc"])) * 1.1))
    ax.set_xlabel(r"$|k|$ [$\mathrm{\AA}^{-1}$]")
    ax.set_ylabel(r"$S_{cc}(|k|)$")
    ax.set_title(rf"Fe-B, {case['pressure']:.1f} GPa, $x_{{{case['species_a']}}}={case['x_a']:.2f}$")
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return a, b, out_png, out_pdf


def plot_pressure_group(composition, cases):
    PRESSURE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.8, 3.5))
    xprod = cases[0]["x_a"] * cases[0]["x_b"]
    ax.axhline(xprod, color="0.55", lw=1.3, ls="--", label=rf"$x_{{Fe}}x_{{B}}={xprod:.2f}$")
    colors = ["#4C88C7", "#d95f02", "#7570b3", "#1b9e77", "#e7298a"]
    for idx, case in enumerate(sorted(cases, key=lambda x: x["pressure"])):
        a, b = fit_oz(case["k"], case["scc"])
        color = colors[idx % len(colors)]
        ax.scatter(
            case["k"],
            case["scc"],
            s=28,
            color=color,
            alpha=0.7,
            edgecolor="none",
            label=f"{case['pressure']:.1f} GPa",
        )
        if np.isfinite(a):
            kfit = np.linspace(0.0, KMAX_FIT, 200)
            ax.plot(kfit, oz(kfit, a, b), color=color, lw=1.7, ls="--")
            ax.scatter([0.0], [a], marker="x", s=70, color=color, lw=1.8)
    ax.axvspan(0, KMAX_FIT, color="0.8", alpha=0.12, lw=0)
    ax.set_xlim(-0.05, K_PLOT_MAX)
    ax.set_ylim(0.0, max(0.6, max(float(np.nanmax(c["scc"])) for c in cases) * 1.1))
    ax.set_xlabel(r"$|k|$ [$\mathrm{\AA}^{-1}$]")
    ax.set_ylabel(r"$S_{cc}(|k|)$")
    ax.set_title(rf"Fe-B, $x_{{Fe}}={cases[0]['x_a']:.2f}, x_{{B}}={cases[0]['x_b']:.2f}$")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    out_png = PRESSURE_OUT_DIR / f"single_Fe_B_scc_{composition:.1f}_pressure_compare.png"
    out_pdf = PRESSURE_OUT_DIR / f"single_Fe_B_scc_{composition:.1f}_pressure_compare.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_png, out_pdf


def write_summary(rows):
    with open(FIT_SUMMARY, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "csv",
                "composition",
                "pressure_GPa",
                "x_Fe",
                "x_B",
                "n_frames",
                "kmax_fit",
                "oz_a_Scc0_fit",
                "oz_b",
                "S_cc0_csv",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def resolve_input_paths(names):
    paths = []
    for name in names:
        if not name or "/" in name or "\\" in name or name in {".", ".."}:
            raise ValueError(f"Input must be a CSV filename inside {ROOT}: {name!r}")
        path = ROOT / "data" / name
        if path.suffix.lower() != ".csv" or not path.is_file():
            raise ValueError(f"Input CSV does not exist: {path}")
        paths.append(path)
    return paths


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv:
        paths = resolve_input_paths(argv)
    else:
        paths = sorted(ROOT.glob("data-*/csv/single_Fe_B_scc_*.csv"))
    cases = [load_case(path) for path in paths]
    summary = []
    grouped = defaultdict(list)
    for case in cases:
        a, b, out_png, _ = plot_case(case)
        print(f"Wrote {out_png}")
        grouped[case["composition"]].append(case)
        summary.append({
            "csv": str(case["path"].relative_to(ROOT)),
            "composition": f"{case['composition']:.8g}",
            "pressure_GPa": f"{case['pressure']:.8g}",
            "x_Fe": f"{case['x_a']:.8g}",
            "x_B": f"{case['x_b']:.8g}",
            "n_frames": case["n_frames"],
            "kmax_fit": KMAX_FIT,
            "oz_a_Scc0_fit": f"{a:.10g}",
            "oz_b": f"{b:.10g}",
            "S_cc0_csv": "" if not np.isfinite(case["scc0"]) else f"{case['scc0']:.10g}",
        })
    for composition, group in sorted(grouped.items()):
        out_png, _ = plot_pressure_group(composition, group)
        print(f"Wrote {out_png}")
    write_summary(summary)
    print(f"Wrote {FIT_SUMMARY}")


if __name__ == "__main__":
    main()
