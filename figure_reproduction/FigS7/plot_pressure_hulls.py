import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- User-editable plot settings ----------------
PLOT_TITLE_0K = "Fe-B convex hull at 0 K"
PLOT_TITLE_1800K = "Fe-B convex hull at 1800 K"
PLOT_TITLE_COMBINED = "Fe-B convex hulls under pressure"
X_LABEL = "B fraction"
Y_LABEL_0K = "Formation energy (eV/atom)"
Y_LABEL_1800K = "Formation energy (eV/atom)"
FIGSIZE_SINGLE = (16, 4.8)
FIGSIZE_COMBINED = (16, 8)
DPI = 300
TITLE_FONTSIZE = 20
SUBPLOT_TITLE_FONTSIZE = 16
AXIS_LABEL_FONTSIZE = 16
TICK_FONTSIZE = 16
ANNOTATION_FONTSIZE = 16
LEGEND_FONTSIZE = 16
HULL_LINEWIDTH = 1.4
POINT_SIZE = 60
GRID_ALPHA = 0.2
LABEL_OFFSET = 7
ENDPOINT_LABEL_OFFSET = -12
# -------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "data"
OUT_DIR = ROOT / "output"
PRESSURES = [0, 5, 10]
PHASES = ["Fe", "Fe2B", "Fe3B", "Fe23B6", "FeB", "FeB2", "FeB4", "boron"]
RES_NAME = {"Fe": "Fe", "Fe2B": "Fe2B", "Fe3B": "Fe3B", "Fe23B6": "Fe23B6", "FeB": "FeB", "FeB2": "FeB2", "FeB4": "FeB4", "boron": "Boron"}
FORMULA = {"Fe": (1, 0), "Fe2B": (2, 1), "Fe3B": (3, 1), "Fe23B6": (23, 6), "FeB": (1, 1), "FeB2": (1, 2), "FeB4": (1, 4), "boron": (0, 1)}
PHASE_COLOR = {"Fe": "#4c78a8", "Fe2B": "#f58518", "Fe3B": "#54a24b", "Fe23B6": "#b279a2", "FeB": "#e45756", "FeB2": "#72b7b2", "FeB4": "#ff9da6", "boron": "#9d755d"}
PV_COEF = 0.006241509074
KJMOL_TO_EV = 96.4853321233
T_HULL = 1800
TOL = 1e-10


def formula_label(phase):
    text = "B" if phase == "boron" else phase
    return re.sub(r"(\d+)", r"$_{\1}$", text)


def parse_res(path):
    parts = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0].split()
    return {"pressure_file": float(parts[2]), "volume": float(parts[3]), "energy_raw_cell_eV": float(parts[4]), "natom": int(parts[7])}


def parse_phonon(path):
    natom = None
    current_t = None
    free_energy = None
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = re.match(r"^natom:\s+(\d+)\s*$", raw)
        if m:
            natom = int(m.group(1))
            continue
        m = re.match(r"^- temperature:\s+([-+0-9.eE]+)\s*$", raw)
        if m:
            current_t = float(m.group(1))
            continue
        m = re.match(r"^\s+free_energy:\s+([-+0-9.eE]+)\s*$", raw)
        if m and current_t is not None:
            if int(round(current_t)) == T_HULL:
                free_energy = float(m.group(1))
                break
            current_t = None
    if natom is None or free_energy is None:
        raise ValueError(f"Bad phonon file: {path}")
    return {"natom": natom, "free_energy_kjmol": free_energy}


def lower_hull(points):
    pts = sorted(points, key=lambda p: (p["x_B"], p["formation_eV_atom"], p["phase"]))
    hull = []
    for point in pts:
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            cross = (b["x_B"] - a["x_B"]) * (point["formation_eV_atom"] - a["formation_eV_atom"]) - (b["formation_eV_atom"] - a["formation_eV_atom"]) * (point["x_B"] - a["x_B"])
            if cross <= TOL:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def interp_hull(hull, x):
    if x <= hull[0]["x_B"]:
        return hull[0]["formation_eV_atom"]
    if x >= hull[-1]["x_B"]:
        return hull[-1]["formation_eV_atom"]
    for i in range(len(hull) - 1):
        x1, y1 = hull[i]["x_B"], hull[i]["formation_eV_atom"]
        x2, y2 = hull[i + 1]["x_B"], hull[i + 1]["formation_eV_atom"]
        if x1 <= x <= x2:
            t = (x - x1) / (x2 - x1)
            return y1 * (1 - t) + y2 * t
    return hull[-1]["formation_eV_atom"]


def collect_data():
    raw = {}
    for pressure in PRESSURES:
        raw[pressure] = {}
        for phase in PHASES:
            phase_dir = DATA_ROOT / f"{pressure}GPa" / phase
            res = parse_res(phase_dir / "H" / f"{RES_NAME[phase]}.res")
            phonon = parse_phonon(phase_dir / "thermal_properties.yaml")
            if res["natom"] != phonon["natom"]:
                raise ValueError(f"natom mismatch for {phase_dir}")
            fe, b = FORMULA[phase]
            x_b = b / (fe + b)
            pv_cell = pressure * res["volume"] * PV_COEF
            h_atom = (res["energy_raw_cell_eV"] + pv_cell) / res["natom"]
            vib_atom = phonon["free_energy_kjmol"] / KJMOL_TO_EV / phonon["natom"]
            raw[pressure][phase] = {
                "phase": phase,
                "x_B": x_b,
                "raw_energy_cell_eV": res["energy_raw_cell_eV"],
                "pv_cell_eV": pv_cell,
                "corrected_energy_atom_eV": h_atom,
                "vib_energy_atom_eV": vib_atom,
                "g1800_atom_eV": h_atom + vib_atom,
            }
    return raw


def build_rows(raw):
    rows = []
    plot_data = {}
    for pressure in PRESSURES:
        plot_data[pressure] = {}
        for temp, key in [(0, "corrected_energy_atom_eV"), (1800, "g1800_atom_eV")]:
            fe_ref = raw[pressure]["Fe"][key]
            b_ref = raw[pressure]["boron"][key]
            points = []
            for phase in PHASES:
                item = raw[pressure][phase]
                formation = item[key] - ((1 - item["x_B"]) * fe_ref + item["x_B"] * b_ref)
                points.append({"phase": phase, "x_B": item["x_B"], "formation_eV_atom": formation})
            hull = lower_hull(points)
            for point in points:
                item = raw[pressure][point["phase"]]
                above = point["formation_eV_atom"] - interp_hull(hull, point["x_B"])
                rows.append(
                    {
                        "pressure_GPa": pressure,
                        "temperature_K": temp,
                        "phase": point["phase"],
                        "x_B": point["x_B"],
                        "formation_eV_atom": point["formation_eV_atom"],
                        "energy_above_hull_eV_atom": above,
                        "stable": abs(above) <= 1e-6,
                        "raw_energy_cell_eV": item["raw_energy_cell_eV"],
                        "pv_cell_eV": item["pv_cell_eV"],
                        "corrected_energy_atom_eV": item["corrected_energy_atom_eV"],
                        "vib_energy_atom_eV": 0.0 if temp == 0 else item["vib_energy_atom_eV"],
                        "total_energy_atom_eV": item[key],
                    }
                )
            plot_data[pressure][temp] = {"points": points, "hull": hull}
    return rows, plot_data


def write_csv(rows):
    with (OUT_DIR / "hull_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_one(temp, plot_data, out_png, out_pdf):
    title = PLOT_TITLE_0K if temp == 0 else PLOT_TITLE_1800K
    ylabel = Y_LABEL_0K if temp == 0 else Y_LABEL_1800K
    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE_SINGLE, sharey=True)
    all_y = [p["formation_eV_atom"] for pressure in PRESSURES for p in plot_data[pressure][temp]["points"]]
    ymin, ymax = min(all_y), max(all_y)
    pad = max(0.03, 0.15 * (ymax - ymin if ymax > ymin else 1.0))
    for ax, pressure in zip(axes, PRESSURES):
        points = plot_data[pressure][temp]["points"]
        hull = plot_data[pressure][temp]["hull"]
        ax.axhline(0, color="0.75", lw=0.8)
        ax.scatter([p["x_B"] for p in points], [p["formation_eV_atom"] for p in points], s=POINT_SIZE, c=[PHASE_COLOR[p["phase"]] for p in points], edgecolors="k", linewidths=0.4, zorder=3)
        ax.plot([p["x_B"] for p in hull], [p["formation_eV_atom"] for p in hull], color="black", lw=HULL_LINEWIDTH, zorder=2)
        for point in points:
            offset = ENDPOINT_LABEL_OFFSET if point["phase"] in ("Fe", "boron") else LABEL_OFFSET
            ax.annotate(formula_label(point["phase"]), (point["x_B"], point["formation_eV_atom"]), xytext=(0, offset), textcoords="offset points", ha="center", va="bottom" if offset > 0 else "top", fontsize=ANNOTATION_FONTSIZE)
        ax.set_title(f"{pressure} GPa", fontsize=SUBPLOT_TITLE_FONTSIZE)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.grid(True, alpha=GRID_ALPHA)
    axes[0].set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    fig.suptitle(title, fontsize=TITLE_FONTSIZE)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_DIR / out_png, dpi=DPI)
    fig.savefig(OUT_DIR / out_pdf)
    plt.close(fig)


def plot_combined(plot_data):
    fig, axes = plt.subplots(2, 3, figsize=FIGSIZE_COMBINED, sharex=True, sharey="row")
    for row_idx, temp in enumerate([0, 1800]):
        ylabel = Y_LABEL_0K if temp == 0 else Y_LABEL_1800K
        all_y = [p["formation_eV_atom"] for pressure in PRESSURES for p in plot_data[pressure][temp]["points"]]
        ymin, ymax = min(all_y), max(all_y)
        pad = max(0.03, 0.15 * (ymax - ymin if ymax > ymin else 1.0))
        for col_idx, pressure in enumerate(PRESSURES):
            ax = axes[row_idx, col_idx]
            points = plot_data[pressure][temp]["points"]
            hull = plot_data[pressure][temp]["hull"]
            ax.axhline(0, color="0.75", lw=0.8)
            ax.scatter([p["x_B"] for p in points], [p["formation_eV_atom"] for p in points], s=POINT_SIZE, c=[PHASE_COLOR[p["phase"]] for p in points], edgecolors="k", linewidths=0.4, zorder=3)
            ax.plot([p["x_B"] for p in hull], [p["formation_eV_atom"] for p in hull], color="black", lw=HULL_LINEWIDTH, zorder=2)
            for point in points:
                offset = ENDPOINT_LABEL_OFFSET if point["phase"] in ("Fe", "boron") else LABEL_OFFSET
                ax.annotate(formula_label(point["phase"]), (point["x_B"], point["formation_eV_atom"]), xytext=(0, offset), textcoords="offset points", ha="center", va="bottom" if offset > 0 else "top", fontsize=ANNOTATION_FONTSIZE)
            ax.set_title(f"{pressure} GPa", fontsize=SUBPLOT_TITLE_FONTSIZE)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(ymin - pad, ymax + pad)
            ax.tick_params(labelsize=TICK_FONTSIZE)
            ax.grid(True, alpha=GRID_ALPHA)
            if row_idx == 1:
                ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
            if col_idx == 0:
                ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    fig.suptitle(PLOT_TITLE_COMBINED, fontsize=TITLE_FONTSIZE)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / "hulls_0K_1800K.png", dpi=DPI)
    fig.savefig(OUT_DIR / "hulls_0K_1800K.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw = collect_data()
    rows, plot_data = build_rows(raw)
    write_csv(rows)
    plot_one(1800, plot_data, "hulls_1800K.png", "hulls_1800K.pdf")
    print(f"WROTE {OUT_DIR / 'hull_summary.csv'}")
    print(f"WROTE {OUT_DIR / 'hulls_1800K.pdf'}")
    print(f"WROTE {OUT_DIR / 'hulls_1800K.png'}")
