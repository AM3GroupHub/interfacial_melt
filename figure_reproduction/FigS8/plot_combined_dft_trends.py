import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- User-editable plot settings ----------------
FIGURE_TITLE = "DFT calculations at r$^2$SCAN level for thermodynamic stability of FeB$_2$ and FeB$_4$"
LEFT_TITLE = "Formation enthalpy from DFT"
RIGHT_TITLE = "0 GPa formation free energy from DFT"
LEFT_X_LABEL = "Pressure (GPa)"
LEFT_Y_LABEL = "Formation enthalpy (meV/atom)"
RIGHT_X_LABEL = "Temperature (K)"
RIGHT_Y_LABEL = "Formation free energy (meV/atom)"
LEGEND_LABEL_FEB2 = "FeB$_2$"
LEGEND_LABEL_FEB4 = "FeB$_4$"
LEGEND_LABEL_REF = "FeB + B"
FIGSIZE = (12.0, 4.8)
DPI = 300
SUPTITLE_FONTSIZE = 16
TITLE_FONTSIZE = 16
AXIS_LABEL_FONTSIZE = 15
TICK_FONTSIZE = 16
LEGEND_FONTSIZE = 16
LINEWIDTH = 1.8
MARKER_SIZE = 5
GRID_ALPHA = 0.2
ZERO_LINE_COLOR = "0.7"
ZERO_LINE_WIDTH = 0.8
ZERO_LINE_STYLE = "--"
LEFT_COLOR_FEB2 = "#1f77b4"
LEFT_COLOR_FEB4 = "#d62728"
RIGHT_COLOR_FEB2 = "#1f77b4"
RIGHT_COLOR_FEB4 = "#d62728"
T_MIN = 0
T_MAX = 2000
T_STEP = 10
# -------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
PHONON_ROOT = ROOT / "data" / "0GPa" / "phonon"
OUT_DIR = ROOT / "output"
OUT_PNG = OUT_DIR / "combined_dft_trends.png"
OUT_PDF = OUT_DIR / "combined_dft_trends.pdf"
KJMOL_TO_EV = 96.4853321233

# Pressure trend data from DFT res relative to FeB + B (meV/atom)
PRESSURE_GPA = [0, 5, 10]
FEB2_PRESSURE = [39.787713, 18.848041, -1.101370]
FEB4_PRESSURE = [13.788283, 0.012111, -12.717167]

# 0 GPa DFT U+PV from res for temperature trend
H_FEB2 = -118.05136 / 12.0
H_FEB4 = -88.511229 / 10.0
H_FEB = -89.144107 / 8.0
H_B = -264.46235 / 36.0


def parse_thermal(path):
    natom = None
    current_t = None
    free_map = {}
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
            free_map[int(round(current_t))] = float(m.group(1))
            current_t = None
    if natom is None:
        raise ValueError(f"Missing natom in {path}")
    return natom, free_map


def vib_ev_atom(path):
    natom, free_map = parse_thermal(path)
    return {t: free_map[t] / KJMOL_TO_EV / natom for t in free_map}


def build_temperature_curves():
    feb2_vib = vib_ev_atom(PHONON_ROOT / "FeB2" / "thermal_properties.yaml")
    feb4_vib = vib_ev_atom(PHONON_ROOT / "FeB4" / "thermal_properties.yaml")
    feb_vib = vib_ev_atom(PHONON_ROOT / "FeB" / "thermal_properties.yaml")
    b_vib = vib_ev_atom(PHONON_ROOT / "Boron" / "thermal_properties.yaml")
    temperatures = [t for t in range(T_MIN, T_MAX + 1, T_STEP) if t in feb2_vib and t in feb4_vib and t in feb_vib and t in b_vib]
    feb2_curve = []
    feb4_curve = []
    for t in temperatures:
        g_feb2 = H_FEB2 + feb2_vib[t]
        g_feb4 = H_FEB4 + feb4_vib[t]
        g_feb = H_FEB + feb_vib[t]
        g_b = H_B + b_vib[t]
        d_feb2 = g_feb2 - ((2.0 / 3.0) * g_feb + (1.0 / 3.0) * g_b)
        d_feb4 = g_feb4 - ((2.0 / 5.0) * g_feb + (3.0 / 5.0) * g_b)
        feb2_curve.append(d_feb2 * 1000.0)
        feb4_curve.append(d_feb4 * 1000.0)
    return temperatures, feb2_curve, feb4_curve


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    temperatures, feb2_temp, feb4_temp = build_temperature_curves()
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)

    # Left panel: pressure
    ax = axes[0]
    ax.axhline(0, color=ZERO_LINE_COLOR, lw=ZERO_LINE_WIDTH, linestyle=ZERO_LINE_STYLE, label=LEGEND_LABEL_REF)
    ax.plot(PRESSURE_GPA, FEB2_PRESSURE, marker="o", markersize=MARKER_SIZE, linewidth=LINEWIDTH, color=LEFT_COLOR_FEB2, label=LEGEND_LABEL_FEB2)
    ax.plot(PRESSURE_GPA, FEB4_PRESSURE, marker="o", markersize=MARKER_SIZE, linewidth=LINEWIDTH, color=LEFT_COLOR_FEB4, label=LEGEND_LABEL_FEB4)
    ax.set_title(LEFT_TITLE, fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(LEFT_X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(LEFT_Y_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.grid(True, alpha=GRID_ALPHA)
    ax.legend(fontsize=LEGEND_FONTSIZE)

    # Right panel: temperature
    ax = axes[1]
    ax.axhline(0, color=ZERO_LINE_COLOR, lw=ZERO_LINE_WIDTH, linestyle=ZERO_LINE_STYLE, label=LEGEND_LABEL_REF)
    ax.plot(temperatures, feb2_temp, linewidth=LINEWIDTH, color=RIGHT_COLOR_FEB2, label=LEGEND_LABEL_FEB2)
    ax.plot(temperatures, feb4_temp, linewidth=LINEWIDTH, color=RIGHT_COLOR_FEB4, label=LEGEND_LABEL_FEB4)
    ax.set_title(RIGHT_TITLE, fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(RIGHT_X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(RIGHT_Y_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.grid(True, alpha=GRID_ALPHA)
    ax.legend(fontsize=LEGEND_FONTSIZE)

    fig.suptitle(FIGURE_TITLE, fontsize=SUPTITLE_FONTSIZE)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PNG, dpi=DPI)
    fig.savefig(OUT_PDF)
    print(f"WROTE {OUT_PNG}")
    print(f"WROTE {OUT_PDF}")


if __name__ == "__main__":
    main()
