from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parent
PRESSURE_DIRS = {
    "0 GPa": ROOT / "0GPa_liquid",
    "10 GPa": ROOT / "10GPa_liquid",
}
PLOT_TEMPERATURES = [1200.0, 1500.0, 1800.0, 2000.0, 2300.0]
PAIR_COLUMNS = [
    (("g_Fe_B", "g_B_Fe"), "Fe-B"),
    (("g_Fe_Fe",), "Fe-Fe"),
    (("g_B_B",), "B-B"),
]
TEMPERATURE_COLORS = {
    1200.0: "#37a9e6",
    1500.0: "#8ed7df",
    1800.0: "#f28e8e",
    2000.0: "#e25555",
    2300.0: "#b51f2e",
}
FONT = 16
REFERENCE_R = 2.37


def formulas_in_folder(folder):
    return sorted(path.stem.removesuffix("_G") for path in folder.glob("*_G.txt"))


def composition_label(formula):
    import re

    counts = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    total = counts.get("Fe", 0) + counts.get("B", 0)
    if total == 0:
        return formula
    x_fe = counts.get("Fe", 0) / total
    x_b = counts.get("B", 0) / total
    return rf"Fe$_{{{x_fe:.2f}}}$B$_{{{x_b:.2f}}}$"


def composition_tag(formula):
    import re

    counts = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    total = counts.get("Fe", 0) + counts.get("B", 0)
    if total == 0:
        return formula
    x_fe = counts.get("Fe", 0) / total
    x_b = counts.get("B", 0) / total
    return f"Fe{x_fe:.2f}B{x_b:.2f}"


def read_rdf(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().split()
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {name: data[:, index] for index, name in enumerate(header)}


def first_present_column(columns, names):
    for name in names:
        if name in columns:
            return name
    return None


def plot_formula(formula, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.8), sharex="col", sharey="col")
    has_any = False

    for row_index, (pressure_label, folder) in enumerate(PRESSURE_DIRS.items()):
        for temperature in PLOT_TEMPERATURES:
            rdf_path = folder / f"{formula}_dir" / f"{temperature:.1f}_dir" / "rdf.dat"
            if not rdf_path.exists():
                continue
            columns = read_rdf(rdf_path)
            r = columns["r_A"]
            for col_index, (column_names, title) in enumerate(PAIR_COLUMNS):
                ax = axes[row_index, col_index]
                column = first_present_column(columns, column_names)
                if column is None:
                    continue
                ax.plot(r, columns[column], color=TEMPERATURE_COLORS[temperature], linewidth=1.2, label=f"{temperature:.0f} K")
                if row_index == 0:
                    ax.set_title(title, fontsize=FONT)
                ax.tick_params(labelsize=FONT, direction="in", width=0.8, length=4, top=True, right=True)
                ax.grid(True, alpha=0.25)
                has_any = True

    if not has_any:
        plt.close(fig)
        return None

    for row_index, pressure_label in enumerate(PRESSURE_DIRS):
        axes[row_index, 0].set_ylabel(f"{pressure_label}\n$g(r)$", fontsize=FONT)
    for col_index in range(3):
        axes[1, col_index].set_xlabel(r"$r$ ($\AA$)", fontsize=FONT)
    for ax in axes.flat:
        ax.set_xlim(0.0, 5.0)
        ax.set_xticks(np.arange(1.0, 6.0, 1.0))
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
    for ax in axes[:, 2]:
        ax.axvline(REFERENCE_R, color="gray", linestyle="--", linewidth=1.2, alpha=0.8, zorder=1.5)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, fontsize=FONT, frameon=False, loc="center left", bbox_to_anchor=(0.98, 0.5))

    fig.suptitle(composition_label(formula), fontsize=FONT)
    fig.tight_layout(rect=[0, 0, 0.92, 0.95], pad=0.5, w_pad=0.5, h_pad=0.5)
    out_path = out_dir / f"rdf_{composition_tag(formula)}_0GPa_10GPa_combined.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"rdf_{composition_tag(formula)}_0GPa_10GPa_combined.pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    formulas = sorted(set(formulas_in_folder(PRESSURE_DIRS["0 GPa"])) & set(formulas_in_folder(PRESSURE_DIRS["10 GPa"])))
    out_dir = ROOT / "liquid_rdf_by_composition_combined"
    out_dir.mkdir(exist_ok=True)
    outputs = []
    for formula in formulas:
        out = plot_formula(formula, out_dir)
        if out is not None:
            outputs.append(out)
            print(f"Wrote {out}")
    print(f"Generated {len(outputs)} combined RDF figures")


if __name__ == "__main__":
    main()
