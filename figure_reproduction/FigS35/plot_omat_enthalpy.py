"""Replot S35 using the two original pressure-specific analysis modules."""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def load_source(pressure):
    path = ROOT / "source_scripts" / f"plot_{pressure}_1500k.py"
    spec = importlib.util.spec_from_file_location(f"omat_{pressure}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.SOURCE = ROOT / "data" / pressure / "energy_composition_b_1500k_raw.csv"
    return module


def main():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, pressure in zip(axes, ["0GPa", "10GPa"]):
        source = load_source(pressure)
        rows = source.normalize_rows(source.read_rows())
        rows.sort(key=lambda row: row["x_ref"])
        x = np.array([row["x_ref"] for row in rows], dtype=float)
        enthalpy = np.array([row["ent_norm"] for row in rows], dtype=float)
        source.plot_one_panel(
            ax,
            x,
            enthalpy,
            color_points="green",
            title=f"Enthalpy vs x(B) at 1500 K and {source.PRESSURE_LABEL}",
            ylabel="Normalized Enthalpy (eV/atom)",
        )

    fig.tight_layout()
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    pdf_path = out_dir / "omat_energy_reproduced.pdf"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(out_dir / "omat_energy_reproduced.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
