"""Plot derivatives of quintic free-energy fits with pointwise OLS uncertainty."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "output"
KB = 8.617333262145e-5
PHASES = ((r"Fe$_{23}$B$_6$", 6.0 / 29.0), (r"Fe$_3$B", 1.0 / 4.0),
          (r"FeB", 1.0 / 2.0), (r"Fe$_2$B", 1.0 / 3.0),
          (r"FeB$_4$", 4.0 / 5.0))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7,
                         "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True,
                             sharey=True, constrained_layout=True)
    grid = np.linspace(0.1, 0.9, 801)
    # A centered coordinate and SVD avoid inverting the normal equations.
    ugrid = (grid - 0.5) / 0.4
    basis = np.column_stack([np.zeros_like(grid), np.zeros_like(grid)] +
                            [k * (k - 1) * ugrid**(k - 2) / 0.4**2
                             for k in range(2, 6)])
    exported = []
    for row, pressure in enumerate((0, 10)):
        source = ROOT / "data" / f"normalized_energy_{pressure}GPa_liquid.csv"
        with source.open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle))
        for col, temperature in enumerate((1500, 1800)):
            data = sorted((float(r["x_B"]), float(r["Free_norm_eV_per_atom"]))
                          for r in records if float(r["temperature_K"]) == temperature)
            x, free = np.array(data).T
            y = free + KB * temperature * (x * np.log(x) + (1-x) * np.log(1-x))
            design = np.polynomial.polynomial.polyvander((x - 0.5) / 0.4, 5)
            coeff, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
            assert rank == 6 and len(x) > 6
            residual = y - design @ coeff
            variance = residual @ residual / (len(x) - 6)
            inverse = np.linalg.pinv(design)
            curvature = basis @ coeff
            sigma = np.sqrt(variance * np.sum((basis @ inverse)**2, axis=1))
            original = np.polyfit(x, y, 5)
            np.testing.assert_allclose(curvature, np.polyval(np.polyder(original, 2), grid),
                                       atol=1e-8, rtol=1e-8)
            assert np.all(np.isfinite(sigma))
            ax = axes[row, col]
            color = ("#1f6b3b", "#5b2a86")[col]
            ax.fill_between(grid, curvature-sigma, curvature+sigma,
                            color=color, alpha=0.20, label=r"Pointwise $\pm1\sigma$")
            ax.plot(grid, curvature, color=color, lw=1.5, label=r"$G''(x)$ (quintic fit)")
            ax.axhline(0, color="0.3", lw=0.8, ls="--")
            for phase, phase_x in PHASES:
                line_style = ":" if phase == r"FeB$_4$" else "--"
                line_color = "0.25" if phase == r"FeB$_4$" else "0.45"
                line_width = 0.9 if phase == r"FeB$_4$" else 0.65
                ax.axvline(phase_x, color=line_color, lw=line_width,
                           ls=line_style, alpha=0.85)
                ax.text(phase_x, 0.98, phase, transform=ax.get_xaxis_transform(),
                        rotation=90, ha="right", va="top", fontsize=9,
                        color="0.35", clip_on=True)
            ax.set_title(f"({chr(97 + row*2 + col)}) {pressure} GPa, {temperature} K", loc="left")
            ax.set_xlim(0.1, 0.9)
            ax.set_xticks(np.arange(0.1, 1.0, 0.2))
            ax.tick_params(direction="in", top=True, right=True)
            if row == 1:
                ax.set_xlabel(r"Boron fraction $x_{\mathrm{B}}$")
            if col == 0:
                ax.set_ylabel(r"$\partial^2 G/\partial x_{\mathrm{B}}^2$ (eV/atom)")
            if row == col == 0:
                ax.legend(fontsize=8, frameon=False, loc="lower left")
            exported.extend(zip([pressure]*len(grid), [temperature]*len(grid),
                                grid, curvature, sigma, curvature-sigma, curvature+sigma))
            index = np.argmin(abs(grid - 0.8))
            print(f"{pressure} GPa {temperature} K: FeB4 G'' = {curvature[index]:.8f} +/- {sigma[index]:.8f}")
    for extension in ("png", "pdf"):
        suffix = "" if extension == "png" else "_phase_guides"
        target = OUT_DIR / f"curvature_fit_5th_order{suffix}.{extension}"
        fig.savefig(target, dpi=600)
        print(target)
    plt.close(fig)
    with (OUT_DIR / "curvature_with_uncertainty.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["P_GPa", "T_K", "x_B", "G_second_eV_per_atom",
                         "sigma_second_eV_per_atom", "lower_1sigma", "upper_1sigma"])
        writer.writerows(exported)


if __name__ == "__main__":
    main()
