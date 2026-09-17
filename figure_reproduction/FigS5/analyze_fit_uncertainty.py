import csv
from dataclasses import dataclass
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from plot_s_test_structure_factors import (
    KMAX_FIT,
    K_PLOT_MAX,
    ROOT,
    load_case,
    output_stem,
    oz,
    resolve_input_paths,
)


OUT_DIR = ROOT / "output"
PRESSURE_OUT_DIR = ROOT / "output"
SUMMARY_PATH = OUT_DIR / "fit_uncertainty.csv"
BOOTSTRAP_PATH = OUT_DIR / "fit_uncertainty_bootstrap.csv"
N_BOOTSTRAP = 5000
RNG_SEED = 20260817


@dataclass
class UncertaintyResult:
    k: np.ndarray
    s: np.ndarray
    popt: np.ndarray
    pcov: np.ndarray
    pair_boot: np.ndarray
    pair_failures: int
    resid_boot: np.ndarray
    resid_failures: int
    loo_rows: list


def initial_guess(k, s):
    a0 = float(s[np.argmin(k)])
    with np.errstate(divide="ignore", invalid="ignore"):
        bg = np.nanmedian((a0 / s - 1.0) / (k * k))
    b0 = bg if np.isfinite(bg) and bg > 0 else 1.0
    return [a0, b0]


def fit_oz_with_cov(k, s):
    popt, pcov = curve_fit(
        oz,
        k,
        s,
        p0=initial_guess(k, s),
        bounds=([0.0, 0.0], [np.inf, np.inf]),
        maxfev=10000,
    )
    return np.asarray(popt, float), np.asarray(pcov, float)


def fit_points(case):
    k = np.asarray(case["k"], float)
    s = np.asarray(case["scc"], float)
    mask = np.isfinite(k) & np.isfinite(s) & (s > 0) & (k <= KMAX_FIT)
    return k[mask], s[mask]


def parameter_errors(pcov):
    if pcov.shape != (2, 2) or not np.all(np.isfinite(pcov)):
        return np.nan, np.nan, np.nan
    variances = np.diag(pcov)
    if np.any(variances < 0):
        return np.nan, np.nan, np.nan
    se_a, se_b = np.sqrt(variances)
    corr = pcov[0, 1] / (se_a * se_b) if se_a > 0 and se_b > 0 else np.nan
    return float(se_a), float(se_b), float(corr)


def rmse(k, s, a, b):
    dof = max(len(k) - 2, 1)
    residuals = s - oz(k, a, b)
    return float(np.sqrt(np.sum(residuals * residuals) / dof))


def bootstrap_fits(k, s, rng):
    fits = []
    failures = 0
    n = len(k)
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, n)
        try:
            popt, _ = fit_oz_with_cov(k[idx], s[idx])
        except (RuntimeError, ValueError, FloatingPointError):
            failures += 1
            continue
        if np.all(np.isfinite(popt)):
            fits.append(popt)
        else:
            failures += 1
    return as_fit_array(fits), failures


def residual_bootstrap_fits(k, s, popt, rng):
    fits = []
    failures = 0
    fitted = oz(k, *popt)
    residuals = s - fitted
    residuals = residuals - np.mean(residuals)
    n = len(k)
    for _ in range(N_BOOTSTRAP):
        sample_s = fitted + rng.choice(residuals, n, replace=True)
        try:
            sample_popt, _ = fit_oz_with_cov(k, sample_s)
        except (RuntimeError, ValueError, FloatingPointError):
            failures += 1
            continue
        if np.all(np.isfinite(sample_popt)):
            fits.append(sample_popt)
        else:
            failures += 1
    return as_fit_array(fits), failures


def as_fit_array(fits):
    if not fits:
        return np.empty((0, 2), float)
    return np.asarray(fits, float)


def leave_one_out(k, s):
    rows = []
    for idx in range(len(k)):
        mask = np.ones(len(k), dtype=bool)
        mask[idx] = False
        try:
            popt, _ = fit_oz_with_cov(k[mask], s[mask])
        except (RuntimeError, ValueError, FloatingPointError):
            continue
        if np.all(np.isfinite(popt)):
            rows.append({
                "k_removed": float(k[idx]),
                "s_removed": float(s[idx]),
                "a": float(popt[0]),
                "b": float(popt[1]),
            })
    return rows


def quantiles(values):
    return np.percentile(values, [2.5, 50.0, 97.5])


def fmt(value):
    return f"{value:.10g}"


def interval_fields(prefix, values):
    if len(values) == 0:
        return {
            f"{prefix}_p025": "nan",
            f"{prefix}_p50": "nan",
            f"{prefix}_p975": "nan",
        }
    q025, q50, q975 = quantiles(values)
    return {
        f"{prefix}_p025": fmt(q025),
        f"{prefix}_p50": fmt(q50),
        f"{prefix}_p975": fmt(q975),
    }


def summarize_case(case, result):
    a, b = result.popt
    se_a, se_b, corr = parameter_errors(result.pcov)
    xi = float(np.sqrt(b))
    loo_a = np.asarray([row["a"] for row in result.loo_rows], float)
    loo_b = np.asarray([row["b"] for row in result.loo_rows], float)
    influence_a = max(result.loo_rows, key=lambda row: abs(row["a"] - a))
    influence_b = max(result.loo_rows, key=lambda row: abs(row["b"] - b))

    row = {
        "csv": str(case["path"].relative_to(ROOT)),
        "composition": f"{case['composition']:.8g}",
        "pressure_GPa": f"{case['pressure']:.8g}",
        "n_fit_points": len(result.k),
        "kmax_fit": KMAX_FIT,
        "a_full": fmt(a),
        "a_cov_se": fmt(se_a),
        "a_loo_min": fmt(np.min(loo_a)),
        "a_loo_max": fmt(np.max(loo_a)),
        "b_full": fmt(b),
        "b_cov_se": fmt(se_b),
        "b_loo_min": fmt(np.min(loo_b)),
        "b_loo_max": fmt(np.max(loo_b)),
        "xi_full_A": fmt(xi),
        "rmse": fmt(rmse(result.k, result.s, a, b)),
        "param_corr_cov": fmt(corr),
        "pair_bootstrap_success": len(result.pair_boot),
        "pair_bootstrap_failures": result.pair_failures,
        "resid_bootstrap_success": len(result.resid_boot),
        "resid_bootstrap_failures": result.resid_failures,
        "most_influential_k_for_a": fmt(influence_a["k_removed"]),
        "a_without_influential": fmt(influence_a["a"]),
        "most_influential_k_for_b": fmt(influence_b["k_removed"]),
        "b_without_influential": fmt(influence_b["b"]),
    }
    row.update(interval_fields("a_pair", result.pair_boot[:, 0]))
    row.update(interval_fields("a_resid", result.resid_boot[:, 0]))
    row.update(interval_fields("b_pair", result.pair_boot[:, 1]))
    row.update(interval_fields("b_resid", result.resid_boot[:, 1]))
    row.update(interval_fields("xi_pair_A", np.sqrt(result.pair_boot[:, 1])))
    row.update(interval_fields("xi_resid_A", np.sqrt(result.resid_boot[:, 1])))
    return row


def plot_case_uncertainty(case, k, s, popt, resid_boot):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    kfit = np.linspace(0.0, KMAX_FIT, 200)
    yfit = oz(kfit, *popt)
    ymax_values = [0.6, float(np.nanmax(case["scc"])), float(np.nanmax(yfit))]

    fig, ax = plt.subplots(figsize=(4.6, 3.5))
    ax.scatter(case["k"], case["scc"], s=22, color="0.78", edgecolor="none", label="all reciprocal")
    ax.scatter(k, s, s=34, color="#4C88C7", alpha=0.78, edgecolor="none", label="fit window")
    if len(resid_boot) > 0:
        yboot = oz(kfit[:, None], resid_boot[:, 0], resid_boot[:, 1])
        ylo, yhi = np.percentile(yboot, [2.5, 97.5], axis=1)
        ax.fill_between(
            kfit,
            ylo,
            yhi,
            color="#d95f02",
            alpha=0.20,
            lw=0,
            label="residual bootstrap 95% CI",
        )
        ax.plot(kfit, ylo, color="#d95f02", lw=0.7, alpha=0.65)
        ax.plot(kfit, yhi, color="#d95f02", lw=0.7, alpha=0.65)
        ymax_values.append(float(np.nanmax(yhi)))
    ax.plot(kfit, yfit, color="#d95f02", lw=1.8, ls="--", label="OZ fit")
    ax.axvspan(0, KMAX_FIT, color="0.8", alpha=0.10, lw=0)
    ax.set_xlim(-0.05, K_PLOT_MAX)
    ax.set_ylim(0.0, max(ymax_values) * 1.08)
    ax.set_xlabel(r"$|k|$ [$\mathrm{\AA}^{-1}$]")
    ax.set_ylabel(r"$S_{cc}(|k|)$")
    ax.set_title(rf"Fe-B, {case['pressure']:.1f} GPa, OZ uncertainty")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    out_png = OUT_DIR / f"{output_stem(case)}_oz_uncertainty.png"
    out_pdf = OUT_DIR / f"{output_stem(case)}_oz_uncertainty.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_png, out_pdf


def plot_pressure_result(ax, kfit, color, case, result):
    pressure = case["pressure"]
    ymax_values = [float(np.nanmax(case["scc"]))]
    ax.scatter(
        case["k"],
        case["scc"],
        s=22,
        color=color,
        alpha=0.42,
        edgecolor="none",
        label=f"{pressure:.1f} GPa data",
    )
    yfit = oz(kfit, *result.popt)
    ymax_values.append(float(np.nanmax(yfit)))
    if len(result.resid_boot) > 0:
        yboot = oz(kfit[:, None], result.resid_boot[:, 0], result.resid_boot[:, 1])
        ylo, yhi = np.percentile(yboot, [2.5, 97.5], axis=1)
        ax.fill_between(kfit, ylo, yhi, color=color, alpha=0.18, lw=0, label=f"{pressure:.1f} GPa 95% CI")
        ax.plot(kfit, ylo, color=color, lw=0.7, alpha=0.65)
        ax.plot(kfit, yhi, color=color, lw=0.7, alpha=0.65)
        ymax_values.append(float(np.nanmax(yhi)))
    ax.plot(kfit, yfit, color=color, lw=1.9, ls="--", label=f"{pressure:.1f} GPa OZ fit")
    return ymax_values


def plot_pressure_uncertainty(case_results):
    PRESSURE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.9, 3.6))
    colors = ["#4C88C7", "#d95f02"]
    kfit = np.linspace(0.0, KMAX_FIT, 200)
    ymax_values = [0.6]
    xprod = case_results[0][0]["x_a"] * case_results[0][0]["x_b"]
    ax.axhline(xprod, color="0.55", lw=1.2, ls=":", label=rf"$x_{{Fe}}x_{{B}}={xprod:.2f}$")

    sorted_results = sorted(case_results, key=lambda item: item[0]["pressure"])
    for color, (case, result) in zip(colors, sorted_results):
        ymax_values.extend(plot_pressure_result(ax, kfit, color, case, result))

    ax.axvspan(0, KMAX_FIT, color="0.8", alpha=0.08, lw=0)
    ax.set_xlim(-0.05, K_PLOT_MAX)
    ax.set_ylim(0.0, max(ymax_values) * 1.08)
    ax.set_xlabel(r"$|k|$ [$\mathrm{\AA}^{-1}$]")
    ax.set_ylabel(r"$S_{cc}(|k|)$")
    composition = case_results[0][0]["composition"]
    ax.set_title(rf"Fe-B, $x_B={composition:.2f}$, OZ fits with 95% CI")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    out_png = PRESSURE_OUT_DIR / f"single_Fe_B_scc_{composition:.1f}_pressure_compare_with_uncertainty.png"
    out_pdf = PRESSURE_OUT_DIR / f"single_Fe_B_scc_{composition:.1f}_pressure_compare_with_uncertainty.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_png, out_pdf


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def resolve_paths(argv):
    if argv:
        return resolve_input_paths(argv)
    # Preserve the original 0 GPa -> 10 GPa order and shared RNG stream.
    return [
        ROOT / "data" / "single_Fe_B_scc_0.8_0GPa.csv",
        ROOT / "data" / "single_Fe_B_scc_0.8_10GPa.csv",
    ]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    rng = np.random.default_rng(RNG_SEED)
    summary_rows = []
    bootstrap_rows = []
    case_results = []
    paths = resolve_paths(argv)
    for case in [load_case(path) for path in paths]:
        k, s = fit_points(case)
        popt, pcov = fit_oz_with_cov(k, s)
        pair_boot, pair_failures = bootstrap_fits(k, s, rng)
        resid_boot, resid_failures = residual_bootstrap_fits(k, s, popt, rng)
        loo_rows = leave_one_out(k, s)
        result = UncertaintyResult(
            k=k,
            s=s,
            popt=popt,
            pcov=pcov,
            pair_boot=pair_boot,
            pair_failures=pair_failures,
            resid_boot=resid_boot,
            resid_failures=resid_failures,
            loo_rows=loo_rows,
        )
        case_results.append((case, result))
        summary_rows.append(summarize_case(case, result))
        for method, boot in [("paired", pair_boot), ("residual", resid_boot)]:
            for idx, sample in enumerate(boot):
                bootstrap_rows.append({
                    "method": method,
                    "case": output_stem(case),
                    "sample": idx,
                    "a": f"{sample[0]:.10g}",
                    "b": f"{sample[1]:.10g}",
                    "xi_A": f"{np.sqrt(sample[1]):.10g}",
                })
        print(f"Fitted {case['pressure']:g} GPa: S_cc(0) = {popt[0]:.10g}", flush=True)

    if len(case_results) > 1:
        out_png, _ = plot_pressure_uncertainty(case_results)
        print(f"Wrote {out_png}")

    write_csv(SUMMARY_PATH, summary_rows, list(summary_rows[0]))
    write_csv(BOOTSTRAP_PATH, bootstrap_rows, ["method", "case", "sample", "a", "b", "xi_A"])
    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote {BOOTSTRAP_PATH}")


if __name__ == "__main__":
    main()
