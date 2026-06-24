#!/usr/bin/env python
"""Export per-figure plotted data to data/csv/ from the reciprocal/KBI caches.
Reruns the M2 OZ fit so each CSV carries the exact S_cc(0)/(a,b) used in the figures.
last_frac comes from config.yaml (0.5 -> the *_last0.5 caches)."""
from __future__ import annotations
import csv
import numpy as np
from src import common, fitting, config
from src import conventions as cv

CSV_DIR = cv.ROOT / "data" / "csv"
_, TAG = common.last_frac_to_window(config.CONFIG["last_frac"])
COMPS_2x2 = [0.6, 0.7, 0.8, 0.9]                       # Fig.3 / SM-KBI panels
PANELS_GEX = [(1500.0, 0.0), (1800.0, 0.0), (1500.0, 10.0), (1800.0, 10.0)]
SM_TP = [(1500, 0), (1500, 10), (1800, 0), (1800, 10)]


def _write(name, header, rows):
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    with open(CSV_DIR / name, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print(f"[csv] {CSV_DIR/name}  ({len(rows)} rows)")


def export_fig3(T=1500):
    rows = []
    for P in (0, 10):
        recs = {round(r["x_B"], 2): r for r in
                common.load_records(common.cache_path("scc_reciprocal", T, P, TAG))}
        for xb in COMPS_2x2:
            r = recs[xb]; s0 = fitting.ozfit_direct(r)
            s0s = f"{s0:.6e}" if np.isfinite(s0) else ""
            for k, scc in zip(np.asarray(r["uk"]), np.asarray(r["scc"])):
                rows.append((f"{xb:.1f}", P, f"{k:.4f}", f"{scc:.6e}", s0s))
    _write("fig3_scc_k.csv", ["x_B", "P_GPa", "k", "S_cc", "S_cc0_M2"], rows)


def export_fig1_gex():
    rows = []
    for T, P in PANELS_GEX:
        for r in sorted(common.load_records(common.cache_path("scc_reciprocal", T, P, TAG)),
                        key=lambda r: r["x_B"]):
            s0 = fitting.ozfit_direct(r)
            if not (np.isfinite(s0) and s0 > 0):
                continue
            xfe, xb = r["x_Fe"], r["x_B"]; G = xfe * xb / s0
            rows.append((int(T), int(P), f"{xb:.4f}", f"{G:.6e}",
                         f"{(G - 1.0) / (xfe * xb):.6e}"))
    _write("fig1_gex_curvature.csv", ["T_K", "P_GPa", "x_B", "Gamma", "d2gex_over_kBT"], rows)


def export_sm_state_diagnostics():
    rows = []
    for T, P in SM_TP:
        for r in sorted(common.load_records(common.cache_path("scc_reciprocal", T, P, TAG)),
                        key=lambda r: r["x_B"]):
            a, b = fitting.oz_nonlin(r["uk"], r["scc"], r.get("w"))
            af = f"{a:.6e}" if np.isfinite(a) else ""
            bf = f"{b:.6e}" if np.isfinite(b) else ""
            for k, scc in zip(np.asarray(r["uk"]), np.asarray(r["scc"])):
                rows.append((int(T), int(P), f"{r['x_B']:.4f}", f"{k:.4f}",
                             f"{scc:.6e}", af, bf, af))
    _write("sm_state_diagnostics.csv",
           ["T_K", "P_GPa", "x_B", "k", "S_cc", "a", "b", "S_cc0"], rows)


def export_sm_kbi(T=1500):
    rows = []
    for P in (0, 10):
        recs = {round(r["x_B"], 2): r for r in
                common.load_records(common.cache_path("scc_kbi", T, P, TAG))}
        for xb in COMPS_2x2:
            r = recs[xb]; s0 = float(np.asarray(r["scc"])[0])     # uk[0]==0
            for k, scc in zip(np.asarray(r["uk"]), np.asarray(r["scc"])):
                rows.append((f"{xb:.1f}", P, f"{k:.4f}", f"{scc:.6e}", f"{s0:.6e}"))
    _write("sm_kbi_scc_2x2.csv", ["x_B", "P_GPa", "k", "S_cc_KBI", "S_cc0_KBI"], rows)


def main():
    export_fig3(); export_fig1_gex(); export_sm_state_diagnostics(); export_sm_kbi()


if __name__ == "__main__":
    main()
