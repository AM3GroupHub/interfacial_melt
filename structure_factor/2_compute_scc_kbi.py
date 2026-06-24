#!/usr/bin/env python
"""Item 1 — real-space S_cc(k) via the windowed sinc-Fourier transform of RDFs.

    Ŝ_αβ(k) = ρ_α δ_αβ + ρ_α ρ_β G_αβ(k),
    G_αβ(k) = plateau over R∈[(1-frac)R_cut, R_cut] of
              4π ∫₀^R r² h_αβ(r) · sinc(kr) · w_KV(r/R_cut) dr,   h=g-1.
At k=0 sinc→1, so G_αβ(0) is exactly kbi.integrate_pair().G_KV; S_cc(0) is the
same as the validated KV route. RDFs are computed (and cached) by this script;
--last-frac selects the trajectory window.
"""
from __future__ import annotations
import argparse
import numpy as np

from src import common
from src import config
from src import conventions as cv
from src.kbi import kv_weight, integrate_pair
from src.rdf import compute_partial_rdf

RDF_DIR = common.RESULTS_DIR / "rdfs"
PLATEAU_FRAC = 0.20


def kbi_sinc_plateau(r, g, k: float, R_cut: float, plateau_frac=PLATEAU_FRAC) -> float:
    """Plateau of the running KV-windowed sinc integral of h=g-1 at wavevector k.
    Reduces exactly to integrate_pair(r,g,R_cut).G_KV when k==0."""
    r = np.asarray(r, float); h = np.asarray(g, float) - 1.0
    if k == 0.0:
        sinc = np.ones_like(r)
    else:
        kr = k * r; sinc = np.ones_like(r); nz = kr != 0.0
        sinc[nz] = np.sin(kr[nz]) / kr[nz]
    integrand = 4.0 * np.pi * r ** 2 * h * sinc * kv_weight(r, R_cut)
    dr = r[1] - r[0]
    run = np.cumsum(integrand) * dr - 0.5 * dr * (integrand[0] + integrand)
    idx = int(np.searchsorted(r, R_cut, side="right") - 1)
    idx = max(0, min(idx, len(r) - 1))
    i_lo = max(0, int(np.searchsorted(r, R_cut * (1.0 - plateau_frac))))
    return float(run[i_lo:idx + 1].mean())


def _get_rdf(state, frame_range, dr, stride, reuse_rdf):
    comp = str(state["comp"]); P = int(float(state["P_GPa"])); T = int(float(state["T_K"]))
    _, tag = common.last_frac_to_window(1.0 if frame_range == (0.0, 1.0)
                                        else 1.0 - frame_range[0])
    suffix = "" if tag == "full" else f"_{tag}"
    RDF_DIR.mkdir(parents=True, exist_ok=True)
    cache = RDF_DIR / f"{comp}_P{P}_T{T}{suffix}.npz"
    if reuse_rdf and cache.exists():
        d = np.load(cache, allow_pickle=False)
        return (d["r"], d["g_FeFe"], d["g_FeB"], d["g_BB"],
                float(d["rho_Fe"]), float(d["rho_B"]), np.asarray(d["L_mean"], float))
    L = np.asarray(state["L_mean"], float); r_max = float(L.min() / 2.0)
    rdf = compute_partial_rdf(common.resolve_dump(state), common.species_map_for(state),
                              state["rho_Fe"], state["rho_B"], r_max=r_max, dr=dr,
                              stride=stride, frame_range=frame_range)
    np.savez_compressed(cache, r=rdf.r, g_FeFe=rdf.g_FeFe, g_FeB=rdf.g_FeB,
                        g_BB=rdf.g_BB, rho_Fe=np.float64(rdf.rho_Fe),
                        rho_B=np.float64(rdf.rho_B), L_mean=L, n_frames=rdf.n_frames)
    return (rdf.r, rdf.g_FeFe, rdf.g_FeB, rdf.g_BB,
            float(rdf.rho_Fe), float(rdf.rho_B), L)


def kernel(state, frame_range, dr=0.02, stride=10, reuse_rdf=True, nk=150, kmax=3.0) -> dict:
    r, gff, gfb, gbb, rho_Fe, rho_B, L = _get_rdf(state, frame_range, dr, stride, reuse_rdf)
    rho = rho_Fe + rho_B; x_Fe = rho_Fe / rho
    R_cut = min(float(L.min() / 2.0), float(r[-1]))
    uk = np.linspace(0.0, kmax, nk)
    S_FeFe = np.empty(nk); S_FeB = np.empty(nk); S_BB = np.empty(nk); scc = np.empty(nk)
    for i, k in enumerate(uk):
        Gff = kbi_sinc_plateau(r, gff, k, R_cut)
        Gfb = kbi_sinc_plateau(r, gfb, k, R_cut)
        Gbb = kbi_sinc_plateau(r, gbb, k, R_cut)
        Shat = np.array([[rho_Fe + rho_Fe ** 2 * Gff, rho_Fe * rho_B * Gfb],
                         [rho_Fe * rho_B * Gfb, rho_B + rho_B ** 2 * Gbb]])
        S_FeFe[i] = Shat[0, 0]; S_FeB[i] = Shat[0, 1]; S_BB[i] = Shat[1, 1]
        scc[i] = common.scc_from_S(Shat, x_Fe, rho)
    return dict(comp=str(state["comp"]), x_Fe=x_Fe, x_B=1.0 - x_Fe, rho_tot=rho,
                rho_Fe=rho_Fe, rho_B=rho_B, uk=uk, scc=scc,
                S_FeFe=S_FeFe, S_FeB=S_FeB, S_BB=S_BB, w=np.ones(nk), n_frames=0)


def main():
    cfg = config.CONFIG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--last-frac", type=float, default=cfg["last_frac"])
    ap.add_argument("--stride", type=int, default=cfg["kbi"]["stride"])
    ap.add_argument("--dr", type=float, default=cfg["kbi"]["dr"])
    ap.add_argument("--kmax", type=float, default=cfg["kbi"]["kmax"])
    ap.add_argument("--nk", type=int, default=cfg["kbi"]["nk"])
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--P", type=float, default=None)
    ap.add_argument("--recompute", action="store_true")
    args = ap.parse_args()
    Ts = [args.T] if args.T is not None else list(cfg["temperatures"])
    Ps = [args.P] if args.P is not None else list(cfg["pressures"])
    for T in Ts:
        for P in Ps:
            print(f"\n=== KBI S_cc(k)  T={int(T)} P={int(P)} ===")
            common.compute_grid(kernel, "scc_kbi", T, P, args.last_frac, args.recompute,
                                dr=args.dr, stride=args.stride, reuse_rdf=True,
                                nk=args.nk, kmax=args.kmax)


if __name__ == "__main__":
    main()
