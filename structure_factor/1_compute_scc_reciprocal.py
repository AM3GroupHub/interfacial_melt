#!/usr/bin/env python
"""Item 2 — reciprocal-space S_cc(k) at commensurate shells k = 2πn/L.

Ŝ_αβ(k) = (1/V)⟨ρ_α(k) ρ_β*(k)⟩ from scaled coordinates with integer modes
(NPT-box-safe), shell-averaged, projected to S_cc(k). --last-frac windows the
strided frames. Writes the unified record schema to results/scc_reciprocal/.
"""
from __future__ import annotations
import argparse
from collections import Counter
import numpy as np

from src import common
from src import config
from src import conventions as cv
from src.mdio import iter_frames


def build_k_grid(L, n_shells=4):
    ns = np.arange(-n_shells, n_shells + 1)
    nx, ny, nz = np.meshgrid(ns, ns, ns, indexing="ij")
    mask = (nx ** 2 + ny ** 2 + nz ** 2) > 0
    n_vec = np.stack([nx[mask], ny[mask], nz[mask]], axis=-1).astype(np.float64)
    k_mag = np.linalg.norm(2.0 * np.pi * n_vec / np.asarray(L, float), axis=1)
    return n_vec, k_mag


def compute_S_at_k(dump_path, type_to_species, n_vec, stride, max_frames,
                   frame_range=(0.0, 1.0)):
    type_Fe = [int(t) for t, s in type_to_species.items() if s == "Fe"]
    type_B = [int(t) for t, s in type_to_species.items() if s == "B"]
    n_k = len(n_vec); two_pi = 2.0 * np.pi
    lo, hi = frame_range; windowed = lo > 0.0 or hi < 1.0
    if windowed:
        frames = list(iter_frames(dump_path, stride=stride, want_vel=False))
        n = len(frames); frames = frames[int(round(lo * n)):int(round(hi * n))]
        if max_frames is not None:
            frames = frames[:max_frames]
    else:
        frames = iter_frames(dump_path, stride=stride, want_vel=False)
    S_sum = np.zeros((n_k, 2, 2), dtype=complex); V_sum = 0.0; n_frames = 0
    for fr in frames:
        if not windowed and max_frames is not None and n_frames >= max_frames:
            break
        types, pos, L = fr.atom_type, fr.pos, np.asarray(fr.L, float)
        V = float(L.prod()); s_all = pos / L
        rhos = np.zeros((n_k, 2), dtype=complex)
        for spec_idx, tlist in enumerate([type_Fe, type_B]):
            m = np.isin(types, tlist)
            if not m.any():
                continue
            phase = np.exp(-1j * two_pi * (n_vec @ s_all[m].T))
            rhos[:, spec_idx] = phase.sum(axis=1)
        S_sum += (rhos[:, :, None] * rhos.conj()[:, None, :]) / V
        V_sum += V; n_frames += 1
    return S_sum / max(n_frames, 1), n_frames


def kernel(state, frame_range, stride=6, max_frames=2500, n_shells=4) -> dict:
    rho_Fe, rho_B = float(state["rho_Fe"]), float(state["rho_B"])
    rho = rho_Fe + rho_B; x_Fe = rho_Fe / rho
    L = np.asarray(state["L_mean"], float)
    n_vec, k_mag = build_k_grid(L, n_shells=n_shells)
    S_mat, n_frames = compute_S_at_k(common.resolve_dump(state), state["type_to_species"], n_vec,
                                     stride=stride, max_frames=max_frames,
                                     frame_range=frame_range)
    S = S_mat.real
    scc_k = np.array([common.scc_from_S(S[m], x_Fe, rho) for m in range(len(n_vec))])
    kr = np.round(k_mag, 2); uk = np.array(sorted(set(kr.tolist())))
    avg = lambda v: np.array([np.nanmean(v[kr == ku]) for ku in uk])
    cnt = Counter(kr.tolist())
    w = np.array([np.sqrt(cnt.get(round(float(k), 2), 1)) for k in uk])
    return dict(comp=str(state["comp"]), x_Fe=x_Fe, x_B=1.0 - x_Fe, rho_tot=rho,
                rho_Fe=rho_Fe, rho_B=rho_B, uk=uk, scc=avg(scc_k),
                S_FeFe=avg(S[:, 0, 0]), S_FeB=avg(S[:, 0, 1]), S_BB=avg(S[:, 1, 1]),
                w=w, n_frames=n_frames)


def main():
    cfg = config.CONFIG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--last-frac", type=float, default=cfg["last_frac"])
    ap.add_argument("--stride", type=int, default=cfg["reciprocal"]["stride"])
    ap.add_argument("--max-frames", type=int, default=cfg["reciprocal"]["max_frames"])
    ap.add_argument("--n-shells", type=int, default=cfg["reciprocal"]["n_shells"])
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--P", type=float, default=None)
    ap.add_argument("--recompute", action="store_true")
    args = ap.parse_args()
    Ts = [args.T] if args.T is not None else list(cfg["temperatures"])
    Ps = [args.P] if args.P is not None else list(cfg["pressures"])
    for T in Ts:
        for P in Ps:
            print(f"\n=== reciprocal S_cc(k)  T={int(T)} P={int(P)} ===")
            common.compute_grid(kernel, "scc_reciprocal", T, P, args.last_frac,
                                args.recompute, stride=args.stride,
                                max_frames=args.max_frames, n_shells=args.n_shells)


if __name__ == "__main__":
    main()
