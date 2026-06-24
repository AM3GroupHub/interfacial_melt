"""Partial radial distribution functions g_αβ(r) under PBC (vendored).

g_FeFe, g_FeB, g_BB, each → 1 at large r. `frame_range=(lo,hi)` selects a
sub-window of the strided trajectory (buffered); default (0,1) streams.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .mdio import iter_frames
from .species import SpeciesMap


@dataclass
class PartialRDF:
    r: np.ndarray
    g_FeFe: np.ndarray
    g_FeB: np.ndarray
    g_BB: np.ndarray
    n_frames: int
    r_max: float
    rho_Fe: float
    rho_B: float


def _pair_hist_AA(pos, L, edges):
    n = len(pos)
    if n < 2:
        return np.zeros(len(edges) - 1, dtype=np.float64)
    delta = pos[:, None, :] - pos[None, :, :]
    delta -= np.round(delta / L) * L
    d = np.sqrt((delta * delta).sum(axis=-1))
    iu, ju = np.triu_indices(n, k=1)
    return np.histogram(d[iu, ju], edges)[0].astype(np.float64)


def _pair_hist_AB(posA, posB, L, edges):
    if len(posA) == 0 or len(posB) == 0:
        return np.zeros(len(edges) - 1, dtype=np.float64)
    delta = posA[:, None, :] - posB[None, :, :]
    delta -= np.round(delta / L) * L
    d = np.sqrt((delta * delta).sum(axis=-1)).ravel()
    return np.histogram(d, edges)[0].astype(np.float64)


def compute_partial_rdf(dump_path, species_map: SpeciesMap, rho_Fe, rho_B, r_max,
                        dr=0.02, stride=1, max_frames=None,
                        frame_range=(0.0, 1.0)) -> PartialRDF:
    fe_t = species_map.fe_type()
    b_t = species_map.b_type()
    if fe_t is None and b_t is None:
        raise ValueError("Neither Fe nor B detected in species_map")

    edges = np.arange(0.0, r_max + dr, dr)
    centers = 0.5 * (edges[:-1] + edges[1:])
    nbins = len(centers)
    h_FeFe = np.zeros(nbins); h_FeB = np.zeros(nbins); h_BB = np.zeros(nbins)
    n_used = 0; L_sum = np.zeros(3); N_Fe_sum = N_B_sum = 0

    lo, hi = frame_range
    if lo > 0.0 or hi < 1.0:                        # buffer strided frames, slice window
        frames = list(iter_frames(dump_path, stride=stride, want_vel=False))
        n = len(frames)
        frames = frames[int(round(lo * n)):int(round(hi * n))]
        if max_frames is not None:
            frames = frames[:max_frames]
    else:
        frames = iter_frames(dump_path, stride=stride, max_frames=max_frames,
                             want_vel=False)

    for fr in frames:
        L = fr.L.astype(np.float64); L_sum += L
        pos = fr.pos.astype(np.float64)
        mask_fe = fr.atom_type == fe_t if fe_t is not None else np.zeros(fr.natoms, bool)
        mask_b = fr.atom_type == b_t if b_t is not None else np.zeros(fr.natoms, bool)
        pos_fe = pos[mask_fe]; pos_b = pos[mask_b]
        if fe_t is not None:
            h_FeFe += _pair_hist_AA(pos_fe, L, edges)
        if b_t is not None:
            h_BB += _pair_hist_AA(pos_b, L, edges)
        if fe_t is not None and b_t is not None:
            h_FeB += _pair_hist_AB(pos_fe, pos_b, L, edges)
        N_Fe_sum += int(mask_fe.sum()); N_B_sum += int(mask_b.sum()); n_used += 1

    if n_used == 0:
        raise RuntimeError(f"No frames read from {dump_path}")

    L_mean = L_sum / n_used
    N_Fe_mean = N_Fe_sum / n_used; N_B_mean = N_B_sum / n_used
    shell_vol = 4.0 * np.pi * centers ** 2 * (edges[1:] - edges[:-1])
    eps = 1e-30
    g_FeFe = (2.0 * h_FeFe / max(N_Fe_mean, eps)) / (n_used * shell_vol * max(rho_Fe, eps))
    g_BB = (2.0 * h_BB / max(N_B_mean, eps)) / (n_used * shell_vol * max(rho_B, eps))
    g_FeB = (h_FeB / max(N_Fe_mean, eps)) / (n_used * shell_vol * max(rho_B, eps))
    return PartialRDF(r=centers, g_FeFe=g_FeFe, g_FeB=g_FeB, g_BB=g_BB,
                      n_frames=n_used, r_max=float(r_max), rho_Fe=rho_Fe, rho_B=rho_B)
