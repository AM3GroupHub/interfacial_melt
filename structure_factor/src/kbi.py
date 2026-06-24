"""Kirkwood–Buff integrals with the Krüger–Vlugt finite-size window (vendored).

G^KV(R) = 4π ∫₀^R r² [g(r)-1] · w(r/R) dr,   w(x)=1-(3/2)x+(1/2)x³.
Production estimate = plateau (mean over the last 20% of R).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def kv_weight(r: np.ndarray, R_cut: float) -> np.ndarray:
    """Krüger–Vlugt window: w(r)=1-(3/2)(r/R_cut)+(1/2)(r/R_cut)³, w(R_cut)=0."""
    x = np.clip(r / R_cut, 0.0, 1.0)
    return 1.0 - 1.5 * x + 0.5 * x ** 3


def kbi_running(r: np.ndarray, g: np.ndarray,
                weight: np.ndarray | None = None) -> np.ndarray:
    """Running KB integral G(R) = 4π ∫₀^R r² [g-1] w dr (cumulative trapezoid)."""
    h = g - 1.0
    w = weight if weight is not None else 1.0
    integrand = 4.0 * np.pi * r ** 2 * h * w
    dr = r[1] - r[0]
    return np.cumsum(integrand) * dr - 0.5 * dr * (integrand[0] + integrand)


@dataclass
class PairKBI:
    G_raw: float
    G_KV: float
    R_plateau: float
    disagreement: float
    R_grid: np.ndarray
    G_raw_run: np.ndarray
    G_KV_run: np.ndarray


def integrate_pair(r: np.ndarray, g: np.ndarray, R_cut: float | None = None,
                   plateau_window_frac: float = 0.20) -> PairKBI:
    """Raw and KV KBI for one channel; production value = plateau mean."""
    if R_cut is None:
        R_cut = float(r[-1])
    G_raw_run = kbi_running(r, g, weight=None)
    G_KV_run = kbi_running(r, g, weight=kv_weight(r, R_cut))
    idx = int(np.searchsorted(r, R_cut, side="right") - 1)
    idx = max(0, min(idx, len(r) - 1))
    R_lo = R_cut * (1.0 - plateau_window_frac)
    i_lo = max(0, int(np.searchsorted(r, R_lo)))
    G_raw = float(G_raw_run[i_lo:idx + 1].mean())
    G_KV = float(G_KV_run[i_lo:idx + 1].mean())
    eps = 1e-12
    disag = abs(G_raw - G_KV) / max(abs(G_KV), eps)
    return PairKBI(G_raw=G_raw, G_KV=G_KV, R_plateau=float(r[idx]),
                   disagreement=disag, R_grid=r,
                   G_raw_run=G_raw_run, G_KV_run=G_KV_run)
