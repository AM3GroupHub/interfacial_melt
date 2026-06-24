"""Extrapolation of S_cc(k) -> S_cc(0): the M2 'OZ-like' route (the only method kept).

M2 (ozfit_direct): single OZ Lorentzian  S_cc(k) = a/(1 + b k^2),  a,b >= 0,  fit
NONLINEARLY in S-space (not 1/S) over |k| < KMAX_FIT, weighted by shell multiplicity
(record['w']).  S_cc(0) = a.  Ported from the updated Gamma_route make_state_diagnostics.py.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import curve_fit

from . import config
KMAX_FIT = config.CONFIG["fit"]["kmax_fit"]     # fit window |k| < KMAX_FIT (Å^-1)


def _window_mask(uk):
    return np.asarray(uk) < KMAX_FIT


def window_kmax(uk):
    """Largest |k| in the fit window (for drawing the fit curve)."""
    uk = np.asarray(uk)
    m = _window_mask(uk)
    return float(np.max(uk[m])) if np.any(m) else KMAX_FIT


def oz_nonlin(uk, S, w=None):
    """Nonlinear LSQ fit of S(k) = a/(1 + b k^2), a,b >= 0, fitting S DIRECTLY over
    the window |k| < KMAX_FIT (weighted). Returns (a, b) with S(0)=a; (nan,nan) if
    fewer than 3 usable points or the fit fails."""
    uk = np.asarray(uk); S = np.asarray(S)
    good = np.isfinite(S) & (S > 0) & _window_mask(uk)
    u, s = uk[good], S[good]
    if len(u) < 3:
        return np.nan, np.nan
    ww = np.asarray(w)[good] if w is not None else np.ones_like(u)
    a0 = float(s[np.argmin(u)])                          # S at smallest |k| as S(0) guess
    with np.errstate(divide="ignore", invalid="ignore"):
        bg = np.nanmedian((a0 / s - 1.0) / (u ** 2))     # rough OZ slope for b
    b0 = bg if (np.isfinite(bg) and bg > 0) else 1.0
    try:
        popt, _ = curve_fit(lambda k, a, b: a / (1.0 + b * k ** 2), u, s,
                            p0=[a0, b0], sigma=1.0 / ww,
                            bounds=([0.0, 0.0], [np.inf, np.inf]), maxfev=10000)
    except Exception:
        return np.nan, np.nan
    return float(popt[0]), float(popt[1])


def ozfit_direct(record):
    """M2 S_cc(0): the intercept a of the OZ Lorentzian a/(1+b k^2) fit to S_cc."""
    a, _ = oz_nonlin(record["uk"], record["scc"], record.get("w"))
    return a if (np.isfinite(a) and a > 0) else np.nan


def ozfit_direct_curve(record, kgrid):
    """The M2 fit curve a/(1 + b k^2) on kgrid; None if the fit failed."""
    a, b = oz_nonlin(record["uk"], record["scc"], record.get("w"))
    if not (np.isfinite(a) and np.isfinite(b)):
        return None
    k = np.asarray(kgrid, float)
    return a / (1.0 + b * k ** 2)


# only M2 is kept
METHODS = {"ozfit_direct": ozfit_direct}
