"""Shared backbone for scc_pipeline (standalone; imports only vendored modules).

Central identity: S_cc(0) = (1/ρ) zᵀ Ŝ(0) z, z=(x_B,−x_Fe);
Γ = x_Fe x_B / S_cc(0); g'' = k_B T / S_cc(0). Species order (Fe,B)=(0,1).
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from . import conventions as cv
from . import config
from .species import SpeciesMap

RESULTS_DIR = cv.ROOT / "results"
FIG_DIR = cv.ROOT / "figures"


def last_frac_to_window(F: float) -> tuple[tuple[float, float], str]:
    """'Use the last fraction F of the trajectory' -> (frame_range, tag)."""
    if not (0.0 < F <= 1.0):
        raise ValueError(f"last_frac must be in (0, 1], got {F}")
    if F == 1.0:
        return (0.0, 1.0), "full"
    return (1.0 - F, 1.0), f"last{F:g}"


def scc_from_S(S, x_Fe: float, rho_tot: float) -> float:
    """Bhatia–Thornton S_cc(0) = (1/ρ) zᵀ Ŝ z, z=(x_B,−x_Fe)."""
    x_B = 1.0 - x_Fe
    z = np.array([x_B, -x_Fe], dtype=np.float64)
    return float(z @ np.asarray(S, dtype=np.float64) @ z) / rho_tot


def gamma_from_scc(x_Fe: float, S_cc: float) -> float:
    return x_Fe * (1.0 - x_Fe) / S_cc


def g2_from_scc(T_K: float, S_cc: float) -> float:
    return cv.KB_EV_PER_K * T_K / S_cc


def load_state_index() -> list[dict]:
    rel = config.CONFIG.get("state_index", "data/state_index.json")
    return json.loads((cv.ROOT / rel).read_text())


def binary_states(T_K: float, P_GPa: float) -> list[dict]:
    """Binary (both densities > 0) states at (T,P), sorted by x_B ascending."""
    sel = [s for s in load_state_index()
           if float(s["rho_Fe"]) > 0 and float(s["rho_B"]) > 0
           and float(s["T_K"]) == T_K and float(s["P_GPa"]) == P_GPa]
    sel.sort(key=lambda s: 1.0 - float(s["rho_Fe"]) / (float(s["rho_Fe"]) + float(s["rho_B"])))
    return sel


def species_map_for(state: dict) -> SpeciesMap:
    return SpeciesMap(
        type_to_species={int(k): v for k, v in state["type_to_species"].items()},
        mean_v2={}, v2_ratio_B_over_Fe=state.get("v2_ratio_B_over_Fe"),
        expected_v2_ratio=state["expected_v2_ratio"], method=state["detection_method"])


def cache_path(method: str, T_K: float, P_GPa: float, tag: str) -> Path:
    d = RESULTS_DIR / method
    d.mkdir(parents=True, exist_ok=True)
    return d / f"T{int(T_K)}_P{int(P_GPa)}_{tag}.npz"


def save_records(path: Path, records: list[dict]) -> None:
    np.savez(path, records=np.array(records, dtype=object))


def load_records(path: Path) -> list[dict]:
    return list(np.load(path, allow_pickle=True)["records"])


def compute_grid(kernel, method: str, T_K: float, P_GPa: float, last_frac: float,
                 recompute: bool, **kernel_kw) -> list[dict]:
    """Run kernel(state, frame_range, **kw) over binary states at (T,P), cache, return."""
    frame_range, tag = last_frac_to_window(last_frac)
    cache = cache_path(method, T_K, P_GPa, tag)
    if cache.exists() and not recompute:
        print(f"[cache] {cache}")
        return load_records(cache)
    recs = []
    for s in binary_states(T_K, P_GPa):
        rec = kernel(s, frame_range, **kernel_kw)
        rec["frac_tag"] = tag
        recs.append(rec)
        print(f"  {rec['comp']:8s} x_B={rec['x_B']:.2f}  S_cc(0)~{rec['scc'][0]:+.4f}")
    save_records(cache, recs)
    print(f"[wrote] {cache}  ({len(recs)} comps)")
    return recs


def resolve_dump(state) -> Path:
    """Absolute path to a state's dump: data_root / <relative dump>. Used only by 1_/2_."""
    p = Path(state["dump"])
    if p.is_absolute():
        return p
    root = config.CONFIG.get("data_root") or ""
    return Path(root) / p
