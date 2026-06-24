"""Vendored project constants (LAMMPS metal units) and local data paths.

Self-contained: no feb_cdft import. Values copied from the validated upstream
feb_cdft.conventions on 2026-06-16.
"""
from __future__ import annotations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
STATE_INDEX = DATA_DIR / "state_index.json"

# Physical constants (positions Å, velocity Å/ps, mass amu, energy eV)
KB_EV_PER_K = 8.617333262e-5     # eV/K

# State grid
TEMPERATURES = [900.0, 1200.0, 1500.0, 1800.0, 2000.0, 2300.0, 2600.0]
PRESSURES_GPA = [0.0, 10.0]
