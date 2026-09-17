# Figures S12-S16

S12, S13, S14, S15, and S16 correspond to 1200, 1500, 1800, 2000, and 2300 K, respectively. Each figure compares 0 GPa (upper row) and 10 GPa (lower row): potential energy, enthalpy, and `H - T S_vib` from left to right.

- `data/`: the two original normalized-energy CSVs, including uncertainties (55 rows per pressure).
- `plot_liquid_energy_by_temperature_pressure.py`: shared script for all five figures.
- `reference/`: the five original PDFs and unmodified source script.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_liquid_energy_by_temperature_pressure.py
```

The script writes `output/energy_0GPa_10GPa_<T>K_combined.pdf` and PNG for all five temperatures. Each panel uses 11 compositions, the supplied error columns, and a fifth-order polynomial fit. Only input/output paths were adapted; data, fitting, and plotting are unchanged.

Sources: `final/plot_liquid_energy_by_temperature_pressure.py`, `final/0GPa_liquid/normalized_energy_0GPa_liquid.csv`, `final/10GPa_liquid/normalized_energy_10GPa_liquid.csv`, and `final/liquid_energy_by_temperature_pressure/`.
