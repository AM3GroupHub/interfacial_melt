# Figure S7

Fe-B free-energy convex hulls at 1800 K and 0/5/10 GPa, using the MACE-MH-1 OMAT-PBE head. The original PDF is `reference.pdf`.

- `data/`: 24 RES files and 24 phonon `thermal_properties.yaml` files, covering eight phases at each pressure.
- `plot_pressure_hulls.py`: reads these files and builds the hulls.
- `reference/`: original script and `hull_summary.csv` for verification.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_pressure_hulls.py
```

Outputs: `output/hulls_1800K.pdf`, PNG, and the 0/1800 K numerical summary. Only local paths and figure-output selection were adapted. The original nominal-pressure PV correction and vibrational free-energy normalization are unchanged. The numerical summary matches the original exactly; minor rendering/spacing differences remain. No MACE or phonon calculation is needed to replot.

Source: `response/DFT_MACE/Fe-B-omat/`. These inputs are distinct from Fig. S9.
