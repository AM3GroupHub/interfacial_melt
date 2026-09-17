# Figure S34

Endpoint-normalized PV term at 10 GPa for 1200, 1500, 1800, 2000, and 2300 K. Each curve is a fifth-order fit to 11 compositions; error bars include liquid and elemental-reference uncertainties combined in quadrature.

- `data/pv_composition_10GPa_liquid.csv`: the original 55-row plotting table.
- `data/source/`: 11 liquid energy tables and Fe250/B144 elemental-reference tables.
- `reference.pdf` and `reference/`: original figure and unmodified scripts.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_pv_composition_10gpa_liquid.py
```

PDF/PNG outputs go to `output/`. Optionally run `python recompute_pv_data.py` to regenerate the identical 55-row CSV in `output/`. Normalization subtracts the composition-weighted elemental PV values per atom. Only paths and summary-only execution were adapted; calculations and plotting are unchanged.

Historical input: `data/source/liquid/Fe216B144_G.txt` was copied from `final/10GPa_liquid/Fe216B144_G - 副本.txt`; the current file without this suffix does not reproduce the original data. Other source tables are from `final/10GPa_liquid/` and `final/10GPa_crystal/`; the CSV-only script is from `final/10GPa_pdf/pv/pv_composition_10gpa_liquid/`.
