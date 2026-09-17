# Figure S29

Sixfold-B fractions at 0/10 GPa and 0 GPa potential energy, at 1500 K. The fifth-order fit uses all 11 energy points; five B-rich compositions are displayed. The original plot does not show error bars.

- `data/`: original fraction summary (10 rows) and potential-energy table (11 rows).
- `reference.pdf` and `reference/`: original figure and unmodified source script.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_sixfold_fraction_energy_overlay_1500k.py
```

PDF/PNG outputs go to `output/`. Only input/output paths and export formats were adapted; numerical methods are unchanged. The replot matches the original data and fit, with minor spacing differences.

The included historical energy table matches the original PDF. Do not substitute the current main energy CSV: its `x_B=0.4` value differs and changes the fit.

Source: `final/3b/sixfold_fraction_energy_overlay_1500k/`; original PDF: `final/3b/sixfold_fraction_energy_overlay_1500K_3b.pdf`.
