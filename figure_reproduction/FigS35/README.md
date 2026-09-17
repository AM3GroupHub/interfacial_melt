# Figure S35

The final figure combines **enthalpy at 0 GPa (left) and 10 GPa (right)**, both at 1500 K. It is not the potential-energy/enthalpy pair produced by either original script alone.

- `data/{0GPa,10GPa}/`: each original 11-row CSV and its 11 source TXT files.
- `source_scripts/`: unmodified pressure-specific analysis/plotting modules, used by the combined driver.
- `reference.pdf` and `reference/`: original final PDF, matching Illustrator file, and original TXT-based scripts.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_omat_enthalpy.py
```

PDF/PNG outputs go to `output/`. Original elemental-reference constants and sixth-order fits are retained; colors follow signed curvature (negative red, positive blue). The y-axis unit is clarified as eV/atom. No error bars are supplied.

The original PDF was assembled in Illustrator. Its leftmost 0 GPa segment is blue, whereas the original algorithm gives negative curvature and draws it red. The replot preserves the algorithm; this color and layout differ from the final PDF/AI.

Sources: `omat/0-l/` and `omat/10-l/`, including each `energy_composition_b_1500k/` package.
