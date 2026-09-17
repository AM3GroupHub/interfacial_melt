# Figure S2

Reference-versus-predicted energy and force comparisons:

- `data/energy_per_atom_comparison.csv`: energy per atom, in eV/atom.
- `data/forces_comparison.csv`: individual force components, in eV/angstrom.
- `reference/`: original `Energy_plot.pdf` and `Forces_plot.pdf`.

Run from this folder with Python 3.9+:

```bash
python -m pip install -r requirements.txt
python plot_comparison.py
```

The script calculates RMSE/MAE and writes both figures as PDF/PNG to `output/`. Arial is used for the original appearance. The original 6-by-6-inch PDF size is preserved.

Sources: CSVs and original PDFs from `test-mixed/`; the existing portable plotting script from `test-mixed/energy_force_comparison/`. Only local paths and PDF output margins were adapted.
