# Reproducing the figures

For the 0 GPa sixfold-fraction and potential-energy overlay figures, these read the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd sixfold_fraction_energy_overlay_0gpa
python plot_sixfold_fraction_energy_overlay_0gpa.py
```

The raw source tables used by the script are:

- `data/csv/sixfold_fraction_summary.csv`
- `data/csv/normalized_energy_0GPa_liquid.csv`

Important: the fit uses the full 0 GPa potential-energy series from `normalized_energy_0GPa_liquid.csv`, while the markers are only shown for the selected six compositions.

For the other components, see their subfolder READMEs.
