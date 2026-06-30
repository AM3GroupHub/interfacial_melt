# Reproducing the figures

For the 0 GPa 1800 K sixfold-fraction versus potential overlay figure (`3c`) — this reads the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd sixfold_fraction_energy_overlay_1800k
python plot_sixfold_fraction_energy_overlay_1800k.py
```

The raw source tables used by the script are:

- `data/csv/sixfold_fraction_summary_1800k.csv`
- `data/csv/potential_energy_0gpa_1800k.csv`

For the other components, see their subfolder READMEs.
