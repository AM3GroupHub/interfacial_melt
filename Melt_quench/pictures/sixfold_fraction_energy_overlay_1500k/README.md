# Reproducing the figures

For the 0 GPa 1500 K sixfold-fraction versus potential overlay figure (`3b`) — this reads the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd sixfold_fraction_energy_overlay_1500k
python plot_sixfold_fraction_energy_overlay_1500k.py
```

The raw source tables used by the script are:

- `data/csv/sixfold_fraction_summary_1500k.csv`
- `data/csv/potential_energy_0gpa_1500k.csv`

For the other components, see their subfolder READMEs.
