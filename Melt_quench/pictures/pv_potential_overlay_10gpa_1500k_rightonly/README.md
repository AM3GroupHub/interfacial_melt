# Reproducing the figures

For the 10 GPa 1500 K right-only potential versus potential+PV overlay figure (`4b`, right panel) — this reads the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd pv_potential_overlay_10gpa_1500k_rightonly
python plot_pv_potential_overlay_10gpa_1500k_rightonly.py
```

The raw source tables used by the script are:

- `data/csv/pv_composition_10gpa_1500k.csv`
- `data/csv/potential_energy_10gpa_1500k.csv`

For the other components, see their subfolder READMEs.
