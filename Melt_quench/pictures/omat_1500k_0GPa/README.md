# Reproducing the figures

For the 0 GPa `energy_composition_B_1500K` figure, this reads the shipped `data/csv/` table, no raw MD trajectories needed:

```bash
cd energy_composition_b_1500k
python plot_energy_composition_b_1500k.py
```

The raw source table used by the script is:

- `data/csv/energy_composition_b_1500k_raw.csv`

For the other components, see their subfolder READMEs.
