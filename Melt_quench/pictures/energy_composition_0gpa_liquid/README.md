# Reproducing the figures

For the 0 GPa energy-composition figures, these read the shipped `data/csv/` table, no raw MD trajectories needed:

```bash
cd energy_composition_0gpa_liquid
python plot_energy_composition_0gpa_liquid.py
```

The raw source table used by the script is:

- `data/csv/normalized_energy_0GPa_liquid.csv`

For the other components, see their subfolder READMEs.
