# Reproducing the figures

For the energy and force comparison figures, these read the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd energy_force_comparison
python plot_comparison.py
```

The raw source tables used by the script are:

- `data/csv/energy_per_atom_comparison.csv`
- `data/csv/forces_comparison.csv`

For the other components, see their subfolder READMEs.
