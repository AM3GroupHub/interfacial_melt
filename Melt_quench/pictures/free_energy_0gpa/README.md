# Reproducing the figures

For the 0 GPa free-energy figures (`a`, with and without configurational entropy) — these read the shipped `data/csv/` tables, no raw MD trajectories needed:

```bash
cd free_energy_0gpa
python plot_free_energy_0gpa_1500_1800.py
python plot_free_energy_0gpa_1500_1800_with_config_entropy.py
```

The raw source table used by both scripts is:

- `data/csv/free_energy_0gpa_1500_1800.csv`

For the other components, see their subfolder READMEs.
