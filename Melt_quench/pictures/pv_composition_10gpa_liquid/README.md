# Reproducing the figures

For the 10 GPa PV-composition figure, this reads the shipped `data/csv/` table, no raw MD trajectories needed:

```bash
cd pv_composition_10gpa_liquid
python plot_pv_composition_10gpa_liquid.py
```

The raw source table used by the script is:

- `data/csv/pv_composition_10GPa_liquid.csv`

For the other components, see their subfolder READMEs.
