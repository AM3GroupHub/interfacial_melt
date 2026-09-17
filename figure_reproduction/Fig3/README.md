# Figure 3

- **a:** `panel_a/`, sixfold-B fractions at 0/10 GPa and potential energy at 0 GPa, 1800 K. The fit uses all 11 composition points; five B-rich points are displayed.
- **b:** `panel_b/`, the structure screenshot, original and relabeled POSCARs, matching VESTA file, coordination table, and analysis/labeling scripts. `C` labels sixfold-coordinated B (B-B cutoff 2.37 angstrom), **not carbon**. There are 54 Fe and 486 B atoms; 270 B atoms carry the `C` label.
- **c:** `panel_c/`, potential energy with/without the PV contribution at 10 GPa and 1500 K.

Replot a and c from this folder (Python 3.9+):

```bash
python -m pip install -r requirements.txt
python panel_a/plot_sixfold_fraction_energy_overlay_1800k.py
python panel_c/plot_10gpa_potential_rightonly.py
```

PDF/PNG outputs go to each panel's `output/`; original PDFs are retained as `reference.pdf`. The input tables match the original PDF curves (other table versions differ at `x_B=0.4`). Numerical methods are unchanged; a's axis labels match the reference. Replots have minor spacing differences.

For b, open `POSCAR-0-1800_C6.vesta` in VESTA. To rebuild the labels from the included coordination table, run inside `panel_b/`:

```bash
python substitute_b_by_coordination.py --poscar POSCAR-0-1800 --coordination-csv b_network_per_atom_POSCAR-0-1800.csv --output POSCAR_C6_rebuilt --mapping 6:C
```

The coordination table can be recalculated with `python analyze_b_network.py --poscar POSCAR-0-1800 --cutoff 2.37 --label rebuilt`.

Sources: a from `final/3c/` and its `sixfold_fraction_energy_overlay_1800k/` package; b from `final/coordinate/`; c script/PDF from `final/4b/`, with matching tables from `data/Melt_quench/pictures/pv_potential_overlay_10gpa_1500k_rightonly/`.
