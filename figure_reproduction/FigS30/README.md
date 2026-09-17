# Figure S30

- **a:** `panel_a/`, sixfold-B fractions at 0/10 GPa and 0 GPa potential energy at 1800 K; same source as main-text Fig. 3a. Includes the matching historical CSVs, plotting script, and original PDF. The fifth-order fit uses all 11 energy points; five B-rich points are shown.
- **b:** `panel_b/POSCAR-0-1800-Fe15B85_C6` and matching PNG.
- **c:** `panel_c/POSCAR-0-1800_C6` and matching PNG.
- **d:** `panel_d/POSCAR-10-1800_C6` and matching PNG.

The extensionless files are POSCARs and can be opened in VESTA. In b-d, **C labels sixfold-coordinated B, not carbon**. Screenshots and structures are copied unchanged.

Replot a from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python panel_a/plot_sixfold_fraction_energy_overlay_1800k.py
```

PDF/PNG outputs go to `panel_a/output/`. Numerical methods are unchanged; paths, export formats, and axis-label text were adapted to match the reference. Minor spacing differences remain. Keep the supplied energy table: other versions differ at `x_B=0.4`.

Sources: a from `final/3c/`; b-d from `final/coordinate/`.
