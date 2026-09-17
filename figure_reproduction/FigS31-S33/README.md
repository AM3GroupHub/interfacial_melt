# Figures S31-S33

S31, S32, and S33 correspond to 1200, 1500, and 1800 K. Bars show sixfold-B fractions at 0/10 GPa (mean +/- SEM); dots show 0 GPa potential energy with its supplied errors. The fifth-order fit uses all 11 energy compositions; six B-rich compositions are displayed.

- `data/`: original fraction summary, matching energy table, and 360 per-frame coordination CSVs (10 frames per state).
- `reference/`: three original PDFs, original scripts, and the original overlay summary.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_sixfold_fraction_with_energy_overlay.py
```

PDF/PNG plots and an overlay summary go to `output/`. Optionally run `python recompute_sixfold_fraction_summary.py` to regenerate `output/sixfold_fraction_summary.csv` from the coordination files. It uses the per-frame fraction of B with CN=6, then the mean and sample standard deviation divided by sqrt(10); all 36 states match the supplied summary to its recorded precision.

The original PDFs match `normalized_energy_0GPa_liquid_updated.csv`, not the current non-`updated` table. Keep the included version. Only paths, selected temperatures, and summary-only execution were adapted; numerical methods and plotting are unchanged.

Sources: `final/plot_sixfold_fraction_with_energy_overlay.py`, `final/plot_sixfold_fraction_bar.py`, `final/sixfold_fraction_bar/`, `final/0GPa_liquid/normalized_energy_0GPa_liquid_updated.csv`, and `final/{0,10}GPa_coor/coor/`.
