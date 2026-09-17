# Figure 2

- **a:** `panel_a/`, 0 GPa at 1500 and 1800 K.
- **b:** `panel_b/`, 10 GPa at the same temperatures.
- Each panel includes its plotting script, two input CSVs in `data/`, and the original PDF (`reference.pdf`).

The energy CSV supplies the normalized vibrational free energy; the script adds ideal mixing entropy and fits a fifth-order polynomial. The `dgmix_recon_*` CSV supplies the structure-factor-derived excess curvature used in the lower plots.

From this folder, using Python 3.9 or newer:

```bash
python -m pip install -r requirements.txt
python panel_a/plot_free_energy_d2gdx2_0gpa_with_config_entropy.py
python panel_b/plot_free_energy_d2gdx2_10gpa_with_config_entropy-final.py
```

Each script writes PDF and PNG files to its own `output/` folder. Data paths are local; fitting and plotting formulas are unchanged. Panel a now exports PDF instead of EPS.

Both scripts have been run successfully. The regenerated PDFs have minor spacing differences from the preserved reference PDFs.

Sources within the FeB project: scripts/reference PDFs from `final/a_new/` and `final/b_new/`; input CSVs from `final/0GPa_liquid/`, `final/10GPa_liquid/`, and `final/d2g_dx2/`.
