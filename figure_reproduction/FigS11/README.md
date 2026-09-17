# Figure S11

10 GPa free-energy curves at 1500 and 1800 K: `upper/` shows `H - T S_vib`; `lower/` additionally includes ideal configurational entropy.

Each panel contains its plotting script, input CSV in `data/`, and original PDF as `reference.pdf`.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python upper/plot_free_energy_10gpa_1500_1800.py
python lower/plot_free_energy_10gpa_1500_1800_with_config_entropy.py
```

Each script writes PDF/PNG to its own `output/`. Fits use 11 compositions per temperature and fifth-order polynomials; the lower panels add `k_B T [x ln x + (1-x) ln(1-x)]` before fitting.

Sources: scripts/PDFs from `final/b/`. **The upper and lower PDFs use different historical data at `x_B = 0.4`; do not interchange the CSVs.** Upper uses `final/10GPa_liquid/normalized_energy_10GPa_liquid.csv`; lower uses `final/b/free_energy_10gpa/data/csv/free_energy_10gpa_1500_1800.csv`. Both were verified against their PDF points and curves.

Only paths, PDF/PNG export, and explicit font selection (Arial above, DejaVu Sans below) were adapted; numerical methods are unchanged. The lower script's duplicated EPS export was corrected to PDF. Replots can differ in spacing from the original PDFs.
