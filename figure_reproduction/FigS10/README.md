# Figure S10

0 GPa free-energy curves at 1500 and 1800 K: `upper/` shows `H - T S_vib`; `lower/` additionally includes ideal configurational entropy.

Each panel contains its plotting script, input CSV in `data/`, and original PDF as `reference.pdf`.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python upper/plot_free_energy_0gpa_1500_1800.py
python lower/plot_free_energy_0gpa_1500_1800_with_config_entropy.py
```

Each script writes PDF/PNG to its own `output/`. Both select 11 compositions at each of the two temperatures and use fifth-order polynomial fits. The lower panels add `k_B T [x ln x + (1-x) ln(1-x)]` before fitting.

Sources: scripts/PDFs from `final/a/`; both input copies are `final/0GPa_liquid/normalized_energy_0GPa_liquid.csv`, verified against the original PDF points. The older portable CSV in `a/free_energy_0gpa/` does not match those points. Only paths, PDF/PNG export, and explicit Arial font selection were adapted; numerical methods are unchanged. Replots can differ in spacing from the original PDFs.
