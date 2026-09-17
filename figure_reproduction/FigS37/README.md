# Figure S37

Fe-B (upper row) and Cr-B (lower row) free energies at 0 GPa, with 1500 K on the left and 1800 K on the right. Each panel adds ideal configurational entropy to the endpoint-normalized vibrational free energy, then fits all 11 compositions with a fifth-order polynomial. No error bars are plotted.

- `data/`: original Fe-B and Cr-B normalized-energy CSVs; `source_crb/` contains the 11 Cr-B energy tables and two elemental-reference tables.
- `reference.pdf` and `reference/`: original figure and unmodified scripts.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_combined_feb_crb_free_energy_0gpa_1500k_1800k.py
```

PDF/PNG outputs go to `output/`. Optionally run `python recompute_crb_data.py` to regenerate the identical 22-row Cr-B CSV in `output/`. Only paths, export formats, and summary-only execution were adapted; fitting and plotting are unchanged.

Sources: `final/Cr-B-2/` and `final/0GPa_liquid/normalized_energy_0GPa_liquid.csv`. Use this Fe-B table, not its `updated` variant, to reproduce this figure.
