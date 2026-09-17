# Figures S18-S28

Shared RDF plots: 0 GPa (upper row), 10 GPa (lower row); Fe-B, Fe-Fe, and B-B (left to right), at 1200, 1500, 1800, 2000, and 2300 K.

| Figure | Composition | Data directory |
|---|---|---|
| S18 | Fe0.90B0.10 | Fe324B36_dir |
| S19 | Fe0.80B0.20 | Fe128B32_dir |
| S20 | Fe0.70B0.30 | Fe252B108_dir |
| S21 | Fe0.60B0.40 | Fe216B144_dir |
| S22 | Fe0.50B0.50 | Fe200B200_dir |
| S23 | Fe0.40B0.60 | Fe160B240_dir |
| S24 | Fe0.30B0.70 | Fe120B280_dir |
| S25 | Fe0.25B0.75 | Fe135B405_dir |
| S26 | Fe0.20B0.80 | Fe72B288_dir |
| S27 | Fe0.15B0.85 | Fe81B459_dir |
| S28 | Fe0.10B0.90 | Fe54B486_dir |

`data/` contains all 110 original `rdf.dat` files (400 radial bins each); trajectories are not needed for replotting. `reference/` contains the original PDFs and unmodified script.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_liquid_rdf_by_composition_combined.py
```

The script writes all 11 PDFs and PNGs to `output/`. Only paths and composition discovery were adapted; RDF values and plotting are unchanged. Columns are read by name, including the Fe-B/B-Fe naming variants.

Sources: `final/plot_liquid_rdf_by_composition_combined.py`, `final/{0GPa,10GPa}_liquid/*_dir/*_dir/rdf.dat`, and `final/liquid_rdf_by_composition_combined/`.
