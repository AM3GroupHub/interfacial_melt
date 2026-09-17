# Figure S6

Kirkwood-Buff real-space structure factors at 1500 K for `x_B = 0.6, 0.7, 0.8, 0.9`, comparing 0 and 10 GPa.

- `data/csv/sm_kbi_scc_2x2.csv`: downloaded curves and zero-wavevector values; the only plotting input.
- `plot_sm_kbi_2x2.py`: downloaded plotting script, with only the output directory changed.
- `reference/plot_sm_kbi_2x2.py`: unmodified upstream script.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_sm_kbi_2x2.py
```

PDF and PNG are written to `output/`. Solid/dashed curves show 0/10 GPa; red x/+ markers show the tabulated `S_cc(0)`. No refitting or MD trajectories are required.

Source: [AM3GroupHub/interfacial_melt, structure_factor](https://github.com/AM3GroupHub/interfacial_melt/tree/276b1655406150db0a026500f2f2eb20ae3f79da/structure_factor), commit `276b1655406150db0a026500f2f2eb20ae3f79da`, downloaded 2026-09-17. The script/export workflow specifies 1500 K; the CSV has no temperature column.
