# Figure 1

- **a, b:** schematics; no numerical source data.
- **c:** `panel_c/FeB-0-1500.png`, the 800 ps frame at 1500 K and nominal 0 GPa. The trajectory is omitted because of its size; rerun the inputs in [two_phase_md](../../two_phase_md/README.md) to reproduce the MD calculation.
- **d:** `panel_d/` contains the plotting script, six layer-MSD CSV files, atom selections, MSD-calculation scripts, and the original plot (`reference.pdf`). The plot uses five groups; `FeB_interface_B` is retained as additional data. CSVs cover 0-844.2 ps at 0.1 ps intervals; each group contains 48 atoms.

## Replot d

From `panel_d/`, using Python 3.9 or newer:

```bash
python -m pip install -r requirements.txt
python plot_layer_mean_msd_compact.py
```

Outputs are saved under `panel_d/output/`. Local paths replace the original machine-specific defaults; plotting and MSD formulas are unchanged. Arial is used when available, otherwise DejaVu Sans.

To recalculate MSD from a newly simulated trajectory (no extra Python packages required):

```bash
python calc_msd_atom_layers.py --dump /path/to/SUPERCELL.dump --max-frames 8443
```

This uses the archived two-phase starting structure and writes `regenerated_data/`. Atom groups are fixed from the initial layers; MSD is relative to the first frame with whole-system center-of-mass drift removed.

Sources: c from `picture-8.27/FIG1/gif/FeB-0-1500.png` in the FeB project; d scripts from `D:/group_meeting/gif/`, data and reference PDF from its `1500_0/layer_msd/` folder.
