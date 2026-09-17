# Figure S4

Concentration structure factors at 1800 K: 11 compositions at 0 and 10 GPa, using the last 50% of each trajectory. The original PDF is `reference.pdf`.

- `plot_github_scc_last0p5_3x4.py`: plotting and weighted Ornstein-Zernike fitting.
- `data/last0p5_records_flat.csv`: structure-factor values, wavevectors, weights, and frame counts.
- `data/last0p5_weighted_oz_bootstrap_summary.csv`: original 44-state fit summary for reference.

Run from this folder with Python 3.9:

```bash
python -m pip install -r requirements.txt
python plot_github_scc_last0p5_3x4.py
```

Outputs are the 1800 K PDF/PNG and a refitted summary in `output/`. The script retains all 44 fits at 1500/1800 K to preserve the original seed ordering: `k < 0.8 angstrom^-1`, 500 residual-bootstrap replicates, seed `20260824`. Shading shows pointwise 95% fit confidence bands, not frame-level errors. No NPZ files or trajectories are needed.

Source: `response/S_test/plot_github_scc_last0p5_3x4.py` and `response/S_test/github_scc_last0p5_pressure_compare/`. Only paths, CSV-loading setup, and output-temperature selection were adapted; fitting and plotting methods are unchanged.
