# Figure 5

Concentration structure factors at 1500 K for `x_B = 0.6, 0.7, 0.8, 0.9`, comparing 0 and 10 GPa. The original PDF is `reference.pdf`; historical `fig3` filenames are retained.

- `plot_fig3_scc_k_with_fit_uncertainty.py`: main plotting script.
- `plot_github_scc_last0p5_3x4.py`: weighted Ornstein-Zernike fitting and bootstrap routines.
- `data/fig3_scc_k.csv`: plotted structure-factor points.
- `data/last0p5_records_flat.csv`: structure-factor values and weights for fitting, from the last 50% of each trajectory.
- `data/last0p5_weighted_oz_bootstrap_summary.csv`: original 44-state fit summary for reference.

Run in a Python 3.9 environment from this folder:

```bash
python -m pip install -r requirements.txt
python plot_fig3_scc_k_with_fit_uncertainty.py
```

The script refits eight displayed states and writes PDF/PNG files to `output/`. It uses `k < 0.8 angstrom^-1`, nonnegative OZ parameters, weights equal to the square root of shell multiplicity, and 500 residual-bootstrap replicates for pointwise 95% fit confidence bands. Seed `20260824` and the original 44-state seed ordering are preserved.

Only paths and CSV-loading setup were adapted; fitting and plotting methods are unchanged. The packaged CSVs are sufficient; no NPZ files or MD trajectories are required.

Sources: `response/S_test/` scripts, `github_interfacial_melt/structure_factor/data/csv/fig3_scc_k.csv`, and `github_scc_last0p5_pressure_compare/` data/reference PDF.
