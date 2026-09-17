# Figure S17

Second derivatives of fifth-order free-energy fits at 0/10 GPa and 1500/1800 K, with phase-composition guides. The original PDF is `reference.pdf`.

- `data/`: the two normalized-energy source CSVs.
- `plot_curvature_fit_5th_order.py`: fits 11 compositions per state after adding ideal configurational entropy.
- `reference/`: original script and exported `curvature_with_uncertainty.csv` for numerical checks.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_curvature_fit_5th_order.py
```

The PDF, PNG, and a 3204-row curvature CSV are written to `output/`. Shading is pointwise +/-1 sigma from ordinary least-squares residual variance propagated to the second derivative, not a 95% interval. The script uses a centered fifth-order basis and retains the original phase guides. Only input/output paths were adapted.

Sources: `final/curvature_fit_5th_order/` and `final/{0,10}GPa_liquid/normalized_energy_{0,10}GPa_liquid.csv`. No MD trajectories are required to replot.
