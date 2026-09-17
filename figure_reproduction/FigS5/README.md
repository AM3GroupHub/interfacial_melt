# Figure S5

Large-cell Fe864B3456 structure factors at 1500 K and 0/10 GPa. The submitted composite is `reference/large_Fe_B_scc_0.8_1500K_snapshot.pdf`.

- `data/`: original reciprocal-space CSVs for both pressures.
- `analyze_fit_uncertainty.py` and `plot_s_test_structure_factors.py`: fitting and curve-plotting scripts.
- `reference/`: original composite, curve-only PDF, fit summary, and bootstrap samples.
- `snapshot/`: original PNG and POSCAR-format `.vasp`; Fe/B/C counts are 864/1881/1575. **C labels six-coordinated B for visualization, not carbon** (B-B cutoff 2.37 angstrom). The optional relabeling script works on an original Fe/B LAMMPS dump, not the POSCAR.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python analyze_fit_uncertainty.py
```

Curve-only PDF/PNG, refitted summary, and bootstrap samples go to `output/`. Fits are unweighted, constrained OZ fits with `k <= 1.2 angstrom^-1`, 5000 paired then 5000 residual-bootstrap samples per pressure, seed `20260817`, in 0 then 10 GPa order. Shading is the pointwise 95% residual-bootstrap fit interval.

Only paths, default input selection, and output selection were adapted. Refitted values can differ in the last digits across library versions; the original results are retained in `reference/`. The snapshot inset was assembled separately; the original composite and its source PNG/POSCAR are retained. Large trajectories are omitted; see [MD inputs](../../large_cell_md/README.md).

Sources: `response/S_test/` scripts and CSVs, `analysis_plots_pressure_compare/` PDFs/PNG, and `gif/long/` POSCAR/relabeling script. The initially supplied `scc_last0p5_T1800_pressure_compare_3x4.pdf` path does not exist in this folder; the composite above is the S5 file referenced by `submit/final-SI.tex`.
