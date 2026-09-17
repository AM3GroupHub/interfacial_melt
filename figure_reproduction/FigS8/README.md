# Figure S8

DFT r2SCAN formation enthalpies versus pressure and 0 GPa free energies versus temperature, relative to FeB + B. The original PDF is `reference.pdf`.

- `data/`: 12 source RES files at nominal 0/5/10 GPa and four 0 GPa phonon `thermal_properties.yaml` files.
- `plot_combined_dft_trends.py`: retains the original six pressure values and static enthalpies; reads the YAML files for 0-2000 K curves.
- `reference/`: original script and temperature-curve CSV for verification.

Run from this folder with Python 3.12:

```bash
python -m pip install -r requirements.txt
python plot_combined_dft_trends.py
```

PDF/PNG are written to `output/`. Only local paths were adapted. The original PDF was subsequently cropped and given (a)/(b) labels; the script retains its original title/layout. Numerical curves are unchanged.

Source: `response/DFT_MACE/DFT/`. Pressure values are tabulated from the RES headers as enthalpies (do not add PV again); the nominal 5 GPa headers report 4.9 GPa. Original filenames, including `FeB2.res.res`, are preserved.
