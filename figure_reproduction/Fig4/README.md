# Figure 4

Cr-B liquid free energies at 0 GPa for 1500 K (left) and 1800 K (right). The horizontal coordinate is the B atomic fraction in `Cr_(1-x)B_x`. The dashed line marks the CrB4 composition (`x_B = 0.8`).

## Method summary

For each of 11 liquid compositions, the script normalizes the calculated free energy against elemental Cr and B endpoints. It then adds the ideal configurational contribution

`k_B T [x ln(x) + (1-x) ln(1-x)]`

and fits all compositions with a fifth-order polynomial. Blue and red curve segments show the convex and concave regions separated at the fitted inflection point nearest `x_B = 0.6`. No error bars are plotted.

The page width, panel width and height, and inter-panel spacing are fixed explicitly so this main-text figure matches the corresponding Fe-B panel geometry. The PDF page is 239.494 x 135.341 pt; each plotting axis is 89.951 x 99.347 pt, with an 11.920 pt horizontal gap.

## Contents

- `plot_free_energy_0gpa_1500k_1800k_with_config_entropy.py`: self-contained normalization and plotting script.
- `data/source_crb/`: 11 alloy `_G.txt` tables plus the elemental Cr and B endpoint tables.
- `data/normalized_energy_CrB_0GPa_1500_1800.csv`: preserved 22-row normalized table for traceability.
- `output/`: regenerated normalized CSV and publication-ready PDF/PNG.
- `reference.pdf`: the accepted main-text rendering.
- `reference/plot_free_energy_0gpa_1500k_1800k_with_config_entropy.original.py`: unadapted source script from the working directory.

## Reproduce

Use Python 3.9 or newer from this folder:

```bash
python -m pip install -r requirements.txt
python plot_free_energy_0gpa_1500k_1800k_with_config_entropy.py
```

The script reads only local files and writes the normalized CSV, PDF, and PNG to `output/`. A successful run reports 22 records.

Only data/output paths and package export behavior were adapted. The thermodynamic normalization, configurational-entropy expression, fit, colors, annotations, and final fixed layout are unchanged from the accepted source figure.

Source: `final/Cr-B-2/`. The curated raw tables are identical to those retained for `FigS37/data/source_crb/`.
