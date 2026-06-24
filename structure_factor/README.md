# structure_factor — concentration–concentration structure factor of Fe–B melts

Our contribution to *Interfacial-melt stability as a thermodynamic prerequisite for
solid-state synthesis* (arXiv:2606.22885). Computes the Bhatia–Thornton structure factor
`S_cc(k)`, the thermodynamic factor `Γ = x_Fe x_B / S_cc(0)`, and the excess free-energy
curvature `(1/k_BT) ∂²g_ex/∂x² = (Γ−1)/(x_Fe x_B)` for liquid Fe–B, along two routes
(reciprocal-space S(k) and real-space Krüger–Vlugt KBI).

## What this folder reproduces

| paper figure | script | output |
|---|---|---|
| **Fig. 3** `S_cc(k)`, 1500 K | `plot_fig3_scc_k.py` | `figures/fig3_scc_k_T1500.pdf` |
| **Fig. 1(c,d) bottom** `(1/k_BT)∂²g_ex/∂x²` | `plot_fig1_gex_curvature.py` | `figures/fig1_gex_curvature.pdf` |
| **SM** per-composition `S_cc(k)` | `plot_sm_state_diagnostics.py` | `figures/state_diagnostics/sm_state_diagnostics_T{T}_P{P}.pdf` |
| **SM** KBI cross-check 2×2 | `plot_sm_kbi_2x2.py` | `figures/sm_kbi_scc_2x2_T1500.pdf` |

It does **not** produce the Fig. 1 top energy panels, Fig. 2 (CN=6 / PV), the MD generation,
or the MACE potential — those are other parts of the project.

## Pipeline (three stages)

```
1_compute_scc_reciprocal.py   dumps  → results/ caches   # needs data_root (raw dumps; internal)
2_compute_scc_kbi.py
3_export_csv.py               caches → data/csv/*.csv     # the published plotted data
plot_*.py                     data/csv → figures/         # portable: CSV + matplotlib only
```

This repo ships **code + data (the CSV tables)**; rendered figures are not included — they
regenerate from `data/csv/`. The plot scripts read **only** the CSVs (no `results/` caches,
no dumps, no `src/` — just numpy + matplotlib):

    python plot_fig3_scc_k.py
    python plot_fig1_gex_curvature.py
    python plot_sm_state_diagnostics.py
    python plot_sm_kbi_2x2.py

The raw MD dumps are not distributed. To rebuild the CSVs from the cached structure factors,
run `python 3_export_csv.py` (reads `results/`); recomputing the caches from scratch needs the
dumps (`1_`/`2_`, with `data_root` set in `config.yaml`). All steps use the last half of each
trajectory (`last_frac: 0.5`).

## Method (see SM §I.A / §I.B)

- **Reciprocal route:** partial `S_αβ(k)` from MD positions at commensurate shells
  `k = 2πn/L` (scaled coords, NPT-safe), projected to `S_cc(k)`, extrapolated to `k=0` with
  a weighted OZ Lorentzian `S_cc(k)=a/(1+b k²)`, `|k|<0.8 Å⁻¹` (M2). `S_cc(0)=a`.
- **KBI route:** windowed sinc-Fourier transform of the partial RDFs (Krüger–Vlugt window);
  at `k=0` it is exactly the static Kirkwood–Buff integral (no extrapolation).
