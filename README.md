# Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis

Code and data for the paper:

> **Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis**
> Zihan Zhang, Mengyi Chen, Qianxiao Li, Peichen Zhong — arXiv:[2606.22885](https://arxiv.org/abs/2606.22885) (2026).

## Overview

Thermodynamic stability on the formation-energy convex hull is the usual screen for
synthesizability, yet many predicted-stable phases resist synthesis. This work proposes an
additional thermodynamic condition for solid-state, interfacial-melt-mediated routes: the
**interfacial melt at the target composition must itself remain locally stable against
spinodal decomposition**. We demonstrate this for the Fe–B system — where thermodynamically
stable FeB₄ is reported under high-pressure but not low-pressure synthesis — using melt-quench
molecular dynamics driven by a fine-tuned machine-learning interatomic potential. At ambient
pressure the B-rich melt near FeB₄ develops a concave free-energy landscape (a demixing
instability); applied pressure restores stability, matching the experimental synthesis boundary.

## The structure-factor descriptor

A homogeneous binary melt is locally stable only where its per-atom free energy `g` is convex
in composition. That curvature is encoded in the **concentration–concentration structure
factor** `S_cc(k)` — the spectrum of composition fluctuations, read directly from the MD melt.
Its `k → 0` limit gives the thermodynamic factor

```
Γ = x_Fe x_B / S_cc(0),        ∂²g/∂x²|_{T,P} = k_B T / S_cc(0)
```

so `Γ = 1` is ideal mixing, `Γ > 1` ordering, `0 < Γ < 1` demixing, and `Γ → 0`
(`S_cc(0) → ∞`) the spinodal. A negative excess curvature
`(1/k_BT) ∂²g_ex/∂x² = (Γ − 1)/(x_Fe x_B)` therefore flags a demixing tendency — the
instability that renders the ambient-pressure FeB₄ interfacial melt unsynthesizable. This
descriptor is computed in [`structure_factor/`](structure_factor/).

## Repository structure

```
interfacial_melt/
├── structure_factor/                 # S_cc(k), Γ, excess curvature — Fig. 3, Fig. 1(c,d) bottom, SM
│   ├── 1_compute_scc_reciprocal.py   #   reciprocal-space S(k)  (dumps → cache)
│   ├── 2_compute_scc_kbi.py          #   real-space KBI route   (dumps → cache)
│   ├── 3_export_csv.py               #   caches → data/csv/  (published plotted data)
│   ├── plot_fig3_scc_k.py            #   the four plot scripts (CSV → figure)
│   ├── plot_fig1_gex_curvature.py
│   ├── plot_sm_state_diagnostics.py
│   ├── plot_sm_kbi_2x2.py
│   ├── src/                          #   shared library (vendored numerics core)
│   ├── data/{state_index.json, csv/} #   state list + per-figure CSV data tables
│   ├── results/                      #   cached structure factors (npz)
│   ├── config.yaml   README.md
│
├── [TBD]/                            # [TBD — Zihan Zhang] energy landscape H−T(S_vib+S_mix),
│                                     #   convex/concave classification — Fig. 1(c,d) top
├── [TBD]/                            # [TBD — Zihan Zhang] CN=6 B fraction + PV analysis — Fig. 2
├── [TBD]/                            # [TBD — Zihan Zhang] melt-quench MD + MACE-MH-1 potential
│
└── README.md
```

### Components

| component | folder | author | paper |
|---|---|---|---|
| Structure-factor melt stability (`S_cc`, `Γ`, curvature) | [`structure_factor/`](structure_factor/) | Mengyi Chen, Peichen Zhong | **Fig. 3**, **Fig. 1(c,d) bottom**, SM |
| Energy-landscape free-energy curvature | `[TBD]` | Zihan Zhang | Fig. 1(c,d) top |
| Six-coordinated-B (CN=6) fraction + `PV` analysis | `[TBD]` | Zihan Zhang | Fig. 2 |
| Melt-quench MD workflow + MACE-MH-1 potential | `[TBD]` | Zihan Zhang | Methods / SM |

> Each subfolder is self-contained with its own README.

## Requirements

Python 3.11+ with `numpy`, `scipy`, `matplotlib`, `pyyaml`. Per-component requirements are in
each subfolder's README (e.g. `structure_factor/`'s plot scripts need only `numpy` + `matplotlib`).

## Reproducing the figures

For the structure-factor figures (Fig. 3, Fig. 1 c,d bottom, SM) — these read the shipped
`data/csv/` tables, no raw MD trajectories needed:

```bash
cd structure_factor
python plot_fig3_scc_k.py
python plot_fig1_gex_curvature.py
python plot_sm_state_diagnostics.py
python plot_sm_kbi_2x2.py
```

For the other components, see their subfolder READMEs (`[TBD]`).

## Citation

```bibtex
@article{Zhang2026InterfacialMelt,
  title   = {Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis},
  author  = {Zhang, Zihan and Chen, Mengyi and Li, Qianxiao and Zhong, Peichen},
  journal = {arXiv preprint arXiv:2606.22885},
  year    = {2026},
}
```
<!-- TODO: update to the journal reference / DOI once published. -->

## Authors

- Zihan Zhang, Mengyi Chen (equal contribution)
- Qianxiao Li
- Peichen Zhong (corresponding author — zhongpc@nus.edu.sg)

National University of Singapore.

## License

**[TBD]** — to be set by the project PI.
