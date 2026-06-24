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
pressure the B-rich interfacial melt near FeB₄ develops a concave free-energy landscape (a
demixing instability), corroborated by the concentration–concentration structure factor;
applied pressure restores melt stability, consistent with the experimental synthesis boundary.

## Repository structure

| folder | contents | paper |
|---|---|---|
| [`structure_factor/`](structure_factor/) | Concentration–concentration structure factor `S_cc(k)`, thermodynamic factor `Γ`, and excess free-energy curvature of the Fe–B melt (reciprocal-space + Kirkwood–Buff routes), with the plotted data as CSV. | **Fig. 3**, **Fig. 1(c,d) bottom row**, SM structure-factor figures |
| _(to be added)_ | Free-energy landscape from melt energetics `H − T(S_vib + S_mix)` and its convex/concave classification. | Fig. 1(c,d) top row |
| _(to be added)_ | Six-coordinated-B (CN=6) fraction and the `PV` analysis of the B-rich melt. | Fig. 2 |
| _(to be added)_ | Melt-quench MD workflow and the fine-tuned MACE-MH-1 potential. | Methods / SM |

> Each subfolder is self-contained and has its own README. `structure_factor/` is contributed
> by Mengyi Chen & Peichen Zhong; the remaining components are added by the respective authors.

## Requirements

Python 3.11+ with `numpy`, `scipy`, `matplotlib`, `pyyaml`. Per-component requirements are in
each subfolder's README (e.g. `structure_factor/` plot scripts need only `numpy` + `matplotlib`).

## Reproducing the figures

See each subfolder. For the structure-factor figures (Fig. 3, Fig. 1 c,d bottom, SM):

```bash
cd structure_factor
python plot_fig3_scc_k.py
python plot_fig1_gex_curvature.py
python plot_sm_state_diagnostics.py
python plot_sm_kbi_2x2.py
```

(These read the shipped `data/csv/` tables — no raw MD trajectories needed.)

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

**TBD** — to be set by the project PI.
