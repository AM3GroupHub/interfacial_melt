# Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis

Code and data for the paper:

> **Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis**
> Zihan Zhang, Mengyi Chen, Qianxiao Li, Peichen Zhong — arXiv:[2606.22885](https://arxiv.org/abs/2606.22885) (2026).

## Subfolders for MD simulations

- `Melt_quench/`: melt-quench molecular dynamics workflow, post-processing scripts, and figure packages
- `MLIP/`: machine-learned interatomic potential training, testing, and comparison plotting
- `Structure_Search/`: AIRSS-based structure search workflow with MACE relaxation
- `snapshot/`: VESTA snapshot files used for structure visualization in the manuscript and SI
- `two_phase_md/`: FeB(010)/alpha-boron(100) interface structures and an example NPT input for 1440 atoms (288 Fe + 1152 B). Supports the interfacial-melting calculation in main-text Fig. 1(c,d) and the upper panel of SI Fig. S1. |
- `large_cell_md/`:  A 4320-atom FeB4-composition structure (864 Fe + 3456 B) and an example NPT input. Supports the finite-size validation of the structure-factor analysis at 1500 K, comparing nominal 0 and 10 GPa in SI Fig. S5. |

See the `README.md` inside each subfolder for details.




## The structure-factor descriptor

We diagnose melt stability from the **concentration–concentration structure factor** `S_cc(k)`
of the MD melt, see its README for the method and figures.


### Reproducing the figures

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


## License

**[TBD]** 
