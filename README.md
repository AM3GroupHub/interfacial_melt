# Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis

Code and data for the paper:

> **Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis**
> Zihan Zhang, Mengyi Chen, Qianxiao Li, Peichen Zhong — arXiv:[2606.22885](https://arxiv.org/abs/2606.22885) (2026).

## Subfolders for MD simulations

- `Melt_quench/`: melt-quench molecular dynamics workflow, post-processing scripts, and figure packages
- `MLIP/`: machine-learned interatomic potential training, testing, and comparison plotting
- `Structure_Search/`: AIRSS-based structure search workflow with MACE relaxation
- `snapshot/`: VESTA snapshot files used for structure visualization in the manuscript and SI

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
