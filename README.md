# Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis

Code and data for the paper:

> **Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis**
> Zihan Zhang, Mengyi Chen, Qianxiao Li, Peichen Zhong — arXiv:[2606.22885](https://arxiv.org/abs/2606.22885) (2026).


## [TBD]




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