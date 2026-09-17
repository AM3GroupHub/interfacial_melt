# Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis

Code and data for the paper:

> **Interfacial-melt stability as a thermodynamic prerequisite for solid-state synthesis**
> Zihan Zhang, Mengyi Chen, Qianxiao Li, Peichen Zhong — arXiv:[2606.22885](https://arxiv.org/abs/2606.22885) (2026).

## Data and scripts

- [figure_reproduction/](figure_reproduction/README.md): source data, plotting scripts, and structure files for main-text Figures 1-5 and Supporting Figures S1-S38. Figure S38 provides a literature citation.
- [Cr-B-input/](Cr-B-input/README.md): Cr-B MD structures and input-generation scripts, following the Fe-B melt-quench workflow.
- [Melt_quench/](Melt_quench/README.md): Fe-B melt-quench MD and post-processing.
- [MLIP/](MLIP/README.md): potential-model training, testing, and usage.
- [Structure_Search/](Structure_Search/README.md): AIRSS structure search and MACE relaxation.
- [snapshot/](snapshot/README.md): structural visualizations.
- [two_phase_md/](two_phase_md/README.md): two-phase interface MD inputs.
- [large_cell_md/](large_cell_md/README.md): 4320-atom MD inputs for finite-size checks.
- [structure_factor/](structure_factor/README.md): concentration-concentration structure-factor analysis, `S_cc(k)`.

## Reproducing figures

Start with [figure_reproduction/README.md](figure_reproduction/README.md), then follow the short README in each figure folder for dependencies and commands. Related figures share a folder where appropriate. Large MD trajectories and potential files are not duplicated in the figure packages.

For MD calculations, configure the model path, LAMMPS executable, and run settings as described in the corresponding input-folder README.

## License

**[TBD]** 
