# Figure Reproduction

This directory collects the source data and scripts for reproducing each figure in the main text and Supporting Information of **Interfacial melting as a thermodynamic indicator of solid-state synthesizability**.

Figure packages will be added individually after their source scripts and input data have been inspected.

## Organization

- `Fig1/`, `Fig2/`, etc.: main-text figures.
- `FigS1/`, `FigS2/`, etc.: Supporting Information figures.
- Panels belonging to the same figure are kept in the same figure folder, with separate panel subfolders when needed.

Each figure package will contain the plotting/analysis scripts, their required source data, and a local `README.md` describing the source locations, dependencies, execution commands, and expected outputs. Any path changes needed to run the copied scripts will be documented. Original source files will be preserved.

## Current status

- [Fig1](Fig1/README.md): a/b schematic notes, c 800 ps structure image, d MSD source data and scripts. Large MD trajectories are not included.
- [Fig2](Fig2/README.md): a/b free-energy and excess-curvature plots at 0/10 GPa, with input CSVs and plotting scripts.
- [Fig3](Fig3/README.md): a sixfold-B/energy overlay, b labeled structure and screenshot, c potential/PV comparison; includes source data and scripts.
- [Fig4](Fig4/README.md): concentration structure factors at 1500 K, with weighted OZ fits, bootstrap confidence bands, and source CSVs.
- [FigS1](FigS1/README.md): three upper-panel trajectory snapshots; trajectories omitted, lower panel is a conceptual diagram.
- [FigS2](FigS2/README.md): energy and force parity plots, with comparison CSVs and a shared plotting script.
- [FigS3](FigS3/README.md): 1500 K structure factors for 11 compositions at 0/10 GPa, with source data, weighted OZ fits, and bootstrap confidence bands.
- [FigS4](FigS4/README.md): 1800 K structure factors for 11 compositions at 0/10 GPa, with source data, weighted OZ fits, and bootstrap confidence bands.
- [FigS5](FigS5/README.md): 4320-atom structure factors at 1500 K and 0/10 GPa, with fit data/scripts, the original composite, and its labeled snapshot/structure.
- [FigS6](FigS6/README.md): KBI real-space structure factors at 1500 K for four compositions and 0/10 GPa, with downloaded source CSV and plotting script.
- [FigS7](FigS7/README.md): OMAT-PBE-head free-energy hulls at 1800 K and 0/5/10 GPa, with RES/phonon inputs and plotting script.
- [FigS8](FigS8/README.md): DFT r2SCAN pressure and temperature trends, with source RES/phonon data and plotting script.
- [FigS9](FigS9/README.md): fine-tuned r2SCAN-head free-energy hulls at 1800 K and 0/5/10 GPa, with RES/phonon inputs and plotting script.
- [FigS10](FigS10/README.md): 0 GPa free energies at 1500/1800 K, before and after configurational entropy; upper/lower scripts and source CSVs.
- [FigS11](FigS11/README.md): 10 GPa free energies at 1500/1800 K, before and after configurational entropy; panel-specific historical data versions retained.
- [FigS12-S16](FigS12-S16/README.md): shared energy-composition data and script for 1200/1500/1800/2000/2300 K, comparing 0/10 GPa potential energy, enthalpy, and vibrational free energy.
- [FigS17](FigS17/README.md): second derivatives of quintic free-energy fits, with pointwise OLS uncertainty and phase-composition guides.
- [FigS18-S28](FigS18-S28/README.md): shared RDF script and 110 source RDF files for 11 compositions, five temperatures, and 0/10 GPa.
- [FigS29](FigS29/README.md): sixfold-B fractions at 0/10 GPa and 0 GPa potential energy at 1500 K, with the historical input tables matching the original fit.
- [FigS30](FigS30/README.md): a 1800 K sixfold-B/potential-energy overlay with source data/script, and b-d structure screenshots with matching POSCARs; C labels sixfold-coordinated B.
- [FigS31-S33](FigS31-S33/README.md): shared sixfold-B/energy overlays at 1200/1500/1800 K, including error bars, matching historical energy data, 360 per-frame coordination files, and statistics-reproduction script.
- [FigS34](FigS34/README.md): 10 GPa endpoint-normalized PV terms, with plotting CSV/script and original liquid/elemental tables for recomputation.
- [FigS35](FigS35/README.md): OMAT enthalpies at 0/10 GPa and 1500 K, with both source datasets/scripts and original PDF/Illustrator file; one post-edited curve-color difference is documented.
- [FigS36](FigS36/README.md): OMAT structure screenshot and matching POSCAR at 0 GPa and 1500 K; C labels sixfold-coordinated B.
- [FigS37](FigS37/README.md): Fe-B/Cr-B free-energy comparison at 0 GPa and 1500/1800 K, including configurational entropy, source CSVs, plotting script, and Cr-B normalization inputs/script.
- [FigS38](FigS38/README.md): literature citation and DOI for Wang et al., Journal of Materiomics 1, 45-51 (2015); no new data or script.

Other figures will be added as their source locations are provided.
