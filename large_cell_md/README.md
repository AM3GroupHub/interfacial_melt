# Large-cell MD validation

This directory contains a 4320-atom Fe-B configuration and a LAMMPS input deck for the large-cell validation of the concentration-concentration structure factor at the FeB4 composition (`x_B = 0.80`). The revised manuscript compares the large-cell results at 1500 K and 0/10 GPa with the smaller-cell calculations to assess the pressure dependence of long-wavelength concentration fluctuations.

## Packaged files

| File | Description |
| --- | --- |
| `4320_supercell.data` | LAMMPS data file containing the simulation box, atomic masses, 4320 atom records with image flags, and 4320 velocity records. The header records LAMMPS version 10 Sep 2025 and timestep 500000. |
| `4320_supercell.vasp` | POSCAR-format structure exported with OVITO Basic 3.10.6, with Cartesian coordinates and a second Cartesian block containing velocities. Useful for inspecting the supercell; the supplied input deck reads the LAMMPS data format. |
| `in.lammps` | Example constant-temperature NPT run at 1500 K, using the fine-tuned Fe-B MACE potential. Paths must be updated before execution. |

The configuration contains 864 Fe atoms (LAMMPS type 1) and 3456 B atoms (type 2), giving Fe:B = 1:4. The orthorhombic box lengths in the supplied configuration are approximately `34.478795 x 23.797892 x 44.633790` angstrom. The LAMMPS data file uses nonzero lower box bounds; the POSCAR representation has the same cell lengths.

## Connection to the manuscript

- The large-cell comparison is discussed in the main text as Fig. S5 and appears in the SI figure labelled `fig:large_S`, using `large_Fe_B_scc_0.8_1500K_snapshot.pdf`.
- The comparison is at 1500 K, with 0 and 10 GPa results. It supports the conclusion that pressure suppresses concentration fluctuations; it is a finite-size check, not a claim of complete thermodynamic-limit convergence.
- The SI also describes Supplementary Movie 1 for a 4320-atom cell annealed at 1500 K and 0 GPa after a preceding 1800 K stage. This folder does not contain that movie or the preceding thermal-history inputs.

The supplied structure is an exported configuration, not a trajectory or a full LAMMPS restart. Its header alone does not establish its pressure, temperature, or complete preparation history. In particular, the input's original `read_data` filename must not be treated as provenance for the differently named packaged structure.

## Parameters in the supplied input

| Setting | Value |
| --- | --- |
| Units and boundaries | `metal`; periodic in all three directions |
| Potential | `pair_style mace no_domain_decomposition`; element order `Fe B` |
| Temperature | Constant 1500 K |
| Pressure | `iso 1 1 1.0`: constant 1 bar, approximately 0 GPa |
| Thermostat / barostat damping | 0.1 ps / 1.0 ps |
| Timestep | 0.001 ps = 1 fs |
| Run length | 500000 steps = 500 ps |
| Thermodynamic output | Every 100 steps = 0.1 ps: step, potential/kinetic/total energy, temperature, pressure, volume |
| Trajectory output | Every 100 steps = 0.1 ps: `id type xu yu zu vx vy vz` in `SUPERCELL.dump` |
| Final configuration | `final.data`, written after the run |

The example applies one 500 ps NPT stage; it does not implement the complete 2600 K melt and staged cooling schedule described for the composition-wide dataset. The `velocity ... create 1800.0 ...` line is commented out, so the run uses velocities read from the supplied data file. Reading a data file does not restore thermostat/barostat state or the original timestep counter.

## Preparation and execution

1. Use a LAMMPS executable built with the MACE pair style and its compatible runtime dependencies. Match CPU/GPU and Kokkos launch options to that build; a generic LAMMPS installation may not provide this pair style.
2. Reuse the fine-tuned Fe-B model and training instructions in [`../MLIP/README.md`](../MLIP/README.md). The large potential is managed in `MLIP/` and is not duplicated in this directory; the source model is [`../MLIP/mace-mh-1_finetuned_FeB_mixed.model`](../MLIP/mace-mh-1_finetuned_FeB_mixed.model). For execution, prepare its LAMMPS-compatible `*.model-lammps.pt` export using the procedure appropriate to the installed MACE/LAMMPS interface, and record the software versions used. A training checkpoint or a generic compiled model is not automatically interchangeable with this export.
3. Work in a separate run directory for each pressure. Copy the input and the chosen starting data file there, then update these lines:

```lammps
read_data       4320_supercell.data
pair_coeff      * * /path/to/mace-mh-1_finetuned_FeB_mixed.model-lammps.pt Fe B
```

The original input refers to `final-0GPa-1500.data`, which is not included, and to an absolute model path on the original HPC system. Using `4320_supercell.data` explicitly starts a new run from the packaged configuration; reproducing the original trajectory also requires its original initialization and thermal history.

4. Set the desired pressure. With `units metal`, pressure is in bar:

```lammps
# Nominal 0 GPa, matching the supplied input (1 bar):
fix npt_run all npt temp 1500.0 1500.0 0.1 iso 1 1 1.0

# 10 GPa (100000 bar), in a separate run:
fix npt_run all npt temp 1500.0 1500.0 0.1 iso 100000 100000 1.0
```

Use only one of these `fix` lines in each input. Exact zero pressure can instead be requested with `iso 0 0 1.0`. For a pressure comparison, prepare and equilibrate the starting configuration at each target pressure and preserve that preparation history.

5. Run from the prepared directory, replacing `lmp` and adding any build-specific GPU/Kokkos options as needed:

```bash
lmp -in in.lammps -log log.lammps
```

This generates `log.lammps`, `SUPERCELL.dump`, and `final.data`. These outputs are not currently included in this folder.

## Structure-factor analysis and remaining reproducibility files

The revised large-cell figure reports an Ornstein-Zernike fit, `S_cc(k) = a / (1 + b k^2)`, constrained to `a, b >= 0`, with `a = S_cc(0)`. Its caption specifies `0 < |k| < 1.2 angstrom^-1`. The response describes unweighted large-cell fits, 5000 residual-bootstrap replicates, and pointwise 95% confidence bands. The ideal-random-mixing reference for this composition is `x_Fe * x_B = 0.16`.

The SI methods text also mentions a `0.8 angstrom^-1` cutoff, while the response's finite-size discussion refers to a common fitting window and weighting convention. These descriptions should be reconciled with the actual analysis script before figure reproduction; no fitting window or weighting choice is encoded in `in.lammps`.

This folder currently provides the configuration and MD input only. Reproducing the submitted structure-factor figure additionally requires the original trajectories or sampled frames for both pressures, their equilibration/sampling selections, the reciprocal-space analysis and fitting scripts, and the numerical structure-factor/fit results. Those files are not included here. The general melt-quench workflow is documented in [`../Melt_quench/README.md`](../Melt_quench/README.md), and model training/testing is documented in [`../MLIP/README.md`](../MLIP/README.md).
