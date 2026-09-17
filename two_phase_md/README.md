# Two-phase Interfacial MD

This directory contains starting structures and an example LAMMPS input for the FeB/alpha-boron two-phase simulations used to examine interfacial melting. The manuscript constructs the interface by joining the FeB (010) surface to the alpha-boron (100) surface and relaxing the combined structure with the fine-tuned MACE-MH-1 potential.

## Connection to the manuscript

- Main text, **Two-phase MD of solid-solid interfaces** and Fig. 1(c,d): the interfacial configuration after 800 ps at 1500 K and nominal 0 GPa, and the corresponding interface/bulk mean square displacement (MSD) comparison.
- Supporting Information, MD workflow and the upper panel of the workflow figure (`fig:workflow`): the initial structure and configurations at 1500 K / 0 GPa, 1500 K / 10 GPa, and 1800 K / 10 GPa.
- Response to Reviewer 1, interface figure (`fig:interface`): the same three state points support the observation that mobility is enhanced near the interface while the bulk regions remain crystalline.

The supplied input covers only 1500 K at approximately ambient pressure. These simulations start from a two-phase structure; they complement the homogeneous melt-quench calculations by probing interfacial melting, rather than directly measuring the onset of demixing in a homogeneous melt.

## Packaged files

| File | Description |
| --- | --- |
| `optimized_tot.vasp` | Two-phase starting configuration in VASP/POSCAR format, with Cartesian coordinates and element order Fe, B |
| `SUPERCELL.lmp` | Starting configuration in LAMMPS `atomic` data format, read directly by the example input; atom type 1 is Fe and type 2 is B |
| `in.lammps-example` | Constant-temperature NPT example using the LAMMPS MACE pair style |

Both structure files contain **1440 atoms: 288 Fe and 1152 B**, giving the overall composition FeB4 (`x_B = 0.8`). This overall composition describes the combined FeB and boron slabs, not a homogeneous FeB4 crystal. The orthorhombic starting box is approximately **12.221925 x 17.021231 x 52.240925 angstrom**, with periodic boundaries in all three directions. The VASP coordinates include periodic images outside the displayed cell; after wrapping, corresponding coordinates in the two files agree to within 0.00003 angstrom.

## Parameters in the supplied input

| Setting | Value |
| --- | --- |
| Unit system | `metal` (time in ps, length in angstrom, pressure in bar) |
| Atom masses | Fe: 55.845; B: 10.810 |
| Potential | `pair_style mace no_domain_decomposition`; element mapping `Fe B` |
| Initial velocities | Gaussian distribution at 1500 K, seed 8257, zero net momentum |
| Ensemble | NPT, isotropic pressure coupling |
| Temperature target / damping time | 1500 K / 0.1 ps |
| Pressure target / damping time | 1 bar / 1.0 ps |
| Time step | 0.001 ps (1 fs) |
| Run length | 1,000,000 steps = 1000 ps |
| Thermodynamic output | Every 100 steps (0.1 ps): step, potential/kinetic/total energy, temperature, pressure, volume |
| Trajectory output | Every 100 steps (0.1 ps): atom ID, type, unwrapped coordinates `xu yu zu`, velocities `vx vy vz` |

Two details matter when relating the example to the paper:

- The `iso 1 1 1.0` setting specifies **1 bar = 0.0001 GPa**, which is approximately the nominal 0 GPa condition in the manuscript. For an exactly zero-pressure target, use `iso 0 0 1.0` in a working copy and record that change.
- The input runs for **1000 ps**; the **800 ps** configuration discussed in Fig. 1(c) corresponds to **step 800000** of this trajectory. The figure time is not the full run length in the example.

## Dependencies and running the example

Required external components are a LAMMPS executable with the `mace` pair style and support for the supplied `no_domain_decomposition` option, its compatible runtime dependencies, and a LAMMPS-exported version of the fine-tuned Fe-B MACE potential. A single-process command is shown below; use the launch configuration supported by your installed MACE/LAMMPS build.

The `pair_coeff` line currently points to an author-specific absolute path ending in:

```text
mace-mh-1_finetuned_FeB_mixed.model-lammps.pt
```

The potential files are maintained centrally in `../MLIP/` and are not duplicated here because of their size. See the [MLIP workflow](../MLIP/README.md) for the training configuration, training inputs, and launch instructions. The underlying fine-tuned model is available at `../MLIP/mace-mh-1_finetuned_FeB_mixed.model` and in `../MLIP/train/`. To run this example, prepare a LAMMPS-compatible export from that model using a MACE version compatible with your LAMMPS build, then replace the absolute path in a working copy of the input. The specific `.model-lammps.pt` export referenced above is not archived here or in `../MLIP/`; merely renaming the `.model` file does not establish compatibility. Exact MACE/LAMMPS versions and the export command are not recorded in this folder.

From this directory, prepare a working copy:

```bash
cp in.lammps-example in.lammps
```

Edit `pair_coeff` in `in.lammps` to point to the available LAMMPS-exported model, keeping the Fe/B mapping:

```text
pair_coeff * * /path/to/mace-mh-1_finetuned_FeB_mixed.model-lammps.pt Fe B
```

Then run with your LAMMPS executable:

```bash
lmp -in in.lammps -log log.lammps
```

The expected generated outputs are `log.lammps`, `SUPERCELL.dump`, and `final.data`; they are not supplied in this directory. `final.data` is written after the full 1000 ps run. To select the configuration corresponding to Fig. 1(c), use the trajectory frame whose `ITEM: TIMESTEP` is `800000`.

For the additional SI state points, prepare separate working directories and input copies. A 10 GPa target corresponds to `iso 100000 100000 1.0` in `metal` units. For 1800 K, change both the initial velocity temperature and the two `fix npt` temperature targets to 1800 K. These are adaptations of the example; independently archived input files and run outputs for those state points are not present here.

## Reproducibility scope

This folder packages the starting configuration and an MD input example, with potential training and model files organized separately under `../MLIP/`. It does not contain the interface construction/relaxation scripts, production trajectories or logs, final snapshots, interface/bulk atom selections, MSD-analysis scripts, or figure-rendering scripts used for the manuscript. Full reproduction of the published snapshots and MSD curves also requires those analysis details and the software versions. The dump records unwrapped coordinates and atom IDs for subsequent displacement analysis, but the manuscript's region selections cannot be recovered from this input alone.
