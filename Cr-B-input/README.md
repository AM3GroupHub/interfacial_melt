# Cr-B MD Inputs

Cr-B counterpart of [Melt_quench/MD](../Melt_quench/MD/), using the same structure-preparation and NVT/NPT workflow with Cr-B structures and a Cr-compatible MACE model. The resulting Cr-B data are used in [Figure S37](../figure_reproduction/FigS37/README.md).

## Contents

- `1_9/` through `9_1/` cover B fractions from 0.10 to 0.90; folder ratios are **B:Cr**. `7.5/` and `8.5/` correspond to B fractions 0.75 and 0.85. `B/` and `Cr/` provide elemental references.
- Each folder contains CIF/POSCAR structures, `cif2vasp.py`, `gen_sc.py`, `gen_input.py`, `cal_MQ.py`, and a PBS job template.
- `cal_MQ.py` calls `gen_sc.read_SC()` to read the supplied CIF cell, then generates and runs the selected MD stages. Outputs are placed under `<structure>/NVT/` and `<structure>/NPT/`.

## Running

Use a Linux/HPC environment with Python, NumPy, ASE, MPI, and MACE-enabled LAMMPS. Pymatgen is needed only to regenerate the supplied `.vasp` files with `cif2vasp.py`.

Before running, update `p_path`, `run_lmp`, and the temperature/step lists in `cal_MQ.py`. The supplied model path names `mace-mh-1-matpes_r2scan.model-lammps.pt`; the model is not included here. Adapt the environment and resource settings in `job_mace.pbs` if using PBS.

From an individual composition folder, after configuration:

```bash
python cal_MQ.py
```

Alternatively, submit `qsub job_mace.pbs`. The timestep is 1 fs; each configured NPT stage runs 150,000 steps (150 ps) at 1 bar, representing nominal 0 GPa. Temperature lists differ between folders.

**Continuation setup:** all alloy folders currently have `T1 = []` and expect an existing `<structure>/final.data` before the first NPT stage; this file is not supplied. Provide the appropriate equilibrated configuration, or configure an initial NVT stage for a fresh run. The elemental `B/` and `Cr/` folders include a 900 K NVT initialization. `Time3 = []` disables the final relaxation loop.

See [Melt_quench](../Melt_quench/README.md) for the shared analysis workflow, adapting any Fe-specific settings to Cr. No simulation outputs or large trajectories are included here.
