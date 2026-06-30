# AIRSS Structure Search Workflow

This directory contains the AIRSS-based structure-search workflow used for Fe-B compositions across a set of fixed stoichiometries.

## Directory layout

```text
Structure_Search/
├── 9_1/
├── 8.5_1.5/
├── 8_2/
├── 7.5_2.5/
├── 7_3/
├── 6_4/
├── 5_5/
├── 4_6/
├── 3_7/
├── 2_8/
└── 1_9/
```

Each composition folder is self-contained and typically includes:

- `BFe.cell`: AIRSS / `buildcell` input template
- `job_mace.pbs`: PBS batch script for the search job
- `mace_SS_airss.py`: Python driver for AIRSS generation plus MACE-based relaxation

## What the Python script does

The main search driver is repeated in each composition-specific folder:

- `1_9/mace_SS_airss.py`
- `2_8/mace_SS_airss.py`
- `3_7/mace_SS_airss.py`
- `4_6/mace_SS_airss.py`
- `5_5/mace_SS_airss.py`
- `6_4/mace_SS_airss.py`
- `7_3/mace_SS_airss.py`
- `7.5_2.5/mace_SS_airss.py`
- `8_2/mace_SS_airss.py`
- `8.5_1.5/mace_SS_airss.py`
- `9_1/mace_SS_airss.py`

These scripts follow the same workflow:

1. generate random candidate structures with AIRSS `buildcell`
2. convert generated `.cell` files to `.cif` using `cabal`
3. read each candidate with ASE
4. attach a MACE calculator
5. relax atomic positions and lattice vectors with `BFGS + ExpCellFilter`
6. keep successfully relaxed structures
7. write:
   - `<seed>_<index>_optimized.cif`
   - `<seed>_<index>.log`
   - `<seed>_<index>.res` in AIRSS-style result format

## Key script settings

All current `mace_SS_airss.py` files use the same main control parameters:

- `N = 1000`
  Number of AIRSS candidates to attempt in one run.
- `infile = 'BFe.cell'`
  AIRSS structure-generation template.
- `timeout 60s buildcell`
  Limits each AIRSS structure generation call to 60 seconds.
- `MACECalculator(..., device="cuda")`
  Uses the specified MACE model on GPU.
- `BFGS(..., maxstep=0.05)`
  Optimizer setup.
- `dyn.run(fmax=0.05, steps=1000)`
  Relaxation stopping criteria.

## Output products

For each accepted candidate, the search script may generate:

- `*.cif`
- `*_optimized.cif`
- `*.log`
- `*.traj` (removed automatically for converged structures)
- `*.res`

The `.res` file stores:

- pressure derived from the final stress tensor
- volume
- total potential energy
- atom count
- symmetry label from `spglib`
- lattice parameters
- fractional coordinates

## External dependencies

### Python packages

- `numpy`
- `ase`
- `mace-torch` / `mace`
- `spglib`

### External executables

- `buildcell` (AIRSS)
- `cabal` (AIRSS utility for `cell -> cif` conversion)

### HPC / scheduler environment

- PBS is assumed by the supplied `job_mace.pbs` files
- GPU availability is assumed because the scripts request `device="cuda"`

## Recommended usage

Inside one composition folder, for example `4_6/`:

```bash
cd 4_6
python mace_SS_airss.py
```

For batch execution on a cluster, use the accompanying PBS script if it matches your environment:

```bash
qsub job_mace.pbs
```

## Things to edit before running

Before launching a search in a new environment, check these fields in `mace_SS_airss.py`:

- `model_paths=...`
  Update the absolute path to the MACE model file.
- `device="cuda"`
  Change to CPU if no GPU is available.
- `N`
  Adjust the number of random candidates per batch.
- `infile`
  Ensure the correct AIRSS `.cell` template is present.

Also review `job_mace.pbs` for:

- project/account IDs
- queue name
- walltime
- GPU / CPU / memory requests

## Notes

- The Python search driver is duplicated across composition folders so that each stoichiometry can be submitted independently on HPC systems.
- If you update the search logic, update all composition copies consistently, or centralize the script before publication.
