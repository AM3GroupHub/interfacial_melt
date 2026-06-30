# Melt_quench Workflow

This directory contains the full melt-quench molecular dynamics workflow used for the Fe-B liquid and glass-forming datasets, from structure preparation to post-processing and figure reproduction.

## Directory layout

```text
Melt_quench/
├── MD/
│   ├── Fe/
│   ├── B/
│   ├── 9_1/
│   ├── 8.5_1.5/
│   ├── 8_2/
│   ├── 7.5_2.5/
│   ├── 7_3/
│   ├── 6_4/
│   ├── 5_5/
│   ├── 4_6/
│   ├── 3_7/
│   ├── 2_8/
│   └── 1_9/
├── script/
└── pictures/
```

## Workflow overview

### 1. Prepare starting structures

Each subdirectory under `MD/` is a composition-specific working folder. The same set of helper scripts is shipped in every `MD/<composition>/` folder.

1. Convert source CIF files to VASP/POSCAR format:

```bash
python cif2vasp.py
```

2. Build the supercell and generate LAMMPS data files.

3. Create LAMMPS input decks for the desired NVT/NPT stages.

4. Run the melt-quench schedule.

### 2. Post-process MD outputs

The `script/` folder contains the trajectory-analysis workflow:

- average potential energy / enthalpy / vibrational entropy / free energy
- explicit PV-term post-processing for high-pressure runs
- density, Onsager transport coefficients, self-diffusion, and RDFs
- B-network frame sampling and bond-matrix extraction
- standalone VDOS calculation from dump files

### 3. Reproduce analysis figures

The `pictures/` folder contains upload-ready, self-contained plotting packages. Each figure package ships its own script, minimal raw data, and a local README.

## Python script inventory

### `script/`

| Script | Purpose | Typical outputs |
| --- | --- | --- |
| `get_avE.py` | Parse `log.lammps`, `SUPERCELL.dump`, and `in.lammps` across `*_dir` folders; compute block-averaged potential energy, enthalpy, VDOS, vibrational entropy, and free energy | `energy_vdos_results.txt`, per-directory `vdos_eV.dat`, `vdos_fe_results.txt` |
| `tot_avE.py` | Same workflow as `get_avE.py`, but also computes the PV contribution explicitly | `tot_energy_vdos_results.txt`, per-directory VDOS outputs |
| `get_transport.py` | Compute density, full Onsager matrix, self-diffusion, and RDFs for `*_dir` trajectories | transport summary tables, `rdf.dat`, Onsager outputs |
| `sample_b_frames.py` | Sample 10 late-time frames from a parent trajectory folder and compute B-B bond matrices / B coordination | `b_frame_analysis/` subfolders with bond matrices and coordination tables |
| `vdos_from_dump.py` | Standalone VDOS calculator from a LAMMPS dump containing velocities | `vdos_eV.dat` |

### `MD/<composition>/`

The following four scripts are repeated in every `MD/<composition>/` folder.

| Script | Purpose |
| --- | --- |
| `cif2vasp.py` | Convert local CIF files to `.vasp` / POSCAR-format inputs using `pymatgen` |
| `gen_sc.py` | Build supercells, prepare `SUPERCELL.xyz`, `SUPERCELL.lmp`, and `POSCAR` |
| `gen_input.py` | Generate LAMMPS `in.lammps` files for NVT and NPT stages |
| `cal_MQ.py` | Orchestrate the melt-quench run sequence for the current composition |

### `pictures/`

Each subfolder contains one figure package with its own plotting script and packaged raw data.

| Folder | Plot script | Notes |
| --- | --- | --- |
| `energy_composition_0gpa_liquid/` | `plot_energy_composition_0gpa_liquid.py` | 0 GPa energy-composition figures |
| `energy_composition_10gpa_liquid/` | `plot_energy_composition_10gpa_liquid.py` | 10 GPa energy-composition figures |
| `energy_composition_crb_1200_1500_1800/` | `plot_energy_composition_crb_1200_1500_1800.py` | Cr-B selected-temperature figures |
| `energy_composition_0gpa_1200_1500_1800_crb_range/` | `plot_energy_composition_0gpa_1200_1500_1800_crb_range.py` | Fe-B high-B window matching Cr-B x-range |
| `free_energy_0gpa/` | `plot_free_energy_0gpa_1500_1800.py`, `plot_free_energy_0gpa_1500_1800_with_config_entropy.py` | 0 GPa free-energy figures |
| `free_energy_10gpa/` | `plot_free_energy_10gpa_1500_1800.py`, `plot_free_energy_10gpa_1500_1800_with_config_entropy.py` | 10 GPa free-energy figures |
| `omat_1500k_0GPa/` | `plot_energy_composition_b_1500k.py` | OMAT comparison at 1500 K, 0 GPa |
| `omat_1500k_10GPa/` | `plot_energy_composition_b_1500k.py` | OMAT comparison at 1500 K, 10 GPa |
| `pv_composition_10gpa_liquid/` | `plot_pv_composition_10gpa_liquid.py` | 10 GPa PV-composition figure |
| `pv_potential_overlay_10gpa_1500k_rightonly/` | `plot_pv_potential_overlay_10gpa_1500k_rightonly.py` | 10 GPa right-only PV/potential overlay |
| `rdf_by_composition_0gpa/` | `plot_rdf_by_composition_0gpa.py` | 0 GPa RDF-by-composition package |
| `rdf_by_composition_10gpa/` | `plot_rdf_by_composition_10gpa.py` | 10 GPa RDF-by-composition package |
| `sixfold_fraction_energy_overlay_0gpa/` | `plot_sixfold_fraction_energy_overlay_0gpa.py` | 0 GPa sixfold-fraction + potential overlays |
| `sixfold_fraction_energy_overlay_1500k/` | `plot_sixfold_fraction_energy_overlay_1500k.py` | 1500 K compact sixfold overlay |
| `sixfold_fraction_energy_overlay_1800k/` | `plot_sixfold_fraction_energy_overlay_1800k.py` | 1800 K compact sixfold overlay |

## Recommended execution order

### Structure preparation and MD

Inside one `MD/<composition>/` folder:

```bash
python cif2vasp.py
python cal_MQ.py
```

Notes:

- `cal_MQ.py` imports `gen_sc.py` and `gen_input.py`
- update `p_path` to the MACE/LAMMPS model you want to use
- update `run_lmp` to the LAMMPS executable and MPI launch command available on your system
- adjust `T1`, `T2`, `T3`, `Time1`, `Time2`, and `Time3` to match your melt, quench, and relaxation schedule

### Thermodynamic post-processing

Run these from a folder that contains temperature subdirectories such as `1200.0_dir/`, `1500.0_dir/`, etc.

```bash
python ../../script/get_avE.py --input in.lammps --dump SUPERCELL.dump -f log.lammps
python ../../script/tot_avE.py --input in.lammps --dump SUPERCELL.dump -f log.lammps
python ../../script/get_transport.py --input in.lammps --dump SUPERCELL.dump -f log.lammps
python ../../script/sample_b_frames.py .
```

### Figure reproduction

Each packaged figure folder under `pictures/` contains its own README and can be run independently from the folder itself.

## Dependencies

### Python packages

- `numpy`
- `matplotlib`
- `ase`
- `pymatgen` (required by `cif2vasp.py`)

### External programs

- `phonopy`
- `LAMMPS`
- MPI launcher (for example `mpirun`, if used by your environment)
- a MACE-compatible LAMMPS potential model file

## Notes

- The `MD/` helper scripts are duplicated across composition folders so that each composition directory remains self-contained for batch execution on HPC systems.
- The `pictures/` packages are already reduced to upload-ready plotting bundles and can be published to GitHub directly.
