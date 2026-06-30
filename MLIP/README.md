# MLIP Workflow

This directory contains the machine-learned interatomic potential (MLIP) training, inference testing, and comparison-plot workflow used for the Fe-B project.

## Directory layout

```text
MLIP/
├── train/
├── test/
└── plot/
```

## Overview

The workflow is split into three stages:

1. **Training / fine-tuning**
2. **Model-based evaluation on XYZ structures**
3. **Plotting parity-style comparison figures for energy and force**

## Training stage

The `train/` folder does not currently contain a Python training launcher. Instead, training is driven by:

- `finetune.yaml`: MACE fine-tuning configuration
- `job-gpu`: PBS batch script that launches training with `mace_run_train`

### Key training inputs

- `train.xyz`: main training set
- `MatPES-R2SCAN-FeB.xyz`: auxiliary / pretraining reference set
- `mace-mh-1.model`: base foundation model

### Key training outputs

- `mace-mh-1_finetuned_FeB.model`
- `mace-mh-1_finetuned_FeB_mixed.model`
- `mace-mh-1_finetuned_FeB_compiled.model`
- `checkpoints/`
- `logs/`
- `results/`

### Training launch command

From `train/`:

```bash
qsub job-gpu
```

The batch script eventually runs:

```bash
mace_run_train --config finetune.yaml
```

## Python scripts

### `test/test.py`

Purpose:

- load a trained MACE model
- read one or more XYZ files
- evaluate predicted energy, force, and stress
- write paired reference/prediction CSV tables for downstream plotting

Main inputs:

- `MODEL_PATH` (default: `mace-mh-1_finetuned_FeB_mixed.model`)
- `XYZ_DIR` (default: current directory)
- `test.xyz` or other XYZ files in the chosen directory

Main outputs:

- `energy_per_atom_comparison.csv`
- `forces_comparison.csv`
- `stress_comparison.csv`

Notes:

- energy is converted to **per-atom** energy before writing
- force and stress are written as flattened component-wise comparisons
- the script first tries `cuda`, then falls back to `cpu`

Run from `test/`:

```bash
python test.py
```

### `plot/plot_comparison.py`

Purpose:

- read comparison CSV tables
- compute RMSE and MAE
- produce parity-style plots for model vs reference values

Current configuration:

- plots **energy per atom**
- plots **forces**
- does **not** plot stress in the packaged version

Main inputs:

- `data/csv/energy_per_atom_comparison.csv`
- `data/csv/forces_comparison.csv`

Main outputs:

- `Energy_plot.pdf`
- `Energy_plot.png`
- `Forces_plot.pdf`
- `Forces_plot.png`

Run from `plot/`:

```bash
python plot_comparison.py
```

## Recommended execution order

### 1. Train or fine-tune the model

From `train/`:

```bash
qsub job-gpu
```

### 2. Run model evaluation on the test structures

From `test/`:

```bash
python test.py
```

### 3. Copy or package the CSV files for plotting

The plotting package expects comparison CSV files under its local data folder.

### 4. Generate comparison figures

From `plot/`:

```bash
python plot_comparison.py
```

## Dependencies

### Python packages

- `numpy`
- `pandas`
- `matplotlib`
- `ase`
- `mace-torch` (or an environment that provides `mace.calculators.MACECalculator`)

### External environment

- CUDA-enabled GPU is preferred for testing and training, but `test.py` can fall back to CPU
- PBS scheduler is assumed by the provided `job-gpu` scripts

## Notes

- The training stage is configured primarily through YAML and batch scripts rather than a dedicated Python launcher.
- The two Python scripts in this module are focused on **testing** (`test.py`) and **post-test plotting** (`plot_comparison.py`).
