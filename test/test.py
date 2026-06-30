import os
import glob
import numpy as np
import pandas as pd
from ase.io import read
from mace.calculators import MACECalculator

# ================= Configuration =================
MODEL_PATH = "mace-mh-1_finetuned_FeB_mixed.model"  # Update this to your model filename if needed
XYZ_DIR = "./"  # Update this to the directory containing the XYZ files if needed

# Output filenames
OUTPUT_ENERGY = "energy_per_atom_comparison.csv" 
OUTPUT_FORCES = "forces_comparison.csv"
OUTPUT_STRESS = "stress_comparison.csv"

# Reference key configuration
KEY_ENERGY_TRUE = "energy"       # Key in atoms.info (typically total energy)
KEY_FORCES_TRUE = "forces"       
KEY_STRESS_TRUE = "stress"       
# ============================================

def main():
    print(f"Loading model: {MODEL_PATH} ...")
    try:
        calculator = MACECalculator(model_paths=MODEL_PATH, device='cuda', default_dtype="float64")
    except Exception as e:
        print(f"Model loading failed: {e}")
        print("Falling back to CPU...")
        calculator = MACECalculator(model_paths=MODEL_PATH, device='cpu', default_dtype="float64")

    xyz_files = glob.glob(os.path.join(XYZ_DIR, "*.xyz"))
    if not xyz_files:
        print(f"No .xyz files were found in {XYZ_DIR}.")
        return

    print(f"Found {len(xyz_files)} XYZ files. Starting evaluation (energy per atom)...")

    data_energy = []  
    data_forces = []  
    data_stress = []  

    count = 0
    
    for xyz_file in xyz_files:
        try:
            atoms_list = read(xyz_file, index=":")
        except Exception as e:
            print(f"Failed to read {xyz_file}: {e}")
            continue

        for atoms in atoms_list:
            # Number of atoms in the current structure
            n_atoms = len(atoms)
            if n_atoms == 0:
                continue  # Skip empty structures

            # --- 1. Energy handling (convert to per-atom energy) ---
            e_true_per_atom = None
            if KEY_ENERGY_TRUE in atoms.info:
                # Assume the stored value is total energy and convert it to per-atom energy
                e_true_per_atom = atoms.info[KEY_ENERGY_TRUE] / n_atoms
            elif KEY_ENERGY_TRUE == "energy":
                try: 
                    e_true_per_atom = atoms.get_potential_energy() / n_atoms
                except: pass
            
            # --- 2. Force handling (forces are already per-atom quantities) ---
            if KEY_FORCES_TRUE in atoms.arrays:
                f_true = atoms.arrays[KEY_FORCES_TRUE].flatten()
            elif KEY_FORCES_TRUE == "forces":
                 try: f_true = atoms.get_forces().flatten()
                 except: f_true = None
            else:
                f_true = None

            # --- 3. Stress handling ---
            if KEY_STRESS_TRUE in atoms.info:
                s_true = np.array(atoms.info[KEY_STRESS_TRUE]).flatten()
            elif KEY_STRESS_TRUE == "stress":
                try: s_true = atoms.get_stress().flatten()
                except: s_true = None
            else:
                s_true = None

            # --- Predictions ---
            atoms.calc = calculator
            
            # Predict total energy and convert to per-atom energy
            try:
                e_total_pred = atoms.get_potential_energy()
                e_pred_per_atom = e_total_pred / n_atoms
                
                if e_true_per_atom is not None:
                    data_energy.append([e_true_per_atom, e_pred_per_atom])
            except Exception as e:
                print(f"Energy evaluation failed: {e}")

            # Predict forces
            try:
                f_pred = atoms.get_forces().flatten()
                if f_true is not None and len(f_true) == len(f_pred):
                    for ft, fp in zip(f_true, f_pred):
                        data_forces.append([ft, fp])
            except Exception as e:
                pass

            # Predict stress
            try:
                s_pred = atoms.get_stress().flatten()
                if s_true is not None and len(s_true) == len(s_pred):
                    for st, sp in zip(s_true, s_pred):
                        data_stress.append([st, sp])
            except Exception as e:
                pass

            count += 1
            if count % 100 == 0:
                print(f"Processed {count} structures...")

    # Save outputs
    print("Saving outputs...")

    if data_energy:
        df_e = pd.DataFrame(data_energy, columns=["True_Energy_per_atom", "Pred_Energy_per_atom"])
        df_e.to_csv(OUTPUT_ENERGY, index=False)
        print(f"Energy-per-atom comparison written to: {OUTPUT_ENERGY}")
    else:
        print("Warning: no energy data were collected.")

    if data_forces:
        df_f = pd.DataFrame(data_forces, columns=["True_Force", "Pred_Force"])
        df_f.to_csv(OUTPUT_FORCES, index=False)
        print(f"Force comparison written to: {OUTPUT_FORCES}")

    if data_stress:
        df_s = pd.DataFrame(data_stress, columns=["True_Stress", "Pred_Stress"])
        df_s.to_csv(OUTPUT_STRESS, index=False)
        print(f"Stress comparison written to: {OUTPUT_STRESS}")

    print("Done.")

if __name__ == "__main__":
    main()
