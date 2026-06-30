#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

def main():
    try:
        from pymatgen.core import Structure
        from pymatgen.io.vasp import Poscar
    except ImportError:
        print("pymatgen was not found. Please install it first: pip install pymatgen")
        sys.exit(1)

    cwd = Path(".").resolve()
    cif_files = sorted(list(cwd.glob("*.cif")) + list(cwd.glob("*.CIF")))

    if not cif_files:
        print("No .cif files were found in the current directory.")
        return

    success, failed = 0, 0
    for cif_path in cif_files:
        stem = cif_path.stem  # Filename without suffix
        out_path = cif_path.with_suffix(".vasp")  # Write a same-name .vasp file in POSCAR format

        try:
            structure = Structure.from_file(str(cif_path))
            # Write VASP POSCAR format (VASP5 style with element labels)
            Poscar(structure).write_file(str(out_path))
            print(f"[OK] {cif_path.name} -> {out_path.name}")
            success += 1
        except Exception as e:
            print(f"[FAIL] {cif_path.name}: {e}")
            failed += 1

    print(f"\nDone: {success} succeeded, {failed} failed. Outputs are same-name .vasp files in POSCAR format.")

if __name__ == "__main__":
    main()
