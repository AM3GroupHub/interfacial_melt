#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LAMMPS average energy workflow with explicit PV output.

This script mirrors ``get_avE.py`` and adds block-averaged PV energy terms
for high-pressure calculations.
"""

import argparse
import os
import re

import numpy as np

import get_avE as base


def parse_lammps_log(log_file, skip_lines, block_size):
    """Parse PotEng, PV, and H = PotEng + KinEng + PV from a LAMMPS log."""
    pe_values = []
    pv_values = []
    enthalpy_values = []
    in_thermo_section = False
    col_indices = {}

    with open(log_file, "r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        lowered = [base.normalize_header(part) for part in parts]
        if "step" in lowered and any(key in lowered for key in ("poteng", "pe")):
            col_indices = {base.normalize_header(header): index for index, header in enumerate(parts)}
            in_thermo_section = True
            continue

        if line.startswith(("Loop time", "WARNING", "ERROR")):
            in_thermo_section = False
            col_indices = {}
            continue

        if not in_thermo_section or not parts or not base.is_float(parts[0]):
            continue

        try:
            pot_i = base.column_index(col_indices, ("PotEng", "pe"))
            kin_i = base.column_index(col_indices, ("KinEng", "ke"))
            press_i = base.column_index(col_indices, ("Press", "pressure"))
            vol_i = base.column_index(col_indices, ("Volume", "vol"))

            pot_eng = float(parts[pot_i].replace("D", "E")) if pot_i is not None else np.nan
            if pot_i is not None:
                pe_values.append(pot_eng)

            if None not in (press_i, vol_i):
                press = float(parts[press_i].replace("D", "E"))
                volume = float(parts[vol_i].replace("D", "E"))
                pv_energy = press * volume * base.PV_BAR_A3_TO_EV
                pv_values.append(pv_energy)
            else:
                pv_energy = np.nan

            if None not in (pot_i, kin_i, press_i, vol_i):
                kin_eng = float(parts[kin_i].replace("D", "E"))
                enthalpy_values.append(pot_eng + kin_eng + pv_energy)
        except (ValueError, IndexError):
            continue

    avg_pe, pe_err, pe_samples, pe_blocks = stats_or_nan(pe_values, skip_lines, block_size, "potential energy")
    avg_pv, pv_err, pv_samples, pv_blocks = stats_or_nan(pv_values, skip_lines, block_size, "PV term")
    avg_h, h_err, _, h_blocks = stats_or_nan(enthalpy_values, skip_lines, block_size, "enthalpy")

    print(f"  - Read {len(pe_values)} potential-energy rows; skipped the first {skip_lines}")
    print(f"  - Potential-energy block error: {base.number_or_na(pe_err)} eV ({pe_blocks} blocks, block_size={block_size})")
    print(f"  - PV-term block error: {base.number_or_na(pv_err)} eV ({pv_blocks} blocks, block_size={block_size})")
    print(f"  - Enthalpy block error: {base.number_or_na(h_err)} eV ({h_blocks} blocks, block_size={block_size})")
    if pv_samples == 0:
        print("  - Warning: Press/Volume columns were not found, so the PV term cannot be computed")

    return avg_pe, pe_err, avg_pv, pv_err, avg_h, h_err, pe_samples, pe_blocks


def stats_or_nan(values, skip_lines, block_size, label):
    if len(values) <= skip_lines:
        print(f"  - Warning: only {len(values)} {label} rows were read, which is not enough after skipping {skip_lines}")
        return np.nan, np.nan, 0, 0
    return base.block_average_stats(values, skip_lines, block_size)


def process_directories(args):
    results = []
    pattern = re.compile(r"(\d+(?:\.\d+)?)_dir$")

    for item in sorted(os.listdir(os.getcwd())):
        match = pattern.match(item)
        if not match or not os.path.isdir(item):
            continue

        directory = item
        dir_temperature = float(match.group(1))
        input_path = os.path.join(directory, args.input)
        log_path = os.path.join(directory, args.filename)
        print(f"\nProcessing directory: {directory}")

        if os.path.exists(input_path):
            settings = base.parse_lammps_input(input_path, args.dump)
        else:
            print(f"  - Warning: {args.input} was not found, so vibrational entropy cannot be evaluated")
            settings = {
                "temperature": None,
                "dump_file": args.dump,
                "masses": {},
                "units": "metal",
                "timestep": None,
                "dump_every": None,
            }

        temperature = settings["temperature"] if settings["temperature"] is not None else dir_temperature
        print(f"  - Temperature: {temperature:.3f} K")

        if os.path.exists(log_path):
            avg_pe, pe_err, avg_pv, pv_err, avg_h, h_err, n_samples, n_blocks = parse_lammps_log(
                log_path, args.skip, args.block_size
            )
            print(f"  - Average potential energy: {base.number_or_na(avg_pe)} eV")
            print(f"  - Average PV term: {base.number_or_na(avg_pv)} eV")
            print(f"  - PV-term error: {base.number_or_na(pv_err)} eV")
            print(f"  - Average enthalpy: {base.number_or_na(avg_h)} eV")
        else:
            avg_pe = pe_err = avg_pv = pv_err = avg_h = h_err = np.nan
            n_samples, n_blocks = 0, 0
            print(f"  - Warning: {args.filename} was not found")

        try:
            thermo = base.analyze_vibrations(directory, settings, temperature, args)
            s_vib = thermo["S_eV_per_K"]
        except Exception as exc:
            s_vib = np.nan
            print(f"  - Warning: vibrational entropy evaluation failed: {exc}")

        ts_vib = temperature * s_vib if not np.isnan(s_vib) else np.nan
        free_energy = avg_h - ts_vib if not np.isnan(avg_h) and not np.isnan(ts_vib) else np.nan
        free_energy_err = h_err
        print(f"  - T*S_vib: {base.number_or_na(ts_vib)} eV")
        print(f"  - Free energy H-T*S_vib: {base.number_or_na(free_energy)} eV")
        print(f"  - Free-energy error: {base.number_or_na(free_energy_err)} eV")
        results.append(
            (
                temperature,
                avg_pe,
                pe_err,
                avg_pv,
                pv_err,
                avg_h,
                h_err,
                s_vib,
                ts_vib,
                free_energy,
                free_energy_err,
                n_samples,
                n_blocks,
                args.block_size,
            )
        )

    write_results(args.output, results)


def write_results(output, results):
    results.sort(key=lambda row: row[0])
    with open(output, "w", encoding="utf-8") as handle:
        handle.write(
            "# Temperature(K)\tAverage_PotEng(eV)\tAverage_PotEng_err(eV)\t"
            "Average_PV(eV)\tAverage_PV_err(eV)\t"
            "Average_Enthalpy(eV)\tAverage_Enthalpy_err(eV)\tS_vib(eV/K)\t"
            "T*S_vib(eV)\tFree_Energy_H_minus_TS(eV)\tFree_Energy_err(eV)\t"
            "N_samples\tN_blocks\tBlock_size\n"
        )
        for row in results:
            values = [f"{row[0]:.2f}"] + [base.number_or_na(value) for value in row[1:11]]
            values += [str(row[11]), str(row[12]), str(row[13])]
            handle.write("\t".join(values) + "\n")

    print(f"\n{'=' * 60}")
    print("Processing complete.")
    print(f"Results written to: {output}")
    print(f"Processed {len(results)} temperature points")
    print(f"{'=' * 60}")


def main():
    parser = argparse.ArgumentParser(
        description="Process LAMMPS log/dump/in.lammps files and output averaged energies, the PV term, VDOS, vibrational entropy, and free energy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-s", "--skip", type=int, default=0, help="Skip the first N thermo lines in the LAMMPS log")
    parser.add_argument("-f", "--filename", type=str, default="log.lammps", help="LAMMPS log filename")
    parser.add_argument("-o", "--output", type=str, default="tot_energy_vdos_results.txt", help="Summary output filename")
    parser.add_argument("--input", type=str, default="in.lammps", help="LAMMPS input filename")
    parser.add_argument("--dump", type=str, default="SUPERCELL.dump", help="Default trajectory dump filename")
    parser.add_argument("--vdos-output", type=str, default="vdos_eV.dat", help="Per-directory VDOS output filename")
    parser.add_argument("--fe-output", type=str, default="vdos_fe_results.txt", help="Per-directory VDOS thermodynamic summary filename")
    parser.add_argument("--block-size", type=int, default=1000, help="Block size for block averaging, in thermo-output lines")
    parser.add_argument("--no-window", action="store_true", help="Disable the Hann window when computing the VDOS")
    args = parser.parse_args()
    if args.block_size <= 0:
        parser.error("--block-size must be a positive integer")
    process_directories(args)


if __name__ == "__main__":
    main()
