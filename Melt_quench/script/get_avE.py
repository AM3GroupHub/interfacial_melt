#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAMMPS energy + vibrational entropy workflow.

For each ``*_dir`` directory this script reads ``log.lammps``,
``SUPERCELL.dump`` and ``in.lammps`` by default, then writes averaged energy,
VDOS, vibrational thermodynamics and total free energy.
"""

import argparse
import os
import re
import numpy as np

from vdos_from_dump import compute_vdos as compute_unit_vdos
from vdos_from_dump import parse_lammps_dump_velocity, remove_com_drift

KB_EV_PER_K = 8.617333262145e-5
EV_TO_J = 1.602176634e-19
NA = 6.02214076e23
PV_BAR_A3_TO_EV = 6.24150913e-7


def is_float(text):
    try:
        float(text.replace("D", "E"))
        return True
    except ValueError:
        return False


def normalize_header(name):
    return name.strip().lower().replace("_", "")


def column_index(columns, aliases):
    for alias in aliases:
        key = normalize_header(alias)
        if key in columns:
            return columns[key]
    return None


def block_average_stats(values, skip_lines, block_size):
    """Return mean and block-averaged standard error after skipping equilibration."""
    if block_size <= 0:
        raise ValueError("block_size must be a positive integer")

    data = np.asarray(values[skip_lines:], dtype=float)
    if data.size == 0:
        return np.nan, np.nan, 0, 0

    mean = float(np.mean(data))
    n_blocks = data.size // block_size
    if n_blocks < 2:
        return mean, np.nan, int(data.size), int(n_blocks)

    used = n_blocks * block_size
    block_means = data[:used].reshape(n_blocks, block_size).mean(axis=1)
    err = float(np.std(block_means, ddof=1) / np.sqrt(n_blocks))
    return mean, err, int(data.size), int(n_blocks)


def parse_lammps_log(log_file, skip_lines, block_size):
    """Parse PotEng and H = PotEng + KinEng + PV from a LAMMPS log."""
    pe_values = []
    enthalpy_values = []
    in_thermo_section = False
    col_indices = {}

    with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        lowered = [normalize_header(p) for p in parts]
        if "step" in lowered and any(k in lowered for k in ("poteng", "pe")):
            col_indices = {normalize_header(header): i for i, header in enumerate(parts)}
            in_thermo_section = True
            continue

        if line.startswith(("Loop time", "WARNING", "ERROR")):
            in_thermo_section = False
            col_indices = {}
            continue

        if not in_thermo_section or not parts or not is_float(parts[0]):
            continue

        try:
            pot_i = column_index(col_indices, ("PotEng", "pe"))
            kin_i = column_index(col_indices, ("KinEng", "ke"))
            press_i = column_index(col_indices, ("Press", "pressure"))
            vol_i = column_index(col_indices, ("Volume", "vol"))
            if pot_i is not None:
                pe_values.append(float(parts[pot_i].replace("D", "E")))
            if None not in (pot_i, kin_i, press_i, vol_i):
                pot_eng = float(parts[pot_i].replace("D", "E"))
                kin_eng = float(parts[kin_i].replace("D", "E"))
                press = float(parts[press_i].replace("D", "E"))
                volume = float(parts[vol_i].replace("D", "E"))
                enthalpy_values.append(pot_eng + kin_eng + press * volume * PV_BAR_A3_TO_EV)
        except (ValueError, IndexError):
            continue

    if len(pe_values) > skip_lines:
        avg_pe, pe_err, pe_samples, pe_blocks = block_average_stats(pe_values, skip_lines, block_size)
        print(f"  - Read {len(pe_values)} potential-energy rows; skipped the first {skip_lines}")
        print(f"  - Potential-energy block error: {number_or_na(pe_err)} eV ({pe_blocks} blocks, block_size={block_size})")
    else:
        avg_pe, pe_err, pe_samples, pe_blocks = np.nan, np.nan, 0, 0
        print(f"  - Warning: only {len(pe_values)} potential-energy rows were read, which is not enough after skipping {skip_lines}")

    if len(enthalpy_values) > skip_lines:
        avg_enthalpy, enthalpy_err, enthalpy_samples, enthalpy_blocks = block_average_stats(
            enthalpy_values, skip_lines, block_size
        )
        print(
            f"  - Enthalpy block error: {number_or_na(enthalpy_err)} eV "
            f"({enthalpy_blocks} blocks, block_size={block_size})"
        )
    else:
        avg_enthalpy, enthalpy_err, enthalpy_samples, enthalpy_blocks = np.nan, np.nan, 0, 0

    return avg_pe, pe_err, avg_enthalpy, enthalpy_err, pe_samples, pe_blocks


def strip_comment(line):
    return line.split("#", 1)[0].strip()


def commands_from_input(lines):
    commands = []
    current = ""
    for raw in lines:
        line = strip_comment(raw)
        if not line and not current:
            continue
        if line.endswith("&"):
            current += line[:-1].strip() + " "
            continue
        current += line
        if current.strip():
            commands.append(current.strip())
        current = ""
    if current.strip():
        commands.append(current.strip())
    return commands


def format_variable_value(value):
    if float(value).is_integer():
        return str(int(value))
    return f"{value:g}"


def resolve_string_token(token, variables):
    token = token.strip().strip("'\"")

    def repl_braced(match):
        name = match.group(1)
        return format_variable_value(variables[name]) if name in variables else match.group(0)

    token = re.sub(r"\$\{([A-Za-z_]\w*)\}", repl_braced, token)
    token = re.sub(
        r"\$([A-Za-z_]\w*)",
        lambda m: format_variable_value(variables[m.group(1)]) if m.group(1) in variables else m.group(0),
        token,
    )
    return token


def parse_numeric_expr(expr, variables):
    if expr is None:
        return None
    expr = resolve_string_token(str(expr), variables)
    expr = re.sub(
        r"\bv_([A-Za-z_]\w*)\b",
        lambda m: format_variable_value(variables[m.group(1)]) if m.group(1) in variables else m.group(0),
        expr,
    )
    expr = expr.replace("D", "E").replace("d", "e")
    if not re.fullmatch(r"[0-9eE+\-*/().\s]+", expr):
        return None
    try:
        return float(eval(expr, {"__builtins__": {}}, {}))
    except Exception:
        return None


def parse_lammps_input(input_file, default_dump):
    """Read units, timestep, masses, temperature and dump settings from in.lammps."""
    settings = {
        "units": "metal",
        "timestep": None,
        "dump_every": None,
        "dump_file": default_dump,
        "masses": {},
        "temperature": None,
        "variables": {},
    }
    temperatures = []

    with open(input_file, "r", encoding="utf-8", errors="ignore") as f:
        commands = commands_from_input(f.readlines())

    variables = settings["variables"]
    for command in commands:
        parts = command.split()
        if not parts:
            continue
        lower = [p.lower() for p in parts]
        cmd = lower[0]

        if cmd == "variable" and len(parts) >= 4:
            value = parse_numeric_expr(" ".join(parts[3:]), variables)
            if value is not None:
                variables[parts[1]] = value
            continue

        if cmd == "units" and len(parts) >= 2:
            settings["units"] = parts[1].lower()
        elif cmd == "timestep" and len(parts) >= 2:
            settings["timestep"] = parse_numeric_expr(parts[1], variables)
        elif cmd == "mass" and len(parts) >= 3:
            atom_type = parse_numeric_expr(parts[1], variables)
            mass = parse_numeric_expr(parts[2], variables)
            if atom_type is not None and mass is not None:
                settings["masses"][int(atom_type)] = mass
        elif cmd == "dump" and len(parts) >= 6:
            every = parse_numeric_expr(parts[4], variables)
            if every is not None:
                settings["dump_every"] = int(round(every))
            settings["dump_file"] = resolve_string_token(parts[5], variables)
        elif cmd == "fix" and "temp" in lower:
            i = lower.index("temp")
            vals = [parse_numeric_expr(parts[j], variables) for j in (i + 1, i + 2) if j < len(parts)]
            vals = [v for v in vals if v is not None]
            if vals:
                temperatures.append(float(np.mean(vals)))
        elif cmd == "velocity" and "create" in lower:
            i = lower.index("create")
            if i + 1 < len(parts):
                temp = parse_numeric_expr(parts[i + 1], variables)
                if temp is not None:
                    temperatures.append(temp)

    if temperatures:
        settings["temperature"] = temperatures[-1]
    else:
        for key, value in variables.items():
            if key.lower() in ("t", "temp", "temperature"):
                settings["temperature"] = value
                break
    return settings


def timestep_to_fs(units, timestep):
    factors = {
        "real": 1.0,
        "metal": 1000.0,
        "si": 1.0e15,
        "cgs": 1.0e15,
        "electron": 1.0,
    }
    unit = (units or "metal").lower()
    if unit not in factors:
        raise ValueError(f"Unsupported LAMMPS units for timestep conversion: {units}")
    return timestep * factors[unit]


def build_masses(types, mass_map):
    masses = np.ones(len(types), dtype=float)
    if not mass_map:
        return masses
    for i, atom_type in enumerate(types):
        atom_type = int(atom_type)
        if atom_type not in mass_map:
            raise ValueError(f"Atom type {atom_type} missing in in.lammps mass commands")
        masses[i] = mass_map[atom_type]
    return masses


def compute_total_vdos(vel, dt_s, masses=None, use_hann=True):
    e_ev, unit_vdos = compute_unit_vdos(vel, dt_s, masses=masses, use_hann=use_hann)
    return e_ev, unit_vdos * (3.0 * vel.shape[1])


def compute_vib_thermo(e_ev, vdos, temperature, emin_eps=1e-12):
    if temperature <= 0:
        raise ValueError("Temperature must be positive for vibrational entropy")
    integral_g = float(np.trapz(vdos, e_ev))
    e_clip = np.clip(e_ev, emin_eps, None)
    kbt_ev = KB_EV_PER_K * temperature
    xarg = np.clip(e_clip / kbt_ev, 1e-12, 1e3)

    ln_one_minus_exp_negx = np.log(-np.expm1(-xarg))
    bose_n = 1.0 / np.expm1(xarg)
    integrand_ln2sinh = np.log(2.0 * np.sinh(0.5 * xarg))
    f_integrand = 0.5 * e_clip + kbt_ev * ln_one_minus_exp_negx
    u_integrand = 0.5 * e_clip + e_clip * bose_n
    s_integrand = KB_EV_PER_K * (xarg * bose_n - ln_one_minus_exp_negx)

    i_ln2sinh = float(np.trapz(integrand_ln2sinh * vdos, e_ev))
    f_via_sinh = kbt_ev * i_ln2sinh
    f_ev = float(np.trapz(f_integrand * vdos, e_ev))
    u_ev = float(np.trapz(u_integrand * vdos, e_ev))
    s_ev_per_k = float(np.trapz(s_integrand * vdos, e_ev))
    s_check = (u_ev - f_ev) / temperature

    return {
        "integral_g": integral_g,
        "I_ln2sinh": i_ln2sinh,
        "F_via_sinh_eV": f_via_sinh,
        "F_eV": f_ev,
        "U_eV": u_ev,
        "S_eV_per_K": s_ev_per_k,
        "S_kB": s_ev_per_k / KB_EV_PER_K,
        "consistency_eV_per_K": s_ev_per_k - s_check,
    }


def write_vdos_fe_results(path, temperature, unit, thermo):
    f_kj_mol = thermo["F_eV"] * EV_TO_J * NA / 1000.0
    u_kj_mol = thermo["U_eV"] * EV_TO_J * NA / 1000.0
    s_j_mol_k = thermo["S_eV_per_K"] * EV_TO_J * NA
    s_kj_mol_k = s_j_mol_k / 1000.0
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"Temperature: {temperature:.3f} K\n")
        f.write(f"Unit of x-axis: {unit}\n")
        f.write(f"Integral of g(E) dE (raw): {thermo['integral_g']:.6f}\n")
        f.write("Using total-cell VDOS normalized to 3N modes.\n")
        f.write(f"I = ∫ ln[2 sinh(E/(2 k_B T))] g(E) dE = {thermo['I_ln2sinh']:.8f} (dimensionless)\n")
        f.write(f"F (via ln 2sinh) = {thermo['F_via_sinh_eV']:.8f} eV\n")
        f.write(f"F (textbook)     = {thermo['F_eV']:.8f} eV   ({f_kj_mol:.6f} kJ/mol)\n")
        f.write(f"U_vib            = {thermo['U_eV']:.8f} eV   ({u_kj_mol:.6f} kJ/mol)\n")
        f.write(
            f"S_vib            = {thermo['S_eV_per_K']:.8f} eV/K  "
            f"({s_kj_mol_k:.6f} kJ/mol/K, {s_j_mol_k:.6f} J/mol/K)\n"
        )
        f.write(f"S_vib in k_B per total cell = {thermo['S_kB']:.6f} k_B\n")
        f.write(f"Consistency check: S - (U-F)/T = {thermo['consistency_eV_per_K']:.3e} eV/K\n")


def output_path_for_dir(directory, filename):
    return os.path.join(directory, os.path.basename(filename) if os.path.isabs(filename) else filename)


def frame_dt_fs(settings, timesteps):
    if settings["timestep"] is None:
        raise ValueError("Cannot find timestep in in.lammps")
    if settings["dump_every"] is not None:
        frame_steps = settings["dump_every"]
    elif len(timesteps) > 1:
        frame_steps = int(round(float(np.median(np.diff(timesteps)))))
    else:
        raise ValueError("Cannot determine dump frame interval")
    return timestep_to_fs(settings["units"], settings["timestep"]) * frame_steps


def analyze_vibrations(directory, settings, temperature, args):
    dump_path = os.path.join(directory, settings["dump_file"])
    if not os.path.exists(dump_path) and settings["dump_file"] != args.dump:
        dump_path = os.path.join(directory, args.dump)
    if not os.path.exists(dump_path):
        raise FileNotFoundError(f"Trajectory file not found: {dump_path}")

    timesteps, _, types_ref, vel = parse_lammps_dump_velocity(dump_path)
    dt_fs = frame_dt_fs(settings, timesteps)
    masses = build_masses(types_ref, settings["masses"])
    vel_corr = remove_com_drift(vel, masses)
    e_ev, vdos = compute_total_vdos(vel_corr, dt_fs * 1.0e-15, masses=masses, use_hann=not args.no_window)
    thermo = compute_vib_thermo(e_ev, vdos, temperature)

    vdos_path = output_path_for_dir(directory, args.vdos_output)
    fe_path = output_path_for_dir(directory, args.fe_output)
    np.savetxt(vdos_path, np.column_stack([e_ev, vdos]), header="E_eV  VDOS_total_3N", comments="")
    write_vdos_fe_results(fe_path, temperature, "eV", thermo)
    print(f"  - VDOS: {os.path.basename(vdos_path)}")
    print(f"  - VDOS thermodynamic check: {os.path.basename(fe_path)}")
    print(f"  - Trajectory frames: {vel.shape[0]}, atoms: {vel.shape[1]}, dt: {dt_fs:.6g} fs")
    print(f"  - S_vib: {thermo['S_eV_per_K']:.8f} eV/K")
    return thermo


def number_or_na(value):
    return "N/A" if np.isnan(value) else f"{value:.6f}"


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

        settings = None
        if os.path.exists(input_path):
            settings = parse_lammps_input(input_path, args.dump)
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
            avg_pe, pe_err, avg_enthalpy, enthalpy_err, n_samples, n_blocks = parse_lammps_log(
                log_path, args.skip, args.block_size
            )
            print(f"  - Average potential energy: {number_or_na(avg_pe)} eV")
            print(f"  - Potential-energy error: {number_or_na(pe_err)} eV")
            print(f"  - Average enthalpy: {number_or_na(avg_enthalpy)} eV")
            print(f"  - Enthalpy error: {number_or_na(enthalpy_err)} eV")
        else:
            avg_pe, pe_err, avg_enthalpy, enthalpy_err = np.nan, np.nan, np.nan, np.nan
            n_samples, n_blocks = 0, 0
            print(f"  - Warning: {args.filename} was not found")

        try:
            thermo = analyze_vibrations(directory, settings, temperature, args)
            s_vib = thermo["S_eV_per_K"]
        except Exception as exc:
            thermo = None
            s_vib = np.nan
            print(f"  - Warning: vibrational entropy evaluation failed: {exc}")

        ts_vib = temperature * s_vib if not np.isnan(s_vib) else np.nan
        free_energy = avg_enthalpy - ts_vib if not np.isnan(avg_enthalpy) and not np.isnan(ts_vib) else np.nan
        free_energy_err = enthalpy_err
        print(f"  - T*S_vib: {number_or_na(ts_vib)} eV")
        print(f"  - Free energy H-T*S_vib: {number_or_na(free_energy)} eV")
        print(f"  - Free-energy error: {number_or_na(free_energy_err)} eV")
        results.append(
            (
                temperature,
                avg_pe,
                pe_err,
                avg_enthalpy,
                enthalpy_err,
                s_vib,
                ts_vib,
                free_energy,
                free_energy_err,
                n_samples,
                n_blocks,
                args.block_size,
            )
        )

    results.sort(key=lambda x: x[0])
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(
            "# Temperature(K)\tAverage_PotEng(eV)\tAverage_PotEng_err(eV)\t"
            "Average_Enthalpy(eV)\tAverage_Enthalpy_err(eV)\tS_vib(eV/K)\t"
            "T*S_vib(eV)\tFree_Energy_H_minus_TS(eV)\tFree_Energy_err(eV)\t"
            "N_samples\tN_blocks\tBlock_size\n"
        )
        for row in results:
            values = [f"{row[0]:.2f}"] + [number_or_na(v) for v in row[1:9]]
            values += [str(row[9]), str(row[10]), str(row[11])]
            f.write("\t".join(values) + "\n")

    print(f"\n{'=' * 60}")
    print("Processing complete.")
    print(f"Results written to: {args.output}")
    print(f"Processed {len(results)} temperature points")
    print(f"{'=' * 60}")


def main():
    parser = argparse.ArgumentParser(
        description="Process LAMMPS log/dump/in.lammps files and compute averaged energies, VDOS, vibrational entropy, and free energy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-s", "--skip", type=int, default=0, help="Skip the first N thermo lines in the LAMMPS log")
    parser.add_argument("-f", "--filename", type=str, default="log.lammps", help="LAMMPS log filename")
    parser.add_argument("-o", "--output", type=str, default="energy_vdos_results.txt", help="Summary output filename")
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
